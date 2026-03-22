#!/bin/bash
# ==============================================================================
# 全自动串行流水线 (最终整合版)
# 功能：断点续传 + 三重QC(CheckM2/GTDB/GUNC) + 外部Python筛选 + 自动分层 + Panaroo + 最终建库
# ==============================================================================

source $(conda info --base)/etc/profile.d/conda.sh
conda activate pathogen_db
set -e
export NCBI_API_KEY="2ef67002ca93460fd13b64d06c67a0379708"

# === 1. 核心配置 ===
INPUT_LIST="taxid_counts.txt"   
FINAL_DB_DIR="04_final_db/00_Species_Repos"
MERGED_DB_DIR="04_final_db/01_Merged_DB"
INDEX_DIR="04_final_db/02_Indices"
mkdir -p "$FINAL_DB_DIR" "$MERGED_DB_DIR" "$INDEX_DIR"

# --- 数据库路径 (请确保正确) ---
CHECKM2_DB="/data40T/liuhch/OmniRef_Pipeline/Bacteria/databases/CheckM2_database/uniref100.KO.1.dmnd"
DB_GTDB="/data40T/liuhch/OmniRef_Pipeline/Bacteria/databases/gtdbtk/release226"
DB_GUNC="/data40T/liuhch/OmniRef_Pipeline/Bacteria/databases/gunc_db/gunc_db_progenomes2.1.dmnd"
DB_T2T="/data40T/liuhch/OmniRef_Pipeline/Bacteria/databases/chm13.v1.1.fa"

# --- 软件命令 ---
PROKKA_BIN="prokka"
PANAROO_BIN="/data40T/liuhch/miniforge3/envs/panaroo_env/bin/panaroo"
CHECKM2_BIN="/data40T/liuhch/miniforge3/envs/checkm2_env/bin/checkm2"

# --- 阈值设置 (传给Python脚本) ---
MIN_COMPLETENESS=90
MAX_CONTAMINATION=5
MAX_GUNC_CSS=0.45
THREADS=20
PY_PATH="/data40T/liuhch/miniforge3/envs/pathogen_db/bin/python"

# --- 工作目录设置 (使用硬盘临时目录) ---
MY_BIG_TEMP="temp_scratch"
mkdir -p "$MY_BIG_TEMP"

# >>>>> [新增] 强制重定向所有软件的临时目录到大硬盘 <<<<<
# 获取 temp_scratch 的绝对路径
ABS_TEMP_PATH=$(readlink -f "$MY_BIG_TEMP")

# 告诉 Linux 和所有软件(CheckM2, Diamond, SortMeRNA等)使用这个路径
export TMPDIR="$ABS_TEMP_PATH"
export TEMP="$ABS_TEMP_PATH"
export TMP="$ABS_TEMP_PATH"

# ==============================================================================
# 处理单个 TaxID 的函数
# ==============================================================================
process_one_taxid() {
    local taxid=$1
    local comp_count=$2
    
    # 定义工作目录
    local work_dir="$MY_BIG_TEMP/temp_run_${taxid}"
    
    echo "----------------------------------------------------------------"
    echo ">>> [Start/Resume] TaxID: $taxid | WorkDir: $work_dir"

    mkdir -p "$work_dir/raw" "$work_dir/clean" "$work_dir/qc" "$work_dir/passed" "$work_dir/stratified"




    # === [Step 1] 智能下载 (脱水+注水版) ===
    if [ -f "$work_dir/.step1_download.done" ]; then
        echo "    [Step 1] Download already done. Skipping."
    else
        echo "    [Step 1] Downloading (Dehydrated Mode)..."
        rm -rf "$work_dir/raw/"* "$work_dir/raw_extract"
        mkdir -p "$work_dir/raw" "$work_dir/raw_extract"

        local zip_file="$work_dir/data.zip"
        local download_status=0

        # ### 修改点 1: 在基础命令中加入 --dehydrated 参数 ###
        # 这样只会下载几KB的元数据，速度极快，不会断连
        local base_cmd="datasets download genome taxon $taxid --include genome,seq-report --filename $zip_file --dehydrated"

        if [ "$comp_count" -gt 0 ]; then
            echo "        [Strategy] Count > 0 ($comp_count). Mode: COMPLETE."
            $base_cmd --assembly-level complete --assembly-source all
            download_status=$?
        else
            echo "        [Strategy] Count = 0. Mode: REFSEQ (Fallback)."
            $base_cmd --assembly-source refseq
            download_status=$?
        fi

        # 1. 检查下载命令 (下载元数据应该瞬间完成)
        if [ $download_status -ne 0 ]; then
            echo "    [Error] datasets CLI (dehydrated) failed with exit code $download_status."
            rm -rf "$work_dir"
            return 1
        fi

        # 2. 检查 ZIP 文件有效性
        if [ ! -s "$zip_file" ]; then
            echo "    [Error] Zip file is empty or does not exist."
            rm -rf "$work_dir"
            return 1
        fi

        # 3. 解压 (这里解压出来的只是 fetch.txt 和 json 元数据)
        echo "        [Action] Unzipping metadata..."
        unzip -o "$zip_file" -d "$work_dir/raw_extract" > "$work_dir/unzip.log" 2>&1
        if [ $? -ne 0 ]; then
            echo "    [Error] Unzip failed. Check $work_dir/unzip.log"
            cat "$work_dir/unzip.log"
            return 1
        fi

        # ### 修改点 2: 增加 Rehydrate 步骤 (真正下载数据) ###
        echo "        [Action] Rehydrating (Fetching actual sequences)..."
        # --max-workers 4 可以并行下载，提高速度
        datasets rehydrate --directory "$work_dir/raw_extract" --max-workers 4 > "$work_dir/rehydrate.log" 2>&1
        
        # 检查 rehydrate 是否成功
        if [ $? -ne 0 ]; then
             echo "    [Error] Rehydration failed. Check $work_dir/rehydrate.log"
             # 可能是网络问题，可以打印日志看看
             cat "$work_dir/rehydrate.log"
             return 1
        fi

        # 4. 移动文件 (这一步逻辑不变，因为 rehydrate 会把文件放回原来的结构中)
        find "$work_dir/raw_extract" -name "*.fna" -exec mv {} "$work_dir/raw/" \; 2>/dev/null

        local count=$(ls "$work_dir/raw/"*.fna 2>/dev/null | wc -l)
        echo "    Downloaded: $count files."

        if [ "$count" -eq 0 ]; then 
            echo "    [Error] No .fna files found after extraction. Skipping."
            # 保留目录以便排查
            return 1
        fi

        touch "$work_dir/.step1_download.done"
        rm -f "$zip_file" 
        rm -rf "$work_dir/raw_extract"
    fi















    # === [Step 2] 并行去宿主 ===
    if [ -f "$work_dir/.step2_dehost.done" ]; then
        echo "    [Step 2] De-hosting already done. Skipping."
    else
        echo "    [Step 2] De-hosting (Parallel)..."
        if [ ! -f "human_t2t.mmi" ]; then minimap2 -d "human_t2t.mmi" "$DB_T2T"; fi
        local MAX_JOBS=5
        local CPU_PER_JOB=4

        for f in "$work_dir/raw/"*.fna; do
            base=$(basename "$f" .fna)
            if [ ! -s "$work_dir/clean/$base.fna" ]; then
                (
                    minimap2 -a -x asm5 -t $CPU_PER_JOB "human_t2t.mmi" "$f" 2>/dev/null \
                    | samtools view -@ 2 -f 4 -b - \
                    | samtools fasta - \
                    | seqkit seq -m 500 - > "$work_dir/clean/$base.fna"
                ) &
                # 兼容老系统的 wait
                while [[ $(jobs -r -p | wc -l) -ge $MAX_JOBS ]]; do sleep 1; done
            fi
        done
        wait
        touch "$work_dir/.step2_dehost.done"
    fi

    # === [Step 3] 三重质量门控 (Trinity QC) ===
    if [ -f "$work_dir/.step3_qc.done" ]; then
        echo "    [Step 3] QC already done. Skipping."
    else
        echo "    [Step 3] Trinity QC (CheckM2 + GTDB + GUNC)..."
        mkdir -p "$work_dir/qc/checkm2" "$work_dir/qc/gtdb" "$work_dir/qc/gunc"

        # 3.1 CheckM2
        if [ ! -f "$work_dir/qc/checkm2/quality_report.tsv" ]; then
            echo "        Running CheckM2..."
            $CHECKM2_BIN predict -t $THREADS -i "$work_dir/clean" -x fna \
                -o "$work_dir/qc/checkm2" --database_path "$CHECKM2_DB" --quiet --force
        fi

        # 3.2 GTDB-Tk
        if [ ! -f "$work_dir/qc/gtdb/gtdbtk.bac120.summary.tsv" ] && [ ! -f "$work_dir/qc/gtdb/gtdbtk.ar53.summary.tsv" ]; then
            echo "        Running GTDB-Tk..."
            export GTDBTK_DATA_PATH="$DB_GTDB"
            gtdbtk classify_wf --cpus $THREADS --genome_dir "$work_dir/clean" --extension fna \
                --out_dir "$work_dir/qc/gtdb" --skip_ani_screen > "$work_dir/qc/gtdb/gtdbtk.log" 2>&1
        fi

        # 3.3 GUNC
        if [ ! -f "$work_dir/qc/gunc/GUNC.progenomes_2.1.maxCSS_level.tsv" ]; then
            echo "        Running GUNC..."
            gunc run --input_dir "$work_dir/clean" --file_suffix .fna -r "$DB_GUNC" \
                --threads $THREADS --out_dir "$work_dir/qc/gunc" > "$work_dir/qc/gunc/gunc.log" 2>&1
        fi

        # 3.4 综合筛选 (调用 scripts/filter_genomes.py)
        echo "        Filtering Genomes (External Python Script)..."
        
        # 准备 GTDB 汇总文件 (处理细菌/古菌并存的情况)
        GTDB_SUMMARY="$work_dir/qc/gtdb/gtdbtk_summary_merged.tsv"
        
        if [ -f "$work_dir/qc/gtdb/gtdbtk.bac120.summary.tsv" ]; then
            cat "$work_dir/qc/gtdb/gtdbtk.bac120.summary.tsv" > "$GTDB_SUMMARY"
            if [ -f "$work_dir/qc/gtdb/gtdbtk.ar53.summary.tsv" ]; then
                tail -n +2 "$work_dir/qc/gtdb/gtdbtk.ar53.summary.tsv" >> "$GTDB_SUMMARY"
            fi
        elif [ -f "$work_dir/qc/gtdb/gtdbtk.ar53.summary.tsv" ]; then
            cat "$work_dir/qc/gtdb/gtdbtk.ar53.summary.tsv" > "$GTDB_SUMMARY"
        else
            touch "$GTDB_SUMMARY"
        fi

        # 调用你的 Python 脚本
        $PY_PATH scripts/filter_genomes.py \
            --checkm "$work_dir/qc/checkm2/quality_report.tsv" \
            --gtdb "$GTDB_SUMMARY" \
            --gunc "$work_dir/qc/gunc/GUNC.progenomes_2.1.maxCSS_level.tsv" \
            --min_comp $MIN_COMPLETENESS \
            --max_cont $MAX_CONTAMINATION \
            --max_css $MAX_GUNC_CSS \
            --out_dir "$work_dir/passed" \
            --source_dir "$work_dir/clean"

        # 检查结果
        local pass_count=$(ls "$work_dir/passed/"*.fna 2>/dev/null | wc -l)
        echo "    $pass_count genomes passed Trinity QC."
        
        if [ "$pass_count" -eq 0 ]; then
            echo "    [Stop] No genomes passed QC. Skipping."
            rm -rf "$work_dir"
            return
        fi
        
        touch "$work_dir/.step3_qc.done"
    fi

# === [Step 3.5] 物种分层 (符合生物学命名法则的智能版) ===
    if [ -f "$work_dir/.step3_5_strat.done" ]; then
        echo "    [Step 3.5] Stratification already done."
    else
        echo "    [Step 3.5] Stratifying species (Smart Biological Grouping)..."
        local jsonl="$work_dir/raw_extract/ncbi_dataset/data/assembly_data_report.jsonl"
        
        python3 -c "
import json, os, shutil, re
jsonl_path = '$jsonl'
passed_dir = '$work_dir/passed'
strat_dir = '$work_dir/stratified'
acc2species = {}

if os.path.exists(jsonl_path):
    with open(jsonl_path) as f:
        for line in f:
            try:
                data = json.loads(line)
                acc = data['accession']
                full_name = data['organism']['organismName']
                
                # --- 智能命名清洗逻辑 ---
                # 移除方括号 (常见于 [Candaditus] ...)
                clean_name = full_name.replace('[', '').replace(']', '')
                parts = clean_name.split()
                
                # 1. 确定基础长度 (Genus + Species)
                # 如果是 Candidatus 开头，基础是3个词，否则是2个词
                if len(parts) > 0 and parts[0] == 'Candidatus':
                    base_idx = 3
                else:
                    base_idx = 2
                
                # 2. 检查是否有亚种 (subsp. / var.)
                # 如果 基础长度 后面紧跟 subsp. 或 var.，则保留后面一个词作为亚种名
                final_len = base_idx
                if len(parts) > base_idx:
                    indicator = parts[base_idx].lower() # 获取下一个词
                    if indicator in ['subsp.', 'var.'] and len(parts) > base_idx + 1:
                        final_len = base_idx + 2 # 保留 'subsp.' 和 'name'
                
                # 3. 截取有效部分
                # 防止名字本身就很短 (比如 sp. 还没定种)
                final_len = min(final_len, len(parts))
                valid_parts = parts[:final_len]
                
                # 4. 拼接并清理非法字符
                short_name = '_'.join(valid_parts)
                sp = re.sub(r'[^a-zA-Z0-9_]', '', short_name)
                # -----------------------
                
                acc2species[acc] = sp
            except: pass

# 分发文件
for f in os.listdir(passed_dir):
    if not f.endswith('.fna'): continue
    acc = f.replace('.fna', '') 
    
    sp_name = 'Unclassified'
    if acc in acc2species:
        sp_name = acc2species[acc]
    else:
        # 模糊匹配
        for k, v in acc2species.items():
            if k in acc or acc in k:
                 sp_name = v
                 break
    
    target_subdir = os.path.join(strat_dir, sp_name)
    os.makedirs(target_subdir, exist_ok=True)
    shutil.copy2(os.path.join(passed_dir, f), os.path.join(target_subdir, f))
"
        touch "$work_dir/.step3_5_strat.done"
    fi

    # === [Step 4] 遍历分析每个物种 ===
    if [ -f "$work_dir/.step4_analysis.done" ]; then
         echo "    [Step 4] Analysis already done."
    else
        for species_path in "$work_dir/stratified/"*; do
            if [ ! -d "$species_path" ]; then continue; fi
            local sp_name=$(basename "$species_path")
            
            local sp_count=$(ls "$species_path"/*.fna 2>/dev/null | wc -l)
            echo "    [Step 4] Processing: $sp_name (Count: $sp_count)"
            
            local sp_work_dir="$work_dir/analysis_$sp_name"
            mkdir -p "$sp_work_dir/final_selection" "$sp_work_dir/prokka" "$sp_work_dir/panaroo"

            if [ "$sp_count" -lt 5 ]; then
                 cp "$species_path"/*.fna "$sp_work_dir/final_selection/"
            else
                 # 4.1 Prokka
                 local MAX_JOBS=5
                 for genome in "$species_path"/*.fna; do
                     base=$(basename "$genome" .fna)
                     if [ ! -f "$sp_work_dir/prokka/$base/$base.gff" ]; then
                         (
                             $PROKKA_BIN --outdir "$sp_work_dir/prokka/$base" --prefix "$base" --locustag "$base" \
                                 --cpus 4 --fast --force "$genome" > "$sp_work_dir/prokka/$base.log" 2>&1
                         ) &
                         while [[ $(jobs -r -p | wc -l) -ge $MAX_JOBS ]]; do sleep 1; done
                     fi
                 done
                 wait

		 # 4.2 Panaroo (极速版: -a core --aligner none)
                 # [修改点1] 删掉了 export PATH 这一行，防止后面报错
                 if [ ! -f "$sp_work_dir/panaroo/gene_presence_absence.csv" ]; then
                     set +e
                     # [修改点2] 在下面这行命令的最开头，加上了临时环境变量设置
                     PATH="/data40T/liuhch/miniforge3/envs/panaroo_env/bin:$PATH" $PANAROO_BIN -i "$sp_work_dir/prokka/"*/*.gff -o "$sp_work_dir/panaroo" \
                         -t $THREADS -a core --aligner none --clean-mode strict --remove-invalid-genes > "$sp_work_dir/panaroo/panaroo.log" 2>&1
                     set -e
                 fi                 
                 # 4.3 Selector
                 if [ -f "$sp_work_dir/panaroo/gene_presence_absence.csv" ]; then
                     $PY_PATH scripts/greedy_selector.py --panaroo_dir "$sp_work_dir/panaroo" \
                         --genome_dir "$species_path" --out_dir "$sp_work_dir/final_selection"
                 else
                     echo "        [Warning] Panaroo failed for $sp_name, keeping all."
                     cp "$species_path"/*.fna "$sp_work_dir/final_selection/"
                 fi
            fi

            # === [Step 5] 归档当前物种 ===
            local target_repo="$FINAL_DB_DIR/$sp_name"
            mkdir -p "$target_repo"
            cp -u "$sp_work_dir/final_selection/"*.fna "$target_repo/" 2>/dev/null
            if [ -f "$sp_work_dir/panaroo/gene_presence_absence.csv" ]; then
                cp "$sp_work_dir/panaroo/gene_presence_absence.csv" "$target_repo/panaroo_matrix.csv"
            fi
            # 复制 QC 报告备查
            cp "$work_dir/qc/checkm2/quality_report.tsv" "$target_repo/checkm2_report.tsv" 2>/dev/null
        done
        touch "$work_dir/.step4_analysis.done"
    fi


# === [Step 6] 智能清理 (滑动窗口 + 垃圾回收) ===
    echo "  [Step 6] Smart Cleaning..."

    local current_pwd=$(pwd)
    
    # 进到大临时目录
    cd "$MY_BIG_TEMP" || true

    # 1. [新增] 清理 Python 多进程残留垃圾 (pymp-*) 和缓存
    # 这一行非常重要，否则硬盘会爆
    rm -rf pymp-* __pycache__

    # 2. 清理旧的物种文件夹，只保留最近的 2 个
    ls -td temp_run_* 2>/dev/null | tail -n +3 | xargs -r rm -rf

    cd "$current_pwd" || true
    echo "  [Info] Cleaned junk files and old runs."
}


# ==============================================================================
# 主循环 (防抢占版 - 使用文件描述符 3)
# ==============================================================================
DONE_FILE="finished_list.txt"  # 定义账本文件

echo ">>> Pipeline Started. Reading $INPUT_LIST..."

# [修改点] 注意这里的 while read -u 3 ... 3< ...
# 这就是为了防止 CheckM2 或 datasets 偷吃输入文件
while read -u 3 taxid count; do
    taxid=$(echo "$taxid" | tr -d '\r')
    count=$(echo "$count" | tr -d '\r')
    
    if [[ "$taxid" =~ ^[0-9]+$ ]]; then 
        # [检查点] 如果账本里有这个号，直接跳过
        if grep -q "^${taxid}$" "$DONE_FILE" 2>/dev/null; then
            echo ">>> [Skip] TaxID $taxid is found in $DONE_FILE. Skipping."
            continue
        fi

        if [ -z "$count" ]; then count=0; fi
        
        # 执行核心函数
        process_one_taxid "$taxid" "$count"
        
        # [记账] 跑完一个，立马记下来
        echo "$taxid" >> "$DONE_FILE"
    fi
done 3< "$INPUT_LIST"

# ==============================================================================
# [Post-Processing] 最终建库模块
# ==============================================================================
# ⚠️ 暂时注释掉最后这一步！
# 等你的 500 个物种全跑完了，再把下面的注释解开。
# 现在不要让它出来混淆视听。

echo ">>> [Info] Batch processing finished. Please check finished_list.txt."
echo ">>> [Info] To build final database, uncomment Module 5 in the script."

# echo ">>> [Module 5] Building Final Database..."
# MERGED_FNA="$MERGED_DB_DIR/all_bacteria_reps.fna"
# find "$FINAL_DB_DIR" -name "*.fna" -print0 | xargs -0 cat > "$MERGED_FNA"
# if [ ! -f "$INDEX_DIR/ref_db.1.bt2" ]; then
#     echo "    Building Bowtie2 Index..."
#     bowtie2-build --threads $THREADS "$MERGED_FNA" "$INDEX_DIR/ref_db"
# fi
# echo ">>> All Done! System Ready."





