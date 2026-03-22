#!/data40T/liuhch/miniforge3/envs/pathogen_db/bin/python

import argparse
import os
import pandas as pd
import numpy as np
import subprocess
import shutil
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial.distance import squareform



def run_mash_dist(genome_files, threads=16):
    """
    计算 Mash 距离矩阵 (修复 Argument list too long 问题版)
    返回: Pandas DataFrame (index=genome, columns=genome, values=distance)
    """
    print(f"[Algorithm] Calculating Mash distances for {len(genome_files)} genomes...")
    
    # --- [修复] 使用文件列表代替命令行参数 ---
    list_file = "temp_genome_list.txt"
    with open(list_file, "w") as f:
        for path in genome_files:
            f.write(path + "\n")
            
    try:
        # 1. Sketch (使用 -l 参数读取文件列表)
        # 注意: mash sketch -l 能够直接读取包含文件路径的列表文件
        cmd_sketch = f"mash sketch -o temp_sketch -p {threads} -s 1000 -k 21 -l {list_file}"
        subprocess.run(cmd_sketch, shell=True, check=True, stderr=subprocess.DEVNULL)

        # 2. Dist
        cmd_dist = f"mash dist -p {threads} temp_sketch.msh temp_sketch.msh > dist.tsv"
        subprocess.run(cmd_dist, shell=True, check=True)

        # 3. Parse to Matrix
        df = pd.read_csv("dist.tsv", sep='\t', names=['Ref', 'Query', 'Dist', 'P', 'Hash'], header=None)
        
        # 清洗文件名
        # Mash 输出的 Ref/Query 也是全路径，需要清洗
        df['Ref'] = df['Ref'].apply(lambda x: os.path.basename(x).replace('.fna','').replace('.fasta',''))
        df['Query'] = df['Query'].apply(lambda x: os.path.basename(x).replace('.fna','').replace('.fasta',''))

        # 转为矩阵
        matrix = df.pivot(index='Ref', columns='Query', values='Dist')
        matrix = matrix.fillna(0)
        
        return matrix

    finally:
        # --- [清理] 无论成功失败都清理临时文件 ---
        for f in ["temp_sketch.msh", "dist.tsv", list_file]:
            if os.path.exists(f):
                os.remove(f)
def perform_clustering(dist_matrix, ani_threshold=0.99):
    """
    执行层次聚类
    ANI > 99% 意味着 Mash Distance < 0.01
    """
    dist_threshold = 1.0 - ani_threshold
    print(f"[Algorithm] Clustering lineages with ANI > {ani_threshold:.1%} (Dist < {dist_threshold:.3f})...")
    
    clustering = AgglomerativeClustering(
        n_clusters=None,
        distance_threshold=dist_threshold,
        metric='precomputed',
        linkage='complete' # 'complete' 保证簇内任意两点距离都小于阈值（最严谨）
    )
    
    labels = clustering.fit_predict(dist_matrix.values)
    return labels

# 把这一段替换你原来的 greedy_selection 函数
def greedy_selection(pa_matrix_df, cluster_strains, genome_map, target_coverage=0.95):
    # 1. 筛选菌株 (同你原来)
    available = [s for s in cluster_strains if s in pa_matrix_df.columns]
    if not available: return []

    # 2. 转换为 Numpy 矩阵 (性能提升的关键)
    # 这一步把 dataframe 变成纯数字矩阵，丢掉 index 和 column 的包袱
    # matrix shape: (n_genes, n_strains)
    sub_df = pa_matrix_df[available]
    matrix = sub_df.values # 0/1 numpy array
    genes = sub_df.index.values
    strains = np.array(available)
    n_genes, n_strains = matrix.shape

    # 3. 找出 Dispensable Genes (Numpy 极速版)
    # sum(axis=1) 按行求和
    gene_counts = matrix.sum(axis=1)
    # 这里的 mask 是一个布尔数组
    target_mask = (gene_counts > 0) & (gene_counts < n_strains)
    
    # 如果没有 Dispensable，直接选第一个
    if not np.any(target_mask):
        return [strains[0]]
    
    # 只保留目标基因的行
    # work_matrix shape: (n_target_genes, n_strains)
    work_matrix = matrix[target_mask, :]
    n_target_genes = work_matrix.shape[0]
    
    print(f"  -> Target Pool: {n_target_genes} dispensable genes.")

    # 4. 贪婪循环 (位运算加速)
    selected_indices = []
    # current_coverage: 布尔向量，长度为 n_target_genes
    current_coverage = np.zeros(n_target_genes, dtype=bool)
    
    while True:
        # 计算当前覆盖率
        cov_pct = np.count_nonzero(current_coverage) / n_target_genes
        if cov_pct >= target_coverage:
            break

        # 核心加速：一次性计算所有菌株的贡献
        # 逻辑：未覆盖的基因 (NOT covered) AND 菌株拥有的基因
        # work_matrix 是 (genes, strains)，转置后是 (strains, genes)
        # 用矩阵乘法或者广播机制
        
        # 找出尚未覆盖的基因索引 (mask)
        uncovered_mask = ~current_coverage
        
        # 这一步是 O(N*M) 的，但在 C 层面执行，极快
        # 计算每个菌株在 uncovered 基因里有多少个 1
        new_counts = np.sum(work_matrix[uncovered_mask, :], axis=0)
        
        # 把已经选过的菌株贡献设为 -1，防止重复选
        new_counts[selected_indices] = -1
        
        # 找最大值
        best_idx = np.argmax(new_counts)
        max_gain = new_counts[best_idx]
        
        if max_gain <= 0:
            break
            
        selected_indices.append(best_idx)
        best_strain_name = strains[best_idx]
        
        # 更新覆盖状态 (按位或运算)
        current_coverage |= (work_matrix[:, best_idx] == 1)
        
        # 打印日志
        print(f"    + Selected {best_strain_name} (+{max_gain}). Cov: {np.count_nonzero(current_coverage)/n_target_genes:.1%}")

    return list(strains[selected_indices])
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--panaroo_dir', required=True, help="Panaroo result dir containing gene_presence_absence.csv")
    parser.add_argument('--genome_dir', required=True, help="Directory containing .fna files")
    parser.add_argument('--out_dir', required=True, help="Output directory for selected genomes")
    args = parser.parse_args()

    # 1. 准备数据
    genomes = [os.path.join(args.genome_dir, f) for f in os.listdir(args.genome_dir) if f.endswith('.fna')]
    if not genomes:
        print("No genomes found!")
        return
        
    genome_map = {os.path.basename(g).replace('.fna', ''): g for g in genomes}
    genome_names = list(genome_map.keys())
    
    # 2. 读取 Panaroo 矩阵
    # 这里的 csv 很大，可能需要一点时间
    pa_file = os.path.join(args.panaroo_dir, "gene_presence_absence.csv")
    print(f"[Input] Loading Pangenome Matrix: {pa_file}")
    # Panaroo 格式: Gene, ..., Annotation, strain1, strain2...
    # 我们只读取 Gene 列和菌株列，跳过中间的 annotation info
    # 简单起见，读取全部，然后 drop 非菌株列
    df_pa = pd.read_csv(pa_file, index_col=0, low_memory=False)
    
    # 过滤列，只保留存在于 genome_names 里的列
    valid_cols = [c for c in df_pa.columns if c in genome_names]
    # 转换为 0/1 矩阵 (Panaroo 输出可能是基因ID或空)
    # 只要非空且非NaN，就是存在
    pa_matrix = df_pa[valid_cols].notna().astype(int)
    
    # 3. 亚种聚类 (Step 2.1)
    if len(genome_names) < 2:
        labels = [0]
    else:
        dist_matrix = run_mash_dist(genomes)
        # 确保矩阵顺序与 genome_names 一致
        dist_matrix = dist_matrix.reindex(index=genome_names, columns=genome_names).fillna(0)
        labels = perform_clustering(dist_matrix, ani_threshold=0.99)
    
    # 4. 基于 Cluster 的贪婪选拔 (Step 3)
    final_reps = []
    
    unique_labels = set(labels)
    print(f"\n[Selection] Found {len(unique_labels)} Clusters (Lineages). Starting Greedy Selection...")
    
    for label in unique_labels:
        # 获取该簇的所有菌株
        cluster_indices = [i for i, x in enumerate(labels) if x == label]
        cluster_strains = [genome_names[i] for i in cluster_indices]
        
        print(f"\n>>> Processing Cluster {label} (Size: {len(cluster_strains)} strains)")
        
        # 调用贪婪算法
        reps = greedy_selection(pa_matrix, cluster_strains, genome_map, target_coverage=0.95)
        final_reps.extend(reps)
        print(f"    [Result] Cluster {label} selected {len(reps)} representatives.")

    # 5. 输出结果 (Step 4)
    print(f"\n[Output] Copying {len(final_reps)} selected genomes to {args.out_dir}...")
    os.makedirs(args.out_dir, exist_ok=True)
    
    for rep in final_reps:
        src = genome_map[rep]
        dst = os.path.join(args.out_dir, os.path.basename(src))
        shutil.copy(src, dst)
        
    print("Done.")

if __name__ == "__main__":
    main()
