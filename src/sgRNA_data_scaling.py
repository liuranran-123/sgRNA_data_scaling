import sys
import time
import warnings
from multiprocessing import Pool, Manager
from itertools import combinations  # 生成两两组合的工具

# 第三方库
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse  
from statsmodels.stats.multitest import multipletests

# 屏蔽警告
warnings.filterwarnings('ignore')

# 路径处理与本地模块导入
sys.path.append('./shared_scripts/') # 这里的路径指向本项目依赖的算法库（由课题组维护）
from adtest_mjy import anderson, ksamp, interp # 导入高精度AD检验函数

# 1. 读取单细胞数据
adata = sc.read("./perturb-seq/hepg2/GSE264667_hepg2_raw_singlecell_01.h5ad")

# 2. 筛选带有non-阴性对照sgRNA的细胞
# 从gene_transcript列匹配含'non-'的行
neg_ctrl_cells = adata.obs[adata.obs['gene_transcript'].str.contains('non-', na=False)]
# 基于筛选结果更新adata（仅保留阴性对照细胞）
adata_neg_ctrl = adata[neg_ctrl_cells.index, :].copy()
# 计算每个基因的平均UMI计数（按细胞）
gene_mean_umi = adata_neg_ctrl.X.mean(axis=0)
# 若结果是二维数组，用flatten()转为一维
if gene_mean_umi.ndim > 1:
    gene_mean_umi = gene_mean_umi.flatten()
# 筛选平均UMI>1的基因名
high_expr_genes = adata_neg_ctrl.var_names[gene_mean_umi > 1]
# 保留这些基因的表达数据
adata_filtered = adata_neg_ctrl[:, high_expr_genes].copy()

#按基因做z标准化（减去均值再除以标准差）
def z_score_normalize(adata):
    if issparse(adata.X):
        expr = adata.X.toarray()
    else:
        expr = adata.X.copy()
    # 计算均值和标准差
    mean = expr.mean(axis=0)
    std = expr.std(axis=0)
    std[std == 0] = 1e-6  # 避免除以0
    # z标准化
    expr_norm = (expr - mean) / std
    adata.X = expr_norm
    return adata

adata_filtered = z_score_normalize(adata_filtered)


# 3. AD检验
# ad_test_single_gene函数
def ad_test_single_gene(sg1_expr, sg2_expr):
    sg1_expr = sg1_expr[~np.isnan(sg1_expr)]
    sg2_expr = sg2_expr[~np.isnan(sg2_expr)]
    p_val = np.nan
    if len(sg1_expr) > 0 and len(sg2_expr) > 0:
        try:
            # anderson_ksamp_interp函数调用
            p_val = anderson_ksamp_interp([sg1_expr, sg2_expr], midrank=True)
        except Exception as e:
            print(f"检验报错: {e}")
            pass
    return p_val

# 改写ad_test_sg_pair_apply
def ad_test_sg_pair_apply(row, sg_grouped, gene_list):
    sg1, sg2 = row['sg1'], row['sg2']
    
    # 修复：先判断sgRNA是否在分组中，再提取数据
    sg1_expr_list = []
    sg2_expr_list = []
    for g in gene_list:
        # 提取sg1的基因表达
        if sg1 in sg_grouped.groups:
            sg1_expr = sg_grouped.get_group(sg1)[g].values if g in sg_grouped.obj.columns else np.array([])
        else:
            sg1_expr = np.array([])
        # 提取sg2的基因表达
        if sg2 in sg_grouped.groups:
            sg2_expr = sg_grouped.get_group(sg2)[g].values if g in sg_grouped.obj.columns else np.array([])
        else:
            sg2_expr = np.array([])
        
        sg1_expr_list.append(sg1_expr)
        sg2_expr_list.append(sg2_expr)
    
    # 构建基因表达DataFrame
    gene_expr_df = pd.DataFrame({
        'sg1_expr': sg1_expr_list,
        'sg2_expr': sg2_expr_list
    }, index=gene_list)
    
    # 应用检验函数
    p_series = gene_expr_df.apply(
        lambda x: ad_test_single_gene(x['sg1_expr'], x['sg2_expr']),
        axis=1
    )
    # 封装结果
    result_df = pd.DataFrame({
        'p_value': [p_series.tolist()]
    }, index=[0])
    return result_df

# 改写mp_worker：修正参数传递，避免全局变量依赖
def mp_worker(pair, expr_df, gene_list):
    sg1, sg2 = pair
    # 子进程内独立分组，避免进程间数据冲突
    sg_grouped = expr_df.groupby('sgRNA', observed=False)
    
    # 构建单行DataFrame供apply调用
    pair_df = pd.DataFrame({'sg1': [sg1], 'sg2': [sg2]})
    
    # 调用检验函数
    apply_result = pair_df.apply(
        lambda row: ad_test_sg_pair_apply(row, sg_grouped, gene_list),
        axis=1
    )
    
    # 展开嵌套结果
    apply_result = pd.concat(apply_result.tolist(), ignore_index=True)
    p_series = apply_result['p_value'].iloc[0]
    
    return f"{sg1};{sg2}", p_series

# 改写load_and_preprocess：修复sg_list生成逻辑
def load_and_preprocess(adata_filtered):
    # 稀疏矩阵转密集
    expr_mat = adata_filtered.X.toarray() if sparse.issparse(adata_filtered.X) else adata_filtered.X
    expr_df = pd.DataFrame(expr_mat, index=adata_filtered.obs.index, columns=adata_filtered.var_names)
    # 修正sgRNA列赋值（确保取到正确的sgRNA字段）
    expr_df['sgRNA'] = adata_filtered.obs['gene_transcript'].values
    
    # 过滤sgRNA列的空值/无效值
    expr_df = expr_df[expr_df['sgRNA'].notna()]
    expr_df = expr_df[expr_df['sgRNA'] != '']
    
    # 重新分组，获取有效键
    sg_grouped = expr_df.groupby('sgRNA', observed=False)
    valid_sg_list = list(sg_grouped.groups.keys())  # 从分组键生成有效sgRNA列表
    
    # 生成sgRNA对
    sg_pairs = list(combinations(valid_sg_list, 2)) if len(valid_sg_list) >=2 else []
    gene_list = adata_filtered.var_names.tolist()
    
    print(f"有效sgRNA数量: {len(valid_sg_list)}")
    print(f"生成sgRNA对数量: {len(sg_pairs)}")
    print(f"数据形状: {expr_df.shape} | 基因数: {len(gene_list)}")
    
    # 返回expr_df而非分组字典，供子进程重新分组
    return expr_df, gene_list, sg_pairs, expr_df.shape

# 主执行逻辑
if __name__ == '__main__':
    n_processes = 10
    
    # 预处理数据
    expr_df, gene_list, sg_pairs, data_shape = load_and_preprocess(adata_filtered)
    
    # 测试单对sgRNA（绕过多进程，直接调用）
    if sg_pairs:
        test_result = mp_worker(sg_pairs[0], expr_df, gene_list)
        pair_name, p_series = test_result
        print(f"单对测试前5个P值: {p_series[:5]}")
    else:
        print("无有效sgRNA对，无法测试")
        exit()
    
    # 多进程执行后半部分 
    print(f"\n开始多进程处理，共{len(sg_pairs)}个sgRNA对，进程数: {n_processes}")
    start_multi = time.time()
    
    # 构造多进程入参：每个任务传入 (pair, expr_df, gene_list)
    task_args = [(pair, expr_df, gene_list) for pair in sg_pairs]
    
    # 启动进程池
    with Pool(processes=n_processes) as pool:
        # 并行执行任务，返回结果列表
        multi_results = pool.starmap(mp_worker, task_args)
    
    # 统计总耗时
    total_multi_time = time.time() - start_multi
    print(f"\n多进程处理完成，总耗时: {total_multi_time:.2f} 秒")
    
    # 结果整合
    # 提取结果到DataFrame
    result_rows = []
    for pair_name, p_vals in multi_results:
        sg1, sg2 = pair_name.split(';')
        result_rows.append({
            'sg1': sg1,
            'sg2': sg2,
            'p_values': p_vals
        })
    
    result_df = pd.DataFrame(result_rows)
    # 按基因列表匹配P值列名
    p_value_df = pd.DataFrame(result_df['p_values'].tolist(), columns=gene_list)
    final_result = pd.concat([result_df[['sg1', 'sg2']], p_value_df], axis=1)
    
    # 保存结果到CSV
    final_result.to_csv("sgRNA_pair_test_results.csv", index=False)
    print(f"\n最终结果已保存至: sgRNA_pair_test_results.csv")
    print(f"结果形状: {final_result.shape}")





