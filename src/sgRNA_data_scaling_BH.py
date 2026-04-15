import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests

# 1. 读取第三步生成的 p 值矩阵
df = pd.read_csv("/DataM/liuranran/sgRNA_pair_test_results.csv")

# 2. 提取所有 p 值列
p_value_columns = df.columns[2:]
p_values_matrix = df[p_value_columns]

# 3. 对每一行（每个 sgRNA 对）单独执行 Benjamini-Hochberg 校正
fdr_matrix = []
for _, row in p_values_matrix.iterrows():
    # 对当前行的所有 p 值进行 BH 校正
    _, adj_p, _, _ = multipletests(row.values, method="fdr_bh")
    fdr_matrix.append(adj_p)

# 4. 构建校正后的结果数据框
fdr_df = pd.DataFrame(fdr_matrix, columns=p_value_columns)

# 合并 sg1 和 sg2 列
result_df = pd.concat([df[["sg1", "sg2"]], fdr_df], axis=1)

# 5. 保存结果
result_df.to_csv("sgRNA_pair_test_results_with_fdr.csv", index=False)
print("BH校正完成，结果已保存为 sgRNA_pair_test_results_with_fdr.csv")



# 1. 读取上一步的输出结果
result_df = pd.read_csv("sgRNA_pair_test_results_with_fdr.csv")

# 2. 提取所有潜在对照 sgRNA
potential_controls = pd.unique(result_df[['sg1', 'sg2']].values.ravel())
print(f"共有 {len(potential_controls)} 个潜在对照 sgRNA")

# 3. 定义显著性阈值
fdr_threshold = 0.05

# 4. 计算每一列（每对比较）的差异基因数
# 提取所有 FDR 列（排除 sg1, sg2）
fdr_columns = [col for col in result_df.columns if col not in ['sg1', 'sg2']]
# 对每一列统计 FDR < 阈值的基因数量
deg_per_pair = (result_df[fdr_columns] < fdr_threshold).sum(axis=0)

# 5. 按 sgRNA 分组，计算平均差异基因数
avg_deg_counts = {}
for sgRNA in potential_controls:
    # 找到所有包含该 sgRNA 的比较列
    pairs_with_sgRNA = [col for col in fdr_columns if sgRNA in col]
    if not pairs_with_sgRNA:
        avg_deg_counts[sgRNA] = 0
        continue
    # 计算这些比较的平均差异基因数
    avg_deg_counts[sgRNA] = deg_per_pair[pairs_with_sgRNA].mean()

# 6. 转换为 DataFrame 并排序
avg_deg_df = pd.DataFrame.from_dict(avg_deg_counts, orient='index', columns=['avg_deg_count'])
avg_deg_df = avg_deg_df.sort_values('avg_deg_count')
print("\n每个对照 sgRNA 的平均差异基因数（已按稳定性排序）：")
print(avg_deg_df)

# 7. 保存结果
avg_deg_df.to_csv("control_sgRNA_stability.csv")
print("\n稳定性结果已保存到 control_sgRNA_stability.csv")
