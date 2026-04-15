sgRNA_data_scaling

本项目复现了 Perturb-seq 实验数据的标准化与缩放（Scaling）预处理流程，目的是实现单细胞水平上 sgRNA 计数的精确校正。

核心算法流程

遵循以下统计学处理管线：

内部标准化 (Internal Normalization)：基于非靶向对照 (NTC) 细胞构建基准，集成 Anderson-Darling 检验 分析分布一致性，并使用 Benjamini-Hochberg (BH) 方法进行多重假设检验校正。

细胞质量控制 (QC)：计算各 gemgroup 间的缩放因子以校正测序深度差异，支持基于 UMI 总量及线粒体百分比的多指标过滤。

表达矩阵标准化：将 UMI 计数缩放至对照细胞中位数水平，并执行 gemgroup 内部的相对 Z-score 标准化。

架构说明

为了应对大规模数据运算并解决交互式环境下长时任务易中断的问题，本项目采用了计算与分析解耦的设计：

src/: 包含核心 Python 脚本，专为服务器后台离线运行优化，负责高耗时的大规模全量数据计算。

notebooks/: 存放 Jupyter Notebook，用于读取服务器生成的中间结果，执行数据二次过滤及可视化分析。

技术栈

生信分析: Scanpy, AnnData

数据科学: Pandas, NumPy, SciPy (Sparse), Statsmodels
