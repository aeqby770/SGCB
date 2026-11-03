# SGCB 0.5.9 复现说明

服务器实体根为 E:\科研项目\SGCB。当前源代码在 software/SGCB，安装库为 library_factor/SGCB，分发包为 release/SGCB_0.5.9.tar.gz。主计算入口是 experiments/45_run_current_method.R，内部调用已安装的 SGCB::sgcbDE。数值结果写入 results，图形写入 manuscript/figures。

基准和案例计算使用 16 线程、单个拟合进程。完整 API 时间与进程峰值分别记录在对应的 summary.csv 和 process_runtime.csv 中。

原始 CSV/RDS 保留拟合时的包标签、原生构建标识和输入标识；fit_package_version 是保存拟合时记录的包标签，release_version 是分发版本。原生构建标识用于追溯保存该拟合的构建记录。投稿工作簿呈现统计方法、设计、效应、区间及完整可用性，详细构建记录保留在对应 CSV/RDS 中。各端点采用的实现与输入见下文。

软件安装包与完整复现资料均由 [GitHub v0.5.9 发行页](https://github.com/aeqby770/SGCB/releases/tag/v0.5.9)提供。[完整复现归档](https://github.com/aeqby770/SGCB/releases/download/v0.5.9/SGCB_0.5.9_reproducibility.zip)包含代码、固定输入、原生竞品结果及论文结果。日常使用只需安装 R 包。

## 设置与安装

在 PowerShell 中将下列路径改为本机复现包根目录和 R 可执行文件位置。复现依赖与 SGCB 统一安装到 `library_factor`；`R_LIBS_USER` 将该库传给本会话启动的 Rcmd 和 Rscript。中文路径直接作为独立参数传入：

~~~powershell
$sgcbRoot = 'E:\科研项目\SGCB'
$sgcbR = 'D:\R\R-4.5.3\bin\x64\Rscript.exe'
$sgcbRcmd = 'D:\R\R-4.5.3\bin\x64\Rcmd.exe'
[Console]::OutputEncoding = [Text.UTF8Encoding]::new($false)
$OutputEncoding = [Console]::OutputEncoding
$env:SGCB_REVISION_ROOT = $sgcbRoot
$env:SGCB_DATA_ROOT = $sgcbRoot
$env:SGCB_LIBRARY = Join-Path $sgcbRoot 'library_factor'
[IO.Directory]::CreateDirectory($env:SGCB_LIBRARY) | Out-Null
$env:R_LIBS_USER = $env:SGCB_LIBRARY
$env:SGCB_THREADS = '16'
$env:OMP_NUM_THREADS = '16'
$env:OMP_STACKSIZE = '64M'
Remove-Item Env:LANG,Env:LC_ALL,Env:LC_CTYPE,Env:SGCB_RNG_INPUT,Env:SGCB_RNG_OUTPUT,Env:SGCB_MEAN_PARENT_OUTPUT -ErrorAction SilentlyContinue
& $sgcbR --encoding=UTF-8 -e "install.packages(c('Rcpp','RcppArmadillo','RcppEigen','matrixStats','data.table','mixsqp','ps','broom'), repos='https://cloud.r-project.org', lib=Sys.getenv('SGCB_LIBRARY'), type='binary')"
& $sgcbRcmd INSTALL ('--library=' + $env:SGCB_LIBRARY) (Join-Path $sgcbRoot 'software\SGCB')
~~~

公开安装入口为 [GitHub v0.5.9发行版](https://github.com/aeqby770/SGCB/releases/tag/v0.5.9)，[安装README](https://github.com/aeqby770/SGCB/blob/main/README.md)提供依赖、可写库目录和逐步操作。Windows R 4.5.x可直接从公开地址免编译安装：

~~~powershell
& $sgcbR --encoding=UTF-8 -e "install.packages('https://github.com/aeqby770/SGCB/releases/download/v0.5.9/SGCB_0.5.9.zip', repos=NULL, type='win.binary', lib=Sys.getenv('SGCB_LIBRARY'))"
& $sgcbR --encoding=UTF-8 -e "library(SGCB); packageVersion('SGCB'); source(system.file('examples','minimal_workflow.R',package='SGCB'))"
~~~

源码安装使用 `remotes::install_github("aeqby770/SGCB", ref="main", dependencies=NA, upgrade="never", build_vignettes=FALSE)`，或下载发行版的 `SGCB_0.5.9.tar.gz` 后运行 `install.packages(..., repos=NULL, type="source")`；源码编译需要C++17工具链，Windows使用匹配R版本的Rtools。公开Windows安装包已通过全新临时库安装和包内120基因示例。软件包README说明当前/旧接口的区别。已有正式结果可直接用于生成图与结果表。

## 在独立Windows环境中复算

先完成上面的环境设置和依赖安装，并在同一个 PowerShell 会话中运行下面的命令。`ps` 用于记录完整进程资源，`broom` 用于导出生存模型结果；两者已列入上述 CRAN 安装命令。病例和注释分析使用对应脚本列出的 Bioconductor 包，也应安装到 `SGCB_LIBRARY` 指定的库中。

数值结果和原生竞品输出已经随包提供，可先直接检查。若需重新生成图和结果表，在Python 3.12中创建根目录下的 `.venv`，安装以下依赖。这里列出服务器已有版本；未列版本号的包按该名称安装。

~~~powershell
py -3.12 -m venv (Join-Path $sgcbRoot '.venv')
& (Join-Path $sgcbRoot '.venv\Scripts\python.exe') -m pip install numpy==1.26.4 pandas==2.3.3 scipy==1.16.3 matplotlib==3.10.9 seaborn==0.13.2 statsmodels==0.14.6 scikit-learn==1.9.0 openpyxl==3.1.5 pillow==12.2.0 pymupdf==1.28.0
~~~

绘图用的 Inter、Montserrat 字体已随包提供；机制图的 Segoe Print 来自 Windows 字体目录。

条件CV²敏感性的C++实现是 `development/57_cv_moment_regression.cpp`，可执行的协变量示例是 `experiments/sle_classical_monocytes_biology.R`。KIRC与肝脏的既有敏感性拟合及样本设计随包提供。

## 六组固定基准与家族消融

输入保持原有矩阵和顺序。每个设计直接取该输入中已定义的前 n0/n1 列；运行过程使用确定性算法。

| 面板 | 固定输入 | runner 的 TYPE | n0/n1 |
|---|---|---|---|
| Positive | data/channel_simulation.rds | abundance | 10/10、10/50、50/50 |
| Observed zeros | data/zero_calibration.rds | abundance_zero | 10/10、10/50、50/50 |

例如复算含零 10:50：

~~~powershell
& $sgcbR --encoding=UTF-8 (Join-Path $sgcbRoot 'experiments\45_run_current_method.R') (Join-Path $sgcbRoot 'data\zero_calibration.rds') 'sgcb_0.5.9_zero_n10_50' 'abundance_zero' 10 50
& $sgcbR --encoding=UTF-8 (Join-Path $sgcbRoot 'experiments\123_evaluate_current_api.R') zero 10 50
& $sgcbR --encoding=UTF-8 (Join-Path $sgcbRoot 'experiments\217_cv_matched_families.R') zero 10_50
~~~

该固定输入的当前 DV 结果为 216 个真阳性、0 个有限真值误报、0 个未定义真值调用。其他五组沿用表中的输入与组数；Positive 的评价和消融面板参数为 channel。完整结果保留效应、分数、BH、95% 区间、可用性、真值和各场景指标。

## 十个空拆分和三个案例

空拆分使用 data/gtex_readcount_independent_null/split_1.rds 至 split_10.rds，TYPE 为 counts。标签为 sgcb_0.5.9_readcount_null_1 至 _10；第 3 项沿用 sgcb_0.5.9_readcount_null_3_streamed。全部完成后运行 experiments/08_recount_existing_null_results.R 汇总十个固定拆分。竞品沿用同一输入和基因家族的既有原生结果，保留在 tmp/revision_0.5.5_integrated/results，由制表脚本的 comparator_results 读取。DESeq2、edgeR、limma、DiffVar、GAMLSS、Brown–Forsythe 有十个拆分；MDSeq 有第 1 个；其余缺失记录在工作簿中。

| 案例 | 固定输入 | TYPE | 结果标签 |
|---|---|---|---|
| SLE classical monocytes | data/new_biology/sle_classical_monocytes_2021/input.rds | counts | sgcb_0.5.9_sle_classical_monocytes_2021 |
| KIRC | data/kirc_independent_discovery.rds | counts | sgcb_0.5.9_kirc |
| MASLD liver | data/new_biology/masld_liver_2023/input.rds | normalized | sgcb_0.5.9_masld_liver_2023 |

三个主拟合仍调用同一个 runner。随后按原入口处理：

- KIRC：20_kirc_patient_validation.R → 21_external_kirc_validation.R → 28_kirc_score_sensitivity_and_examples.R → 48_kirc_distribution_stability.R → 66_kirc_current_enrichment.R → kirc_variability_score.R。
- MASLD：95_population_biology_comparison.R → 80_new_biology_enrichment.R → 81_new_biology_examples.R，三者参数均为 masld_liver_2023。
- SLE：sle_classical_monocytes_biology.R。

上述脚本均在 experiments 中，读取同一结果标签、原样本映射与既有注释。条件矩模型、均值竞品和临床缓存按已有输入与归一化匹配复用。KIRC 的独立肿瘤分布与生存分析保留各自实际患者集合。SLE 保留七个预先指定基因，肝脏保留 SIAH2、ABCG8、EPHX1，主例图使用 EPHX1。

## 沿用的均值验证

均值代码和输入保持不变，SEQC、四个蛋白 spike-in、六个既有 NB 矩阵、Kang 及其配对敏感性使用保存的 0.5.7 均值结果。qPCR 与蛋白比较中的 SGCB 专属来源字段使用 sgcb_fit_version、sgcb_used_endpoint 和适用的 sgcb_reuse_reason，仅填写 method=SGCB 的记录；竞品保留其原生分数来源和输入说明。原始保存结果及其版本记录保持在 CSV/RDS 中。当前完整 API 资源图使用六组基准、十个空拆分和三个案例的 19 次调用。原始输入、原生竞品与历史诊断保留其实际名称与版本。

## 生成图与结果表

在上述会话中，运行绘图与统计汇总入口。绘图读取已提供的结果，使用 seaborn、statsmodels 和 Matplotlib；输出矢量 PDF/SVG 及 300 dpi PNG。

~~~powershell
$sgcbPython = Join-Path $sgcbRoot '.venv\Scripts\python.exe'
$env:PYTHONIOENCODING = 'utf-8'
$env:MPLBACKEND = 'Agg'
& $sgcbPython -X utf8 (Join-Path $sgcbRoot 'visualization\30_publication_figures.py') $sgcbRoot
& $sgcbPython -X utf8 (Join-Path $sgcbRoot 'visualization\82_new_publication_figures.py') $sgcbRoot
& $sgcbPython -X utf8 (Join-Path $sgcbRoot 'experiments\29_publication_tables.py') $sgcbRoot
~~~

图形写入 manuscript/figures，统计汇总写入 results，工作簿为 manuscript/Supplementary_Tables.xlsx。公开归档包括软件、分析脚本、固定输入、原生竞品结果和数值输出。
