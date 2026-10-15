# 开发与仓库文件约定

本文补充 `AGENTS.md` 中的执行规则，说明本仓库中哪些文件应作为源码维护，哪些文件应视为本地或可再生成产物。

## 文件归类

- 源码：`fastDynamicFvMesh/fastDynamicFvMesh/`、`myInterFoam/`、`myPimpleFoam/`、`myRhoPimpleFoam/` 下的 `.C/.H` 和 `Make/` 文件。
- 文档：`fastDynamicFvMesh/README.md`、`fastDynamicFvMesh/Doc/快速动网格程序手册.md` 以及本文件（`Plan.md` 等计划稿是本地草稿，已被 `.gitignore` 忽略）。
- 图表源文件：`fastDynamicFvMesh/Doc/Figures/*.drawio` 可以作为文档源文件跟踪。
- 图表导出文件：PNG/JPG 是否跟踪应按文档需要逐个决定；可再生成的临时 SVG/debug SVG 不应提交。
- 本地工具配置：`.vscode/mcp.json`、工作区 `.codex`、Python 虚拟环境和缓存不应提交。
- OpenFOAM 构建产物：`Make/linux64GccDPInt32Opt/`、`.cache/`、`lnInclude/` 不应作为源码提交。

## draw.io 图表流程

- 优先编辑 `.drawio` 源文件。
- 导出图片时保留最终文档需要引用的图片，删除中间 SVG/debug 文件。
- 如果保留导出脚本，建议放入 `tools/` 并在脚本头部说明 Python 依赖。

## smoke test 基线

烟囱测试（smoke）走**仓库内自洽**的流程：`case/Damfailure/FDM/run_pipeline.sh`。它自带网格生成与模态阶段，**不需要任何外部基线压缩包**；先 `--dry-run` 预检，再用 `--steps` 跑一个短子集，日志落在 run root 下。

## fastDynamicFvMesh 变更要求

修改 `fastDynamicFvMesh/fastDynamicFvMesh/` 下代码时，必须同步更新：

- `fastDynamicFvMesh/README.md`
- `fastDynamicFvMesh/Doc/快速动网格程序手册.md`

并记录：

- `wmake libso fastDynamicFvMesh` 编译结果。
- 至少一次短 smoke run（`case/Damfailure/FDM/run_pipeline.sh`，如 `--dry-run` + 短 `--steps`）或明确说明未能运行的原因。
