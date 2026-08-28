# DamFailure.cas.h5 算例设置记录

**来源文件:** `/root/Workspace/tmp/DamFailure.cas.h5`（5.7 MB）
**解析方式:** Python `h5py` 直接读取 HDF5 结构；设置位于 `settings/Rampant Variables`、`settings/Thread Variables` 等 Scheme 格式文本块中
**解析日期:** 2026-08-17

## 基本信息

| 项        | 值                                                                                    |
| --------- | ------------------------------------------------------------------------------------- |
| 软件      | ANSYS Fluent **22.1** (Build 10213)                                                   |
| 格式      | HDF5 case（网格内嵌，无 .dat）                                                        |
| 维度/精度 | 3D，双精度                                                                            |
| 运行状态  | flow-time = 0，time-step = 0 — **从未运行的初始算例**                                 |
| 原始路径  | `D:/Special/202311_FastDynamicMesh/`（本仓库 `fastDynamicFvMesh` 的 Fluent 原型算例） |

## 网格

- 单流体域 `fluiddomain`：**181,200 单元**（cellType 4，六面体类）、**200,109 节点**、562,204 面
- 计算域：**0.584 × 0.365 × 0.012 m**，z 向仅 12 mm，准二维水槽
- 经典溃坝冲击弹性障碍 FSI 基准构型（Walhorn 类型）

### 边界面区（包围盒由节点坐标实算确认）

| 区域          | ID  |                              位置 |    面数 | 类型                                                        |
| ------------- | --- | --------------------------------: | ------: | ----------------------------------------------------------- |
| front         | 13  |                  z = 0.012 大侧面 |  15,100 | 对称边界（symmetryPlane）                                   |
| back          | 17  |                      z = 0 大侧面 |  15,100 | 对称边界（symmetryPlane）                                   |
| top           | 14  |                         y = 0.365 |   1,920 | 静止无滑移壁面                                              |
| bottom_fluid  | 18  |                             y = 0 |   1,776 | 静止无滑移壁面                                              |
| left          | 16  |                             x = 0 |   1,164 | 静止无滑移壁面                                              |
| right         | 15  |                         x = 0.584 |   1,164 | 静止无滑移壁面                                              |
| **fsi_fluid** | 19  | x ∈ [0.292, 0.304]，y ∈ [0, 0.08] | **984** | 竖直弹性薄板（12 mm 厚 × 80 mm 高，立于底部中央），静止壁面 |

## 物理模型

| 模型   | 设置                                                                                                                              |
| ------ | --------------------------------------------------------------------------------------------------------------------------------- |
| 多相流 | **VOF 两相**：phase-1 = air，phase-2 = water-liquid；无表面张力、无壁面粘附                                                       |
| 材料   | air（ρ = 1.225 kg/m³，μ = 1.7894e-5 Pa·s）；water-liquid（ρ = 998.2 kg/m³，μ = 1.003e-3 Pa·s）；壁面材料 aluminum（仅热属性占位） |
| 重力   | (0, **-9.81**, 0) m/s²，已启用                                                                                                    |
| 湍流   | **k-ω SST**（`kw/solve? #t`，SST 常数在位，曲率修正关）                                                                           |
| 能量   | 关闭                                                                                                                              |
| 其他   | 无 System Coupling、无 6DOF、无本征结构求解器（intrinsic FSI）                                                                    |

## 结构（弹性挡板）材料属性

结构材料参数（对应论文表 3-2 固有频率 6.02 Hz 的唯一正确版本，来源 `damfailure_fine_symm/mode/system.in`）：

| 参数       | 值                  |
| ---------- | ------------------- |
| 弹性模量 E | **1 MPa**（1e6 Pa） |
| 泊松比 ν   | **0**               |
| 密度 ρ     | **2500 kg/m³**      |

注：论文正文写"弹性模量 10 MPa"为笔误（10 MPa 对应一阶约 19 Hz，与论文表 3-2 的 6.02 Hz 不符）。

## 求解设置

| 项                      | 值                                                                                                                         |
| ----------------------- | -------------------------------------------------------------------------------------------------------------------------- |
| 时间类型                | 瞬态                                                                                                                       |
| 时间步长                | **0.001 s**                                                                                                                |
| 每步最大迭代            | 50                                                                                                                         |
| 压力-速度耦合           | **PISO**（`piso/coupling? #t`，1 次 neighbor 修正）                                                                        |
| 时间格式                | 二阶隐式                                                                                                                   |
| 离散格式（Fluent 编号） | pressure/scheme 14；mom/scheme 1（二阶迎风）；k/omega/scheme 1；mp/scheme 16（编号与格式名称的映射以 Fluent GUI 显示为准） |

- 无 report 定义、无监控器
- 有一个每时间步保存的 `phase-2-vof` 云图动画（`animation-1`，输出到 `D:/.../animation/DamFailure`）
- 初始水柱分布在 .dat 中定义，本 case 文件内无 patch 区域记录

## 动网格 / FSI（UDF 实现）

这是该算例最关键的部分，FSI 完全由编译型 UDF `libudf` 实现：

| 挂钩               | UDF                          |
| ------------------ | ---------------------------- |
| 网格运动（动态区） | `FDM_method::libudf`         |
| execute-at-end     | `Set_next_time_step::libudf` |
| execute-at-exit    | `Finish_process::libudf`     |

- 申请 **31 个 node UDM**：在每个节点上存储模态坐标/速度/位移等状态，与本仓库 `fastDynamicFvMesh` 的模态叠加（modal superposition）方法完全对应
- 动态区条目：deforming 类型，`slide/remesh/layer` 均标记启用

### 两个疑点（以 Fluent 打开后 GUI 显示为准）

1. **动态区挂在线程 17 = `back`（z=0 大侧面），而非 `fsi_fluid` 板面。** 可能是 UDF 内部自行定位 FSI 节点，也可能是该条目配置不完整。
2. **顶层 `dynamesh/method` 的 smoothing / remeshing / layering 三个开关均为 `#f`**，与动态区条目中的标记矛盾。
3. 顶层 `rp-unsteady? #f`、`sg-mphase? #f` 与 `case-config`（unsteady + VOF + dynmesh）矛盾 — 保存时 rpvar 顶层标志未刷新，**实际生效配置以 case-config / GUI 为准**。

## 脚本流水线总览

`case/Damfailure/` 下现有两条从零开始的流水线，互不依赖：**FDM**（fastDynamicFvMesh 内建 FSI）与 **preCICE**（外部分区耦合）。两者共享同源的网格生成脚本（01），但组装、运行、视频导出各自独立。脚本已从本目录移入 `FDM/` 与 `preCICE/` 两个子目录：

```
case/Damfailure/
├── CaseDescription.md                  # 本文档
├── FDM/                                # FDM（fastDynamicFvMesh 内建 FSI）流水线
│   ├── 01_generate_mesh.py              # 网格生成（blockMesh 流体 + CalculiX 固体 + VTK + --render 渲染）
│   ├── 02_run_modal_virtual_elastic.py  # 模态分析 + 导出 mode/
│   ├── 03_assemble_fdm_case.py          # 组装 FDM 算例（默认 8 核）
│   ├── 04_run_fdm.py                    # 运行（串行/并行，不重建）
│   ├── 05_postprocess.py                # 并行重建 + 位移/性能报告
│   ├── 06_export_fdm_video.py           # 视频导出（pvpython）
│   └── run_pipeline.sh                  # 流水线驱动入口
└── preCICE/                             # preCICE 双求解器外部耦合流水线
    ├── 01_generate_mesh.py              # 网格生成（同 FDM 版）
    ├── 03_assemble_precice_case.py      # 组装 preCICE 算例（fluid-openfoam + solid-calculix + precice-config.xml）
    ├── 04_run_precice.py                # 双求解器运行（interFoam + ccx_preCICE）
    ├── 05_export_precice_video.py       # 视频导出
    ├── run_modal_analysis.py            # 独立模态验证（可选）
    └── run_pipeline.sh                  # 流水线驱动入口
```

## FDM 流水线（fastDynamicFvMesh 内建 FSI）

从零生成并运行 FDM（fastDynamicFvMesh）算例的完整流程，**不依赖任何已有算例**。FSI 由 `fastDynamicFvMesh` 内建实现（模态叠加），整个耦合只需单个 `myInterFoam` 进程：

| 脚本                              | 功能                                                                                     |
| --------------------------------- | ---------------------------------------------------------------------------------------- |
| `01_generate_mesh.py`             | 生成 OpenFOAM 流体网格（blockMesh）+ CalculiX 固体网格（dam_gate.inp）+ VTK 可视化产物   |
| `02_run_modal_virtual_elastic.py` | 挡板+虚拟弹性体组合模态分析（CalculiX ccx），导出 FDM 模态文件（mode/）                  |
| `03_assemble_fdm_case.py`         | 自包含组装 FDM 算例（初始/边界条件、仿真参数全部内嵌，只从输入复制 polyMesh + mode/）    |
| `04_run_fdm.py`                   | 运行 FDM（myInterFoam + fastDynamicFvMesh），支持 fresh/restart/continue、并行、时间估算；不重建 |
| `05_postprocess.py`               | 重建并行结果，提取挡板右上角位移并生成性能报告；`--self-test` 运行内置回归测试 |
| `06_export_fdm_video.py`          | FDM 结果导出视频（ParaView pvpython 逐帧 PNG，ffmpeg 合成 MP4）                          |

### 运行命令

推荐从仓库根目录经 `run_pipeline.sh` 驱动（`--steps` 可选步骤，默认全流程；`--dry-run` 先预览）：

```bash
cd /root/OpenFOAM/user-v2412/src
bash case/Damfailure/FDM/run_pipeline.sh \
    --case-root /root/Workspace/<root> --steps 1,2,3,4,5,6 \
    --res 115 --z-layers 12 --end-time 0.5 --nprocs 8 --video
```

或单独调用各脚本（路径前缀 `case/Damfailure/FDM/`）：

```bash
# 1. 网格生成（--res 控制密度，--z-layers 默认 12；--view 导出 VTK，--render 生成网格图）
python3 case/Damfailure/FDM/01_generate_mesh.py --res 115 --z-layers 12 \
    --out-dir /root/Workspace/<mesh> --view --render

# 2. 模态分析 + 导出（--export-mode 导出 mode/ 文件）
python3 case/Damfailure/FDM/02_run_modal_virtual_elastic.py \
    --fluid-dir /root/Workspace/<mesh>/fluid-openfoam \
    --out-dir /root/Workspace/<modal> --export-mode

# 3. 组装算例（默认 endTime 0.5, deltaT 0.001, writeInterval 50, 8 核）
python3 case/Damfailure/FDM/03_assemble_fdm_case.py \
    --mesh-dir /root/Workspace/<mesh>/fluid-openfoam \
    --mode-dir /root/Workspace/<modal>/mode \
    --out-dir /root/Workspace/<fdm>

# 4. 运行（--dry-run 先预览；-f 全新 / -c 续算 / --nprocs 并行）
python3 case/Damfailure/FDM/04_run_fdm.py -f --case /root/Workspace/<fdm> --end-time 0.5

# 5. 并行重建 + 位移/性能报告
python3 case/Damfailure/FDM/05_postprocess.py \
    --case /root/Workspace/<fdm> --out-dir /root/Workspace/<postprocess>

# 6. 导出视频（PNG 帧 → MP4；xvfb 提供无头渲染 DISPLAY）
xvfb-run -a pvpython case/Damfailure/FDM/06_export_fdm_video.py \
    --case /root/Workspace/<fdm> --out /root/Workspace/<frames>
ffmpeg -y -framerate 10 -i /root/Workspace/<frames>/frame_%04d.png \
    -c:v libx264 -pix_fmt yuv420p /root/Workspace/fdm_result.mp4
```

### 关键默认参数（内嵌在脚本中）

| 项                     | 默认值                                                   |
| ---------------------- | -------------------------------------------------------- |
| 时间步长 deltaT        | 0.001 s                                                  |
| 仿真时长 endTime       | 0.5 s（论文/参考算例时长）                               |
| 输出间隔 writeInterval | 50 步（0.05 s）                                          |
| 并行 nprocs            | **8 核**（--nprocs 默认 8）                              |
| 挡板材料               | **E=1 MPa, ν=0, ρ=2500 kg/m³**（论文基准）               |
| 流体物性               | 水 ρ=998.2/ν=1.0048e-6，空气 ρ=1.225/ν=1.4607e-5，σ=0.07 |
| 湍流                   | kOmegaSST                                                |
| FDM 动网格             | theta=1.4, fsiPatches=FSI_FLUID, pRef=101325             |
| 网格质量               | **MAX_EXPANSION=4.0**（相邻单元尺寸比上限，可用 --max-expansion 调） |
| 水柱初始               | box (0,0,0)-(0.146,0.292,0.012) 内 alpha.water=1         |
| 视频                   | 含 t=0 初始帧（--include-t0 默认开）                     |

### 新功能说明

- **网格质量保证**（01）：新增 `--max-expansion`（默认 4.0），限制每个方向相邻单元尺寸比 ≤4，首层厚度在必要时放宽。这保证任何 res 下网格质量（之前 10.2 的扩张比降到 4.0）。
- **网格渲染**（01）：`--render` 生成 `mesh_plot_xy.png` + `mesh_dx_dist.png` 两张网格图（x-y 截面 + x 向尺寸分布）。
- **可选 AMR**（03）：`--amr` 启用 dynamicRefineFvMeshCoeffs（refineInterval 默认 1、lowerRefineLevel 20、maxRefinement 1、maxCells 500000），以 fixedDt_refineEveryStep_band 为蓝本。
- **视频从 t0 开始**（05）：OpenFOAMReader 默认跳过 time=0，05 用 `SkipZeroTime=0` 让视频含初始水柱帧（默认开，`--no-t0` 关）。
- **输出过滤**（run_pipeline.sh）：控制台只显示带 `[FDM` 前缀的主要行，工具原始输出只落盘 `run_step_<n>.log`。
- **GLX 修复**（05）：`os._exit(0)` 绕过 pvpython 在 xvfb 下的 GLXBadContext 干净退出。

### 视频导出要点（06_export_fdm_video.py）

- 运行方式：`xvfb-run -a pvpython`（Ubuntu ParaView 是 X11 版，无头环境需 xvfb 提供 DISPLAY）
- 默认样式：白底 + 蓝色水相云图（alpha.water）+ 黑色网格线（`--no-wireframe` 关闭线框）
- 参数：`--case`、`--out`、`--res`（默认 1280x720）、`--field`（默认 alpha.water）、`--wireframe`/`--no-wireframe`、`--wire-color`、`--include-t0`/`--no-t0`
- 边线实现：默认用 `ExtractSurface` + `Surface With Edges`；若 xvfb/OSMesa 下仍出现 AMR 表面三角化斜纹，使用 `--no-wireframe` 导出纯 `Surface` 画面
- AMR 多面体：导出脚本关闭 OpenFOAMReader 的 `Decomposepolyhedra`，防止粗细过渡多面体被拆成 pyramid/tetra 后出现三角面内部对角线
- 位移后处理（06）：AMR 改变网格点数/编号时，从各时步 `FSI_FLUID` patch 的 z=0 前缘点中定位顶部候选带内 x 最大点（`--node-tol` 为候选带厚度，默认 5e-4 m），不在普通流体点中回退；即使位移提取无有效行，默认仍生成包含 N/A 说明的性能报告
- 默认从 t=0 开始（--include-t0 默认开）；GLXBadContext 已修复（`os._exit(0)` 干净退出）
- 输出过滤：经 run_pipeline.sh 运行时控制台只显 `[FDM` 前缀行
- 前置依赖：`xvfb`、`pvpython`（paraview + python3-paraview）、`ffmpeg` 已安装

### 边界层设计（01_generate_mesh.py）

- 11 块结构化网格，x 方向 5 切分（0/0.146/0.292/0.304/0.584）
- 块 3/7/11（挡板右~RIGHT 两端都是壁面）用 blockMesh **多段 grading** 同时加密两端
- 所有壁面（BOTTOM/TOP/LEFT/RIGHT/挡板三面）**第一层厚度统一**（`--bl-thickness` 控制，默认 0.5×最小平均单元尺寸）

## preCICE 流水线（外部分区耦合）

preCICE 与 FDM 的**本质区别**在于耦合机制：

- **FDM**：`fastDynamicFvMesh` 内建 FSI（模态叠加），单 myInterFoam 进程
- **preCICE**：`dynamicMotionSolverFvMesh` + `displacementLaplacian`，通过 preCICE 3.2 库**外部耦合**两个求解器：
  - 流体 = 原版 `interFoam`（preciceDict 驱动 preCICEAdapter）
  - 固体 = `ccx_preCICE -i dam_gate -precice-participant Solid`（CalculiX）

脚本位于 `preCICE/` 子目录：

| 脚本                          | 功能                                                                                     |
| ----------------------------- | ---------------------------------------------------------------------------------------- |
| `01_generate_mesh.py`         | 网格生成（同 FDM 版，流体 + 固体网格 + VTK）                                              |
| `03_assemble_precice_case.py` | 组装 preCICE 算例（fluid-openfoam + solid-calculix + precice-config.xml）                |
| `04_run_precice.py`           | 双求解器运行（interFoam + ccx_preCICE），支持 fresh/restart/continue、并行、时间估算     |
| `05_export_precice_video.py`  | preCICE 结果导出视频（读 fluid-openfoam/case.foam，其余同 FDM 05）                       |
| `run_modal_analysis.py`       | 独立模态验证（可选）：对 dam_gate.inp 做 ccx 模态分析，与 APDL 参考频率对比              |
| `run_pipeline.sh`             | 流水线驱动入口（步骤 1,3,4,5）                                                           |

### 运行命令

```bash
cd /root/OpenFOAM/user-v2412/src
bash case/Damfailure/preCICE/run_pipeline.sh \
    --case-root /root/Workspace/<root> --steps 1,3,4,5 \
    --res 115 --z-layers 12 --end-time 0.5 --nprocs 8 --video
```

### preCICE 关键参数

| 项                     | 值                                                       |
| ---------------------- | -------------------------------------------------------- |
| 流体系                 | `interFoam`（原版，非 myInterFoam）                      |
| 时间步 deltaT          | **0.0005 s**（=preCICE time-window-size）                |
| 耦合                   | parallel-implicit + IQN-ILS + RBF thin-plate-splines     |
| time-window            | 0.0005 s，max-time 0.5 s，max-iterations 50              |
| 固体                   | ccx_preCICE，`--solid-threads` 默认 8                    |
| 挡板材料               | **E=1 MPa, ν=0, ρ=2500 kg/m³**（与 FDM 统一，原 zip 模板 1.4e6/0.4/10000 已改） |
| 并行                   | `--nprocs` 默认 8                                        |

### 为什么没有步骤 2

preCICE 无独立模态脚本（01 已生成固体网格，run_modal_analysis.py 是可选验证）。步骤为 1→3→4→5。

### 与 FDM 的效率差异

preCICE 因 time-window 减半（1000 窗口 vs 500 步）+ 外部耦合开销 + solid 参与（CalculiX 4 线程），串行时约 5.85x 慢于 FDM；但两者机制不同（preCICE 是更通用的分区耦合，FDM 是专用模态叠加）。

## 验证结果（实测）

- 模态一阶频率 6.045 Hz（vs APDL 基准 6.018，误差 0.45%）
- 小网格（51,792 单元）0.5s 串行耗时 9,664 s（~2h41m）
- 全分辨率（181,200 单元）8 核并行 0.5s 耗时 5,510 s
- 04 脚本内置时间估算：`~0.187 s/单元/步`（串行标定）
- 网格扩张比从 10.2 降到 4.0（MAX_EXPANSION 钳制 grading）
- 全分辨率网格 Max dx ≈ 9.3mm（改后）
- 统一边界层后：全部 7 个壁面第一层厚度一致（实测 +0.0% 偏差）
- checkMesh：Mesh OK，skewness ~2.5e-13

## 解析产物

设置全文（Scheme 格式）在解析时转储至临时目录，查具体条目可直接 grep：

- `/root/Workspace/tmp/opencode/Rampant_Variables.txt`（320 KB，求解器全局变量）
- `/root/Workspace/tmp/opencode/Thread_Variables.txt`（218 KB，各边界/单元区设置）
- `/root/Workspace/tmp/opencode/Domain_Variables.txt`（15 KB，各 domain 求解状态）
- `/root/Workspace/tmp/opencode/Cortex_Variables.txt`（5.7 KB，参考系/图形状态）

对比用算例解压位置：

- `/root/Workspace/tmp/opencode/damfailure_fine_symm/`（fastDynamicFvMesh 模态算例，含 `mode/` 模态文件）
- `/root/Workspace/tmp/opencode/precice_damfailure_fine/`（preCICE 分区耦合算例，含 `solid-calculix/` 固体）
