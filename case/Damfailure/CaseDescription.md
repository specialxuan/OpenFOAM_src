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

## 五脚本流水线（FDM 算例生成、运行与可视化）

从零生成并运行 FDM（fastDynamicFvMesh）算例的完整流程，**不依赖任何已有算例**。脚本位于本目录：

| 脚本                              | 功能                                                                                     |
| --------------------------------- | ---------------------------------------------------------------------------------------- |
| `01_generate_mesh.py`             | 生成 OpenFOAM 流体网格（blockMesh）+ CalculiX 固体网格（dam_gate.inp）+ VTK 可视化产物   |
| `02_run_modal_virtual_elastic.py` | 挡板+虚拟弹性体组合模态分析（CalculiX ccx），导出 FDM 模态文件（mode/）                  |
| `03_assemble_fdm_case.py`         | 自包含组装 FDM 算例（初始/边界条件、仿真参数全部内嵌，只从输入复制 polyMesh + mode/）    |
| `04_run_fdm.py`                   | 运行 FDM（myInterFoam + fastDynamicFvMesh），支持 fresh/restart/continue、并行、时间估算 |
| `05_export_fdm_video.py`          | FDM 结果导出视频（ParaView pvpython 逐帧 PNG，ffmpeg 合成 MP4）                          |

### 完整流程

```bash
cd /root/OpenFOAM/user-v2412/src

# 1. 网格生成（--res 控制密度，--z-layers 默认 12）
python3 case/Damfailure/01_generate_mesh.py --res 115 --z-layers 12 \
    --out-dir /root/Workspace/<mesh> --view

# 2. 模态分析 + 导出（--export-mode 导出 mode/ 文件）
python3 case/Damfailure/02_run_modal_virtual_elastic.py \
    --fluid-dir /root/Workspace/<mesh>/fluid-openfoam \
    --out-dir /root/Workspace/<modal> --export-mode

# 3. 组装算例（默认 endTime 0.5, deltaT 0.001, writeInterval 50）
python3 case/Damfailure/03_assemble_fdm_case.py \
    --mesh-dir /root/Workspace/<mesh>/fluid-openfoam \
    --mode-dir /root/Workspace/<modal>/mode \
    --out-dir /root/Workspace/<fdm>

# 4. 运行（--dry-run 先预览；-f 全新 / -c 续算 / --nprocs 并行）
python3 case/Damfailure/04_run_fdm.py -f --case /root/Workspace/<fdm> --end-time 0.5

# 5. 导出视频（PNG 帧 → MP4；xvfb 提供无头渲染 DISPLAY）
xvfb-run -a pvpython case/Damfailure/05_export_fdm_video.py \
    --case /root/Workspace/<fdm> --out /root/Workspace/<frames>
ffmpeg -y -framerate 10 -i /root/Workspace/<frames>/frame_%04d.png \
    -c:v libx264 -pix_fmt yuv420p /root/Workspace/fdm_result.mp4
```

### 视频导出要点（05_export_fdm_video.py）

- 运行方式：`xvfb-run -a pvpython`（Ubuntu ParaView 是 X11 版，无头环境需 xvfb 提供 DISPLAY）
- 默认样式：白底 + 蓝色水相云图（alpha.water）+ 黑色网格线（`--no-wireframe` 关闭线框）
- 参数：`--case`、`--out`、`--res`（默认 1280x720）、`--field`（默认 alpha.water）、`--wireframe`/`--no-wireframe`、`--wire-color`
- 前置依赖：`xvfb`、`pvpython`（paraview + python3-paraview）、`ffmpeg` 已安装

### 关键默认参数（内嵌在脚本中）

| 项                     | 默认值                                                   |
| ---------------------- | -------------------------------------------------------- |
| 时间步长 deltaT        | 0.001 s                                                  |
| 仿真时长 endTime       | 0.5 s（论文/参考算例时长）                               |
| 输出间隔 writeInterval | 50 步（0.05 s）                                          |
| 挡板材料               | E=1 MPa, ν=0, ρ=2500 kg/m³                               |
| 流体物性               | 水 ρ=998.2/ν=1.0048e-6，空气 ρ=1.225/ν=1.4607e-5，σ=0.07 |
| 湍流                   | kOmegaSST                                                |
| FDM 动网格             | theta=1.4, fsiPatches=FSI_FLUID, pRef=101325             |
| 水柱初始               | box (0,0,0)-(0.146,0.292,0.012) 内 alpha.water=1         |

### 边界层设计（01_generate_mesh.py）

- 11 块结构化网格，x 方向 5 切分（0/0.146/0.292/0.304/0.584）
- 块 3/7/11（挡板右~RIGHT 两端都是壁面）用 blockMesh **多段 grading** 同时加密两端
- 所有壁面（BOTTOM/TOP/LEFT/RIGHT/挡板三面）**第一层厚度统一**（`--bl-thickness` 控制，默认 0.5×最小平均单元尺寸）

### 验证结果（实测）

- 模态一阶频率 6.045 Hz（vs APDL 基准 6.018，误差 0.45%）
- 小网格（51,792 单元）0.5s 串行耗时 9,664 s（~2h41m）
- 全分辨率（181,200 单元）8 核并行 0.5s 耗时 5,510 s
- 04 脚本内置时间估算：`~0.187 s/单元/步`（串行标定）

### 边界层与网格质量说明

- 统一边界层后：全部 7 个壁面第一层厚度一致（实测 +0.0% 偏差）
- RIGHT 端面修复：从 15.35mm（未加密）→ 目标值
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
