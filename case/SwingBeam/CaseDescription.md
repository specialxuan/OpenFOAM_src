# StructTest 封闭水箱悬臂梁基础算例说明

## 计算域与流体几何

- 三维准二维封闭水箱，坐标范围为 `x=[0,0.584] m`、`y=[0,0.365] m`、`z=[0,0.012] m`，外形尺寸与 DamFailure 完全一致。
- 水箱完全封闭，内部充满单相水，无空气、无自由液面：`alpha.water=1` 全域。
- 流体网格拓扑与 DamFailure 相同（`x` 分割点 `0.146/0.292/0.304`），其中 `x=[0.292,0.304] m`、`y=[0,0.120] m` 的块为梁体空腔，不生成流体单元。
- 初始网格由 `blockMesh` 生成，局部边界层采用分段 grading（`--res`、`--z-layers` 控制）。

## 弹性悬臂梁几何与材料

- 悬臂梁位于 `x=[0.292,0.304] m`（宽 `12 mm`）、`y=[0,0.120] m`（高 `120 mm`）、`z=[0,0.012] m`（厚 `12 mm`）。
- 梁底 `y=0` 完全固定（clamped），梁顶 `y=0.120 m` 为自由端。
- 流体侧耦合面 `FSI_FLUID` 由三段组成：梁左面 `x=0.292`、梁右面 `x=0.304`、梁顶面 `y=0.120`。
- 结构为实体六面体网格（`6 x 45 x 3`，1288 节点 / 810 单元），通过 CalculiX 模态分析生成流体虚拟弹性体的模态数据。

结构材料参数：

| 参数 | 值 |
|---|---:|
| 杨氏模量 `E` | `1.0e6 Pa` |
| 泊松比 `nu` | `0` |
| 密度 `rho` | `2500 kg/m3` |

## 前 10 阶固有频率参考

以下频率来自当前算例 `mode/FluidNodeDisp1..10.csv` 首行的 CalculiX 模态结果（`res=30`、`z-layers=4` 的冒烟网格）：

| 阶次 | 频率 `Hz` |
|---:|---:|
| 1 | 2.707922 |
| 2 | 16.372890 |
| 3 | 41.668780 |
| 4 | 43.572710 |
| 5 | 80.118360 |
| 6 | 123.612300 |
| 7 | 125.057100 |
| 8 | 172.079600 |
| 9 | 208.597900 |
| 10 | 224.165900 |

## 周期梁端激励

- 载荷形式：`F(t) = A*sin(2*pi*f*t)`，方向 `+x`，作用于自由端节点 `(0.304, 0.120, 0)`。
- 默认幅值 `A = 0.02 N`，频率 `f = 2.5 Hz`（均可在流水线中通过 `--force-amp` / `--force-freq` 调整）。
- 由 `02_run_modal_virtual_elastic.py --excite-force` 生成 `mode/StructForce1.csv`（时间步长 `1e-4 s`，覆盖 `[0, max(endTime, 0.5)]`）。
- 运行时由 `fastDynamicFvMesh` 的 `structuralForceEnabled` 机制读取：
  `nStructuralForces 1; structuralForceFilePrefix StructForce; structuralTargetTolerance 1e-3;`
  并将目标坐标映射到最近的梁节点后投影为模态力。

## FSI 耦合与稳定性调参

- 结构文件：`mode/StructNodeCoor.csv`、`mode/StructNodeDisp1..N.csv`（梁节点，升序节点号）；流体文件：`mode/FluidNodeCoor.csv`、`mode/FluidNodeDisp1..N.csv`。
- 耦合方式：`fsiCoupling integrated`，模态叠加动网格（`fastDynamicFvMesh`），力残差判据 `fsiResidual < fsiTolerance`。
- **关键调参 `couplingRelaxation = 0.01`**：封闭满水水箱没有自由液面泄压，梁体被水包围，流固迭代呈附加质量主导。实测 `res=30` 网格下 `couplingRelaxation=1/0.2/0.05` 均在第一个时间步内发散（位移每迭代放大数倍），`0.01` 稳定（`endTime=0.5 s` 长程验证，可用 `run_pipeline.sh --coupling-relax` 覆盖）。
- 模态阻尼试验：`modalDamp=2e5`（约 `zeta=6000`）可稳定迭代，但梁响应被冻结（`t=0.05 s` 位移约 `1e-7 m`），因此默认不使用；如需可加 `--modal-zeta`。
- 默认参数下 `t=0.05 s` 梁自由端位移约 `0.25 mm`（见 `postprocess/fdm_gate_topright_disp.csv`），随激励时间继续增长。

## 边界条件

### 流体边界

- `FRONT` / `BACK`：`symmetryPlane`，位于 `z=0.012 m` / `z=0 m`。
- `TOP` / `BOTTOM_FLUID` / `LEFT` / `RIGHT`：固定壁面（封闭水箱）。
- `FSI_FLUID`：梁体移动壁面（`movingWallVelocity`），与梁模态位移同步。

### 结构边界

- 梁底 `y=0` 固定位移约束，`z` 向由对称面约束。
- 梁与流体通过 `FSI_FLUID` 交换压力/剪切载荷与模态位移响应。

## 初始条件

- `alpha.water=1`：全域水体（`0/alpha.water` 内部场 = 1，`setFieldsDict` 盒体覆盖整个水箱）。
- `U=(0,0,0) m/s`：全域初始静止。
- `p_rgh`、`p`：按算例参考点初始化（FDM 保持 `(0.392, 0.36, 0.006)`；preCICE 组装为 `(0.146, 0.18, 0.006)`，并将 `p_rgh` 在 `TOP` 固定为 0）。
- `k`、`omega`、`nut`：按组装脚本中的湍流初值初始化。
- 从 `t=0` 开始，梁在静水压力与周期端部激励下振动。

## preCICE 直接分区耦合（`preCICE/` 目录）

- 求解器：`interFoam`（流体）↔ `ccx_preCICE`（固体），`parallel-implicit` + IQN-ILS，窗口 `0.0005 s`。
- 周期端部激励在 CalculiX deck 内实现：`*AMPLITUDE` 表采样 `A*sin(2*pi*f*t)`，经辅助节点 `NtipAux` 的线性 `*EQUATION` 均匀分配到梁顶面 28 个节点。直接用 `*CLOAD` 写在 `Ninterface` 上不可行：adapter 的 `setNodeForces` 每窗口覆盖接口节点的 `xforc`，且 CalculiX 拒绝同节点同自由度混用不同 amplitude（`forcadd.f`）。
- 封闭水箱无自由液面：`p_rgh` 在 `TOP` 取 `fixedValue 0`，避免 `adjustPhi` 在准静止工况下因相对容差 `1e-8` 触发 FatalError，并给出物理压力基准。
- 收敛判据仅保留 `Force`：preCICE 对多判据取 OR，原 `DisplacementDelta < 5e-7` 会在耦合力收敛前提前结束迭代并注入能量。
- 稳定性：`res=30` 下 `E=1.0e6 Pa` 的分区耦合在约 15 个时间窗口后、`t≈0.008 s` 处因流体侧 SIGFPE 中止（独立验证）；因此 preCICE 冒烟验证使用 10 倍刚度的 `--gate-e 1.0e7`。
- 冒烟命令：

```bash
cd case/StructTest/preCICE
bash run_pipeline.sh --steps 1,3,4,5,6 --end-time 0.02 --res 30 --z-layers 4 --gate-e 1.0e7 --write-interval 0.002 --case-root /root/Workspace/StructTest/precice
```

## 流水线

```bash
cd case/StructTest/FDM
bash run_pipeline.sh --steps 1,2 --res 30 --z-layers 4 --case-root /root/Workspace/StructTest/fdm
bash run_pipeline.sh --steps 3   --res 30 --z-layers 4 --end-time 0.05 --case-root /root/Workspace/StructTest/fdm
bash run_pipeline.sh --steps 4   --res 30 --z-layers 4 --end-time 0.05 --case-root /root/Workspace/StructTest/fdm
bash run_pipeline.sh --steps 5,6 --res 30 --z-layers 4 --end-time 0.05 --case-root /root/Workspace/StructTest/fdm
```

- `--no-excite` 关闭周期激励；`--force-end-time` 指定激励序列时长；`--coupling-relax` 覆盖耦合松弛；`--modal-zeta` 生成 `modalMass`/`modalDamp`。
