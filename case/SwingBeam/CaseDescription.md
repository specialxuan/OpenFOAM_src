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
- **组装默认值**（`03_assemble_fdm_case.py` 生成，`run_pipeline.sh` 同步）：`couplingRelaxation 0.7`、顶层 `fsiTolerance 1e-3`、`modalMass (1..1)` 与 `modalDamp (2*zeta*omega_m)`（`zeta=0.03`）、`nOuterCorrectors 100`；梁材料（`E=1e6, nu=0, rho=2500`）与周期端部激励保持不变。
- **稳定性实测**（`res=30`、`z-layers=4`、`A=0.1 N`、`dt=1e-3`、8 ranks）：
  - 请求值 `couplingRelaxation=0.7`：在第一个时间步内发散——位移增量逐迭代放大（`2.1e-5 → 1.4e-4 → 3.1e-3 m`，×7/×22），随后压力求解 GAMG SIGFPE（`t=0.001`，退出码 136）。
  - 本库常数松弛的稳定窗口约为 `α ≤ 6e-3`：`α=0.01`（100 迭代）在 `t=0.029` 崩溃；`α=0.007/0.008` 残差不收敛；`α=0.003/0.004/0.005` 全程稳定。更小并非更好（`α=0.002` 启动段更差）；把 `nOuterCorrectors` 提到 300（`α=0.005`）反而在 `t=0.021` 崩溃；`dt` 减半只重新采样启动窗口、不能消除。库内没有 Aitken/自适应松弛实现（`couplingRelaxation` 为常数），因此上述常值扫描即为可用诊断手段。
  - **最佳诊断值 `α=0.004`（2.0 s 全程）**：2000/2000 步完成，`ExecutionTime 1782 s` / `ClockTime 1893 s`（8 ranks）。**每一步末端的最大网格位移增量都 < 1e-6 m**（中位 `1.1e-9 m`，最大 `4.2e-8 m`），即物理层面的位移耦合每步都已收敛。
  - 力残差：末端残差中位 `6.7e-4`，`1204/2000` 步（60.2%）`fsiResidual < 1e-3`；796 步打印 “FSI force residual did not converge”。告警集中在启动段（`t≤0.011`）与模态力过零窗口。`writeDiagnostics true` 的诊断运行（`t=0.411`）显示最大相对残差由 mode 4 主导：其力在 `±3e-4` 附近过零（其余模态 `|F|=0.1–2.4`），单迭代绝对变化仅约 `5e-5`，而位移增量已收敛到 `1.2e-9 m`。也就是说，库内 `res = |ΔF|/|F|` 的逐模态相对力残差在近零力模态上病态，该告警并不表示物理未收敛。
  - 结论：在现有常数松弛库下，`fsiTolerance=1e-3` 无法在启动段与过零窗口逐步满足；`α=0.004` 是当前可稳定实现的最佳耦合设置，其物理响应与旧 `α=0.01` 基线一致（梁端 `dx` p2p `24.8 mm` vs `25.4 mm`）。
- `A=0.1 N`、2.0 s 梁端位移（`postprocess/fdm_gate_topright_disp.csv`）：`dx` p2p **24.80 mm**（min `-12.99 mm @ t=0.44 s`，max `+11.81 mm @ t=0.66 s`），每 0.4 s 周期 p2p `20.6/24.8/17.0/16.6/17.9 mm`，持续振荡无衰减；`dy` p2p `1.10 mm`。
- 复现（组装默认 + 诊断松弛值）：

```bash
cd case/SwingBeam/FDM
bash run_pipeline.sh --steps 3,4 --res 30 --z-layers 4 --end-time 2.0 --write-interval 5 \
    --case-root /root/Workspace/StructTest/fdm --force-amp 0.1 --coupling-relax 0.004
# 400 个时间目录的 8 路并行重建（本机 DrvFs 串行重建约 90 min）：
bash /root/Workspace/StructTest/probes/reconstruct_fdm_full.sh
bash run_pipeline.sh --steps 5,6 --res 30 --z-layers 4 --end-time 2.0 \
    --case-root /root/Workspace/StructTest/fdm
```

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

- `--no-excite` 关闭周期激励；`--force-end-time` 指定激励序列时长；`--coupling-relax` 覆盖耦合松弛（默认 `0.7`；本算例实测稳定诊断值 `0.004`）；`--modal-zeta` 生成 `modalMass`/`modalDamp`（默认 `0.03`）。
