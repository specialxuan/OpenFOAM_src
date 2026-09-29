# FDM + HiSA + RANS + AGARD 445.6 —— 颤振结论汇总

日期：2026-09-28
范围：在 AGARD 445.6 半展机翼、M=0.96、q=2936 Pa 工况下，用 **HiSA**（密度基耦合隐式求解器，RANS `kOmegaSST`）驱动 **FDM**（`fastDynamicFvMesh`，本仓库的模态叠加动网格库）做气动弹性（颤振）计算。
状态：**跑通并完成 5 次独立运行 + 3 层验证**。这是与 `HiSA_CONCLUSIONS.md` 所记录的 `displacementLaplacian` 路线**不同**的一条路线。
判决：**该工况下存在模态 1 主导的持续失稳，同时经受住 dt 细化与耦合拧紧；但其速率未被收敛量化。三层研究整体仍为 INCONCLUSIVE，本结论不给出可发表的颤振边界。**

---

## 1. 目标与路线

| 项 | 内容 |
|---|---|
| 目标 | AGARD 445.6 半展机翼，M=0.96，q=2936 Pa 的气动弹性（颤振）计算 |
| 求解器 | **HiSA**：密度基耦合隐式求解器，RANS `kOmegaSST` |
| 动网格 | **FDM**（`fastDynamicFvMesh`）：本仓库的模态叠加动网格库 |
| 与既有路线的关系 | 与 `HiSA_CONCLUSIONS.md` 所记录的 `displacementLaplacian` 路线**不同**的另一条路线 |
| 状态 | **跑通并完成 5 次独立运行 + 3 层验证** |

§2 给出的正面结论说明：**HiSA 并不绑定某一种动网格实现**，FDM 可以整体接上；既有的“尺度失配”判断也不是 FDM 路线的障碍。

---

## 2. 核心结论（正面）

### 2.1 HiSA 能驱动任意 `dynamicFvMesh` 子类

源码层面确认：

- `createDynamicFvMesh.H` 按 `dynamicMeshDict` 里的 `dynamicFvMesh` 关键字做**运行时选择**；
- `hisaModule.C::outerIteration()` 用 `isA<dynamicFvMesh>(mesh)`，随后 `dynamic_cast<dynamicFvMesh&>(mesh).update()`；
- **不** cast 到 `dynamicMotionSolverFvMesh`，**不**查 `motionSolver` / `cellMotionU`。

因此 FDM 作为 `dynamicFvMesh` 的普通子类即可被无差别接管。

另一面：`fsiCoupling` / `maxFsiIter` / `fsiTolerance` 对 HiSA **完全无效**（在 HiSA 源码树里零匹配），这些键仅为 `my*Foam` 家族服务。FDM 与 HiSA 之间的耦合语义必须由 `dynamicMeshDict` 自身给出，不能指望 `fsiCoupling`。

### 2.2 尺度失配不是 FDM 路线的障碍（本次最重要的正向发现）

此前该几何上用 `displacementLaplacian` 做气弹瞬态时：近壁首层 **12.79 µm** 与 O(1 mm) 的模态位移相差约 **100 倍**，导致 Euler 在 **t = 0.0575 s**、RANS 在 **t = 0.022 s** 崩溃。

本次 FDM 路线下：

| 量 | 值 |
|---|---|
| 翼尖物理位移 | **1.597e-3 m = 124.9 个首层厚度** |
| 倒转单元 | **0** |
| 最小体积变化 | **+0.2%** |
| 总体积守恒 | **7 位有效数字** |
| 非正交与斜度 | 平坦 |

**机理**：FDM 的渐变模态位移场把近壁层**近乎刚体**地带走；真正的判据是局部应变量（由 `checkMesh` 监测），而**不是**绝对位移 / 首层厚度之比。绝对位移远超首层厚度，本身并不必然毁网格。

### 2.3 网格与模态映射位级精确

- `polyMesh/points` 与 `FluidNodeCoor.csv` 在全部 **522300 行**上 `max|Δcoord| = 0`。
- 模态 1 的最大位移节点 **`256545`** 正是翼尖（x=1.17566, y=0.76200），位移为纯 z 向。
- 故模态 1 = 一阶弯曲。

### 2.4 结构 ODE 的 ω 正确

- 实测 `(Force_1 − Acc_1)/Disp_1 = 3666.4`。
- 模态文件表头 `ω² = (2π·9.640002)² = 3668.7`。
- 二者之比 **0.9994**。

### 2.5 力投影行为正确

模态 6 的广义力最大（`Force_6 ≈ −2.24`）。原因：它是唯一的**流向（x）主导**模态，而本工况 `α_eff = 0.0°`（升力≈0、阻力有限），稳态气动载荷主要投影到流向模态。**这不是 bug。**

---

## 3. 算例配置

### 3.1 网格

| 项 | 值 |
|---|---|
| 名称 | Fluid8 |
| 规模 | 505000 单元 / 522300 点 |
| patch | `FSI`(wall, 1754 面) / `FRONT` / `BACK` / `IN` / `OUT` / `TOP` / `BOTTOM` |

### 3.2 模态

- 10 阶 **soft** 翼模态（`mode_soft`），**f1 = 9.640002 Hz**，**522300 节点**。
- 选 soft 的依据：实验 M=0.96 点的 `freq ratio = 0.3647` = 13.8942 Hz / 38.1 Hz（soft 二阶），精确吻合。

### 3.3 来流

`ρ = 0.061499 kg/m³`，`T = 257.7964 K`，`p = 4551.0913 Pa`，`U = (309, 0, 0)`。

### 3.4 `constant/dynamicMeshDict`

```
dynamicFvMesh   fastDynamicFvMesh;
fsiPatches      ("FSI");
pressureField   p;
pRef            4551.0913;      // 来流静压，关键
theta           1.4;
mappingTolerance 1e-09;
modalMass       (1 1 1 1 1 1 1 1 1 1);
modalDamp       (0 0 0 0 0 0 0 0 0 0);
```

- `pressureField p`；`pRef 4551.0913`（来流静压）为**关键**：HiSA 的 `p` 是**绝对静压 [Pa]**，故 FDM 的 `pressureScale = 1`，`pRef` 必须取来流静压而非 0。
- `modalMass (1 ×10)`、`modalDamp (0 ×10)`；AMR 关闭。
- 初始场：复用已收敛的定常 RANS 解（**未重跑定常段**）。

### 3.5 气弹响应特征

- 模态 1 响应频率（气弹）：**≈12.3 Hz**；实验 M=0.96 颤振频率 13.8942 Hz（**−11%**）。
- 相位 **F₁ 超前 q₁ 112°**（90° 对应纯负阻尼），是典型颤振特征。
- 实验参考：`/root/Workspace/AGARD/validation/agard4456_flutter_boundary_experimental.csv`。

---

## 4. 证据链：五次运行

五次算例目录均在 `/root/Workspace/AGARD/HiSA/`。ζ 均在窗口 **[0.15, 0.35]** 内取值。

| # | 算例目录 | dt | 步数 / 时长 | ζ(+峰) / ζ(−峰) | 抽取结果与结论 |
|---|---|---|---|---|---|
| 0 | `fluid8_M096_fdm_rans_flutter`（探针，`log.probe_dt1e-4`） | 1e-4 | 50 步 / 0.005 s | — | **机制验证通过**：`Selecting dynamicFvMesh fastDynamicFvMesh`；`Mapped 590654 points out of 590654`；50/50 步 `pseudoTime: converged`；EXIT 0 |
| 1 | `fluid8_M096_fdm_rans_flutter`（正式） | 1e-3 | 430 步 / 0.43 s（4.14 个一阶周期） | **−0.174 / −0.150** | 模态 1 强发散，ζ≈−0.17，f≈12.3 Hz，网格完好 |
| 2 | `fluid8_M096_fdm_rans_flutter_dt5e-4` | 5e-4 | 700 步 / 0.35 s | **−0.113 / −0.064** | 比值 0.43–0.65 → **INCONCLUSIVE** |
| 3 | `fluid8_M096_fdm_rans_flutter_dt2p5e-4` | 2.5e-4 | 1400 步 / 0.35 s | **−0.058 / −0.078** | 三点最小二乘外推 dt→0：ζ_phys = **−0.027 / −0.035**（与 0 不可分）；但总包络 σ 非单调（**13.2→9.0→12.2 s⁻¹**）；模态 2 的 σ 随细化**单调上升**（**0.21→9.01→15.56 s⁻¹**）→ 模态内容不随 dt 收敛 |
| 4 | `fluid8_M096_fdm_rans_flutter_coupling_tight` | 5e-4 | 700 步 ×4 子迭代 / 0.35 s | **−0.053 / −0.094** | 把 `moveMeshOuterCorrectors` 0→1、`nPseudoCorrMin` 2→4（每步内迭代 2→4、每步网格更新 1→4 次，子迭代相对收敛 **4.8e−5**）后：ζ+ 比值 **0.47**、ζ− 比值 **1.45**；**模态 1 增长未塌陷**，**模态 2 增长反升 1.8×** → 耦合误差假说被削弱 |

**运行 #1 于 0.43 s 内完成 4.14 个一阶周期**，是判定模态 1 频率与阻尼的主依据。运行 #2、#3、#4 的时长均为 0.35 s。

---

## 5. 最终判决

### 5.1 Δt 三层研究：INCONCLUSIVE

- ζ+ 比值 **0.65 / 0.51**，ζ− 比值 **0.43 / 1.21**。
- 总包络 σ **非单调**。
- 严格规则下判为 **INCONCLUSIVE**。

### 5.2 耦合拧紧：仍 INCONCLUSIVE，但证据偏向物理

- 严格规则下 0.47 落在“耦合误差”桶、1.45 落空，故整体**仍判 INCONCLUSIVE**。
- **但耦合误差假说被削弱，证据偏向物理**：
  - 没有任何稳健统计量出现 “≤0.5× 或 |ζ|<0.03” 的塌陷；
  - 模态 2 的增长反而**上升**；
  - 0.47 那个数来自**单峰比值与小样本**。

### 5.3 综合

该工况下存在**模态 1 主导的持续失稳**（**σ ≈ +6 ~ +13 s⁻¹**，**ζ ≈ −0.05 ~ −0.13**）。它**同时经受住了 dt 细化与耦合拧紧**；但其**速率未被收敛量化**（响应已非单模态，模态 2 被卷入）。

---

## 6. 这不是什么 / 局限

本结论**不**给出可发表的颤振边界。

本结论**不**能区分“物理失稳”与“残余的空间 / 受力误差”。

本结论**不**支持“纯数值伪像”或“已收敛的物理颤振”**任一**强主张。

换言之：本次能确认的是“存在一个稳健的、模态 1 主导的失稳，且它不因减 dt 或加内迭代而消失”；不能确认的是它的物理归属与收敛速率。

---

## 7. 未解物理缺口与下一步

### 7.1 缺口

实验点 **M=0.96 / q=2936 Pa 正好落在实验颤振边界上**，收敛的求解器应给出 **ζ≈0**；实测却**持续为负**。

### 7.2 一个尚未测试的合理因素

本路线全部使用 `modalDamp (0 …)`，即**结构阻尼为零**；而 AGARD 445.6 实验边界通常按结构阻尼 **g ≈ 0.02** 定义。这是一个**合法**的贡献因素，但**尚未测试**。

### 7.3 指明的下一步（均未执行，属待决）

1. 给 `modalDamp` 赋予 ζ = 0.01–0.02 后重跑，检验边界是否移动；
2. 用更长记录 / 更多周期做单模态隔离（例如只保留模态 1–2）；
3. 网格与格式收敛研究（而**不是**继续减 dt 或继续加内迭代，这两条已被本次证否）。

---

## 8. 教训（工程层面）

| 教训 | 说明 |
|---|---|
| **`modal_diagnostics.csv` 会多行 / 步** | FDM 在拧紧耦合后每步有多行（每子迭代一行）。分析前必须按 `Time` 去重取 step-final，否则寻峰器会把重复行当成大量假峰。本次主会话就因此先算错了比值（**0.01**），去重后才与正确值（**0.47**）吻合 |
| **重启会重置模态状态** | FDM 在非零时刻启动会重跑 2 步 startup 初始化，因此每次运行必须是**单一连续进程**。探针的 `modal_diagnostics.csv` 里因此出现过第二个表头（两次运行拼接），其 `t=0.0053` 行曾被误读为“ODE 不自洽”，实为重启后仅积分 1 步的产物 |
| **不要用绝对位移 / 首层厚度比作为失稳判据** | 见 §2.2：FDM 下近壁层近乎刚体随动，绝对比值没有判据力 |
| **`/dev/shm` 只有 64 MB** | 并行运行前需清空并加 `--mca btl self,tcp` |
| **`mpirun` root 保护** | 需 `--allow-run-as-root` + `OMPI_ALLOW_RUN_AS_ROOT{,_CONFIRM}=1` |
| **并发纪律** | 一次只跑一个 16 rank 作业 |

---

## 9. 产物索引

四份报告均在 `/root/Workspace/AGARD/HiSA/fluid8_M096_fdm_rans_flutter/`：

| 路径 | 内容 |
|---|---|
| `/root/Workspace/AGARD/HiSA/fluid8_M096_fdm_rans_flutter/REPORT_fdm_hisa_flutter.md` | 运行 #1（含 Phase A 三道闸门、列语义、物理位移表） |
| `/root/Workspace/AGARD/HiSA/fluid8_M096_fdm_rans_flutter/REPORT_dt_refinement.md` | 运行 #2（两点对比） |
| `/root/Workspace/AGARD/HiSA/fluid8_M096_fdm_rans_flutter/REPORT_dt_convergence_final.md` | 运行 #3 + 三点 `σ(dt)` 拟合 + §7 综合 + §9 冻结证据 md5 |
| `/root/Workspace/AGARD/HiSA/fluid8_M096_fdm_rans_flutter/REPORT_coupling_resolution.md` | 运行 #4（耦合控制清单 + 拧紧对比 + §8 判决） |

其余产物：

- 各算例目录下的 `log.*` 与 `modal_diagnostics.csv`；
- 三层算例均**冻结保留**，原算例 **md5 未变**；
- 探针：`log.probe_dt1e-4`、`modal_diagnostics.probe.csv`。

各运行目录：

- `/root/Workspace/AGARD/HiSA/fluid8_M096_fdm_rans_flutter/`（运行 #0 探针 + 运行 #1）
- `/root/Workspace/AGARD/HiSA/fluid8_M096_fdm_rans_flutter_dt5e-4/`（运行 #2）
- `/root/Workspace/AGARD/HiSA/fluid8_M096_fdm_rans_flutter_dt2p5e-4/`（运行 #3）
- `/root/Workspace/AGARD/HiSA/fluid8_M096_fdm_rans_flutter_coupling_tight/`（运行 #4）

---

## 10. 复现

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
cd /root/Workspace/AGARD/HiSA/fluid8_M096_fdm_rans_flutter
rm -f /dev/shm/*
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
mpirun --allow-run-as-root --mca btl self,tcp -np 16 hisa -parallel > log.flutter_transient 2>&1
```

- dt 与耦合设置见 `system/controlDict` 与 `system/fvSolution`；
- 耦合拧紧点见 `_coupling_tight/` 兄弟目录。

---

## 11. 待决事项

1. `modalDamp` 赋 ζ = 0.01–0.02 后重跑的边界移动检验（§7.3）。
2. 更长记录 / 更多周期的单模态隔离（例如只保留模态 1–2）。
3. 网格与格式收敛研究；**不再**沿“继续减 dt / 继续加内迭代”方向推进（已被本次证否）。
