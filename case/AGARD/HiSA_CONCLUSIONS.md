# HiSA + AGARD 445.6 —— 结论汇总

日期：2026-09-28
范围：用 HiSA（密度基耦合隐式求解器）做 AGARD 445.6 机翼的**纯流体**与**气动弹性（颤振）**计算。
状态：**纯流体阶段完成；气弹阶段被网格形态阻断**。所有算例已停止。

---

## 1. 工具链

| 组件 | 说明 |
|---|---|
| 求解器 | **HiSA**（CSIR, Heyns/Oxtoby）：`AUSMPlusUp` + `dualTime rPseudoDeltaT backward` + `GMRES(LaxFriedrichs/LUSGS)`，`nPseudoCorr 20`，`maxCo 20000` |
| 气弹模块 | HiSA 自带 `libmodal.so` 的 `modalSolver` function object + `constant/modal/{mass,damping,stiffness}Matrix` |
| OpenFOAM | OpenCFD **v2412**（build `_f00230df-20250528`）；`WM_PROJECT_VERSION=v2412` ⇒ HiSA 走非 FOUNDATION 分支，读 `constant/turbulenceProperties`（旧式 `simulationType RAS; RAS { RASModel ...; }`） |
| 网格 A | **Fluid8**：505,000 hex，首层 12.79 µm，max AR 1.68e5，3,220 个 >70° 面（y+≈1 粘性网格） |
| 网格 B | **官方 cfMesh**：1,258,205 单元（+AMR ~1.6M），max AR 6.39，仅 3 个 >70° 面，`boundaryLayers nLayers 1` |

**性能**：官方 1.26M 单元 / 8 核 —— 定常 500 步 **39 min**；非定常 0.5 s **6.5 h**。
对比本项目早先的 `rhoPimpleFoam` 路线（1.039 s 物理时间需 ~5.8 天 / 64 核）：**HiSA 快约 2 个数量级**。
关键机制：dual-time 把声学 CFL 约束从物理步里剥离，物理 dt 可达 ~1.7e-3（原路线被压在 1e-6）。

---

## 2. Fluid8 定常结果（两条方程都收敛）

工况：M=0.96 / Pd=2936 Pa（`runs.orig` seq 43）；远场 `U=(309 0 0)`、`ρ=0.061499`、`T=257.7964 K`、`p=4551.0913 Pa`。

| 量 | Euler | RANS `kOmegaSST` |
|---|---|---|
| 迭代 | 1000（16 核，899.63 s） | 1000（16 核，1103.79 s，**1.23×**） |
| 收敛判据 | 首次于 it285（<0.5%/100 步），但不持续（~300 步极限环） | 首次于 **it523**（<1%/100 步） |
| **去假象 `Cd_x`** | 0.00073 | **0.008785**（摩擦主导：粘性 0.00725 + 压力 0.00153） |
| **去假象 `Cl_y`** | +4.8e−4 | −2.37e−4 |
| **去假象 `Cm_z`** | +2.9e−3 | −4.33e−3 |
| `y+`（翼面） | — | min 0.1075 / **mean 0.4299** / max 1.6065（99.7% 面 <1） |
| 极限环 | ⚠️ 有 ~1%（不衰减） | 无 |
| 湍流场界限 | — | `nut` 2.2e−11…1.93e−2（≥0），228 次轻微 `bounding k` |

**结论：RANS 是 Fluid8 上的生产质量基线；Euler 次之（极限环，非生产网格）。**

---

## 3. 力输出的「开口壳」假象（教训，必记）

`FSI` 翼面 patch 是**开口壳**（根部端盖不在其上），`forces` 积的是**绝对压力**：

```
raw Fy  =  ∫ p n dA  =  p∞·A_net  +  ∫ Δp n dA
             ≈ −38.3 N（假象，占 99%） + 真实气动力
A_net = (0, −0.0084108, 0) m²   （Fluid8 实测）
```

⇒ **`raw Fy ≈ −38.75 N` 不是升力**。官方算例同样如此。
Fluid8 真实气动力 ≈ **(0.78, 0.52, 0.02) N**（α_eff = 0.0° ± 0.01°，Cp 上下完全对称，max|ΔCp|=1.73e−5，升力物理上为零）。
**任何跨网格比较必须扣除该项，且永远不要把 raw Fy 当升力。**

---

## 4. Cp 分辨率：Fluid8 不足（代码间对比）

| 项 | 结果 |
|---|---|
| 对比对象 | 官方 cfMesh AMR 网格（Euler @500） |
| 平均 \|ΔCp\| | Euler 0.019 / RANS 0.023 |
| 最差局部 \|ΔCp\| | 0.120 / 0.136（压缩脚） |
| 原因 | Fluid8 每排仅 **~61 个弦向点**、中弦有 ~0.10 x/c 空档 → 激波被抹平，给出光滑无激波的重压 |
| 官方 | 分辨出陡峭后掠压缩（几个 % 弦长内 ΔCp≈0.16–0.18） |

**注意**：本地文献（Mihaila-Andres-2015，**计算**值）报告超临界但**无激波**——所以「官方网格 = 真值」本身也不确定。
抽取方法已自证：把官方 patch 压力独立积分 → `(1.34798, −38.75014, −0.38340) N`，与求解器日志 5 位有效数字吻合。

---

## 5. 实验数据：AGARD 445.6 **没有**实验 Cp

- 445.6 的 TDT 试验**只测了颤振**。NASA 2023 (NTRS 20220017734) 明确写「no other data (pressures, flow fields, etc.) exist」；Silva 2014 确认无表面压力测量。
- 文献里的「span station Cp」是**计算值**；真正的测压孔数据属于 **ONERA M6**（AGARD-AR-138）——不同机翼、不同 Mach，易混淆。
- **交付真实验颤振边界** → `/root/Workspace/AGARD/validation/agard4456_flutter_boundary_experimental.csv`（8 点，含 M=0.957）。
  自校验：M=0.96 → 速度指数 **0.3078**、频率比 **0.3647**，与 Mihaila-Andres 引用的实验值 0.307 / 0.365 **完全吻合**。
- **官方 `simulation_43` 的定量特征**（从原始数据重算）：一阶频率 **13.732 Hz**（实验 13.894，差 1.2%），对数衰减 **ζ ≈ −0.006 ± 0.008** → **近中性，恰好压在颤振边界**。

---

## 6. 实验压力验证：改用 ONERA M6（成功）

因 AGARD 445.6 无实验 Cp，用官方 HiSA `oneraM6` 例子对 AGARD-AR-138 测压孔做验证：

| 项 | 结果 |
|---|---|
| 网格 | 294,912 hex（306,577 点） |
| 运行 | 2000 步，16 核 **21.5 min**，exit 0 |
| 来流 | M 0.8395 / α 3.0601° / Re_MAC 1.172e7（与标准 2308 **完全一致**） |
| 实验源 | NASA TMR `case_2308.dat`（Schmitt & Charpin, AGARD-AR-138 1979），271 孔 |
| 对比 | mean \|ΔCp\| **0.049–0.098**；核心区 max 0.066–0.277；**激波 x/c 误差 ≤0.042（外翼 ≤0.019）**；LE 吸力峰吻合 ~0.02 |
| 结论 | **HiSA RANS 对跨音速机翼压力可信**；内翼 λ 激波脚偏弱 |

---

## 7. 【核心负面结论】Fluid8 不适合气动弹性

### 7.1 用的**不是 FDM**

崩掉的两个气弹算例 `constant/dynamicMeshDict` 逐字相同：

```
dynamicFvMesh      dynamicMotionSolverFvMesh;         // OpenFOAM 原生
solver             displacementLaplacian;             // 解 Laplace 方程求点位移场
displacementLaplacianCoeffs { diffusivity inverseFaceDistance (FSI); }
```

- **动网格**：`dynamicMotionSolverFvMesh` + **`displacementLaplacian`**（`inverseFaceDistance` 让扩散系数在翼面附近很大 → 位移被**集中吸收在最薄的头几层**）
- **结构/模态**：HiSA `modalSolver` function object（`libmodal.so`）
- **`libfastDynamicFvMesh`（FDM）完全没参与**

对照：FDM 路线（`/root/Workspace/AGARD/FDM/soft`）用 `dynamicFvMesh fastDynamicFvMesh;` + `fsiCoupling integrated;`，机制是**直接把模态振型场加到点上**（`fastDynamicFvMeshUpdate.C:251-277`）并带 `maxDispChange_` 限幅。

### 7.2 A/B：同一网格，不动不崩

| 运行 | 网格运动 | 结果 |
|---|---|---|
| Euler 定常 | ❌ 无 | ✅ 1000 步收敛 |
| RANS 定常 | ❌ 无 | ✅ 1000 步收敛 |
| 气弹 Euler | ✅ 有 | ❌ **t = 0.0575 s**（T→3.07e9 K，`thermo::T()` 迭代失败） |
| 气弹 RANS | ✅ 有 | ❌ **t = 0.014–0.022 s**（`bounding omega`→1e25、`bounding k`→1e13） |

**同一网格、同一求解器，唯一差别是有没有网格运动。**

### 7.3 直接证据：近壁网格被破坏

RANS 气弹的 `y+`（翼面）：

```
t≤0.015   min 0.107→0.053   max 1.64→2.10    mean ~0.425   （健康）
t=0.016   min 0.0455        max 17391        mean 18.5     ← 首次崩坏
t=0.018   min 0.0232        max 2.16         mean 4.28
末端       min 0.0151        max 22289        mean 17.29
```

**mean 恒为 0.43（附面层整体健康），但 max 尖峰、min 塌陷** → 极少数近壁单元被拉坏/压扁。均值被健康单元主导，**必须看 max/min**。

### 7.4 换方程没用、减小 dt 也没用

| 改动 | 结果 |
|---|---|
| Euler → RANS | 崩得**更早**（0.0575 → 0.022 s）|
| dt 减半（maxDeltaT 2e-3 → 5e-4） | 崩得**更早**（0.022 → 0.013 s）|

若瓶颈在方程或时间推进，这两处修改至少应有帮助；它们都使情况恶化 → **瓶颈在几何**。

### 7.5 根因：尺度失配 ≈ 100 倍

- Fluid8 **首层 12.79 µm**
- 一阶模态翼面位移 **O(1 mm)**
- 比值 **≈ 78–100 倍**

网格运动求解器要把位移分摊到 ≥ 位移量级的距离；Fluid8 近壁只给了 12.79 µm，位移全砸进最薄的头几层 → 单元长宽比爆炸/翻负 → `omegaWallFunction` 的 ω=6ν/(β₁y²) 发散 → 湍流爆 → 能量方程 T 反演失败。

### 7.6 官方网格（同一运动算法 + 同一扩散系数）活得很好

| 网格 | 翼面附近最细单元 | 运动/单元 | 结果 |
|---|---|---|---|
| 官方 cfMesh | 0.4 m ÷ 2⁸ ≈ **1.56 mm** | ≈ **0.6 倍** | ✅ 跑满 4.8 周期 |
| Fluid8 | **12.79 µm** | ≈ **78–100 倍** | ❌ <0.25 周期崩溃 |

官方 `meshDict` 原话：`nLayers 1; // Boundary layer cause stability issues with Euler runs` —— HiSA 作者自己就只用 1 层边界层。

### 7.7 证据边界（未做的一步）

气弹算例与定常算例实际差三件事：① 物理时间推进；② **网格运动**；③ 模态耦合反馈。
证据强烈指向 ②（y+ 直接测到近壁破坏；官方网格用同样的 ①③ 却没事），但**未做**隔离实验。

**确认实验（未执行）**：冻结流场、只动网格（或 `simulationType laminar`），用同一模态运动驱动一个周期。
- 也崩 → 100% 锁定几何/动网格
- 不崩 → 说明湍流耦合也有份

### 7.8 官方 cfMesh 网格上的**定常复现**（独立完成 ✅）

在 `/root/Workspace/AGARD/HiSA/agard445_flutter/` 用官方工作流独立复现 M=0.96 的定常段（16 rank, scotch）：

| 项 | 本次（16 rank scotch） | 官方（8 rank hierarchical） |
|---|---|---|
| AMR 后单元数 | **1,599,770** | ~1.6 M |
| 步骤 | 500 | 500 |
| `ExecutionTime`（求解器） | **1780.99 s** | 2344.55 s |
| `ClockTime`（含 redistribute） | **2705 s ≈ 45.1 min** | 2381 s ≈ 39.7 min |
| `wing` 上最终力 | **(1.34777, −38.74977, −0.38247) N** | **(1.34798, −38.75014, −0.38340) N** |

**与官方吻合到 4–5 位有效数字**（同样注意：raw `Fy` 是开口壳压力假象，不是升力）。
注意：16 rank scotch 的**墙钟反而比官方 8 rank 更长**——AMR 的 `rebalance`/redistributor 不能并行。

**该算例的气弹非定常段未完成**（被本会话的"卡死会话"误判 kill 掉，见 §10 教训），只在
`t=5.7e-3 s`（0.055 T1）留下 6 行历史；其起始段与官方参考逐点吻合（`t=4.2986e-3` 时
`x0` 相对差 2.2e-6 %）。Part 2（M=0.901 / q=4277 Pa）已 staged（`simulation_34/`）但 0 步执行。

---

## 8. Fluid8 的能力边界（已用实验数据锁死）

| 用途 | 结论 |
|---|---|
| 定常 RANS 纯流体 | ✅ **生产质量**（Cd_x=0.008785，y+ 均值 0.43，已验证求解器可信） |
| 定常 Euler | ⚠️ 可算但有 ~1% 极限环，非生产网格 |
| Cp 分辨（vs AMR 网格） | ⚠️ 每排 ~61 点，激波被抹平 |
| **气弹 / 颤振（含网格运动）** | ❌ **物理上不兼容**（首层 12.79 µm vs 位移 O(1 mm)） |

---

## 9. 交付物索引

| 路径 | 内容 |
|---|---|
| `/root/Workspace/AGARD/HiSA/fluid8_M096/` | Fluid8 **Euler 定常**基线（1000 步）+ `REPORT_fluid8_purefluid.md` |
| `/root/Workspace/AGARD/HiSA/fluid8_M096_rans/` | Fluid8 **RANS 定常**基线（生产）+ `REPORT_rans_fluid8.md` |
| `/root/Workspace/AGARD/HiSA/fluid8_M096_flutter/` | Fluid8 **Euler 气弹**（负面结果）：脚手架 + 部分模态历史 |
| `/root/Workspace/AGARD/HiSA/fluid8_M096_flutter_rans/` | Fluid8 **RANS 气弹**（负面结果）+ `REPORT_flutter_fluid8_rans.md`、`tools/*.dat` |
| `/root/Workspace/AGARD/HiSA/oneraM6/` | **ONERA M6 实验压力验证** + `REPORT_oneraM6_validation.md`、`results/*.csv`、`experimental/` |
| `/root/Workspace/AGARD/HiSA/cp_codeto_code/` | Cp **代码间对比** + `REPORT_cp_codeto_code.md`、`cp_*.csv` |
| `/root/Workspace/AGARD/HiSA/agard445/` | **官方参考**算例（含已完成的 `simulation_43` 气弹） |
| `/root/Workspace/AGARD/validation/` | **实验颤振边界** CSV + 出处的 PROVENANCE/README |
| `/root/Workspace/AGARD/HiSA/PLAN_fluid8_purefluid.md` | 原始计划 |
| `/root/OpenFOAM/user-v2412/src/REPORT_oneraM6_validation.md` | ONERA M6 报告（repo 内副本）|

---

## 10. 环境坑（值得记）

| 坑 | 表现 | 处理 |
|---|---|---|
| **`/dev/shm` 只有 64 MB** | OpenMPI vader 共享内存段被占 → `Signal: Bus error (7)` | 清理 `/dev/shm` + `--mca btl self,tcp` |
| **~~物理 dt 硬顶 1e-3~~（已证伪，勿再引用）** | 早期报告称 `adjustTimeStep yes` 时 dt 硬顶 1e-3（`hisaModule::timeStepScaling` → `fvOptions_->maxDeltaT()`）。**实测反驳：官方 dt = 1.667e-3 > 1e-3。源码查明：v2412 下 `MAJ_VER=2412 ≥ 1000` ⇒ `buildDefs.include` **不定义 `FOUNDATION`** ⇒ `#if FOUNDATION >= 10` 的 clamp **根本没编译进去** | 该 clamp 在本构建**不存在**。实测 dt = 1.667e-3 与官方一致，源于 `adjustableRunTime` 把 dt 截到 `writeInterval` 网格（0.5 s / 300 步）。想缩短运行：提高 `maxDeltaT`/`writeInterval`，或 `adjustTimeStep no` + 显式 `deltaT`；**不要假设有 1e-3 上限** |
| **`mpirun` root 保护** | 直接报错 | `--allow-run-as-root` + `OMPI_ALLOW_RUN_AS_ROOT{,_CONFIRM}=1` |
| **HiSA 需要 `cellLevel`/`pointLevel`** | 官方 `runSimFixed` 脚本硬性检查（AMR 需要） | 静态网格路径可绕过；否则补全零字段 |
| **`runSimTransient` 的拷贝步骤易漏** | 该脚本要先把**定常末态网格**（`fixed/<TLAST>/polyMesh`，AMR 细化版）和**末态场**拷进 `transient/constant/polyMesh` 与 `transient/0/`，再做 `getPatchPoints`/`genModeShape`/`decomposePar`。漏掉就报 `Cannot find file "points" in directory "polyMesh"` / `cannot open case directory ".../transient/processor0"` | 逐行照抄官方 `runSimTransient`；`$LASTMESH` 会解析为 `500`（不是 `constant`）|
| **`FRONT` 是普通 `patch`** | 标量场用 `symmetry` BC 会 FATAL | 用 `zeroGradient`(标量) + `slip`(U)，在平面 patch 上等价于对称面 |
| **`decomposeParDict` 的 `n` 与 `PROC` 必须匹配** | 官方 fixed 用 `method hierarchical; n (2 2 2)`=8 域；把 `PROC` 改成 16 后 `decomposePar` 直接 abort（`Wrong number of divisions (2 2 2) for 16 domains`）| 改 `method scotch`（或把 `n` 改成 `(4 2 2)`）|
| **⚠️ `runSimTransient` 的 `cp -r $TLAST/*` 陷阱（危险）** | 若定常段没产出时间目录，`TLAST=$(foamListTimes -latestTime)` 为空 → `cp -r $TLAST/* ../transient/0/` 变成 **`cp -r /* ../transient/0/`，把根文件系统拷进算例**（实测一次 2.6 GB 的 `/etc`、`/proc`…）| 在 `TLAST=` 之后加 `if [ -z "$TLAST" ]; then exit 1; fi`；且 **`runSim` 没有 `set -e`**，fixed 失败后 transient 仍会跑 |
| **并发 agent 互相 kill** | 多个 16 rank 作业互杀，负载 36/20 核，模态历史反复重置 | 一次只跑一个 16 rank 作业；不要写自愈/看门狗重启脚本 |
| **慢 ≠ 卡死（本轮教训）** | 一个仍在推理/读文件的子 agent 有 ~27 min 无文件写入，被主会话误判为"卡死"并 kill 其进程 → **直接毁掉一次气弹瞬态**（t=5.7e-3 s）。反过来，该子 agent 也在报告里把主会话（轮询 + kill）识别为"恶意并发会话" | 判断子 agent 是否卡死要看**是否仍在消耗 CPU / 是否有子进程**，不能只看文件 mtime；kill 之前先确认所有权 |
| **僵尸进程无法回收** | 700+ 个 `[hisa]` 僵尸 | 容器内 PID 1 是 `sleep infinity`，不 reap；无害（不占 CPU），重启容器才清 |

---

## 11. 待决事项

1. **气弹网格路线**
   - **A**：官方 cfMesh 网格（已验证可跑、ζ≈−0.006、一阶频率 13.73 Hz）——最快拿到可用颤振结果
   - **B**：按气弹需求重画 Fluid8 几何的网格（**y+ 30–100 壁面函数思路，首层 ~0.5 mm**，运动/首层 ≈ 2）
   - **C**：两者都做
2. **是否做「冻结流场、只动网格」的确认实验**（把 ② 与 ①③ 分离，见 §7.7）
3. **是否把 FDM（`fastDynamicFvMesh`）与 `displacementLaplacian` 在同网格同初场下对照**——FDM 直接施加全局模态场且带位移限幅，破坏机制不同，但尺度失配的根本矛盾仍在

---

## 12. 如何续跑（容器重启后自包含）

**唯一在跑的作业 = 无。** 官方网格上的**定常段已完成**（`simulation_43/fixed/500`），只差**气弹瞬态**。

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

# 0. 必须确认没有竞争进程；清理 64 MB /dev/shm
ps -eo stat,comm | awk '$2 ~ /hisa|mpirun/ && $1 !~ /Z/'   # 必须为空
rm -f /dev/shm/*

# 1. Part 1 气弹瞬态 M=0.96（定常段已在 simulation_43/fixed/500）
cd /root/Workspace/AGARD/HiSA/agard445_flutter/simulation_43
nohup ./runSimTransient > /root/Workspace/tmp/run_part1_M096_transient.log 2>&1 &

# 2. Part 2：M=0.901 / q=4277（已 staged）—— 等 Part 1 结束再跑，一次只跑一个 16 rank 作业
cd /root/Workspace/AGARD/HiSA/agard445_flutter/simulation_34
nohup ./runSim > /root/Workspace/tmp/run_part2_M0901.log 2>&1 &

# 3. 分析（需 ≥4 个一阶周期；窗口取末段 ≥0.25 s）
cd /root/Workspace/AGARD/HiSA/agard445_flutter
python3 tools/analyze_modal_official.py \
  "M096=simulation_43/transient/postProcessing/modalMotion/0/modalSolver.dat" \
  "M0901=simulation_34/transient/postProcessing/modalMotion/0/modalSolver.dat" --window 0.25
```

**工作副本里已打的 3 个补丁（仅副本，官方 example 未动）**：
1. `settings`: `PROC=8 → 16`
2. `simulation/*/system/decomposeParDict`: `hierarchical (2 2 2)` → `scotch`（16 域必需）
3. `simulation/runSimTransient`: 加 `if [ -z "$TLAST" ]; then exit 1; fi`（防 `cp -r /*` 事故）；`mpirun` 加 `--allow-run-as-root --mca btl self,tcp`

**成本**：定常 45 min + 瞬态 ~2.9 h（16 rank）→ 单个完整工况点 **≈3.5–4 h**；两点 campaign **7–8 h**。
**预期 Part 1 结果**（来自参考）：f₁ ≈ **13.73 Hz**、ζ ≈ **−0.006 ± 0.008**、近中性。

**纪律**：
- **一次只跑一个 16 rank 作业**；
- **绝不再创建自愈/看门狗重启脚本**；
- 官方 example `/root/OpenFOAM/user-v2412/hisa/hisa/examples/agard445/` 与既有已完成算例**只读**；
- 工作根目录 `/root/Workspace/AGARD/HiSA/agard445_flutter/`，scratch 用 `/root/Workspace/tmp/`。
