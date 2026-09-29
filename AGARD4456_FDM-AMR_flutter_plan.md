# AGARD 445.6 机翼 FDM-AMR 颤振分析计划

## 目标

用模态叠加 FDM + AMR 动网格做 AGARD 445.6 半展机翼的颤振分析，输出颤振边界并与 Yates 实验对比。

## 现有资产

- 模态：`/root/Workspace/AGARD/Mode/mode_soft/`、`mode_hard/`（各 23 CSV，10 阶，已过门）
  - `FluidNodeCoor.csv` = 388,465 个远场流侧节点，`StructNodeCoor.csv` = 977,086 个机翼节点
  - 软翼 f1 = 9.62 Hz，硬翼 f1 = 14.52 Hz
- 流水线：`case/AGARD/Mode/`（assemble → generate → export，v4.2.2）
- 远场网格：`case/AGARD/Mode/mesh/farfield_c3d8.zip`（388,465 节点 / 372,760 C3D8，含 `FSI/FRONT/BACK/IN/OUT/TOP/BOTTOM` NSET）
- 求解器：`myRhoPimpleFoam`（可压缩 FSI）+ `fastDynamicFvMesh`（AMR）

## 缺口

- OpenFOAM 流体 case 不存在（Workspace 内无 controlDict）
- v4.2 移除 `scripts/mesh/` 后无 inp→polyMesh 转换器

## 进度

- [x] 第 1 步 转换器 + 网格：`case/AGARD/Mesh/inp_to_polymesh.py` → `/root/Workspace/AGARD/CFD/soft/constant/polyMesh`
  - 388,465 点 / 372,760 单元 / 1,102,739 内面 / 31,082 边界面，7 个 patch（FSI 为 wall）
  - 与 `FluidNodeCoor.csv` 逐点坐标 `max|diff| = 0.0`；单元体积全正；拓扑/上三角序/面顶点检查通过
  - ⚠️ `checkMesh` 报 `Failed 1`：高长宽比（峰值 1.36e6，87,596 单元）；另有 13,818 个面非正交 >70°（峰值 89.998°）
  - ⚠️ 量化：13,316 单元体积 <1e-12 m³（最小 9.5e-17）；18,886 单元长宽比 >1e5；翼尖有 4,726 条边 <1e-6 m（无零长度边，是真实薄单元，非重合节点）
  - 该质量源自 `Assembly3.inp` 的逐行裁剪（PROVENANCE），非转换产物
- [ ] 第 1.5 步（待定）：决定是否修网格后再进第 2 步

## 步骤

1. **远场网格 → polyMesh**：新写 inp→polyMesh 转换脚本（放 `case/AGARD/`），C3D8 转 hex，NSET 归类为边界 patch，`FSI_NSET` 为动网格壁面。
   验收：点数 = 388,465 且与 `FluidNodeCoor.csv` 坐标一一对应，`checkMesh` 全过，体积全正。
   状态：已完成；坐标/计数门全过，`checkMesh` 因源网格固有翼尖薄单元剩 1 项未过。
2. **搭 myRhoPimpleFoam case**：`0/ constant/ system/`，空气自由来流、FSI 壁面 movingWall、k-omega SST。
3. **配 dynamicMeshDict + AMR**：10 阶模态 CSV、Wilson-Theta、耦合松弛；`meshRefinementSupport true`，密度/压力梯度追激波 + 尺寸下限。
   待验证：CalculiX 特征向量归一化与 `modalMass` 一致。
4. **本机 smoke**：极短 endTime，验证模态映射满配、`fsiResidual` 发布、AMR 触发、并行集合通信正常。
5. **同步 FSIDocker 正式算**：rsync case + 代码（远端仓库仍 v4.1.2，需同步并重新 wmake），64 核 decomposePar；先亚临界收敛，再穿颤振边界。
6. **颤振边界提取**：从 `modal_diagnostics.csv` 提各阶阻尼，扫速压/马赫数拟合零阻尼点，对比 Yates 跨音速 dip。

## 待定

- 软翼 or 硬翼（建议先软翼 Model No.3）
- 工况点：先 M=0.499 验证，还是直接跨音速 M≈0.9~1.0
- 结构阻尼比取值（文献常用 0~2%）
