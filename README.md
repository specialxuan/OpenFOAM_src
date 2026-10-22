# OpenFOAM v2412 User Code

本仓库是基于 OpenFOAM v2412 的用户代码树，面向气动弹性/流固耦合（FSI）与动态网格算例的开发与验证。

## 组件

- `fastDynamicFvMesh/`：模态叠加动态网格库（AMR、重启、力投影、模态扩展）
- `fsiCoupling/`：FSI/PIMPLE 收敛契约（header-only）
- `myPimpleFoam/`：不可压 FSI 求解器（含 overset / SRF 变体）
- `myInterFoam/`：两相 VOF FSI 求解器（含 overset 变体与三相子求解器）
- `myRhoPimpleFoam/`：可压 FSI 求解器家族
- `myHisa/`：HiSA v1.12.7 副本，在伪时间循环中接入 FSI 模态力判据
- `case/`：算例生成与流水线（Python + Bash，不走 `wmake`）
  - `Damfailure/`：溃坝 FSI（FDM / preCICE）
  - `AGARD/`：AGARD 445.6 机翼（模态标定、网格工具、LFS 网格包）
  - `SwingBeam/`：封闭水箱悬臂梁 FSI（FDM / preCICE）

## 文档

- `AGENTS.md`：构建、验证、文件约定与协作规则
- `fastDynamicFvMesh/README.md`、`fastDynamicFvMesh/Doc/快速动网格程序手册.md`：网格库文档
- `case/*/CaseDescription.md`、`case/*/AGENTS.md`：各算例的物理说明与流水线契约
