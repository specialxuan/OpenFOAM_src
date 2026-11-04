# AGARD 445.6 FDM 自包含流水线

从零生成并运行一个颤振算例，**不依赖任何外部工作区**（网格、模态、字典全部在本目录内）。
默认档位是一次约 1.5 小时的冒烟验证；完整 4 个周期用 `--end-time 0.286`。

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
. $WM_PROJECT_DIR/bin/tools/change-userdir.sh $HOME/OpenFOAM/user-v2412 -platforms

./run_pipeline.sh --case-root /root/Workspace/AGARD/FDM-run/demo
```

## 阶段

| 阶段 | 脚本 | 输入 | 输出 | 用时（8 核） |
|---|---|---|---|---|
| 01 | `01_generate_mesh.py` | `data/fluid8_mesh.zip` | `constant/polyMesh` | 秒 |
| 02 | `02_prepare_modes.py` | `data/modes_4m.zip` + 网格 | `mode/`（4 阶）+ 节点数门禁 | 秒 |
| 03 | `03_make_steady.py` | 网格 + `templates/steady.zip` | 稳态解 → `0/` | **25 min** |
| 04 | `04_assemble_case.py` | `templates/myHisa.zip` + 参数 | `system/`、`constant/`、`0/` 字典 | 秒 |
| 05 | `05_run_flutter.py` | 组装好的算例 | 时间目录 + `log.run` | **20 步 ≈ 65 min** |
| 06 | `06_postprocess.py` | `modal_diagnostics.csv` | 逐周期 p2p、σ、f | 秒 |

`run_pipeline.sh` 支持 `--steps 1,2,3,4,5,6`、`--dry-run`、`--force`、`--no-amr`。

## 默认值

| 参数 | 默认 | 说明 |
|---|---|---|
| `--mach` | 0.96 | 模板对应 M=0.96 |
| `--pressure` | 2936 Pa | 实验颤振点；p、ρ 按 q 线性缩放，T/U 不变 |
| `--modes` | 2 | 写成 `nActiveModes 2`（`mode/` 固定 4 阶） |
| `--end-time` | 0.02（20 步） | 冒烟档；完整 4 周期用 0.286 |
| `--nprocs` | 8 | 本地只用 8 个 P 核 |
| AMR | 开 | `--no-amr` 关闭（`refineInterval 1e8`、`meshRefinementSupport false`） |

## 目录内容

```
FDM/
├── data/fluid8_mesh.zip    25 MB  Fluid8 网格（5 个 polyMesh 文件 + sets）
├── data/modes_4m.zip       40 MB  FluidNodeCoor + Disp1..4（4 阶模态集）
├── templates/myHisa.zip    12 KB  瞬态模板（system 4 + constant 4 + 0/ 8 个场）
├── templates/steady.zip   5.8 KB  稳态模板（staticFvMesh + dualTime steadyState）
├── tpl.py                         从上面两个 zip 读取模板（脚本不再解包）
└── 01..06_*.py  run_pipeline.sh
```

`0/` 模板只保留边界条件结构（`internalField` 与壁面 `nonuniform List` 已收成 `uniform`），
所以只有 36 KB；真正的内部场由阶段 03 的稳态解填入。

## 三条容易踩的坑（已在脚本里处理）

1. **`modalMass`/`modalDamp` 必须列满 `mode/` 里的全部阶数**（4 个），`nActiveModes` 才决定叠加几阶。
   只写 2 个会直接 fatal：`modalMass size (2) does not match the number of available modes (4)`。
2. **并行稳态解不会写根目录时间文件夹**：`processor*/1000` 需要 `reconstructPar -latestTime` 才能取回。
3. **冒烟档拿不到 σ**：0.02 s 不足一个结构周期（一阶约 0.1 s），阶段 06 会明确报
   `INSUFFICIENT DATA`；要 σ/f 至少 `--end-time 0.3`。

## 实测（2026-11-04，本地 8 核）

```
01  5 s     02  2 s     03  24.9 min     04  10 s     05  65.6 min     06  4 s
合计 1 h 32 min，20/20 步，4 次 AMR 加密，0 FATAL，平均 196 s/步
```

稳态 1000 次伪迭代约 1.5 s/次。
