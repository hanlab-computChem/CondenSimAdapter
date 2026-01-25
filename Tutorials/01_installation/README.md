# 安装教程（01_installation）

本教程面向 `CondenSimAdapter` 的本地安装与验证，覆盖 adapter `cg/backmap/minimize` 的核心依赖，并以“运行并通过 tests”作为最后一步。

## 1) 前置要求

- 建议使用 **Python 3.11**（与 `install.sh`、`requirements.txt` 保持一致）。
- backmap 相关依赖默认使用 **CPU 版** PyTorch/DGL（便于安装与验证）。
- OpenMM 仍由 Conda 安装，保持 GPU 可用性（如系统已安装 CUDA 驱动）。
- 建议使用 Conda 环境。
- 如需完整 `minimize` 流程，请确保系统已安装 **GROMACS**（`gromacswrapper` 仅是 Python 绑定）。

> 注意：仓库规则里提示 `conda activate multiscale2`，而 `install.sh` 使用 `CondenSimAdapter`。建议优先使用 `CondenSimAdapter`，若你已有 `multiscale2`，可将脚本中的环境名统一为现有环境，确保 Python 版本为 3.11。

## 2) 创建并进入环境

```bash
conda create -n CondenSimAdapter python=3.11 -y
conda activate CondenSimAdapter
```

## 3) 安装依赖

### 3.1 使用安装脚本（推荐）

```bash
bash install.sh
```

脚本会安装以下核心依赖（与 adapter 相关）：
- **cg**：`openmm`, `mdtraj`, `mdanalysis`, `numpy`, `pandas`, `scipy`
- **backmap**：`torch==2.1.2`, `dgl==2.1.0`, `e3nn`, `se3-transformer`, `ml-collections`
- **minimize**：`openmm`, `gromacswrapper`（外部 GROMACS 需自行安装）

其中 OpenMM 默认版本为 **8.4**。

### 3.2 备用：requirements.txt

若不使用脚本，可参考 `requirements.txt` 手动安装。注意 `dgl` 与 `se3-transformer` 需要特殊安装源。

## 4) 验证安装（可选）

```bash
python -c "import torch; print(torch.__version__)"
python -c "import openmm; print(openmm.__version__)"
python -c "import MDAnalysis; print(MDAnalysis.__version__)"
```

## 5) 最后一步：运行 tests 并通过

测试主要覆盖 adapter `cg/backmap/minimize` 的轻量烟雾用例。

```bash
pytest -q
```

如仅运行 adapter 相关测试：

```bash
pytest -q tests/test_adapter_*.py
```

说明：
- backmap 模型权重已在仓库内提供：`CondenSimAdapter/extern/ms2_cg2all/model/Martini3-FIX.ckpt`
- tests 只做轻量验证，不会执行长时间模拟或 GPU 计算
