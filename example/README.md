# CALVADOS 多尺度模拟示例

这个文件夹包含了使用 `multiscale2` 软件运行CALVADOS粗粒化模拟的示例。

## 文件夹结构

```
example/
├── README.md                    # 本说明文件
└── IDP/                        # IDP任务示例
    ├── config.yaml             # IDP配置文件
    ├── protein.fasta           # 蛋白质序列文件
    ├── run_example.py          # IDP运行脚本
    └── USAGE.md                # IDP使用说明
```

## 任务类型

### IDP (Intrinsically Disordered Proteins)
- **用途**: 模拟内在无序蛋白质区域
- **特点**: 使用CALVADOS2的residues文件
- **配置**: 自动生成components.yaml文件

### MDP (Membrane-Disrupting Peptides)
- **用途**: 模拟膜破坏肽
- **特点**: 使用CALVADOS3的residues文件
- **配置**: 需要domains.yaml文件定义功能域
- **状态**: 待实现

## 快速开始

### 运行IDP任务
```bash
cd example/IDP
conda activate cg2all
export CUDA_VISIBLE_DEVICES=0
python run_example.py
```

## 重要说明

1. **components.yaml**: 该文件由CALVADOS的 `prepare.py` 脚本自动生成，不需要预先提供
2. **路径处理**: 所有路径都是相对路径，确保在不同电脑上的可移植性
3. **环境要求**: 需要CALVADOS环境和CUDA支持

## 文件说明

### 配置文件
- `config.yaml`: 主配置文件，包含所有模拟参数
- `protein.fasta`: 蛋白质序列文件（可替换）

### 脚本文件
- `run_example.py`: 交互式运行脚本
- `USAGE.md`: 详细使用说明

## 输出文件

运行完成后会生成：
- `prepare.py`: CALVADOS准备脚本
- `run.py`: CALVADOS运行脚本
- `components.yaml`: CALVADOS组件文件（自动生成）
- `*.dcd`: 轨迹文件
- `*.pdb`: 结构文件
- `*.log`: 日志文件

## 注意事项

1. 确保CALVADOS环境正确安装
2. 确保CUDA可用
3. 根据系统大小调整模拟参数
4. 大系统可能需要更长的模拟时间
