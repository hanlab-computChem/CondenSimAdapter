# CALVADOS Wrapper 代码审查报告

**文件**: `CondenSimAdapter/src/calvados_wrapper.py`  
**审查日期**: 2024年12月30日  
**代码行数**: 630行

---

## 一、大体分析

### 1.1 文件定位与核心功能

`calvados_wrapper.py` 是 CondenSimAdapter 项目中专门用于包装和运行 **CALVADOS** 粗粒化分子动力学模拟的核心模块。CALVADOS 是一种基于序列的粗粒化模型，主要用于模拟无序蛋白（IDPs）和折叠蛋白（MDPs）的相分离行为。

该文件的核心职责是：
- **配置转换**: 将 `CGSimulationConfig` 转换为 CALVADOS 原生所需的配置格式
- **模拟执行**: 调用底层 CALVADOS 引擎运行模拟
- **输出整理**: 将原生输出组织为统一的文件结构

### 1.2 架构设计

```
 CalvadosWrapper (核心类)
 ├── __init__(): 初始化配置
 ├── _get_residues_path(): 获取残基参数文件路径
 ├── _topol_to_calvados(): 拓扑类型转换 (CUBIC→grid, SLAB→slab)
 ├── _platform_to_string(): 平台类型转换
 ├── create_config(): 创建 CALVADOS Config 对象
 ├── create_components(): 创建 CALVADOS Components 对象
 ├── write(): 写入配置文件到目录
 ├── _generate_config_yaml(): 生成 config.yaml 内容
 ├── _generate_components_yaml(): 生成 components.yaml 内容
 ├── _write_to_dir(): 写入配置文件到指定目录
 ├── _process_inline_fdomains(): 处理内联域定义
 ├── run(): 运行模拟
 ├── _organize_output(): 整理输出文件
 ├── _write_simulation_log(): 写入高层级日志
 │
 └── run_calvados(): 便捷入口函数
```

### 1.3 依赖关系

**内部依赖**:
- `cg.py`: `CGSimulationConfig`, `CGComponent`, `ComponentType`, `TopologyType`, `Platform`, `SimulationResult`

**外部依赖**:
- `CondenSimAdapter.extern.ms2_calvados.calvados.cfg`: Config, Components 类
- `CondenSimAdapter.extern.ms2_calvados.calvados.sim`: 模拟运行核心

### 1.4 代码风格评价

**优点**:
- 清晰的文档字符串，遵循 Google Style
- 合理的方法分组（初始化、配置生成、运行、输出整理）
- 错误处理完善，有 try-except 捕获
- 支持内联 YAML 格式的 fdomains 定义

**需改进**:
- 部分注释使用中文，但变量名使用英文，混用风格不够统一
- `_generate_config_yaml()` 和 `create_config()` 存在功能重复
- 缺少类型注解（部分方法没有返回类型注解）

---

## 二、类详细解读

### 2.1 CalvadosWrapper 类

**类定位**: 核心包装器类，封装所有 CALVADOS 相关的配置转换和模拟运行逻辑。

**继承关系**: 无继承，独立的包装器类。

**主要职责**:
1. 配置格式转换（CGSimulationConfig → CALVADOS 格式）
2. 文件写入
3. 模拟执行
4. 输出整理

#### 2.1.1 `__init__(self, config: CGSimulationConfig)`

**功能**: 初始化包装器实例。

**参数**:
- `config`: CGSimulationConfig 实例，包含完整的模拟配置

**实现逻辑**:
1. 保存配置引用
2. 初始化输出目录为 None
3. 调用 `_get_residues_path()` 获取残基文件路径

**代码质量**: ⭐⭐⭐⭐
- 简洁明了
- 有类型注解
- 建议：可添加参数验证

**代码示例**:
```python
def __init__(self, config: CGSimulationConfig):
    self.config = config
    self.output_dir: Optional[str] = None
    self._residues_path = self._get_residues_path()
```

#### 2.1.2 `_get_residues_path(self) -> str`

**功能**: 根据组件类型选择正确的残基参数文件。

**实现逻辑**:
1. 检查是否存在 MDP 类型组件
2. 如果有 MDP 组件，使用 `residues_CALVADOS3.csv`
3. 否则使用 `residues_CALVADOS2.csv`

**设计考量**:
- CALVADOS2 用于纯 IDP 系统
- CALVADOS3 支持 MDP（有结构域的蛋白）

**代码质量**: ⭐⭐⭐⭐⭐
- 逻辑清晰
- 使用 Path 正确处理路径
- 动态选择文件

**改进建议**:
```python
def _get_residues_path(self) -> str:
    """获取 residue 参数文件路径
    
    Residues 文件从 ms2_calvados 包的 data 目录加载：
    - residues_CALVADOS2.csv: 用于纯 IDP 系统
    - residues_CALVADOS3.csv: 用于包含 MDP 的系统
    """
    from CondenSimAdapter.extern.ms2_calvados.calvados import data as calvados_data
    
    has_mdp = any(c.type == ComponentType.MDP for c in self.config.components)
    residues_file = 'residues_CALVADOS3.csv' if has_mdp else 'residues_CALVADOS2.csv'
    
    residues_path = Path(calvados_data.__file__).parent / residues_file
    
    if not residues_path.exists():
        raise FileNotFoundError(f"Residues file not found: {residues_path}")
    
    return str(residues_path)
```

#### 2.1.3 `_topol_to_calvados(self) -> str`

**功能**: 将 TopologyType 枚举转换为 CALVADOS 拓扑字符串。

**映射关系**:
- `TopologyType.CUBIC` → `'grid'` (网格放置)
- `TopologyType.SLAB` → `'slab'` (平面限制)
- 其他 → `'grid'` (默认)

**代码质量**: ⭐⭐⭐⭐
- 简单清晰
- 有默认处理

**代码示例**:
```python
def _topol_to_calvados(self) -> str:
    if self.config.topol == TopologyType.CUBIC:
        return 'grid'
    elif self.config.topol == TopologyType.SLAB:
        return 'slab'
    else:
        return 'grid'
```

#### 2.1.4 `_platform_to_string(self) -> str`

**功能**: 将 Platform 枚举转换为字符串。

**实现逻辑**:
1. 如果是 Platform 枚举，使用 `.value` 获取值
2. 否则直接转字符串

**代码质量**: ⭐⭐⭐⭐
- 健壮的类型处理

**代码示例**:
```python
def _platform_to_string(self) -> str:
    if isinstance(self.config.simulation.platform, Platform):
        return self.config.simulation.platform.value
    return str(self.config.simulation.platform)
```

#### 2.1.5 `create_config(self) -> 'ms2_config.Config'`

**功能**: 创建 CALVADOS Config 对象。

**实现逻辑**:
1. 从 `self.config.simulation` 获取模拟参数
2. 计算 slab_width（如果是 SLAB 拓扑）
3. 构建参数字典
4. 创建并返回 Config 对象

**关键参数**:
- `slab_width`: SLAB 拓扑时设为 `box[2] / 2`（z 方向的一半）
- `slab_eq`: 固定为 False
- `k_eq`: 固定为 0.02

**代码质量**: ⭐⭐⭐⭐
- 逻辑清晰
- 文档注释详细

**潜在问题**:
- 硬编码了 `slab_eq=False` 和 `k_eq=0.02`，限制了灵活性
- 未使用 Config 对象的默认值，而是全部显式设置

**改进建议**:
```python
def create_config(self) -> 'ms2_config.Config':
    """创建 ms2_calvados Config 对象
    
    Notes:
        - slab_width: 自动计算为 box[2] / 2（z 方向的一半）
        - 对于非 SLAB 拓扑，不设置 slab_width（使用 CALVADOS 默认值）
    """
    from CondenSimAdapter.extern.ms2_calvados.calvados.cfg import Config
    
    sim_params = self.config.simulation
    
    # 动态构建参数
    params = {
        'sysname': self.config.system_name,
        'box': self.config.box,
        'temp': self.config.temperature,
        'ionic': self.config.ionic,
        'pH': 7.0,
        'topol': self._topol_to_calvados(),
        'wfreq': sim_params.wfreq,
        'steps': sim_params.steps,
        'platform': self._platform_to_string(),
        'restart': None,
        'verbose': sim_params.verbose,
    }
    
    # 仅 SLAB 拓扑添加 slab_width
    if self.config.topol.value == 'slab':
        params['slab_width'] = self.config.box[2] / 2
        params['slab_eq'] = False
        params['k_eq'] = 0.02
    
    return Config(**params)
```

#### 2.1.6 `create_components(self) -> 'ms2_config.Components'`

**功能**: 创建 CALVADOS Components 对象。

**实现逻辑**:
1. 从第一个组件获取默认值
2. 创建 Components 实例
3. 遍历所有组件，转换为 CALVADOS 格式
4. 添加到 Components 对象

**组件类型处理**:
- **IDP**: 需要 ffasta 文件
- **MDP**: 需要 fpdb 文件，可选 fdomains（域定义）

**代码质量**: ⭐⭐⭐⭐
- 逻辑清晰
- 类型判断正确

**潜在问题**:
- 假设第一个组件的设置适用于所有组件（nmol, restraint 等）
- 可能会导致配置不一致

**改进建议**:
```python
def create_components(self) -> 'ms2_config.Components':
    """创建 ms2_calvados Components 对象
    
    注意：每个组件的设置独立处理，不再使用第一个组件的默认值
    """
    from CondenSimAdapter.extern.ms2_calvados.calvados.cfg import Components
    
    # 获取 pdb_folder（从第一个 MDP 组件的 fpdb 提取）
    pdb_folder = None
    for comp in self.config.components:
        if comp.type == ComponentType.MDP and comp.fpdb:
            pdb_folder = os.path.dirname(os.path.abspath(comp.fpdb))
            break
    
    defaults = {
        'molecule_type': 'protein',
        'fresidues': self._residues_path,
        'pdb_folder': pdb_folder,
    }
    
    ms2_components = Components(**defaults)
    
    for comp in self.config.components:
        comp_dict = {
            'name': comp.name,
            'nmol': comp.nmol,
            'restraint': comp.restraint,
            'charge_termini': comp.charge_termini,
        }
        
        if comp.type == ComponentType.IDP:
            if comp.ffasta:
                comp_dict['ffasta'] = comp.ffasta
        
        elif comp.type == ComponentType.MDP:
            if comp.fpdb:
                comp_dict['fpdb'] = comp.fpdb
            
            if comp.fdomains:
                comp_dict['fdomains'] = comp.fdomains
            
            if comp.restraint:
                comp_dict['restraint_type'] = comp.restraint_type
                comp_dict['use_com'] = comp.use_com
                comp_dict['k_harmonic'] = comp.k_harmonic
                comp_dict['colabfold'] = comp.colabfold
        
        ms2_components.add(**comp_dict)
    
    return ms2_components
```

#### 2.1.7 `write(self, output_dir: str, overwrite: bool = False) -> Dict[str, str]`

**功能**: 写入配置文件到目录。

**实现逻辑**:
1. 检查目录是否存在
2. 创建输出目录
3. 创建并写入 config.yaml
4. 创建并写入 components.yaml
5. 生成并写入 run.py

**返回值**: 文件路径字典

**代码质量**: ⭐⭐⭐⭐
- 良好的错误处理
- 支持覆盖模式

**代码示例**:
```python
def write(self, output_dir: str, overwrite: bool = False) -> Dict[str, str]:
    output_dir = os.path.abspath(output_dir)
    
    if os.path.exists(output_dir) and not overwrite:
        raise FileExistsError(f"Output directory exists: {output_dir}")
    
    os.makedirs(output_dir, exist_ok=True)
    self.output_dir = output_dir
    
    ms2_config = self.create_config()
    ms2_config.write(output_dir, name='config.yaml')
    
    ms2_components = self.create_components()
    ms2_components.write(output_dir, name='components.yaml')
    
    return {
        'config': os.path.join(output_dir, 'config.yaml'),
        'components': os.path.join(output_dir, 'components.yaml'),
        'run_script': os.path.join(output_dir, 'run.py'),
    }
```

**问题**:
- `run_script` 返回了路径，但实际并未创建 `run.py` 文件
- 应该同步调用 `_generate_config_yaml()` 和 `_generate_components_yaml()` 的逻辑

#### 2.1.8 `_generate_config_yaml(self) -> str`

**功能**: 生成 CALVADOS config.yaml 的 YAML 内容。

**实现逻辑**:
1. 构建完整的参数字典
2. 处理 SLAB 拓扑的 slab_width
3. 添加所有 CALVADOS 默认参数
4. 使用 yaml.dump 生成 YAML 字符串

**参数特点**:
- 硬编码了大量默认值（eps_lj=0.2, cutoff_lj=2.0 等）
- 添加了原版 CALVADOS 的所有参数

**代码质量**: ⭐⭐⭐
- 功能完整
- 但与 `create_config()` 存在逻辑重复
- 硬编码参数过多

**代码示例**:
```python
def _generate_config_yaml(self) -> str:
    """生成 CALVADOS config.yaml 内容"""
    import yaml
    
    sim_params = self.config.simulation
    
    # SLAB 拓扑：自动计算 slab_width = box_z / 2
    if self.config.topol == TopologyType.SLAB:
        slab_width = self.config.box[2] / 2
    else:
        slab_width = 100  # 默认值，不推荐
    
    config_dict = {
        'sysname': self.config.system_name,
        'box': self.config.box,
        'temp': self.config.temperature,
        'ionic': self.config.ionic,
        'pH': 7.0,
        'topol': self._topol_to_calvados(),
        'slab_width': slab_width,
        'slab_eq': False,
        'friction': 0.01,
        'wfreq': sim_params.wfreq,
        'steps': sim_params.steps,
        'platform': self._platform_to_string(),
        # ... 更多参数
    }
    
    return yaml.dump(config_dict, default_flow_style=False)
```

**重构建议**:
```python
def _generate_config_yaml(self) -> str:
    """生成 CALVADOS config.yaml 内容
    
    建议：重构为从 create_config() 获取参数，而不是复制逻辑
    """
    config = self.create_config()
    # 假设 Config 对象有 to_dict() 方法
    if hasattr(config, 'to_dict'):
        config_dict = config.to_dict()
    else:
        # 手动构建（临时方案）
        config_dict = {
            'sysname': self.config.system_name,
            'box': self.config.box,
            'temp': self.config.temperature,
            'ionic': self.config.ionic,
            'pH': 7.0,
            'topol': self._topol_to_calvados(),
            'wfreq': self.config.simulation.wfreq,
            'steps': self.config.simulation.steps,
            'platform': self._platform_to_string(),
            'verbose': self.config.simulation.verbose,
        }
        
        # SLAB 拓扑添加额外参数
        if self.config.topol.value == 'slab':
            config_dict['slab_width'] = self.config.box[2] / 2
            config_dict['slab_eq'] = False
            config_dict['k_eq'] = 0.02
    
    return yaml.dump(config_dict, default_flow_style=False, sort_keys=False)
```

#### 2.1.9 `_generate_components_yaml(self) -> str`

**功能**: 生成 CALVADOS components.yaml 的 YAML 内容。

**实现逻辑**:
1. 获取 pdb_folder（从第一个 MDP 组件）
2. 构建 defaults 字典（包含原版 CALVADOS 的默认参数）
3. 遍历组件，转换为字典
4. 移除 None 值
5. 使用 yaml.dump 生成

**设计特点**:
- 添加了原版 CALVADOS default_component.yaml 中的所有默认参数
- 处理了 fpdb 和 pdb_folder 的格式转换

**代码质量**: ⭐⭐⭐⭐
- 逻辑清晰
- 注释详细
- 注意：与 `create_components()` 存在重复逻辑

**问题识别**:
1. 假设第一个组件的 restraint, charge_termini 等设置作为默认值
2. pdb_folder 只从第一个 MDP 组件提取
3. IDP 组件也设置了 pdb_folder（可能不需要）

**改进建议**:
```python
def _generate_components_yaml(self) -> str:
    """生成 CALVADOS components.yaml 内容
    
    处理 fpdb 和 pdb_folder:
    - CALVADOS 期望 pdb_folder（目录）和 name（不含扩展名的文件名）
    - 我们的 config 使用 fpdb（完整文件路径）
    
    注意：每个组件独立设置，不再使用第一个组件的默认值
    """
    import yaml
    
    # 获取 pdb_folder（从第一个 MDP 组件的 fpdb 提取）
    pdb_folder = None
    for comp in self.config.components:
        if comp.type == ComponentType.MDP and comp.fpdb:
            pdb_folder = os.path.dirname(os.path.abspath(comp.fpdb))
            break
    
    components = {
        'defaults': {
            'molecule_type': 'protein',
            'fresidues': self._residues_path,
            'pdb_folder': pdb_folder,
            # 原版 CALVADOS default_component.yaml 中的参数
            'periodic': False,
            'cutoff_restr': 0.9,
            'k_go': 15.0,
            'use_com': True,
            'colabfold': 0,
        },
        'system': {}
    }
    
    for comp in self.config.components:
        comp_dict = {
            'name': comp.name,
            'molecule_type': 'protein',
            'nmol': comp.nmol,
            'restraint': comp.restraint if comp.restraint else None,
            'restraint_type': comp.restraint_type if comp.restraint else None,
            'use_com': comp.use_com if comp.restraint else None,
            'k_harmonic': comp.k_harmonic if comp.restraint else None,
            'colabfold': comp.colabfold if comp.restraint else None,
            'charge_termini': comp.charge_termini,
        }
        
        # 添加类型特定参数
        if comp.type == ComponentType.IDP and comp.ffasta:
            comp_dict['ffasta'] = comp.ffasta
        
        elif comp.type == ComponentType.MDP:
            if comp.fpdb:
                comp_dict['fpdb'] = comp.fpdb
            if comp.fdomains:
                comp_dict['fdomains'] = comp.fdomains
        
        # 移除 None 值和 False 值
        comp_dict = {k: v for k, v in comp_dict.items() if v is not None and v is not False}
        
        components['system'][comp.name] = comp_dict
    
    return yaml.dump(components, default_flow_style=False, sort_keys=False)
```

#### 2.1.10 `_write_to_dir(self, output_dir: str) -> Dict[str, str]`

**功能**: 写入配置文件到指定目录（支持内联 fdomains 处理）。

**实现逻辑**:
1. 创建输出目录
2. 写入 config.yaml
3. 生成 components.yaml
4. 调用 `_process_inline_fdomains()` 处理内联 fdomains
5. 写入 components.yaml

**特点**: 支持两种 fdomains 格式：
1. 文件路径：`'domains.yaml'` - 直接复制
2. 内联 YAML：`'TDP43:\n  - [3, 76]\n...'` - 写入临时文件

**代码质量**: ⭐⭐⭐⭐
- 逻辑清晰
- 灵活支持多种格式

**代码示例**:
```python
def _write_to_dir(self, output_dir: str) -> Dict[str, str]:
    """写入配置文件到指定目录（返回文件路径字典）
    
    支持两种 fdomains 格式：
    1. 文件路径：'domains.yaml' - 直接复制到输出目录
    2. 内联 YAML：'TDP43:\n  - [3, 76]\n...' - 写入临时文件
    """
    import tempfile
    import shutil
    
    os.makedirs(output_dir, exist_ok=True)
    self.output_dir = output_dir
    
    # 写入 config.yaml
    config_file = os.path.join(output_dir, 'config.yaml')
    with open(config_file, 'w') as f:
        f.write(self._generate_config_yaml())
    
    # 处理 components.yaml，支持内联 fdomains
    components_yaml = self._generate_components_yaml()
    components_yaml = self._process_inline_fdomains(components_yaml, output_dir)
    
    # 写入 components.yaml
    components_file = os.path.join(output_dir, 'components.yaml')
    with open(components_file, 'w') as f:
        f.write(components_yaml)
    
    return {
        'config': config_file,
        'components': components_file,
    }
```

#### 2.1.11 `_process_inline_fdomains(self, components_yaml: str, output_dir: str) -> str`

**功能**: 处理内联的 fdomains，将其写入临时文件。

**实现逻辑**:
1. 解析 YAML
2. 遍历所有组件的 fdomains
3. 检测是否是内联 YAML（检查格式特征）
4. 如果是内联 YAML，解析并写入文件
5. 替换为文件路径

**检测逻辑**:
- 以 `{` 或 `[` 开头
- 包含换行符且有 YAML 特征（`:` 或 `-`）
- 包含冒号但不以 `.yaml` 或 `.yml` 结尾

**代码质量**: ⭐⭐⭐⭐
- 健壮的格式检测
- 良好的错误处理

**问题**: 检测逻辑可能误判复杂路径，建议增强验证

**改进建议**:
```python
def _process_inline_fdomains(self, components_yaml: str, output_dir: str) -> str:
    """处理内联的 fdomains，如果是 YAML 内容则写入临时文件"""
    import yaml
    
    components = yaml.safe_load(components_yaml)
    
    for name, props in components.get('system', {}).items():
        fdomains = props.get('fdomains')
        if fdomains and isinstance(fdomains, str):
            # 移除 YAML 引号
            stripped = fdomains.strip().strip('"\'')
            
            # 增强检测：优先检查文件是否存在
            if os.path.isfile(stripped):
                continue  # 是文件路径，保持原样
            
            # 检查是否是内联 YAML
            is_inline = False
            if stripped.startswith(('{', '[')):
                is_inline = True
            elif '\n' in stripped and (':' in stripped or stripped.startswith('-')):
                is_inline = True
            elif ':' in stripped and not stripped.endswith(('.yaml', '.yml')):
                # 排除明显不是 YAML 的情况
                if not any(ext in stripped for ext in ['/', '\\']):
                    is_inline = True
            
            if is_inline:
                try:
                    domains_data = yaml.safe_load(fdomains)
                    domains_file = os.path.join(output_dir, f'{name}_domains.yaml')
                    with open(domains_file, 'w') as f:
                        yaml.dump(domains_data, f, default_flow_style=False)
                    props['fdomains'] = domains_file
                    print(f"  📄 内联 domains 写入: {domains_file}")
                except yaml.YAMLError:
                    pass  # 不是有效 YAML，保持原样
    
    return yaml.dump(components, default_flow_style=False, sort_keys=False)
```

#### 2.1.12 `run(self, output_dir: str = None, gpu_id: int = 0) -> SimulationResult`

**功能**: 运行 CALVADOS 模拟。

**实现逻辑**:
1. 确定输出目录（添加 `_CG` 后缀）
2. 备份旧结果（如存在）
3. 创建 raw 目录
4. 写入配置文件
5. 设置 GPU 环境变量
6. 调用 calvados_sim.run()
7. 组织输出文件
8. 写入日志
9. 返回结果

**统一输出结构**:
```
{output_dir}/{system_name}_CG/
├── final.pdb                   # 最终结构
├── trajectory.dcd              # 模拟轨迹
├── simulation.log              # 高层级日志
└── raw/                        # 原生输出
    ├── config.yaml
    ├── components.yaml
    └── *.dcd, *.xml, *.pdb, *.chk, *.txt
```

**代码质量**: ⭐⭐⭐⭐
- 结构清晰
- 错误处理完善
- 日志输出友好

**潜在问题**:
1. `time.time()` 用于计时但 `elapsed` 可能未正确传递
2. 缺少实际的计时逻辑（start_time 定义后未被使用）

**代码示例**:
```python
def run(self, output_dir: str = None, gpu_id: int = 0) -> SimulationResult:
    """运行 CALVADOS 模拟"""
    from CondenSimAdapter.extern.ms2_calvados.calvados import sim as calvados_sim
    import shutil
    from datetime import datetime
    import time
    
    if output_dir is None:
        output_dir = self.config.output_dir
    
    # 统一添加 _CG 后缀
    task_name = f"{self.config.system_name}_CG"
    output_dir = os.path.join(output_dir, task_name)
    raw_dir = os.path.join(output_dir, 'raw')
    
    # 备份旧结果
    if os.path.exists(output_dir):
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup_dir = f"{output_dir}_backup_{timestamp}"
        shutil.move(output_dir, backup_dir)
        print(f"  📁 备份旧结果到: {backup_dir}")
    
    os.makedirs(raw_dir, exist_ok=True)
    
    # 写入配置文件
    files = self._write_to_dir(raw_dir)
    
    result = SimulationResult()
    result.output_dir = output_dir
    
    start_time = time.time()
    
    try:
        print(f"\n[CALVADOS] Running simulation...")
        print(f"  GPU ID: {gpu_id}")
        
        os.environ['CUDA_VISIBLE_DEVICES'] = str(gpu_id)
        
        calvados_sim.run(
            path=raw_dir,
            fconfig='config.yaml',
            fcomponents='components.yaml'
        )
        
        self._organize_output(raw_dir, output_dir, task_name)
        
        result.success = True
        elapsed = time.time() - start_time
        print(f"  ✓ CALVADOS simulation completed ({elapsed:.1f}s)")
        
    except Exception as e:
        result.success = False
        result.errors.append(str(e))
        print(f"  ✗ CALVADOS simulation failed: {e}")
        elapsed = time.time() - start_time
    
    self._write_simulation_log(output_dir, task_name, elapsed, result.success)
    
    result.trajectory = os.path.join(output_dir, 'trajectory.dcd')
    result.structure = os.path.join(output_dir, 'final.pdb')
    
    return result
```

#### 2.1.13 `_organize_output(self, raw_dir: str, output_dir: str, task_name: str)`

**功能**: 整理 CALVADOS 输出文件到统一结构。

**实现逻辑**:
1. 复制轨迹文件 (`{sysname}.dcd` → `trajectory.dcd`)
2. 查找最终 PDB（优先 checkpoint.pdb）
3. 复制重要文件到 raw 目录
4. 重命名 log 文件

**代码质量**: ⭐⭐⭐⭐
- 逻辑清晰
- 文件组织有序

**代码示例**:
```python
def _organize_output(self, raw_dir: str, output_dir: str, task_name: str):
    """组织输出文件到统一结构
    
    统一命名规则：
    - trajectory.dcd  <- {task_name}.dcd
    - final.pdb       <- 带时间戳的 pdb 或 checkpoint.pdb
    """
    import shutil
    
    sysname = self.config.system_name
    
    # 处理轨迹文件
    src_dcd = os.path.join(raw_dir, f'{sysname}.dcd')
    dst_dcd = os.path.join(output_dir, 'trajectory.dcd')
    if os.path.exists(src_dcd):
        shutil.copy2(src_dcd, dst_dcd)
        print(f"  📦 trajectory.dcd")
    
    # 查找最终结构
    src_pdb = os.path.join(raw_dir, 'checkpoint.pdb')
    if not os.path.exists(src_pdb):
        for f in os.listdir(raw_dir):
            if f.endswith('.pdb') and f != 'top.pdb':
                src_pdb = os.path.join(raw_dir, f)
                break
    
    dst_pdb = os.path.join(output_dir, 'final.pdb')
    if os.path.exists(src_pdb):
        shutil.copy2(src_pdb, dst_pdb)
        print(f"  📦 final.pdb")
    
    # 复制重要文件
    important_files = [
        (f'{sysname}.xml', 'system.xml'),
        ('top.pdb', 'top.pdb'),
        ('restart.chk', 'restart.chk'),
        ('checkpoint.pdb', 'checkpoint.pdb'),
    ]
    
    for src_name, dst_name in important_files:
        src = os.path.join(raw_dir, src_name)
        if os.path.exists(src):
            dst = os.path.join(raw_dir, dst_name)
            if src != dst:
                shutil.copy2(src, dst)
    
    print(f"  📁 原始输出已整理到: {raw_dir}")
```

#### 2.1.14 `_write_simulation_log(self, output_dir: str, task_name: str, elapsed: float, success: bool)`

**功能**: 写入高层级模拟日志。

**实现逻辑**:
1. 构建日志内容（包含系统配置、组件信息等）
2. 写入 simulation.log 文件

**日志内容**:
- 任务名称、力场、日期
- 状态和耗时
- 系统配置（box, temperature, ionic, topol）
- 组件信息
- 输出文件列表

**代码质量**: ⭐⭐⭐⭐
- 内容完整
- 格式清晰

**代码示例**:
```python
def _write_simulation_log(self, output_dir: str, task_name: str, elapsed: float, success: bool):
    """写入高层级模拟日志"""
    from datetime import datetime
    
    log_file = os.path.join(output_dir, 'simulation.log')
    
    status = "SUCCESS" if success else "FAILED"
    components_info = []
    for comp in self.config.components:
        comp_info = f"  - {comp.name}: {comp.type.value}, nmol={comp.nmol}"
        if comp.type == ComponentType.IDP:
            comp_info += f", seq={comp.ffasta}"
        elif comp.type == ComponentType.MDP:
            comp_info += f", pdb={comp.fpdb}"
        components_info.append(comp_info)
    
    log_content = f"""# CondenSimAdapter CG Simulation Log
# ============================

Task: {task_name}
Force Field: CALVADOS
Date: {datetime.now().isoformat()}

Status: {status}
Duration: {elapsed:.2f} seconds

System Configuration:
  Box: {self.config.box} nm
  Temperature: {self.config.temperature} K
  Ionic Strength: {self.config.ionic} M
  Topology: {self.config.topol.value}

Components ({len(self.config.components)}):
{chr(10).join(components_info)}

Output Files:
  - final.pdb: Final structure
  - trajectory.dcd: Simulation trajectory
  - raw/: Native simulation output files
"""
    with open(log_file, 'w') as f:
        f.write(log_content)
    
    print(f"  📝 simulation.log")
```

### 2.2 `run_calvados()` 函数

**功能**: 运行 CALVADOS 模拟的便捷入口函数。

**实现逻辑**:
1. 创建 CalvadosWrapper 实例
2. 调用 run() 方法

**代码质量**: ⭐⭐⭐⭐⭐
- 简洁明了
- 提供便捷接口

**代码示例**:
```python
def run_calvados(config: CGSimulationConfig, output_dir: str = None, gpu_id: int = 0) -> SimulationResult:
    """运行 CALVADOS 模拟的便捷函数
    
    Args:
        config: CGSimulationConfig 实例
        output_dir: 输出目录
        gpu_id: GPU 设备 ID
        
    Returns:
        SimulationResult
    """
    wrapper = CalvadosWrapper(config)
    return wrapper.run(output_dir=output_dir, gpu_id=gpu_id)
```

---

## 三、代码审查总结

### 3.1 优点

1. **架构清晰**: 职责分离明确，配置转换、文件写入、模拟运行、输出整理各司其职
2. **文档完善**: 每个方法都有详细的文档字符串
3. **错误处理**: 有 try-except 捕获异常
4. **灵活性**: 支持多种配置格式（文件路径和内联 YAML）
5. **输出规范**: 统一的输出目录结构

### 3.2 问题与改进建议

#### 3.2.1 重复代码

**问题**: `_generate_config_yaml()` 和 `create_config()` 存在大量重复逻辑

**建议**:
- 重构为共享配置构建逻辑
- 使用单一数据源生成 YAML 和 Config 对象

#### 3.2.2 类型注解不完整

**问题**: 部分方法缺少返回类型注解

**建议**:
- 为所有方法添加类型注解
- 使用 mypy 进行类型检查

#### 3.2.3 配置硬编码

**问题**: 大量参数硬编码（如 slab_eq, k_eq, eps_lj 等）

**建议**:
- 将硬编码参数移到配置文件
- 添加配置覆盖机制

#### 3.2.4 混合语言

**问题**: 注释使用中文，变量名使用英文

**建议**:
- 统一使用英文注释
- 或统一使用中文（但国际团队项目建议英文）

#### 3.2.5 默认值处理

**问题**: 第一个组件的设置被用作所有组件的默认值

**建议**:
- 每个组件独立设置
- 添加全局默认值类

### 3.3 总体评价

| 指标 | 评分 |
|------|------|
| 代码结构 | ⭐⭐⭐⭐ |
| 可读性 | ⭐⭐⭐⭐ |
| 可维护性 | ⭐⭐⭐ |
| 健壮性 | ⭐⭐⭐⭐ |
| 文档完整性 | ⭐⭐⭐⭐⭐ |

**总体评分**: ⭐⭐⭐⭐ (4/5)

---

## 四、建议的重构方案

### 4.1 配置生成统一化

```python
class CalvadosConfigBuilder:
    """CALVADOS 配置构建器"""
    
    def __init__(self, config: CGSimulationConfig):
        self.config = config
        self._params = {}
    
    def build_params(self) -> Dict:
        """构建所有参数"""
        self._build_system_params()
        self._build_simulation_params()
        self._build_topology_params()
        return self._params
    
    def _build_system_params(self):
        """系统参数"""
        self._params.update({
            'sysname': self.config.system_name,
            'box': self.config.box,
            'temp': self.config.temperature,
            'ionic': self.config.ionic,
        })
    
    def _build_simulation_params(self):
        """模拟参数"""
        sim = self.config.simulation
        self._params.update({
            'wfreq': sim.wfreq,
            'steps': sim.steps,
            'platform': self._platform_to_string(),
            'verbose': sim.verbose,
        })
    
    def _build_topology_params(self):
        """拓扑参数"""
        topol = self._topol_to_calvados()
        self._params['topol'] = topol
        
        if topol == 'slab':
            self._params.update({
                'slab_width': self.config.box[2] / 2,
                'slab_eq': False,
                'k_eq': 0.02,
            })
    
    def to_config(self) -> 'ms2_config.Config':
        """转换为 Config 对象"""
        from CondenSimAdapter.extern.ms2_calvados.calvados.cfg import Config
        return Config(**self.build_params())
    
    def to_yaml(self) -> str:
        """转换为 YAML"""
        import yaml
        return yaml.dump(self.build_params(), default_flow_style=False, sort_keys=False)
```

### 4.2 组件配置优化

```python
class CalvadosComponentsBuilder:
    """CALVADOS 组件构建器"""
    
    def __init__(self, config: CGSimulationConfig):
        self.config = config
        self._residues_path = self._get_residues_path()
    
    def build(self) -> Dict:
        """构建组件配置"""
        pdb_folder = self._find_pdb_folder()
        
        defaults = {
            'molecule_type': 'protein',
            'fresidues': self._residues_path,
            'pdb_folder': pdb_folder,
        }
        
        system = {}
        for comp in self.config.components:
            system[comp.name] = self._build_component(comp)
        
        return {'defaults': defaults, 'system': system}
    
    def _build_component(self, comp: CGComponent) -> Dict:
        """构建单个组件配置"""
        result = {
            'name': comp.name,
            'nmol': comp.nmol,
            'charge_termini': comp.charge_termini,
        }
        
        if comp.type == ComponentType.IDP and comp.ffasta:
            result['ffasta'] = comp.ffasta
        
        elif comp.type == ComponentType.MDP:
            if comp.fpdb:
                result['fpdb'] = comp.fpdb
            if comp.fdomains:
                result['fdomains'] = comp.fdomains
            if comp.restraint:
                result.update({
                    'restraint': True,
                    'restraint_type': comp.restraint_type,
                    'use_com': comp.use_com,
                    'k_harmonic': comp.k_harmonic,
                    'colabfold': comp.colabfold,
                })
        
        return {k: v for k, v in result.items() if v is not None}
```

---

## 五、测试建议

### 5.1 单元测试

```python
import pytest
from CondenSimAdapter.src.cg import CGSimulationConfig, CGComponent, ComponentType
from CondenSimAdapter.src.calvados_wrapper import CalvadosWrapper

class TestCalvadosWrapper:
    """CalvadosWrapper 单元测试"""
    
    @pytest.fixture
    def sample_config(self):
        """示例配置"""
        from CondenSimAdapter.src.cg import CGSimulationConfig, CGComponent, ComponentType
        
        config = CGSimulationConfig(
            system_name="test_simulation",
            box=[25.0, 25.0, 30.0],
            topol=TopologyType.SLAB,
        )
        
        config.add_component(CGComponent(
            name="test_protein",
            type=ComponentType.IDP,
            nmol=10,
            ffasta="test.fasta",
        ))
        
        return config
    
    def test_topol_to_calvados_cubic(self, sample_config):
        """测试 CUBIC 拓扑转换"""
        sample_config.topol = TopologyType.CUBIC
        wrapper = CalvadosWrapper(sample_config)
        assert wrapper._topol_to_calvados() == 'grid'
    
    def test_topol_to_calvados_slab(self, sample_config):
        """测试 SLAB 拓扑转换"""
        sample_config.topol = TopologyType.SLAB
        wrapper = CalvadosWrapper(sample_config)
        assert wrapper._topol_to_calvados() == 'slab'
    
    def test_generate_config_yaml(self, sample_config):
        """测试配置 YAML 生成"""
        wrapper = CalvadosWrapper(sample_config)
        yaml_content = wrapper._generate_config_yaml()
        
        import yaml
        config = yaml.safe_load(yaml_content)
        
        assert config['sysname'] == 'test_simulation'
        assert config['topol'] == 'slab'
        assert config['slab_width'] == 15.0  # box[2] / 2
    
    def test_process_inline_fdomains(self, sample_config):
        """测试内联 fdomains 处理"""
        wrapper = CalvadosWrapper(sample_config)
        
        # 创建临时目录
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            yaml_content = """
defaults:
  molecule_type: protein
system:
  test:
    name: test
    fdomains: |
      test:
        - [1, 10]
        - [20, 30]
"""
            result = wrapper._process_inline_fdomains(yaml_content, tmpdir)
            
            import yaml
            data = yaml.safe_load(result)
            assert 'fdomains' in data['system']['test']
            assert data['system']['test']['fdomains'].endswith('_domains.yaml')
```

---

**报告完成**

*审查人: AI Assistant*  
*审查日期: 2024年12月30日*
