#!/usr/bin/env python3
"""
CALVADOS Wrapper

将 CGSimulationConfig 转换为 ms2_calvados 格式并运行模拟。

拓扑映射：
- CUBIC → CALVADOS 'grid'
- SLAB → CALVADOS 'slab'

Usage:
    from CondenSimAdapter.src import CGSimulationConfig
    from CondenSimAdapter.src.calvados_wrapper import run_calvados
    
    config = CGSimulationConfig.from_yaml("config.yaml")
    result = run_calvados(config, output_dir="output/", gpu_id=0)
"""

import os
from pathlib import Path
from typing import Dict, Optional

from .cg import (
    CGSimulationConfig,
    CGComponent,
    ComponentType,
    TopologyType,
    ComputePlatform,
    SimulationResult,
)


class CalvadosWrapper:
    """
    CALVADOS 模拟包装器
    
    将 CGSimulationConfig 转换为 ms2_calvados 格式并运行。
    
    Attributes:
        config: CGSimulationConfig 实例
        output_dir: 输出目录
        ms2_config: ms2_calvados Config 对象
        ms2_components: ms2_calvados Components 对象
    """
    
    def __init__(self, config: CGSimulationConfig):
        """
        初始化 wrapper
        
        Args:
            config: CGSimulationConfig 实例
        """
        self.config = config
        self.output_dir: Optional[str] = None
        
        # 获取 residue 文件路径
        self._residues_path = self._get_residues_path()
    
    def _get_residues_path(self) -> str:
        """获取 residue 参数文件路径
        
        Residues 文件从 ms2_calvados 包的 data 目录加载：
        - residues_CALVADOS2.csv: 用于纯 IDP 系统
        - residues_CALVADOS3.csv: 用于包含 MDP 的系统
        """
        from CondenSimAdapter.extern.ms2_calvados.calvados import data as calvados_data
        
        has_mdp = any(c.type == ComponentType.MDP for c in self.config.components)
        residues_file = 'residues_CALVADOS3.csv' if has_mdp else 'residues_CALVADOS2.csv'
        
        # calvados_data 是命名空间包，使用 __path__ 而不是 __file__
        data_path = calvados_data.__path__[0]
        residues_path = Path(data_path) / residues_file
        
        if not residues_path.exists():
            raise FileNotFoundError(f"Residues file not found: {residues_path}")
        
        return str(residues_path)
    
    def _topol_to_calvados(self) -> str:
        """将 TopologyType 转换为 CALVADOS 拓扑字符串"""
        if self.config.topol == TopologyType.CUBIC:
            return 'grid'
        elif self.config.topol == TopologyType.SLAB:
            return 'slab'
        else:
            return 'grid'  # 默认
    
    def _platform_to_string(self) -> str:
        """将 ComputePlatform 转换为字符串"""
        if isinstance(self.config.simulation.platform, ComputePlatform):
            return self.config.simulation.platform.value
        return str(self.config.simulation.platform)
    
    def create_config(self) -> 'ms2_config.Config':
        """创建 ms2_calvados Config 对象
        
        Notes:
            - 只传递用户需要修改的参数，让 Config 类使用 default_config.yaml 的默认值
            - CALVADOS 的物理常量（eps_lj, cutoff_lj, friction_coeff 等）保持不变
            - slab_width: SLAB 拓扑时自动计算为 box[2] / 2
        """
        from CondenSimAdapter.extern.ms2_calvados.calvados.cfg import Config
        
        sim_params = self.config.simulation
        
        # SLAB 拓扑：自动计算 slab_width = box_z / 2
        # 其他拓扑：使用 CALVADOS 默认值（100）
        if self.config.topol.value == 'slab':
            slab_width = self.config.box[2] / 2
        else:
            slab_width = None  # 使用 CALVADOS 默认值
        
        # 只传递用户实际配置的参数
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
            'verbose': sim_params.verbose,
        }
        
        # SLAB 拓扑需要指定 slab_width
        if slab_width is not None:
            params['slab_width'] = slab_width
        
        return Config(**params)
    
    def create_components(self) -> 'ms2_config.Components':
        """创建 ms2_calvados Components 对象"""
        from CondenSimAdapter.extern.ms2_calvados.calvados.cfg import Components
        
        first_comp = self.config.components[0] if self.config.components else None
        
        defaults = {
            'molecule_type': 'protein',
            'nmol': first_comp.nmol if first_comp else 1,
            'restraint': first_comp.restraint if first_comp else False,
            'charge_termini': first_comp.charge_termini if first_comp else 'both',
            'fresidues': self._residues_path,
        }
        
        ms2_components = Components(**defaults)
        
        for comp in self.config.components:
            comp_dict = {
                'name': comp.name,
                'nmol': comp.nmol,
                'restraint': comp.restraint,
                'charge_termini': comp.charge_termini,
                'fresidues': self._residues_path,
            }
            
            if comp.type == ComponentType.IDP:
                if comp.ffasta:
                    comp_dict['ffasta'] = comp.ffasta
            
            elif comp.type == ComponentType.MDP:
                if comp.fpdb:
                    comp_dict['fpdb'] = comp.fpdb
                    comp_dict['pdb_folder'] = os.path.dirname(os.path.abspath(comp.fpdb))
                
                if comp.fdomains:
                    comp_dict['fdomains'] = comp.fdomains
                
                if comp.restraint:
                    comp_dict['restraint_type'] = comp.restraint_type
                    comp_dict['use_com'] = comp.use_com
                    comp_dict['k_harmonic'] = comp.k_harmonic
                    comp_dict['colabfold'] = comp.colabfold
            
            ms2_components.add(**comp_dict)
        
        return ms2_components
    
    def write(self, output_dir: str, overwrite: bool = False) -> Dict[str, str]:
        """
        写入配置文件
        
        Args:
            output_dir: 输出目录
            overwrite: 是否覆盖
            
        Returns:
            生成的文件路径字典
        """
        output_dir = os.path.abspath(output_dir)
        
        if os.path.exists(output_dir) and not overwrite:
            raise FileExistsError(f"Output directory exists: {output_dir}")
        
        os.makedirs(output_dir, exist_ok=True)
        self.output_dir = output_dir
        
        # 创建并写入 config
        ms2_config = self.create_config()
        ms2_config.write(output_dir, name='config.yaml')
        
        # 创建并写入 components
        ms2_components = self.create_components()
        ms2_components.write(output_dir, name='components.yaml')
        
        return {
            'config': os.path.join(output_dir, 'config.yaml'),
            'components': os.path.join(output_dir, 'components.yaml'),
            'run_script': os.path.join(output_dir, 'run.py'),
        }
    
    def _generate_config_yaml(self, gpu_id: int = 0, verbose: bool = False, continue_from: str = None) -> str:
        """生成 CALVADOS config.yaml 内容

        策略：
        1. 加载 CALVADOS 的 default_config.yaml 作为基础配置
        2. 只覆盖用户实际配置的参数
        3. 这样避免了硬编码所有物理常量（eps_lj, cutoff_lj 等）

        这种方式与原始 CALVADOS 的 Config 类保持一致的设计理念。

        Args:
            gpu_id: GPU 设备 ID（用户指定的 GPU）
            verbose: 是否输出详细日志
            continue_from: 继续模拟的坐标文件路径（PDB格式）
        """
        import yaml
        from CondenSimAdapter.extern.ms2_calvados.calvados.cfg import Config
        
        sim_params = self.config.simulation
        
        # SLAB 拓扑：自动计算 slab_width = box_z / 2
        # 其他拓扑：使用 CALVADOS 默认值
        if self.config.topol == TopologyType.SLAB:
            slab_width = self.config.box[2] / 2
        else:
            slab_width = None
        
        # 使用 Config 类加载默认配置
        config_obj = self.create_config()
        config_dict = config_obj.config.copy()
        
        # 覆盖用户配置的参数（包括 gpu_id 和 verbose）
        config_dict.update({
            'sysname': self.config.system_name,
            'box': self.config.box,
            'temp': self.config.temperature,
            'ionic': self.config.ionic,
            'pH': 7.0,
            'topol': self._topol_to_calvados(),
            'wfreq': sim_params.wfreq,
            'steps': sim_params.steps,
            'platform': self._platform_to_string(),
            'verbose': verbose,  # 控制 CALVADOS 详细输出
            'gpu_id': gpu_id,  # 用户指定的 GPU ID
        })
        
        # 如果提供了 continue_from，设置 restart='pdb' 和 frestart
        if continue_from:
            config_dict.update({
                'restart': 'pdb',
                'frestart': continue_from,
            })
        
        # SLAB 拓扑需要指定 slab_width
        if slab_width is not None:
            config_dict['slab_width'] = slab_width
        
        return yaml.dump(config_dict, default_flow_style=False, sort_keys=False)
    
    def _generate_components_yaml(self) -> str:
        """生成 CALVADOS components.yaml 内容
        
        处理 fpdb 和 pdb_folder:
        - CALVADOS 期望 pdb_folder（目录）和 name（不含扩展名的文件名）
        - 我们的 config 使用 fpdb（完整文件路径）
        
        Notes:
            添加了原版 CALVADOS default_component.yaml 中的所有默认参数：
            - periodic: false
            - cutoff_restr: 0.9
            - k_go: 15.
            - use_com: true
            - colabfold: 0
        """
        import yaml
        
        first_comp = self.config.components[0] if self.config.components else None
        
        # 计算 pdb_folder（从第一个 MDP 组件的 fpdb 提取）
        pdb_folder = None
        for comp in self.config.components:
            if comp.type.value == 'mdp' and comp.fpdb:
                pdb_folder = os.path.dirname(os.path.abspath(comp.fpdb))
                break
        
        components = {
            'defaults': {
                'molecule_type': 'protein',
                'nmol': first_comp.nmol if first_comp else 1,
                'restraint': first_comp.restraint if first_comp else False,
                'charge_termini': first_comp.charge_termini if first_comp else 'both',
                'fresidues': self._residues_path,
                'alpha': 0,
                'kb': 8033.0,
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
            # 对于 MDP，验证 fpdb 文件是否存在
            # CALVADOS 期望文件路径为 {pdb_folder}/{name}.pdb
            if comp.type == ComponentType.MDP and comp.fpdb:
                pdb_folder = os.path.dirname(os.path.abspath(comp.fpdb))
                expected_pdb = os.path.join(pdb_folder, f"{comp.name}.pdb")
                actual_pdb = os.path.abspath(comp.fpdb)
                
                # 验证文件存在
                if not os.path.exists(actual_pdb):
                    raise FileNotFoundError(
                        f"PDB file not found: {actual_pdb}\n"
                        f"  Component: {comp.name}\n"
                        f"  Expected by CALVADOS: {expected_pdb}\n"
                        f"\n"
                        f"解决方案（二选一）：\n"
                        f"  1. 重命名 PDB 文件: mv h1.pdb H1.pdb\n"
                        f"  2. 修改 component name: name: h1  (小写)"
                    )
                
                # 如果实际文件路径与 CALVADOS 期望的不匹配，给出警告
                if actual_pdb != expected_pdb:
                    print(f"\n  ⚠️  Warning: Component '{comp.name}' has fpdb='{comp.fpdb}'")
                    print(f"      CALVADOS expects file named: {comp.name}.pdb")
                    print(f"      This mismatch WILL cause errors!")
                    print(f"\n  解决方案：")
                    print(f"    mv {comp.fpdb} {os.path.join(pdb_folder, comp.name)}.pdb")
                    print(f"")
            
            comp_dict = {
                'name': comp.name,
                'molecule_type': 'protein',
                'nmol': comp.nmol,
                'ffasta': comp.ffasta,
                'fdomains': comp.fdomains,
                'fpdb': comp.fpdb,
                'restraint': comp.restraint,
                'restraint_type': comp.restraint_type,
                'use_com': comp.use_com,
                'k_harmonic': comp.k_harmonic,
                'colabfold': comp.colabfold,
                'charge_termini': comp.charge_termini,
            }
            # 移除 None 值
            comp_dict = {k: v for k, v in comp_dict.items() if v is not None}
            components['system'][comp.name] = comp_dict
        
        return yaml.dump(components, default_flow_style=False, sort_keys=False)
    
    def _write_to_dir(self, output_dir: str, gpu_id: int = 0, verbose: bool = False, continue_from: str = None) -> Dict[str, str]:
        """写入配置文件到指定目录（返回文件路径字典）

        支持两种 fdomains 格式：
        1. 文件路径：'TDP43_domains.yaml' - 直接复制到输出目录
        2. 内联 YAML：'TDP43:\n  - [3, 76]\n...' - 写入临时文件

        Args:
            output_dir: 输出目录
            gpu_id: GPU 设备 ID（用于写入 config.yaml）
            verbose: 是否输出详细日志
            continue_from: 继续模拟的坐标文件路径（PDB格式，会被复制到输出目录）
        """
        import tempfile
        import shutil

        os.makedirs(output_dir, exist_ok=True)
        self.output_dir = output_dir

        # 写入 config.yaml（传入 gpu_id、verbose 和 continue_from）
        config_file = os.path.join(output_dir, 'config.yaml')
        with open(config_file, 'w') as f:
            f.write(self._generate_config_yaml(gpu_id=gpu_id, verbose=verbose, continue_from=continue_from))
        
        # 处理 components.yaml，支持内联 fdomains
        components_yaml = self._generate_components_yaml()
        
        # 检查是否有内联 fdomains 需要处理
        components_yaml = self._process_inline_fdomains(components_yaml, output_dir)
        
        # 写入 components.yaml
        components_file = os.path.join(output_dir, 'components.yaml')
        with open(components_file, 'w') as f:
            f.write(components_yaml)
        
        return {
            'config': config_file,
            'components': components_file,
        }
    
    def _process_inline_fdomains(self, components_yaml: str, output_dir: str) -> str:
        """处理内联的 fdomains，如果是 YAML 内容则写入临时文件"""
        import yaml

        # 解析 YAML
        components = yaml.safe_load(components_yaml)

        for name, props in components.get('system', {}).items():
            fdomains = props.get('fdomains')
            if fdomains and isinstance(fdomains, str):
                # 移除 YAML 引号（单引号或双引号）
                stripped = fdomains.strip()
                if stripped.startswith('"') and stripped.endswith('"'):
                    stripped = stripped[1:-1]
                elif stripped.startswith("'") and stripped.endswith("'"):
                    stripped = stripped[1:-1]

                # 检查是否是内联 YAML（不是文件路径）
                is_inline = False
                if stripped.startswith('{') or stripped.startswith('['):
                    is_inline = True
                elif '\n' in stripped and (':' in stripped or stripped.startswith('-')):
                    # 多行内容且包含 YAML 特征
                    is_inline = True
                elif ':' in stripped and not stripped.endswith('.yaml') and not stripped.endswith('.yml'):
                    # 包含冒号但不像是文件路径
                    is_inline = True

                if is_inline:
                    try:
                        # 尝试解析为 YAML
                        domains_data = yaml.safe_load(stripped)

                        # 确保解析结果是字典
                        if isinstance(domains_data, dict):
                            # 只写入当前蛋白的域数据，使用蛋白名称作为 key
                            protein_domains = {name: domains_data.get(name, [])}
                        elif isinstance(domains_data, list):
                            # 直接是域列表 [[3, 76], ...]
                            protein_domains = {name: domains_data}
                        else:
                            continue

                        # 写入临时文件
                        domains_file = os.path.join(output_dir, f'{name}_domains.yaml')
                        with open(domains_file, 'w') as f:
                            yaml.dump(protein_domains, f, default_flow_style=False)

                        # 替换为文件路径
                        props['fdomains'] = domains_file

                    except yaml.YAMLError:
                        # 不是有效的 YAML，保持原样（可能是文件路径）
                        pass

        return yaml.dump(components, default_flow_style=False, sort_keys=False)
    
    def run(self, output_dir: str = None, gpu_id: int = 0, verbose: bool = False) -> SimulationResult:
        """
        运行 CALVADOS 模拟
        
        统一输出结构：
        {system_name}_CG/
        ├── final.pdb                   # 最终结构
        ├── trajectory.xtc              # 模拟轨迹
        ├── simulation.log              # 高层级日志
        └── raw/                        # 原生输出
            ├── config.yaml
            ├── components.yaml
            ├── *.xtc, *.xml, *.pdb, *.chk, *.txt
        
        Args:
            output_dir: 输出目录（默认使用 config 中的 output_dir，如果传入则直接使用）
            gpu_id: GPU 设备 ID
            
        Returns:
            SimulationResult
        """
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
        
        # 如果目录已存在，备份后重建
        if os.path.exists(output_dir):
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            backup_dir = f"{output_dir}_backup_{timestamp}"
            shutil.move(output_dir, backup_dir)
            print(f"  📁 备份旧结果到: {backup_dir}")
        
        os.makedirs(raw_dir, exist_ok=True)

        # 写入配置文件到 raw 目录（传入 gpu_id 和 verbose）
        files = self._write_to_dir(raw_dir, gpu_id=gpu_id, verbose=verbose)
        
        result = SimulationResult()
        result.output_dir = output_dir
        
        start_time = time.time()
        
        try:
            print(f"\n[CALVADOS] Running simulation...")
            print(f"  GPU ID: {gpu_id}")
            print(f"  Task: {task_name}")
            print(f"  Raw output: {raw_dir}")
            print(f"  Topology: {self._topol_to_calvados()}")
            
            # 运行模拟（输出到 raw 目录）
            calvados_sim.run(
                path=raw_dir,
                fconfig='config.yaml',
                fcomponents='components.yaml'
            )
            
            # 组织输出文件
            self._organize_output(raw_dir, output_dir, task_name)
            
            result.success = True
            elapsed = time.time() - start_time
            print(f"  ✓ CALVADOS simulation completed ({elapsed:.1f}s)")
            
        except Exception as e:
            result.success = False
            result.errors.append(str(e))
            print(f"  ✗ CALVADOS simulation failed: {e}")
            elapsed = time.time() - start_time
        
        # 写入高层级日志
        self._write_simulation_log(output_dir, task_name, elapsed, result.success)
        
        # 设置结果文件路径
        result.trajectory = os.path.join(output_dir, 'trajectory.xtc')
        result.structure = os.path.join(output_dir, 'final.pdb')
        
        for key in ['trajectory', 'structure']:
            path = getattr(result, key)
            if path and not os.path.exists(path):
                setattr(result, key, None)
        
        return result
    
    def _organize_output(self, raw_dir: str, output_dir: str, task_name: str):
        """
        组织输出文件到统一结构
        
        统一命名规则：
        - trajectory.xtc  <- {task_name}.xtc
        - final.pdb       <- 带时间戳的 pdb 或 checkpoint.pdb
        """
        import shutil
        
        sysname = self.config.system_name
        
        # 1. 处理轨迹文件
        src_xtc = os.path.join(raw_dir, f'{sysname}.xtc')
        dst_xtc = os.path.join(output_dir, 'trajectory.xtc')
        if os.path.exists(src_xtc):
            shutil.copy2(src_xtc, dst_xtc)
            print(f"  📦 trajectory.xtc")
        
        # 2. 查找并复制最终结构（优先使用 checkpoint.pdb，否则找时间戳 PDB）
        src_pdb = os.path.join(raw_dir, 'checkpoint.pdb')
        if not os.path.exists(src_pdb):
            # 找带时间戳的 PDB
            for f in os.listdir(raw_dir):
                if f.endswith('.pdb') and f != 'top.pdb':
                    src_pdb = os.path.join(raw_dir, f)
                    break
        
        dst_pdb = os.path.join(output_dir, 'final.pdb')
        if os.path.exists(src_pdb):
            shutil.copy2(src_pdb, dst_pdb)
            print(f"  📦 final.pdb")
        
        # 3. 复制重要文件到 raw 目录（如果不在那里）
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
        
        # 4. 重命名 log 文件
        for f in os.listdir(raw_dir):
            if f.endswith('.log') or f.endswith('.txt'):
                pass  # 保留原样
        
        print(f"  📁 原始输出已整理到: {raw_dir}")
    
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
  - trajectory.xtc: Simulation trajectory
  - raw/: Native simulation output files
"""
        with open(log_file, 'w') as f:
            f.write(log_content)
        
        print(f"  📝 simulation.log")


def run_calvados(config: CGSimulationConfig, output_dir: str = None, gpu_id: int = 0) -> SimulationResult:
    """
    运行 CALVADOS 模拟的便捷函数
    
    Args:
        config: CGSimulationConfig 实例
        output_dir: 输出目录
        gpu_id: GPU 设备 ID
        
    Returns:
        SimulationResult
    """
    wrapper = CalvadosWrapper(config)
    return wrapper.run(output_dir=output_dir, gpu_id=gpu_id)


__all__ = [
    'CalvadosWrapper',
    'run_calvados',
]

