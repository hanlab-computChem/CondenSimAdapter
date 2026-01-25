#!/usr/bin/env python3
"""
Backmap Module

Handles conversion from coarse-grained (CG) structures to all-atom (AA) representations.
Supports two input modes:
1. adapter cg output: Reads final.pdb from CG simulation output
2. user provided: Standardizes PDB using calvados and then backmaps
"""

import os
import tempfile
from pathlib import Path
from typing import Optional, List, Dict
from enum import Enum
from dataclasses import dataclass

from openmm import Vec3
import openmm as mm
import openmm.unit as unit
from openmm.app import PDBFile

from .cg import CGSimulationConfig, ComponentType, BackmapConfig
from .calvados_wrapper import CalvadosWrapper


class SourceType(Enum):
    """Input source type."""
    MS2_CG = "ms2_cg"
    USER_PROVIDED = "user_provided"


@dataclass
class PreparedInput:
    """Prepared input data."""
    pdb_path: str
    model_type: str  # "ResidueBasedModel" or "CalphaBasedModel"
    config: Optional[CGSimulationConfig] = None  # Original config (if available)


@dataclass
class BackmapResult:
    """Backmap result."""
    success: bool
    output_pdb: str
    input_pdb: str
    model_type: str
    errors: List[str] = None
    
    def __post_init__(self):
        if self.errors is None:
            self.errors = []


def standardize_pdb_with_calvados(pdb_path: str, config: CGSimulationConfig, output_pdb: str) -> str:
    """
    Build a calvados system and standardize PDB format (CA + 3-letter codes).

    Extracted from cg.py _convert_mpipi_pdb_to_calvados_format and can be reused
    to standardize any PDB file.
    
    Args:
        pdb_path: Input PDB path
        config: CGSimulationConfig (contains component info)
        output_pdb: Output PDB path
        
    Returns:
        Output PDB path
    """
    temp_dir = tempfile.mkdtemp(prefix='calvados_topology_')
    try:
        # 1. Build system via CalvadosWrapper (generates top.pdb)
        wrapper = CalvadosWrapper(config)
        
        # Write config files to temp directory
        wrapper._write_to_dir(temp_dir, gpu_id=0, verbose=False)
        
        # Create Calvados Sim and build system
        from CondenSimAdapter.extern.ms2_calvados.calvados import sim as calvados_sim
        from yaml import safe_load
        
        with open(f'{temp_dir}/config.yaml', 'r') as stream:
            calvados_config = safe_load(stream)
        with open(f'{temp_dir}/components.yaml', 'r') as stream:
            components = safe_load(stream)
        
        # Create Sim object (build system only, no simulation)
        calvados_sim_obj = calvados_sim.Sim(temp_dir, calvados_config, components)
        calvados_sim_obj.build_system()
        
        # 2. Read topology from calvados top.pdb
        calvados_top_pdb = os.path.join(temp_dir, 'top.pdb')
        if not os.path.exists(calvados_top_pdb):
            raise FileNotFoundError(f"Calvados top.pdb not found: {calvados_top_pdb}")
        
        calvados_pdb_file = PDBFile(calvados_top_pdb)
        calvados_topology = calvados_pdb_file.topology
        
        # 3. Read coordinates from input PDB
        input_pdb_file = PDBFile(pdb_path)
        input_positions = input_pdb_file.positions
        
        # Validate atom count match
        n_calvados = calvados_topology.getNumAtoms()
        n_input = len(input_positions)
        
        if n_calvados != n_input:
            raise ValueError(
                f"Atom count mismatch: calvados={n_calvados}, input={n_input}. "
                f"Check component definitions for consistency."
            )
        
        # 4. Convert units: Angstrom -> nm (PDB uses Angstrom, OpenMM uses nm)
        positions_nm_quantity = input_positions.in_units_of(unit.nanometer)
        
        # 5. Get box vectors (from input topology or config)
        box_vectors = None
        try:
            box_vectors = input_pdb_file.topology.getPeriodicBoxVectors()
            if box_vectors is None:
                raise ValueError("No box vectors in input topology")
        except:
            # Fall back to config box
            pass
        
        # Normalize box vectors to a single Quantity to avoid OpenMM write errors
        def _normalize_box_vectors(vectors):
            if vectors is None:
                return None
            if unit.is_quantity(vectors):
                return vectors
            if isinstance(vectors, (list, tuple)) and all(unit.is_quantity(v) for v in vectors):
                vecs_nm = [v.value_in_unit(unit.nanometer) for v in vectors]
                return unit.Quantity(vecs_nm, unit.nanometer)
            return vectors

        # If none from PDB, build from config box
        if box_vectors is None:
            box = config.box
            box_vectors = unit.Quantity(
                [
                    mm.Vec3(box[0], 0, 0),
                    mm.Vec3(0, box[1], 0),
                    mm.Vec3(0, 0, box[2]),
                ],
                unit.nanometer,
            )
        else:
            box_vectors = _normalize_box_vectors(box_vectors)
        
        # Set box vectors in topology
        calvados_topology.setPeriodicBoxVectors(box_vectors)
        
        # 6. Save PDB using calvados topology + input coordinates
        with open(output_pdb, 'w') as f:
            PDBFile.writeFile(
                calvados_topology,
                positions_nm_quantity,
                f,
                keepIds=False  # Do not keep original IDs; let OpenMM renumber
            )
        
        return output_pdb
        
    finally:
        # Clean up temp directory
        import shutil
        shutil.rmtree(temp_dir, ignore_errors=True)


class BackmapSimulator:
    """Backmap simulator: CG to AA conversion."""
    
    def __init__(self, config: Optional[CGSimulationConfig] = None, backmap_config: Optional[BackmapConfig] = None):
        """
        Initialize backmap configuration.
        
        Args:
            config: CGSimulationConfig (optional, for user-provided mode)
            backmap_config: BackmapConfig (optional, from config.backmap or CLI)
        """
        self.config = config
        self.backmap_config = backmap_config or BackmapConfig()
        
        # Merge backmap config from config (CLI overrides)
        if config and config.backmap:
            if not backmap_config or not backmap_config.model_type:
                if config.backmap.model_type:
                    self.backmap_config.model_type = config.backmap.model_type
            if not backmap_config or not backmap_config.output_dir:
                if config.backmap.output_dir:
                    self.backmap_config.output_dir = config.backmap.output_dir
        # Backmap always uses CPU; ignore device/gpu settings
        self.backmap_config.device = "cpu"
    
    def center_slab_in_z(self, pdb_path: str, config: CGSimulationConfig):
        """
        SLAB topology: center condensate along z.
        
        Args:
            pdb_path: PDB file path
            config: CGSimulationConfig instance
        """
        try:
            import MDAnalysis as mda
            import numpy as np
        except ImportError:
            print("    ⚠ MDAnalysis not available, skipping z-centering")
            return
        
        # Use MDAnalysis to read PDB (handles large systems well)
        u = mda.Universe(pdb_path)
        
        # Compute protein COM (mass-weighted, Angstrom)
        # Select protein atoms only (exclude solvent)
        protein = u.select_atoms("protein")
        if len(protein) == 0:
            # If no protein detected, use all atoms
            protein = u.atoms
        
        # COM z coordinate
        com = protein.center_of_mass()  # 返回 [x, y, z]
        z_com = com[2]  # Angstrom
        
        # Box size (MDAnalysis uses Angstrom)
        if u.dimensions is not None:
            box_z = u.dimensions[2]  # Angstrom
        else:
            # From config (nm -> Angstrom)
            box_z = config.box[2] * 10.0
        
        # Compute box center and offset
        box_center_z = box_z / 2.0
        offset_z = box_center_z - z_com
        
        print(f"    Protein COM z: {z_com/10.0:.2f} nm")
        print(f"    Box z center: {box_center_z/10.0:.2f} nm")
        print(f"    Offset (z only): {offset_z/10.0:+.2f} nm")
        
        # Apply translation: shift only along z; x/y unchanged
        # Note: modify positions and assign back
        positions = u.atoms.positions.copy()
        positions[:, 2] += offset_z
        u.atoms.positions = positions
        
        # Save to original file
        u.atoms.write(pdb_path)
        self._fix_c_terminus_atom_names(pdb_path)
        self._insert_ter_after_oxt(pdb_path)
        
        print(f"    ✓ Condensate centered in z direction")

    def _fix_c_terminus_atom_names(self, pdb_path: str):
        """
        Fix C-terminal atom naming from MDAnalysis (OT1/OT2 -> O/OXT).
        """
        rename_map = {
            "OT1": " O  ",
            "OT2": " OXT",
        }
        try:
            with open(pdb_path, "r", encoding="utf-8") as handle:
                lines = handle.readlines()
        except OSError:
            return

        updated = False
        fixed_lines = []
        for line in lines:
            if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 16:
                atom_name = line[12:16]
                key = atom_name.strip()
                if key in rename_map:
                    line = f"{line[:12]}{rename_map[key]}{line[16:]}"
                    updated = True
            fixed_lines.append(line)

        if updated:
            with open(pdb_path, "w", encoding="utf-8") as handle:
                handle.writelines(fixed_lines)

    def _insert_ter_after_oxt(self, pdb_path: str):
        """
        Insert TER after OXT for each chain (OXT is treated as chain end).
        """
        try:
            with open(pdb_path, "r", encoding="utf-8") as handle:
                lines = handle.readlines()
        except OSError:
            return

        updated = False
        out_lines = []
        total = len(lines)
        for idx, line in enumerate(lines):
            out_lines.append(line)
            if not line.startswith("ATOM  "):
                continue
            if len(line) < 16:
                continue
            atom_name = line[12:16].strip()
            if atom_name != "OXT":
                continue
            next_line = lines[idx + 1] if idx + 1 < total else ""
            if next_line.startswith("TER"):
                continue
            ter_line = self._format_ter_from_atom_line(line)
            out_lines.append(ter_line)
            updated = True

        if updated:
            with open(pdb_path, "w", encoding="utf-8") as handle:
                handle.writelines(out_lines)

    def _format_ter_from_atom_line(self, atom_line: str) -> str:
        """
        基于 ATOM 行生成 TER 行，保持链/残基信息。
        """
        try:
            serial = int(atom_line[6:11]) + 1
        except ValueError:
            serial = None
        resname = atom_line[17:20]
        chain = atom_line[21:22]
        resseq = atom_line[22:26]
        icode = atom_line[26:27]
        if serial is None:
            return "TER\n"
        return f"TER   {serial:>5}      {resname:>3} {chain}{resseq}{icode}\n"
    
    def detect_source_type(self, input_path: str) -> SourceType:
        """
        检测输入是 ms2 cg 输出还是 user provided
        
        Args:
            input_path: 输入路径（目录或文件）
            
        Returns:
            SourceType.MS2_CG 或 SourceType.USER_PROVIDED
        """
        path = Path(input_path)
        
        # 检查是否是目录（ms2 cg 输出）
        if path.is_dir():
            # 检查目录结构特征
            has_final_pdb = (path / "final.pdb").exists()
            has_simulation_log = (path / "simulation.log").exists()
            
            if has_final_pdb or has_simulation_log:
                return SourceType.MS2_CG
        
        # 检查是否是 PDB 文件（user provided）
        if path.is_file() and path.suffix.lower() == '.pdb':
            return SourceType.USER_PROVIDED
        
        raise ValueError(f"Cannot determine source type for: {input_path}")
    
    def find_config_yaml(self, cg_output_dir: str, system_name: str, explicit_config: Optional[str] = None) -> Optional[str]:
        """
        查找 config.yaml 文件
        
        Args:
            cg_output_dir: CG 输出目录
            system_name: 系统名称
            explicit_config: 用户显式提供的 config 路径
            
        Returns:
            config.yaml 路径，如果找不到返回 None
        """
        # 优先使用用户显式提供的
        if explicit_config:
            config_path = Path(explicit_config)
            if config_path.exists():
                return str(config_path.resolve())
            else:
                raise FileNotFoundError(f"Config file not found: {explicit_config}")
        
        # 备选：查找当前工作目录
        pwd_config = Path.cwd() / f"{system_name}.yaml"
        if pwd_config.exists():
            return str(pwd_config)
        
        # 备选：查找父目录
        parent_config = Path(cg_output_dir).parent / f"{system_name}.yaml"
        if parent_config.exists():
            return str(parent_config)
        
        return None
    
    def select_model_type(self, config: CGSimulationConfig, force_field: Optional[str] = None) -> str:
        """
        选择 CG model 类型
        
        规则：
        - 如果用户指定了 model_type，使用用户指定的
        - 否则：Calvados + MDP → ResidueBasedModel，其他 → CalalphaBasedModel
        
        Args:
            config: CGSimulationConfig
            force_field: 力场名称（可选，用于检测）
            
        Returns:
            "ResidueBasedModel" 或 "CalphaBasedModel"
        """
        # 如果用户指定了 model_type，使用用户指定的
        if self.backmap_config.model_type and self.backmap_config.model_type != 'auto':
            return self.backmap_config.model_type
        
        # 自动选择逻辑
        has_mdp = any(c.type == ComponentType.MDP for c in config.components)
        
        # 检测力场（如果未提供）
        if force_field is None:
            force_field = self._detect_force_field_from_dir(config)
        
        is_calvados = force_field == 'calvados'
        
        if is_calvados and has_mdp:
            return "ResidueBasedModel"
        else:
            return "CalphaBasedModel"
    
    def _detect_force_field_from_dir(self, config: CGSimulationConfig) -> Optional[str]:
        """
        从配置或目录结构检测力场
        
        这是一个简单的启发式方法，可以通过检查输出目录结构来推断力场
        """
        # 如果 config 有相关信息，可以在这里添加检测逻辑
        # 目前返回 None，让调用者提供
        return None
    
    def prepare_ms2_cg_input(self, cg_output_dir: str, config_path: Optional[str] = None) -> PreparedInput:
        """
        准备 ms2 cg 输出用于 backmap
        
        Args:
            cg_output_dir: CG 输出目录路径
            config_path: 显式提供的 config.yaml 路径（可选）
            
        Returns:
            PreparedInput 对象
        """
        cg_output_path = Path(cg_output_dir)
        
        # 1. 读取 final.pdb
        final_pdb = cg_output_path / "final.pdb"
        if not final_pdb.exists():
            raise FileNotFoundError(f"final.pdb not found in {cg_output_dir}")
        
        # 2. 查找并加载 config.yaml
        # 从目录名推断 system_name（假设格式为 {system_name}_CG）
        system_name = cg_output_path.name.replace('_CG', '')
        
        config_yaml_path = self.find_config_yaml(str(cg_output_dir), system_name, config_path)
        
        if config_yaml_path:
            config = CGSimulationConfig.from_yaml(config_yaml_path)
            # 更新 simulator 的 config（用于后续使用）
            self.config = config
            # 合并 backmap 配置（CLI 参数优先）
            if config.backmap:
                if not self.backmap_config.model_type and config.backmap.model_type:
                    self.backmap_config.model_type = config.backmap.model_type
                if not self.backmap_config.device and config.backmap.device:
                    self.backmap_config.device = config.backmap.device
                if not self.backmap_config.output_dir and config.backmap.output_dir:
                    self.backmap_config.output_dir = config.backmap.output_dir
        else:
            # 如果没有找到 config，创建一个最小配置（仅用于 model type 选择）
            # 这种情况下，我们无法准确判断是否有 MDP，默认使用 CalalphaBasedModel
            config = None
        
        # 3. 检测力场（通过目录结构）
        force_field = self._detect_force_field_from_directory(cg_output_path)
        
        # 4. 选择 CG model
        if config:
            model_type = self.select_model_type(config, force_field)
        else:
            # 如果没有 config，默认使用 CalalphaBasedModel
            model_type = "CalphaBasedModel"
        
        return PreparedInput(
            pdb_path=str(final_pdb),
            model_type=model_type,
            config=config
        )
    
    def _detect_force_field_from_directory(self, cg_output_path: Path) -> Optional[str]:
        """从目录结构检测力场类型"""
        # 检查子目录名称
        subdirs = [d.name for d in cg_output_path.iterdir() if d.is_dir()]
        
        if 'Mpipi-Recharged' in subdirs or 'mpipi_recharged' in subdirs:
            return 'mpipi_recharged'
        elif 'HPS' in subdirs or 'hps' in subdirs:
            return 'hps_urry'
        elif 'COCOMO' in subdirs or 'cocomo' in subdirs:
            return 'cocomo'
        elif 'raw' in subdirs:
            # raw 目录通常表示 calvados
            return 'calvados'
        
        return None
    
    def prepare_user_provided_input(self, pdb_path: str, config: CGSimulationConfig) -> PreparedInput:
        """
        准备 user provided PDB 用于 backmap
        
        Args:
            pdb_path: 用户提供的 PDB 文件路径
            config: CGSimulationConfig（必须包含 components 信息）
            
        Returns:
            PreparedInput 对象
        """
        if not config or not config.components:
            raise ValueError("User provided mode requires config with components")
        
        # 1. 使用 calvados 构建标准化 PDB
        temp_standardized = tempfile.mktemp(suffix='.pdb', prefix='standardized_')
        standardized_pdb = standardize_pdb_with_calvados(pdb_path, config, temp_standardized)
        
        # 2. 选择 model type
        model_type = self.select_model_type(config, force_field='calvados')
        
        return PreparedInput(
            pdb_path=standardized_pdb,
            model_type=model_type,
            config=config
        )
    
    def run(self, input_path: str, config_path: Optional[str] = None, output_dir: Optional[str] = None) -> BackmapResult:
        """
        执行 backmap
        
        Args:
            input_path: 输入路径（CG 输出目录或 PDB 文件）
            config_path: 配置文件路径（可选）
            output_dir: 输出目录（可选，默认 {system_name}_backmap）
            
        Returns:
            BackmapResult 对象
        """
        result = BackmapResult(
            success=False,
            output_pdb="",
            input_pdb=input_path,
            model_type="",
            errors=[]
        )
        
        try:
            # 1. 检测输入类型
            source_type = self.detect_source_type(input_path)
            
            # 2. 准备输入
            if source_type == SourceType.MS2_CG:
                prepared = self.prepare_ms2_cg_input(input_path, config_path)
            else:  # USER_PROVIDED
                if not config_path:
                    raise ValueError("User provided mode requires config.yaml via -f option")
                config = CGSimulationConfig.from_yaml(config_path)
                # 更新 simulator 的 config
                self.config = config
                # 合并 backmap 配置（CLI 参数优先）
                if config.backmap:
                    if not self.backmap_config.model_type and config.backmap.model_type:
                        self.backmap_config.model_type = config.backmap.model_type
                    if not self.backmap_config.output_dir and config.backmap.output_dir:
                        self.backmap_config.output_dir = config.backmap.output_dir
                # Backmap 强制使用 CPU，忽略配置中的 device/gpu 相关设置
                self.backmap_config.device = "cpu"
                prepared = self.prepare_user_provided_input(input_path, config)
            
            result.model_type = prepared.model_type
            
            # 3. 确定输出目录（优先级：CLI 参数 > config.backmap.output_dir > 默认）
            if output_dir is None:
                if self.backmap_config.output_dir:
                    output_dir = self.backmap_config.output_dir
                elif prepared.config:
                    output_dir = f"{prepared.config.system_name}_backmap"
                else:
                    # 从输入路径推断
                    input_name = Path(input_path).stem
                    if input_name.endswith('_CG'):
                        input_name = input_name[:-3]
                    output_dir = f"{input_name}_backmap"
            
            os.makedirs(output_dir, exist_ok=True)
            
            # 4. 确定输出文件名（统一为 final.aa.pdb）
            output_pdb = os.path.join(output_dir, "final.aa.pdb")
            
            # 5. 执行 backmap
            from CondenSimAdapter.extern.ms2_cg2all import convert_cg2all
            
            device = "cpu"
            convert_cg2all(
                in_pdb_fn=prepared.pdb_path,
                out_fn=output_pdb,
                model_type=prepared.model_type,
                fix_atom=False,
                device=device,
                write_ssbond=False  # 默认不写入二硫键记录
            )
            
            # 6. SLAB topology: 在 z 方向居中 condensate
            if prepared.config and prepared.config.topol.value == 'slab':
                print(f"  SLAB topology detected: centering condensate in z direction...")
                self.center_slab_in_z(output_pdb, prepared.config)
            
            result.success = True
            result.output_pdb = output_pdb
            
        except Exception as e:
            result.errors.append(str(e))
            result.success = False
        
        return result
