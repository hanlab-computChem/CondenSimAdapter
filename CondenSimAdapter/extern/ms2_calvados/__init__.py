"""
ms2_calvados - Internalized CALVADOS coarse-grained simulation package
保留核心功能：体系构建 + OpenMM 模拟
"""

# 核心模块
from .calvados import cfg
from .calvados import components
from .calvados import build
from .calvados import sim
from .calvados import sequence
from .calvados import interactions

# 暴露关键类
from .calvados.cfg import Config, Job, Components
from .calvados.components import Component, Protein, RNA, Lipid, Crowder

__all__ = [
    'cfg',
    'components',
    'build',
    'sim',
    'sequence',
    'interactions',
    'Config',
    'Job',
    'Components',
    'Component',
    'Protein',
    'RNA',
    'Lipid',
    'Crowder',
]

