"""
CG -> all-atom backmapping.

Wraps the cg2all neural-network model with a minimal clean interface.
"""

from .backmapper import Backmapper, BackmapResult

__all__ = ["Backmapper", "BackmapResult"]
