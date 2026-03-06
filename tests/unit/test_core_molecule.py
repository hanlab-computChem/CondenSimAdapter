"""
Unit tests for core/molecule.py

Tests cover:
- IDP chain building (spiral/linear/compact)
- MDP chain extraction from PDB
- Multi-chain placement (grid/slab/random)
- build_all_chains complete workflow
- Edge cases: empty components, clash handling
"""

from __future__ import annotations

import numpy as np
import pytest
from pathlib import Path
from typing import List, Tuple

from CondenSimAdapter.core.molecule import (
    build_idp_chain,
    build_mdp_chain,
    place_chains_grid,
    place_chains_slab,
    place_chains_random,
    build_all_chains,
    _build_spiral,
    _build_linear,
    _build_compact,
    _build_xyzgrid,
    _assemble_chains,
)
from CondenSimAdapter.core.config import (
    Component,
    ComponentType,
    CGConfig,
    TopologyType,
)


# =============================================================================
# IDP Chain Building Tests
# =============================================================================

class TestBuildIDPChain:
    """Tests for build_idp_chain function."""

    def test_spiral_method_basic(self):
        """Test spiral method creates correct shape."""
        coords = build_idp_chain(n_beads=10, method="spiral", spacing=0.38)
        assert coords.shape == (10, 3)
        assert coords.dtype == np.float64
        # Should be centered at origin
        np.testing.assert_allclose(coords.mean(axis=0), [0, 0, 0], atol=1e-10)

    def test_spiral_bond_length(self):
        """Test spiral method has reasonable bond lengths."""
        coords = build_idp_chain(n_beads=20, method="spiral", spacing=0.38)
        # Calculate bond lengths
        bonds = np.linalg.norm(coords[1:] - coords[:-1], axis=1)
        # Spiral method doesn't have fixed bond length, but should be reasonable
        assert np.all(bonds > 0.1) and np.all(bonds < 1.0)

    def test_linear_method(self):
        """Test linear method creates straight chain along z."""
        coords = build_idp_chain(n_beads=10, method="linear", spacing=0.38)
        assert coords.shape == (10, 3)
        # Linear chain along z means x and y should be near zero
        np.testing.assert_allclose(coords[:, 0], 0, atol=1e-10)
        np.testing.assert_allclose(coords[:, 1], 0, atol=1e-10)
        # z should be evenly spaced
        z_diff = np.diff(coords[:, 2])
        np.testing.assert_allclose(z_diff, 0.38, rtol=0.1)

    def test_compact_method(self):
        """Test compact method creates cubic lattice."""
        coords = build_idp_chain(n_beads=27, method="compact", spacing=0.38)
        assert coords.shape == (27, 3)
        # Should be centered at origin
        np.testing.assert_allclose(coords.mean(axis=0), [0, 0, 0], atol=1e-10)

    def test_invalid_method_raises(self):
        """Test that invalid method raises ValueError."""
        with pytest.raises(ValueError, match="Unknown IDP build method"):
            build_idp_chain(n_beads=10, method="invalid")

    def test_single_bead(self):
        """Test chain with single bead."""
        coords = build_idp_chain(n_beads=1, method="spiral")
        assert coords.shape == (1, 3)
        np.testing.assert_allclose(coords[0], [0, 0, 0], atol=1e-10)


class TestBuildSpiralInternal:
    """Tests for _build_spiral helper."""

    def test_spiral_centered(self):
        """Test spiral is centered at origin."""
        coords = _build_spiral(n=50, d=0.38)
        com = coords.mean(axis=0)
        np.testing.assert_allclose(com, [0, 0, 0], atol=1e-10)

    def test_spiral_xy_increases(self):
        """Test that xy radius increases with index."""
        coords = _build_spiral(n=50, d=0.38)
        radii = np.linalg.norm(coords[:, :2], axis=1)
        # Later beads should generally be further from origin in xy
        assert radii[-1] > radii[0]

    def test_spiral_z_range(self):
        """Test z span is reasonable."""
        n = 50
        d = 0.38
        coords = _build_spiral(n=n, d=d)
        z_span = coords[:, 2].max() - coords[:, 2].min()
        expected_span = n * d
        assert z_span < expected_span * 1.5  # Allow some tolerance


class TestBuildLinearInternal:
    """Tests for _build_linear helper."""

    def test_linear_z_spacing(self):
        """Test linear chain has correct z spacing."""
        n = 10
        d = 0.38
        coords = _build_linear(n=n, d=d)
        z_coords = coords[:, 2]
        # Should be centered
        assert z_coords.mean() == pytest.approx(0, abs=1e-10)
        # Spacing should be d
        dz = np.diff(z_coords)
        np.testing.assert_allclose(dz, d, rtol=1e-10)


class TestBuildCompactInternal:
    """Tests for _build_compact helper."""

    def test_compact_lattice(self):
        """Test compact creates cubic lattice points."""
        coords = _build_compact(n=27, d=1.0)  # 3x3x3 = 27
        # All coordinates should be multiples of d
        unique_x = np.unique(coords[:, 0])
        unique_y = np.unique(coords[:, 1])
        unique_z = np.unique(coords[:, 2])
        assert len(unique_x) == 3
        assert len(unique_y) == 3
        assert len(unique_z) == 3

    def test_compact_centered(self):
        """Test compact is centered at origin."""
        coords = _build_compact(n=27, d=1.0)
        com = coords.mean(axis=0)
        np.testing.assert_allclose(com, [0, 0, 0], atol=1e-10)


# =============================================================================
# MDP Chain Building Tests
# =============================================================================

class TestBuildMDPChain:
    """Tests for build_mdp_chain function (requires MDAnalysis)."""

    @pytest.fixture
    def sample_pdb(self, tmp_path: Path) -> str:
        """Create a sample PDB file with 5 residues."""
        pdb_content = """\
ATOM      1  N   ALA A   1      11.104  14.773  14.040  1.00  0.00           N
ATOM      2  CA  ALA A   1      11.381  13.323  14.114  1.00  0.00           C
ATOM      3  C   ALA A   1      10.254  12.440  14.678  1.00  0.00           C
ATOM      4  O   ALA A   1       9.171  12.908  14.966  1.00  0.00           O
ATOM      5  CB  ALA A   1      12.600  13.090  14.983  1.00  0.00           C
ATOM      6  N   GLY A   2      10.551  11.156  14.781  1.00  0.00           N
ATOM      7  CA  GLY A   2      11.705  10.576  14.121  1.00  0.00           C
ATOM      8  C   GLY A   2      12.978  11.292  14.511  1.00  0.00           C
ATOM      9  O   GLY A   2      13.138  12.515  14.436  1.00  0.00           O
ATOM     10  N   SER A   3      13.878  10.561  14.982  1.00  0.00           N
ATOM     11  CA  SER A   3      15.150  11.089  15.459  1.00  0.00           C
ATOM     12  C   SER A   3      15.126  12.486  16.025  1.00  0.00           C
ATOM     13  O   SER A   3      14.103  13.133  16.020  1.00  0.00           O
ATOM     14  CB  SER A   3      16.276  10.938  14.424  1.00  0.00           C
ATOM     15  OG  SER A   3      16.384   9.614  13.938  1.00  0.00           O
ATOM     16  N   VAL A   4      16.220  12.967  16.497  1.00  0.00           N
ATOM     17  CA  VAL A   4      16.341  14.292  17.087  1.00  0.00           C
ATOM     18  C   VAL A   4      15.088  14.785  17.797  1.00  0.00           C
ATOM     19  O   VAL A   4      14.695  15.922  17.710  1.00  0.00           O
ATOM     20  CB  VAL A   4      17.519  14.312  18.082  1.00  0.00           C
ATOM     21  CG1 VAL A   4      17.796  15.705  18.629  1.00  0.00           C
ATOM     22  CG2 VAL A   4      18.777  13.693  17.539  1.00  0.00           C
ATOM     23  N   LEU A   5      14.476  13.952  18.614  1.00  0.00           N
ATOM     24  CA  LEU A   5      13.288  14.323  19.385  1.00  0.00           C
ATOM     25  C   LEU A   5      13.658  15.170  20.595  1.00  0.00           C
ATOM     26  O   LEU A   5      14.784  15.599  20.716  1.00  0.00           O
ATOM     27  CB  LEU A   5      12.323  15.088  18.465  1.00  0.00           C
ATOM     28  CG  LEU A   5      11.067  14.383  17.920  1.00  0.00           C
ATOM     29  CD1 LEU A   5      11.397  13.062  17.259  1.00  0.00           C
ATOM     30  CD2 LEU A   5       9.993  14.250  18.953  1.00  0.00           C
END
"""
        pdb_path = tmp_path / "sample.pdb"
        pdb_path.write_text(pdb_content)
        return str(pdb_path)

    def test_mdp_chain_ca_coords(self, sample_pdb: str):
        """Test extracting CA coordinates from PDB."""
        coords, seq = build_mdp_chain(sample_pdb, use_com=False)
        assert coords.shape == (5, 3)  # 5 residues
        assert len(seq) == 5
        assert seq == "AGSVL"  # ALA-GLY-SER-VAL-LEU
        # Should be centered
        np.testing.assert_allclose(coords.mean(axis=0), [0, 0, 0], atol=1e-10)

    def test_mdp_chain_com_coords(self, sample_pdb: str):
        """Test extracting COM coordinates from PDB."""
        coords, seq = build_mdp_chain(sample_pdb, use_com=True)
        assert coords.shape == (5, 3)
        assert seq == "AGSVL"
        np.testing.assert_allclose(coords.mean(axis=0), [0, 0, 0], atol=1e-10)

    def test_mdp_chain_units_in_nm(self, sample_pdb: str):
        """Test that coordinates are in nanometers (not Angstroms)."""
        coords, _ = build_mdp_chain(sample_pdb, use_com=False)
        # PDB coordinates are in Angstroms (10s of nm)
        # Typical CA-CA distance is ~3.8 Angstroms = 0.38 nm
        bonds = np.linalg.norm(coords[1:] - coords[:-1], axis=1)
        # Should be around 0.38 nm, not 3.8
        assert np.all(bonds < 1.0), f"Bonds seem to be in Angstroms: {bonds}"


# =============================================================================
# Grid Placement Tests
# =============================================================================

class TestBuildXYZGrid:
    """Tests for _build_xyzgrid helper."""

    def test_grid_count_matches(self):
        """Test that grid produces exactly N points."""
        for n in [1, 5, 10, 27, 64]:
            grid = _build_xyzgrid(n, [10.0, 10.0, 10.0])
            assert len(grid) == n, f"Expected {n} points, got {len(grid)}"

    def test_grid_within_box(self):
        """Test that all grid points are within box bounds."""
        box = [20.0, 15.0, 10.0]
        grid = _build_xyzgrid(50, box)
        assert np.all(grid >= 0)
        assert np.all(grid[:, 0] <= box[0])
        assert np.all(grid[:, 1] <= box[1])
        assert np.all(grid[:, 2] <= box[2])

    def test_grid_staggered_pattern(self):
        """Test that grid has staggered pattern (alternating offsets)."""
        grid = _build_xyzgrid(27, [3.0, 3.0, 3.0])
        # In staggered grid, adjacent z-planes should have different xy positions
        # Check that not all xy positions are the same
        xy = grid[:, :2]
        unique_xy = np.unique(xy, axis=0)
        # Should have more than just regular grid positions
        assert len(unique_xy) > 3


class TestPlaceChainsGrid:
    """Tests for place_chains_grid function."""

    def test_single_chain(self):
        """Test placing single chain."""
        chain = np.array([[0, 0, 0], [0, 0, 0.38], [0, 0, 0.76]])
        positions = place_chains_grid([chain], [10.0, 10.0, 10.0])
        assert positions.shape == (3, 3)
        # Single chain should be translated to first grid point
        assert positions[0, 0] >= 0

    def test_multiple_chains(self):
        """Test placing multiple chains."""
        chains = [
            np.array([[0, 0, 0], [0, 0, 0.38]]),
            np.array([[0, 0, 0], [0, 0, 0.38]]),
            np.array([[0, 0, 0], [0, 0, 0.38]]),
        ]
        positions = place_chains_grid(chains, [10.0, 10.0, 10.0])
        assert positions.shape == (6, 3)  # 3 chains * 2 beads each
        # Chains should be at different positions
        chain1_com = positions[:2].mean(axis=0)
        chain2_com = positions[2:4].mean(axis=0)
        chain3_com = positions[4:6].mean(axis=0)
        # At least one coordinate should differ significantly
        assert not np.allclose(chain1_com, chain2_com, atol=0.1)
        assert not np.allclose(chain2_com, chain3_com, atol=0.1)


class TestPlaceChainsSlab:
    """Tests for place_chains_slab function."""

    def test_slab_z_centered(self):
        """Test that slab is centered in z."""
        chains = [np.array([[0, 0, 0], [0, 0, 0.38]]) for _ in range(10)]
        box = [20.0, 20.0, 100.0]
        positions = place_chains_slab(chains, box, slab_width=20.0)
        # Check z range is centered around box[2]/2 = 50
        z_mean = positions[:, 2].mean()
        assert 45 < z_mean < 55

    def test_slab_width_default(self):
        """Test slab with default width."""
        chains = [np.array([[0, 0, 0]]) for _ in range(5)]
        box = [20.0, 20.0, 50.0]
        positions = place_chains_slab(chains, box, slab_width=None)
        # Default width is min(box[0], box[1]) / 2 = 10
        # Just verify it runs without error
        assert positions.shape[0] == 5


class TestPlaceChainsRandom:
    """Tests for place_chains_random function."""

    def test_random_within_box(self):
        """Test that random placement keeps chains within box."""
        np.random.seed(42)  # For reproducibility
        chains = [np.array([[0, 0, 0], [0.5, 0, 0]]) for _ in range(5)]
        box = [10.0, 10.0, 10.0]
        positions = place_chains_random(chains, box, clash_cutoff=0.1)
        # All positions should be within box
        assert np.all(positions >= 0)
        assert np.all(positions <= np.array(box))

    def test_random_no_clash(self):
        """Test that random placement avoids clashes."""
        np.random.seed(42)
        chains = [np.array([[0, 0, 0]]) for _ in range(10)]
        box = [20.0, 20.0, 20.0]
        positions = place_chains_random(chains, box, clash_cutoff=1.0)
        # Check pairwise distances
        for i in range(len(positions)):
            for j in range(i + 1, len(positions)):
                dist = np.linalg.norm(positions[i] - positions[j])
                assert dist >= 1.0, f"Clash detected: distance {dist}"

    def test_random_impossible_raises(self):
        """Test that impossible placement raises ValueError."""
        np.random.seed(42)
        # Try to place many large chains in small box
        chains = [np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]]) for _ in range(100)]
        box = [2.0, 2.0, 2.0]  # Too small
        with pytest.raises(ValueError, match="Giving up"):
            place_chains_random(chains, box, clash_cutoff=1.0, max_tries=100)


# =============================================================================
# High-level Builder Tests
# =============================================================================

class TestBuildAllChains:
    """Tests for build_all_chains function."""

    def test_single_idp_cubic(self):
        """Test building single IDP in cubic box."""
        comp = Component(
            name="FUS",
            comp_type=ComponentType.IDP,
            sequence="GSMASAS",
            nmol=1,
        )
        config = CGConfig(
            system_name="test",
            force_field="calvados2",
            components=[comp],
            box=[20.0, 20.0, 20.0],
            topology=TopologyType.CUBIC,
        )
        positions, meta = build_all_chains(config)
        assert positions.shape == (7, 3)  # 7 residues
        assert len(meta) == 1
        assert meta[0]["name"] == "FUS"
        assert meta[0]["sequence"] == "GSMASAS"
        assert meta[0]["start"] == 0
        assert meta[0]["end"] == 7

    def test_multiple_idps_slab(self):
        """Test building multiple IDPs in slab."""
        comp = Component(
            name="ProteinA",
            comp_type=ComponentType.IDP,
            sequence="AAAA",
            nmol=3,
        )
        config = CGConfig(
            system_name="test",
            force_field="hps",
            components=[comp],
            box=[30.0, 30.0, 100.0],
            topology=TopologyType.SLAB,
            slab_width=20.0,
        )
        positions, meta = build_all_chains(config)
        assert positions.shape == (12, 3)  # 3 chains * 4 beads
        assert len(meta) == 3
        # Check chain metadata
        assert meta[0]["start"] == 0
        assert meta[0]["end"] == 4
        assert meta[1]["start"] == 4
        assert meta[1]["end"] == 8
        assert meta[2]["start"] == 8
        assert meta[2]["end"] == 12

    def test_mdp_component(self, tmp_path: Path):
        """Test building with MDP component."""
        # Create a simple PDB
        pdb_content = """\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  CA  ALA A   2       0.000   0.000   3.800  1.00  0.00           C
ATOM      3  CA  ALA A   3       0.000   0.000   7.600  1.00  0.00           C
END
"""
        pdb_path = tmp_path / "mdp.pdb"
        pdb_path.write_text(pdb_content)
        
        comp = Component(
            name="MDP1",
            comp_type=ComponentType.MDP,
            pdb_path=str(pdb_path),
            folded_domains=[(1, 3)],
            nmol=1,
        )
        config = CGConfig(
            system_name="test",
            force_field="calvados3",  # Use calvados3 for MDP
            components=[comp],
            box=[20.0, 20.0, 20.0],
            topology=TopologyType.CUBIC,
        )
        positions, meta = build_all_chains(config)
        assert positions.shape == (3, 3)
        assert len(meta) == 1
        assert meta[0]["comp_type"] == ComponentType.MDP
        assert meta[0]["folded_domains"] == [(1, 3)]

    def test_droplet_topology(self):
        """Test building in droplet topology."""
        np.random.seed(42)  # For reproducible placement
        comp = Component(
            name="IDP1",
            comp_type=ComponentType.IDP,
            sequence="AAA",  # Shorter chains for faster/more stable test
            nmol=3,
        )
        config = CGConfig(
            system_name="test",
            force_field="hps",
            components=[comp],
            box=[20.0, 20.0, 20.0],
            topology=TopologyType.DROPLET,
            droplet_radius=5.0,
        )
        positions, meta = build_all_chains(config)
        assert positions.shape == (9, 3)  # 3 chains * 3 beads
        assert len(meta) == 3
        # All positions should be within the droplet region (centered in box)
        # With droplet_radius=5 and box=[20,20,20], droplet box is [10,10,10]
        # centered at [10,10,10], so positions should be in [5,15] range roughly
        com = positions.mean(axis=0)
        box_center = np.array([10.0, 10.0, 10.0])
        # Check that COM is reasonably close to center (allowing for random spread)
        distance_from_center = np.linalg.norm(com - box_center)
        assert distance_from_center < 10.0, f"COM too far from center: {com}"

    def test_multi_component_system(self):
        """Test building system with multiple different components."""
        comp_a = Component(name="A", comp_type=ComponentType.IDP, sequence="AAA", nmol=2)
        comp_b = Component(name="B", comp_type=ComponentType.IDP, sequence="BBBB", nmol=3)
        config = CGConfig(
            system_name="multi",
            force_field="cocomo",
            components=[comp_a, comp_b],
            box=[25.0, 25.0, 25.0],
            topology=TopologyType.CUBIC,
        )
        positions, meta = build_all_chains(config)
        # 2*3 + 3*4 = 18 beads
        assert positions.shape == (18, 3)
        assert len(meta) == 5  # 2 + 3 chains
        # Check first two are component A
        assert meta[0]["name"] == "A"
        assert meta[1]["name"] == "A"
        # Check next three are component B
        assert meta[2]["name"] == "B"
        assert meta[3]["name"] == "B"
        assert meta[4]["name"] == "B"
