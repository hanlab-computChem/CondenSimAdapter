"""Tests for SASA calculation module."""

import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent / "CondenSimAdapter" / "core"))

from sasa import (
    _generate_sphere_points,
    _guess_element,
    calc_sasa_from_pdb,
    calc_sasa_shrake_rupley,
)


class TestSpherePoints:
    """Test sphere point generation."""

    def test_unit_length(self):
        """Generated points should be on unit sphere."""
        points = _generate_sphere_points(1920)
        norms = np.linalg.norm(points, axis=1)
        assert np.allclose(norms, 1.0, atol=1e-10)

    def test_shape(self):
        """Check output shape."""
        points = _generate_sphere_points(100)
        assert points.shape == (100, 3)


class TestElementGuessing:
    """Test element guessing from atom names."""

    def test_common_atoms(self):
        """Test common atom names."""
        # CA is alpha carbon in PDB, not calcium
        assert _guess_element("CA") == "C"
        assert _guess_element("N") == "N"
        assert _guess_element("O") == "O"
        assert _guess_element("CB") == "C"
        assert _guess_element("1HB") == "H"

    def test_two_letter_elements(self):
        """Test two-letter element names."""
        assert _guess_element("CL") == "CL"
        assert _guess_element("BR") == "BR"
        assert _guess_element("FE") == "FE"
        # Note: 'CA' is alpha carbon, not calcium (would be 'CA  ' with spaces in PDB)


class TestSASACalculation:
    """Test SASA calculation."""

    def test_single_atom(self):
        """Single atom should have full surface area."""
        coords = np.array([[0.0, 0.0, 0.0]])
        radii = np.array([1.7])  # Carbon
        sasa = calc_sasa_shrake_rupley(coords, radii, n_sphere_points=960)

        expected = 4 * np.pi * (1.7 + 1.4) ** 2  # 4*pi*(R+probe)^2
        assert np.isclose(sasa[0], expected, rtol=0.02)  # 2% tolerance

    def test_two_separate_atoms(self):
        """Two far apart atoms should have independent surfaces."""
        coords = np.array([[0.0, 0.0, 0.0], [100.0, 0.0, 0.0]])
        radii = np.array([1.7, 1.7])
        sasa = calc_sasa_shrake_rupley(coords, radii, n_sphere_points=960)

        # Each should have approximately full surface
        expected_single = 4 * np.pi * (1.7 + 1.4) ** 2
        assert np.allclose(sasa, expected_single, rtol=0.02)

    def test_two_overlapping_atoms(self):
        """Two close atoms should have reduced SASA."""
        coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
        radii = np.array([1.7, 1.7])
        sasa = calc_sasa_shrake_rupley(coords, radii, n_sphere_points=960)

        expected_single = 4 * np.pi * (1.7 + 1.4) ** 2
        # Each should have reduced SASA
        assert sasa[0] < expected_single
        assert sasa[1] < expected_single
        # Total should be less than 2*single
        assert np.sum(sasa) < 2 * expected_single

    def test_empty(self):
        """Empty input should return empty array."""
        coords = np.array([]).reshape(0, 3)
        radii = np.array([])
        sasa = calc_sasa_shrake_rupley(coords, radii)
        assert len(sasa) == 0


class TestSASAFromPDB:
    """Test SASA calculation from PDB files."""

    @pytest.fixture
    def simple_pdb(self):
        """Create a simple PDB file with 3 residues."""
        pdb_content = """ATOM      1  N   ALA A   1      27.210  12.870  -1.620  1.00  0.00           N
ATOM      2  CA  ALA A   1      27.980  12.660  -2.850  1.00  0.00           C
ATOM      3  C   ALA A   1      29.470  12.450  -2.550  1.00  0.00           C
ATOM      4  O   ALA A   1      29.940  12.230  -1.420  1.00  0.00           O
ATOM      5  CB  ALA A   1      27.760  13.890  -3.730  1.00  0.00           C
ATOM      6  N   GLY A   2      30.180  12.550  -3.650  1.00  0.00           N
ATOM      7  CA  GLY A   2      31.600  12.350  -3.550  1.00  0.00           C
ATOM      8  C   GLY A   2      32.310  13.650  -3.200  1.00  0.00           C
ATOM      9  O   GLY A   2      32.110  14.690  -3.820  1.00  0.00           O
ATOM     10  N   VAL A   3      33.190  13.520  -2.210  1.00  0.00           N
ATOM     11  CA  VAL A   3      33.950  14.710  -1.830  1.00  0.00           C
ATOM     12  C   VAL A   3      35.410  14.350  -1.550  1.00  0.00           C
ATOM     13  O   VAL A   3      35.730  13.210  -1.180  1.00  0.00           O
ATOM     14  CB  VAL A   3      33.780  15.820  -2.880  1.00  0.00           C
ATOM     15  CG1 VAL A   3      32.350  16.220  -3.110  1.00  0.00           C
ATOM     16  CG2 VAL A   3      34.570  17.080  -2.520  1.00  0.00           C
END
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
            f.write(pdb_content)
            return f.name

    @pytest.fixture
    def ca_only_pdb(self):
        """Create a CA-only PDB file."""
        pdb_content = """ATOM      1  CA  ALA A   1      27.980  12.660  -2.850  1.00  0.00           C
ATOM      2  CA  GLY A   2      31.600  12.350  -3.550  1.00  0.00           C
ATOM      3  CA  VAL A   3      33.950  14.710  -1.830  1.00  0.00           C
END
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
            f.write(pdb_content)
            return f.name

    def test_simple_pdb(self, simple_pdb):
        """Test SASA calculation from full-atom PDB."""
        sasa = calc_sasa_from_pdb(simple_pdb, n_sphere_points=960)

        assert sasa is not None
        assert len(sasa) == 3  # 3 residues
        # All values should be positive
        assert np.all(sasa > 0)
        # Values should be in reasonable range (nm²)
        assert np.all(sasa < 10)  # Very generous upper bound

        os.unlink(simple_pdb)

    def test_ca_only_returns_none(self, ca_only_pdb):
        """CA-only PDB should return None."""
        sasa = calc_sasa_from_pdb(ca_only_pdb, n_sphere_points=960)
        assert sasa is None
        os.unlink(ca_only_pdb)

    def test_nonexistent_file(self):
        """Non-existent file should return None."""
        sasa = calc_sasa_from_pdb("/nonexistent/path.pdb")
        assert sasa is None

    def test_units(self, simple_pdb):
        """Verify output is in nm²."""
        sasa = calc_sasa_from_pdb(simple_pdb, n_sphere_points=960)

        # Typical residue SASA is 0.5-3.0 nm²
        # (50-300 Å²)
        assert np.all(sasa > 0.1)  # nm²
        assert np.all(sasa < 5.0)  # nm²

        os.unlink(simple_pdb)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
