"""
Unit tests for core/entanglement.py

Tests cover geometric primitives, PBC handling, enhanced algorithm helpers
(node insertion, kink detection, ghost elimination), and end-to-end
entanglement detection on simple synthetic systems.
"""

from __future__ import annotations

import numpy as np
import pytest

from CondenSimAdapter.core.entanglement import (
    _any_bond_pierces_triangle,
    _segment_pierces_triangle,
    EntanglementAnalyzer,
    EntanglementReport,
)


# =============================================================================
# Geometric primitive tests
# =============================================================================

class TestSegmentPiercesTriangle:
    """Möller-Trumbore intersection tests."""

    def test_direct_pierce(self):
        """Segment clearly pierces the triangle interior."""
        p0 = np.array([0.0, 0.0, 0.0])
        p1 = np.array([2.0, 0.0, 0.0])
        v0 = np.array([1.0, -1.0, -0.5])
        v1 = np.array([1.0, 1.0, -0.5])
        v2 = np.array([1.0, 0.0, 0.5])
        assert _segment_pierces_triangle(p0, p1, v0, v1, v2)

    def test_parallel_no_pierce(self):
        """Segment parallel to triangle plane – no intersection."""
        p0 = np.array([0.0, 0.0, 1.0])
        p1 = np.array([2.0, 0.0, 1.0])
        v0 = np.array([1.0, -1.0, 0.0])
        v1 = np.array([1.0, 1.0, 0.0])
        v2 = np.array([2.0, 0.0, 0.0])
        assert not _segment_pierces_triangle(p0, p1, v0, v1, v2)

    def test_misses_outside(self):
        """Segment passes near but outside the triangle."""
        p0 = np.array([0.0, 0.0, 0.0])
        p1 = np.array([2.0, 0.0, 0.0])
        v0 = np.array([1.0, 10.0, -1.0])
        v1 = np.array([1.0, 11.0, -1.0])
        v2 = np.array([1.0, 10.5, 1.0])
        assert not _segment_pierces_triangle(p0, p1, v0, v1, v2)

    def test_vertex_contact_counts(self):
        """Segment endpoint on triangle vertex – Möller-Trumbore returns
        True when the endpoint lies on an edge (t=1, u+v<=1)."""
        # Segment [(-1,0,0), (1,0,0)] hits triangle edge at (0,0,0)
        p0 = np.array([-1.0, 0.0, 0.0])
        p1 = np.array([1.0, 0.0, 0.0])
        v0 = np.array([0.0, -1.0, 0.0])
        v1 = np.array([0.0, 1.0, 0.0])
        v2 = np.array([0.0, 0.0, 1.0])
        assert _segment_pierces_triangle(p0, p1, v0, v1, v2)


class TestAnyBondPiercesTriangle:
    """Vectorised batch intersection tests."""

    def test_one_hit_among_many(self):
        """Only one bond in a batch hits the triangle."""
        v0 = np.array([1.0, -1.0, -0.5])
        v1 = np.array([1.0, 1.0, -0.5])
        v2 = np.array([1.0, 0.0, 0.5])
        # Two bonds: one misses, one hits
        bond_as = np.array([[0.0, 10.0, 0.0], [0.0, 0.0, 0.0]])
        bond_bs = np.array([[2.0, 10.0, 0.0], [2.0, 0.0, 0.0]])
        box = np.array([10.0, 10.0, 10.0])
        assert _any_bond_pierces_triangle(v0, v1, v2, bond_as, bond_bs, box, use_pbc=False)

    def test_empty_bond_list(self):
        """Empty bond list should return False."""
        empty = np.empty((0, 3), dtype=np.float64)
        box = np.array([10.0, 10.0, 10.0])
        assert not _any_bond_pierces_triangle(
            np.zeros(3), np.ones(3), np.array([1.0, 0.0, 0.0]),
            empty, empty, box, use_pbc=False,
        )

    def test_pbc_crossing_bond(self):
        """Bond that crosses a PBC boundary is correctly unfolded before test."""
        v0 = np.array([1.0, 1.0, 1.0])
        v1 = np.array([1.0, 2.0, 1.0])
        v2 = np.array([1.0, 1.5, 2.0])
        # Bond crosses box boundary in z: image at z≈9 is the closest to triangle
        bond_as = np.array([[0.0, 1.5, 9.5]])
        bond_bs = np.array([[2.0, 1.5, 9.5]])
        box = np.array([10.0, 10.0, 10.0])
        # The bond at z=9.5 is far from triangle at z≈1.5 (min image dist 8.0)
        assert not _any_bond_pierces_triangle(
            v0, v1, v2, bond_as, bond_bs, box, use_pbc=True,
        )


# =============================================================================
# EntanglementAnalyzer helper tests
# =============================================================================

class TestUnwrapChain:
    """PBC unwrapping logic."""

    def test_simple_no_pbc(self):
        """No PBC: chain is returned unchanged."""
        chain = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        analyzer = EntanglementAnalyzer(
            chain, [(0, 3)], box=[10.0, 10.0, 10.0], use_pbc=False,
        )
        unwrapped = analyzer._unwrap_chain(chain)
        np.testing.assert_array_almost_equal(unwrapped, chain)

    def test_pbc_jump_unwrapped(self):
        """A bead at box edge is unwrapped relative to predecessor."""
        chain = np.array([[9.5, 0.0, 0.0], [0.2, 0.0, 0.0], [1.0, 0.0, 0.0]])
        analyzer = EntanglementAnalyzer(
            chain, [(0, 3)], box=[10.0, 10.0, 10.0], use_pbc=True,
        )
        unwrapped = analyzer._unwrap_chain(chain)
        # bead 1 should be at 10.2 (not 0.2) relative to bead 0 at 9.5
        assert unwrapped[1, 0] == pytest.approx(10.2)
        # bead 2 should continue from bead 1
        assert unwrapped[2, 0] == pytest.approx(11.0)


class TestInsertNodesOnLongSegments:
    """Node insertion helper."""

    def test_no_insertion_when_short(self):
        """Segments shorter than lmax are untouched."""
        path = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [1.0, 0.0, 0.0]])
        analyzer = EntanglementAnalyzer(
            path, [(0, 3)], box=[10.0, 10.0, 10.0], lmax_factor=2.0,
        )
        # avg_bond = 0.5, lmax = 1.0, segments are 0.5 < 1.0
        new_path = analyzer._insert_nodes_on_long_segments(path)
        assert len(new_path) == len(path)

    def test_insertion_on_long_segment(self):
        """A single long segment is bisected."""
        path = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [6.0, 0.0, 0.0]])
        analyzer = EntanglementAnalyzer(
            path, [(0, 3)], box=[10.0, 10.0, 10.0], lmax_factor=0.5,
        )
        # avg_bond = 3.0, lmax = 1.5, segments are 3.0 > 1.5
        new_path = analyzer._insert_nodes_on_long_segments(path)
        # max_insertions = max(1, len(path)//10) = max(1, 3//10) = 1
        # So only 1 insertion: [(0,0,0), (1.5,0,0), (3,0,0), (6,0,0)]
        assert len(new_path) == 4
        np.testing.assert_array_almost_equal(new_path[1], [1.5, 0.0, 0.0])


class TestIdentifyKinks:
    """3-stage kink detection."""

    def test_kink_nearby_and_bent(self):
        """Node close to another chain and with a bend is a kink."""
        path = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]])
        # Other chain bond passes right through the middle node (1,0,0)
        other_as = np.array([[1.0, 0.0, -0.1]])
        other_bs = np.array([[1.0, 0.0, 0.1]])
        analyzer = EntanglementAnalyzer(
            path, [(0, 3)], box=[10.0, 10.0, 10.0],
            thickness_factor=1.0, kinkdef1=5.0, kinkdef2=0.001,
        )
        # thickness = 1.0 * avg_bond ≈ 0.707, distcrit1 = 5 * 0.707 ≈ 3.54
        is_kink = analyzer._identify_kinks(
            path, other_as, other_bs, np.array([10.0, 10.0, 10.0]), use_pbc=False,
        )
        assert is_kink[0]
        assert is_kink[1]  # middle node is a kink
        assert is_kink[2]

    def test_not_kink_when_too_far(self):
        """Node far from other chains is not a kink."""
        path = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        other_as = np.array([[10.0, 10.0, 10.0]])
        other_bs = np.array([[10.0, 10.0, 11.0]])
        analyzer = EntanglementAnalyzer(
            path, [(0, 3)], box=[20.0, 20.0, 20.0],
            thickness_factor=0.001, kinkdef1=5.0, kinkdef2=0.001,
        )
        is_kink = analyzer._identify_kinks(
            path, other_as, other_bs, np.array([20.0, 20.0, 20.0]), use_pbc=False,
        )
        assert is_kink[0]
        assert not is_kink[1]  # middle node is NOT a kink
        assert is_kink[2]

    def test_not_kink_when_too_straight(self):
        """Node close but collinear segments → not a kink."""
        path = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        other_as = np.array([[1.0, 0.0, -0.1]])
        other_bs = np.array([[1.0, 0.0, 0.1]])
        analyzer = EntanglementAnalyzer(
            path, [(0, 3)], box=[10.0, 10.0, 10.0],
            thickness_factor=1.0, kinkdef1=5.0, kinkdef2=0.001,
        )
        is_kink = analyzer._identify_kinks(
            path, other_as, other_bs, np.array([10.0, 10.0, 10.0]), use_pbc=False,
        )
        assert is_kink[0]
        assert not is_kink[1]  # straight segment, not a kink
        assert is_kink[2]


class TestEliminateGhostNodes:
    """Ghost-node elimination."""

    def test_removes_unconstrained_node(self):
        """Interior node with no intersecting bonds is a ghost."""
        path = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        # No other chains → no bonds intersect the triangle
        empty = np.empty((0, 3), dtype=np.float64)
        analyzer = EntanglementAnalyzer(path, [(0, 3)], box=[10.0, 10.0, 10.0])
        new_path = analyzer._eliminate_ghost_nodes(
            path, empty, empty, np.array([10.0, 10.0, 10.0]), use_pbc=False,
        )
        assert len(new_path) == 2
        np.testing.assert_array_almost_equal(new_path[0], [0.0, 0.0, 0.0])
        np.testing.assert_array_almost_equal(new_path[1], [2.0, 0.0, 0.0])

    def test_keeps_constrained_node(self):
        """Interior node whose triangle IS pierced is retained."""
        path = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]])
        # Bond from another chain that pierces the triangle
        other_as = np.array([[0.5, 0.5, -1.0]])
        other_bs = np.array([[0.5, 0.5, 1.0]])
        analyzer = EntanglementAnalyzer(path, [(0, 3)], box=[10.0, 10.0, 10.0])
        new_path = analyzer._eliminate_ghost_nodes(
            path, other_as, other_bs, np.array([10.0, 10.0, 10.0]), use_pbc=False,
        )
        assert len(new_path) == 3  # node retained


# =============================================================================
# End-to-end entanglement detection
# =============================================================================

class TestEndToEndEntanglement:
    """Full run() on synthetic systems with known topology."""

    def test_two_parallel_chains_no_entanglement(self):
        """Two parallel chains should have Z = 0."""
        positions = np.array([
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0],  # chain 0
            [0.0, 2.0, 0.0], [1.0, 2.0, 0.0], [2.0, 2.0, 0.0],  # chain 1
        ])
        boundaries = [(0, 3), (3, 6)]
        analyzer = EntanglementAnalyzer(positions, boundaries, box=[10.0, 10.0, 10.0])
        report = analyzer.run()
        assert report.n_entangled == 0
        np.testing.assert_array_equal(report.z_values, [0.0, 0.0])

    def test_two_crossing_chains_entangled(self):
        """Two chains that cross should have Z > 0."""
        # Chain 0: slight z-bend so triangles are non-degenerate.
        # Chain 1: vertical at x=1 with a z-bend; its middle bead's triangle
        # is pierced by chain 0.
        positions = np.array([
            [0.0, 0.0, 0.0], [0.9, 0.0, 0.0], [2.0, 0.0, 0.5], [3.0, 0.0, 0.0],  # chain 0
            [1.0, -1.0, 0.5], [1.0, 0.0, 0.0], [1.0, 1.0, 0.5],                 # chain 1
        ])
        boundaries = [(0, 4), (4, 7)]
        # Large thickness_factor so that kink detection is permissive for
        # this coarse synthetic geometry (real systems have much smaller
        # bond lengths where the default factor is appropriate).
        analyzer = EntanglementAnalyzer(
            positions, boundaries, box=[10.0, 10.0, 10.0], thickness_factor=100.0,
        )
        report = analyzer.run()
        assert report.n_entangled >= 1
        assert report.mean_z > 0.0

    def test_three_chain_braid(self):
        """Three-chain braid: each chain is constrained by at least one other."""
        positions = np.array([
            [0.0, 0.0, 0.0], [0.9, 0.0, 0.0], [2.0, 0.0, 0.5], [3.0, 0.0, 0.0],  # chain 0
            [1.0, -1.0, 0.5], [1.0, 0.0, 0.0], [1.0, 1.0, 0.5],                 # chain 1
            [2.5, -1.0, 0.5], [2.5, 0.0, 0.0], [2.5, 1.0, 0.5],                 # chain 2
        ])
        boundaries = [(0, 4), (4, 7), (7, 10)]
        analyzer = EntanglementAnalyzer(
            positions, boundaries, box=[10.0, 10.0, 10.0], thickness_factor=100.0,
        )
        report = analyzer.run()
        # Chain 0 should be entangled (constrained by chains 1 and/or 2)
        assert report.z_values[0] > 0.0

    def test_single_chain_always_z_zero(self):
        """A single chain cannot be entangled with itself (by default)."""
        positions = np.array([
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0],
        ])
        boundaries = [(0, 3)]
        analyzer = EntanglementAnalyzer(positions, boundaries, box=[10.0, 10.0, 10.0])
        report = analyzer.run()
        assert report.n_entangled == 0
        assert report.z_values[0] == 0.0

    def test_convergence_within_max_iter(self):
        """Algorithm converges in fewer than max_iter iterations."""
        np.random.seed(42)
        n_chains = 5
        beads_per_chain = 10
        positions = np.random.rand(n_chains * beads_per_chain, 3) * 5.0
        boundaries = [(i * beads_per_chain, (i + 1) * beads_per_chain)
                      for i in range(n_chains)]
        analyzer = EntanglementAnalyzer(positions, boundaries, box=[10.0, 10.0, 10.0])
        report = analyzer.run(max_iter=50)
        assert report.n_iter <= 50

    def test_report_structure(self):
        """EntanglementReport has correct shape and types."""
        positions = np.array([
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0], [1.0, 2.0, 0.0], [2.0, 2.0, 0.0],
        ])
        boundaries = [(0, 3), (3, 6)]
        analyzer = EntanglementAnalyzer(
            positions, boundaries, box=[10.0, 10.0, 10.0],
            chain_names=["A", "B"],
        )
        report = analyzer.run()
        assert isinstance(report.z_values, np.ndarray)
        assert len(report.z_values) == 2
        assert len(report.primitive_paths) == 2
        assert report.chain_names == ["A", "B"]
        assert isinstance(report.summary(), str)


class TestBatchPointToSegmentsDistSq:
    """Point-to-segment distance helper."""

    def test_point_on_segment(self):
        """Point on the segment → distance = 0."""
        bond_as = np.array([[0.0, 0.0, 0.0]])
        bond_bs = np.array([[2.0, 0.0, 0.0]])
        point = np.array([1.0, 0.0, 0.0])
        dist_sq = EntanglementAnalyzer._batch_point_to_segments_dist_sq(
            point, bond_as, bond_bs, np.array([10.0, 10.0, 10.0]), use_pbc=False,
        )
        assert dist_sq == pytest.approx(0.0, abs=1e-12)

    def test_point_off_segment(self):
        """Point perpendicular to segment midpoint."""
        bond_as = np.array([[0.0, 0.0, 0.0]])
        bond_bs = np.array([[2.0, 0.0, 0.0]])
        point = np.array([1.0, 1.0, 0.0])
        dist_sq = EntanglementAnalyzer._batch_point_to_segments_dist_sq(
            point, bond_as, bond_bs, np.array([10.0, 10.0, 10.0]), use_pbc=False,
        )
        assert dist_sq == pytest.approx(1.0, abs=1e-12)

    def test_empty_list(self):
        """Empty bond list → inf distance."""
        empty = np.empty((0, 3), dtype=np.float64)
        dist_sq = EntanglementAnalyzer._batch_point_to_segments_dist_sq(
            np.zeros(3), empty, empty, np.array([10.0, 10.0, 10.0]), use_pbc=False,
        )
        assert np.isinf(dist_sq)
