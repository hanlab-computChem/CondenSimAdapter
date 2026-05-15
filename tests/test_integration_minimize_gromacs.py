import shutil
from pathlib import Path

import pytest
from CondenSimAdapter.src.pdb2gmx_utils import (
    run_pdb2gmx_for_structure,
    run_pdb2gmx_for_topology,
)
from CondenSimAdapter.src.top_to_softcore_system import (
    NONBONDED_STANDARD,
    GromacsTopFileWithSoftcore,
)


@pytest.mark.slow
def test_gromacs_pdb2gmx_generates_topology_and_structure(tmp_path: Path) -> None:
    if shutil.which("gmx") is None:
        pytest.skip("GROMACS (gmx) not available")

    repo_root = Path(__file__).resolve().parents[1]
    input_pdb = repo_root / "tests/data/FUS_LC_backmap/final.aa.pdb"

    topology_dir = tmp_path / "topology"
    structure_dir = tmp_path / "structure"

    top_path = run_pdb2gmx_for_topology(
        input_pdb=input_pdb,
        output_dir=topology_dir,
        ff_name="amber99sb-ildn",
        molecule_name="Protein_chain_A",
        water_model="none",
        disable_disulfide=True,
        his_type=0,
        his_repeat_count=500,
    )
    gro_path = run_pdb2gmx_for_structure(
        input_pdb=input_pdb,
        output_dir=structure_dir,
        ff_name="amber99sb-ildn",
        water_model="none",
        disable_disulfide=True,
        his_type=0,
        his_repeat_count=500,
    )

    assert top_path.exists()
    assert gro_path.exists()


@pytest.mark.slow
def test_openmm_builds_system_from_gromacs_topology(tmp_path: Path) -> None:
    if shutil.which("gmx") is None:
        pytest.skip("GROMACS (gmx) not available")

    pytest.importorskip("openmm")

    import openmm.unit as unit
    from openmm.app import GromacsGroFile
    from openmm.app import forcefield as ff

    repo_root = Path(__file__).resolve().parents[1]
    input_pdb = repo_root / "tests/data/FUS_LC_backmap/final.aa.pdb"

    topology_dir = tmp_path / "topology"
    structure_dir = tmp_path / "structure"

    top_path = run_pdb2gmx_for_topology(
        input_pdb=input_pdb,
        output_dir=topology_dir,
        ff_name="amber99sb-ildn",
        molecule_name="Protein_chain_A",
        water_model="none",
        disable_disulfide=True,
        his_type=0,
        his_repeat_count=500,
    )
    gro_path = run_pdb2gmx_for_structure(
        input_pdb=input_pdb,
        output_dir=structure_dir,
        ff_name="amber99sb-ildn",
        water_model="none",
        disable_disulfide=True,
        his_type=0,
        his_repeat_count=500,
    )

    gro = GromacsGroFile(str(gro_path))
    top = GromacsTopFileWithSoftcore(
        str(top_path),
        periodicBoxVectors=gro.getPeriodicBoxVectors(),
        forcefield_type="AMBER",
    )
    system = top.createSystem(
        nonbondedCutoff=1.1 * unit.nanometer,
        nonbondedMethod=ff.CutoffPeriodic,
        nonbonded_type=NONBONDED_STANDARD,
        add_implicit_solvent=False,
    )

    assert system.getNumParticles() > 0
