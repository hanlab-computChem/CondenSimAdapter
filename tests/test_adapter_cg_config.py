from pathlib import Path

import yaml

from CondenSimAdapter.src.cg import CGSimulationConfig


def test_cg_config_from_yaml_resolves_relative_paths(tmp_path: Path) -> None:
    fasta_path = tmp_path / "seq.fasta"
    fasta_path.write_text(">A\nAAAA\n", encoding="utf-8")

    config_data = {
        "system_name": "demo_system",
        "box": [10.0, 10.0, 12.0],
        "topol": "grid",
        "components": [
            {
                "name": "ProteinA",
                "type": "idp",
                "nmol": 1,
                "ffasta": "seq.fasta",
            }
        ],
    }

    yaml_path = tmp_path / "config.yaml"
    yaml.safe_dump(config_data, yaml_path.open("w", encoding="utf-8"), sort_keys=False)

    config = CGSimulationConfig.from_yaml(str(yaml_path))

    assert config.components[0].ffasta == str(fasta_path)
    assert config.total_molecules() == 1
    assert config.validate() == []
