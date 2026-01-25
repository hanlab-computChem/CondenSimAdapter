from CondenSimAdapter.src.minimize import MinimizeConfig, MinimizeSimulator


def test_minimize_defaults_create_simulator() -> None:
    config = MinimizeConfig()
    simulator = MinimizeSimulator(config)

    assert simulator.pdb2gmx_name
    assert simulator.forcefield_type in {"AMBER", "CHARMM"}


def test_minimize_optimization_mode_updates() -> None:
    config = MinimizeConfig()
    config.set_optimization_mode("low")

    assert config.optimization_mode == "low"
    assert len(config.softcore_lambda_values) == 2
