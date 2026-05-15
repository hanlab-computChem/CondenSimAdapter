"""Internal module: softcore force construction functions."""

from __future__ import absolute_import

import openmm as mm
import openmm.unit as unit


def _create_gaussian_nb_force(nonbonded_cutoff, gaussian_width=0.11):
    """Create Gaussian repulsion force for initial untangling (Stage 1).

    Formula: E = ga_k * ga_h * exp(-(r/ga_w)^2)

    This provides pure geometric repulsion without chemical information,
    useful for quickly removing atomic overlaps and untangling knots.

    Parameters
    ----------
    nonbonded_cutoff : Quantity
        Cutoff distance for nonbonded interactions
    gaussian_width : float=0.11
        Width parameter (ga_w) in nm. Default 0.11 nm (1.1 Å).

    Returns
    -------
    CustomNonbondedForce
        Gaussian repulsion force object
    """
    nb_force = mm.CustomNonbondedForce("ga_k * ga_h * exp(-(r/ga_w)^2);")
    nb_force.addGlobalParameter("ga_k", 1.0)  # ON/OFF switch
    nb_force.addGlobalParameter("ga_h", 800.0)  # Height (kJ/mol)
    nb_force.addGlobalParameter("ga_w", gaussian_width)  # Width (nm)
    nb_force.addPerParticleParameter("dummy_type")
    nb_force.addPerParticleParameter("dummy_q")
    nb_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    nb_force.setCutoffDistance(nonbonded_cutoff.value_in_unit(unit.nanometer))
    return nb_force


def _create_standard_custom_nb_force(nonbonded_cutoff, switchDistance=None):
    """Create standard LJ+Coulomb CustomNonbondedForce for debugging/comparison.

    This implements the exact same formula as NonbondedForce but using
    CustomNonbondedForce to match the architecture of softcore mode.

    Parameters
    ----------
    nonbonded_cutoff : Quantity
        Cutoff distance (e.g. 1.0 * unit.nanometer)
    switchDistance : Quantity or None
        Switching distance for smoothing

    Returns
    -------
    nb_force : openmm.CustomNonbondedForce
        CustomNonbondedForce implementing standard LJ+Coulomb potential
    """
    # Standard LJ + Coulomb formula
    # LJ: 4*epsilon*[(sigma/r)^12 - (sigma/r)^6]
    # Coulomb: (1/4πε0) * q1*q2/r = 138.935456 * q1*q2/r (kJ/mol/nm/e^2)

    energy_expr = """
4.0*sqrt(epsilon1*epsilon2)*((0.5*(sigma1+sigma2)/r)^12 - (0.5*(sigma1+sigma2)/r)^6)
+ ONE_4PI_EPS0*q1*q2/r
"""

    nb_force = mm.CustomNonbondedForce(energy_expr)

    # Global parameters
    nb_force.addGlobalParameter("ONE_4PI_EPS0", 138.935456)  # kJ/mol/nm/e^2

    # Per-particle parameters
    nb_force.addPerParticleParameter("q")  # Charge
    nb_force.addPerParticleParameter("sigma")  # LJ sigma
    nb_force.addPerParticleParameter("epsilon")  # LJ epsilon

    # Set cutoff
    nb_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    nb_force.setCutoffDistance(nonbonded_cutoff.value_in_unit(unit.nanometer))

    # Use switching function if provided
    if switchDistance is not None:
        nb_force.setUseSwitchingFunction(True)
        if unit.is_quantity(switchDistance):
            nb_force.setSwitchingDistance(switchDistance.value_in_unit(unit.nanometer))
        else:
            nb_force.setSwitchingDistance(switchDistance)
    else:
        nb_force.setUseSwitchingFunction(False)

    # Disable LRC (to match standard mode behavior with CutoffPeriodic)
    nb_force.setUseLongRangeCorrection(False)

    return nb_force


def _create_standard_pair_force_expression():
    """Create energy expression for standard 1-4 pair interactions.

    This implements the standard LJ+Coulomb formula for CustomBondForce,
    matching what NonbondedForce does for exceptions.

    For CustomBondForce, we use per-bond parameters:
    - charge_prod: q1 * q2 (already scaled by fudgeQQ)
    - sigma: 0.5 * (sigma1 + sigma2)
    - epsilon: sqrt(epsilon1 * epsilon2) (already scaled by fudgeLJ)

    Returns
    -------
    str
        Energy expression string for CustomBondForce
    """
    # Standard LJ + Coulomb (parameters are already combined and scaled)
    energy_expr = """
4.0*epsilon*((sigma/r)^12 - (sigma/r)^6) + ONE_4PI_EPS0*charge_prod/r
"""
    return energy_expr


def _create_gapsys_pair_force_expression(has_nbfix_terms=False):
    """Create energy expression for Gapsys softcore pair interactions (1-4 pairs).

    This function generates the same Gapsys linearized formula used in
    _create_gapsys_linearized_nb_force, but adapted for CustomBondForce.

    For CustomBondForce, we use per-bond parameters:
    - charge_prod: q1 * q2 (already scaled by fudgeQQ)
    - sigma: 0.5 * (sigma1 + sigma2) OR derived from C6/C12 for NBFIX
    - epsilon: sqrt(epsilon1 * epsilon2) (already scaled by fudgeLJ) OR derived from C6/C12 for NBFIX

    Parameters
    ----------
    has_nbfix_terms : bool=False
        If True, use C6/C12-based formula compatible with NBFIX parameters

    Returns
    -------
    str
        Energy expression string for CustomBondForce
    """
    if has_nbfix_terms:
        # NBFIX mode: Parameters are [charge_prod, c6, c12]
        # Formula: E_LJ = C12/r^12 - C6/r^6
        # Need to compute sigma_eff and epsilon_eff for softcore switching

        # For NBFIX pairs, we store C6 and C12 directly (scaled by fudge factors)
        # sigma_eff = (C12/C6)^(1/6)
        # epsilon_eff = C6^2 / (4*C12)

        q_prod = "charge_prod"
        soft = "(1.0 - lambda_val)^0.1666667"
        c6_expr = "c6"
        c12_expr = "c12"

        # Compute sigma_eff and epsilon_eff
        sig_eff = f"(({c12_expr})/({c6_expr}))^0.1666667"
        eps_eff = f"({c6_expr})*({c6_expr})/(4.0*({c12_expr}))"

        # LJ switching parameters
        rsw_lj_inline = f"{sig_eff} * alpha_lj * 1.1224 * {soft}"
        rsw_lj_safe_inline = f"max({rsw_lj_inline}, 1.0e-6)"
        rsw_lj_sq_inline = f"{rsw_lj_safe_inline}^2"
        u_lj_inline = f"{sig_eff} / {rsw_lj_safe_inline}"
        V_lj_sw_inline = f"4.0 * {eps_eff} * ({u_lj_inline}^12 - {u_lj_inline}^6)"
        F_lj_sw_inline = f"(24.0 * {eps_eff} / {rsw_lj_safe_inline}) * (2.0 * {u_lj_inline}^12 - {u_lj_inline}^6)"
        dF_lj_sw_inline = f"(24.0 * {eps_eff} / {rsw_lj_sq_inline}) * (-26.0 * {u_lj_inline}^12 + 7.0 * {u_lj_inline}^6)"

        # Standard LJ from C6/C12
        lj_standard = f"({c12_expr})/(r^12) - ({c6_expr})/(r^6)"

        # Coulomb parameters
        rsw_q_inline = f"alpha_coul * (1.0 + sigma_coul * abs({q_prod})) * {soft}"
        rsw_q_safe_inline = f"max({rsw_q_inline}, 1.0e-6)"
        rsw_q_sq_inline = f"{rsw_q_safe_inline}^2"
        C_pre_inline = f"ONE_4PI_EPS0 * {q_prod}"
        V_q_sw_inline = f"{C_pre_inline} / {rsw_q_safe_inline}"
        F_q_sw_inline = f"{C_pre_inline} / {rsw_q_sq_inline}"
        dF_q_sw_inline = f"-2.0 * {C_pre_inline} / ({rsw_q_sq_inline} * {rsw_q_safe_inline})"

        energy_expr = f"""
lambda_val * (
    select(step(r - {rsw_lj_safe_inline}),
           {lj_standard},
           {V_lj_sw_inline} + {F_lj_sw_inline} * ({rsw_lj_inline} - r) - 0.5 * {dF_lj_sw_inline} * ({rsw_lj_inline} - r)^2)
    +
    select(step(r - {rsw_q_safe_inline}),
           ONE_4PI_EPS0 * {q_prod} / r,
           {V_q_sw_inline} + {F_q_sw_inline} * ({rsw_q_inline} - r) - 0.5 * {dF_q_sw_inline} * ({rsw_q_inline} - r)^2)
);
"""
    else:
        # Standard mode: Parameters are [charge_prod, sigma, epsilon]
        # For CustomBondForce, we have per-bond parameters: charge_prod, sigma, epsilon
        # These are already combined and scaled appropriately

        # Helper expressions for inline use
        e = "epsilon"  # Already combined and scaled
        sig = "sigma"  # Already combined
        q_prod = "charge_prod"  # Already scaled by fudgeQQ
        soft = "(1.0 - lambda_val)^0.1666667"

        # LJ parameters (inline)
        rsw_lj_inline = f"{sig} * alpha_lj * 1.1224 * {soft}"
        rsw_lj_safe_inline = f"max({rsw_lj_inline}, 1.0e-6)"
        rsw_lj_sq_inline = f"{rsw_lj_safe_inline}^2"
        u_lj_inline = f"{sig} / {rsw_lj_safe_inline}"
        V_lj_sw_inline = f"4.0 * {e} * ({u_lj_inline}^12 - {u_lj_inline}^6)"
        F_lj_sw_inline = (
            f"(24.0 * {e} / {rsw_lj_safe_inline}) * (2.0 * {u_lj_inline}^12 - {u_lj_inline}^6)"
        )
        dF_lj_sw_inline = f"(24.0 * {e} / {rsw_lj_sq_inline}) * (-26.0 * {u_lj_inline}^12 + 7.0 * {u_lj_inline}^6)"

        # Coulomb parameters (inline)
        rsw_q_inline = f"alpha_coul * (1.0 + sigma_coul * abs({q_prod})) * {soft}"
        rsw_q_safe_inline = f"max({rsw_q_inline}, 1.0e-6)"
        rsw_q_sq_inline = f"{rsw_q_safe_inline}^2"
        C_pre_inline = f"ONE_4PI_EPS0 * {q_prod}"
        V_q_sw_inline = f"{C_pre_inline} / {rsw_q_safe_inline}"
        F_q_sw_inline = f"{C_pre_inline} / {rsw_q_sq_inline}"
        dF_q_sw_inline = f"-2.0 * {C_pre_inline} / ({rsw_q_sq_inline} * {rsw_q_safe_inline})"

        energy_expr = f"""
lambda_val * (
    select(step(r - {rsw_lj_safe_inline}),
           4.0 * {e} * ({sig}/r)^12 - 4.0 * {e} * ({sig}/r)^6,
           {V_lj_sw_inline} + {F_lj_sw_inline} * ({rsw_lj_inline} - r) - 0.5 * {dF_lj_sw_inline} * ({rsw_lj_inline} - r)^2)
    +
    select(step(r - {rsw_q_safe_inline}),
           ONE_4PI_EPS0 * {q_prod} / r,
           {V_q_sw_inline} + {F_q_sw_inline} * ({rsw_q_inline} - r) - 0.5 * {dF_q_sw_inline} * ({rsw_q_inline} - r)^2)
);
"""
    return energy_expr


def _create_gapsys_linearized_nb_force(
    nonbonded_cutoff,
    current_lambda=1.0,
    alpha_lj=0.85,
    alpha_coul=0.3,
    sigma_coul=1.0,
    use_implicit_solvent=False,
    switchDistance=None,
    has_nbfix_terms=False,
):
    """Create Gapsys "New Soft-Core" (Linearized Force) CustomNonbondedForce.

    Based on Gapsys et al., J. Chem. Theory Comput. 2015, 11, 11, 5920–5930
    "Interaction of Legolane with the Outer Membrane of Gram-Negative Bacteria"

    This implementation linearizes both LJ and Coulomb forces near r=0,
    eliminating singularities and enabling stable Energy Minimization.

    Formula Overview:
    - For r > r_sw (switching distance): Use standard LJ/Coulomb
    - For r <= r_sw: Use Taylor-expanded linearized force

    Parameters
    ----------
    nonbonded_cutoff : Quantity
        Cutoff distance (e.g. 1.0 * unit.nanometer)
    current_lambda : float=1.0
        Decoupling parameter (1.0 = fully coupled, 0.0 = decoupled).
        NOTE: For softcore EM, we use a "softness" lambda (default 0.85)
        to control the softening, NOT the decoupling.
    alpha_lj : float=0.85
        LJ soft-core control parameter (Gapsys paper recommends 0.85)
    alpha_coul : float=0.3
        Coulomb soft-core control parameter (Gapsys paper recommends 0.3)
    sigma_coul : float=1.0
        Coulomb soft-core charge scaling parameter (Gapsys paper recommends 1.0)
    use_implicit_solvent : bool=False
        If True, use CutoffPeriodic with GBSA forces
    has_nbfix_terms : bool=False
        If True, use NBFIX-compatible formula with type parameters

    Returns
    -------
    nb_force : openmm.CustomNonbondedForce
        CustomNonbondedForce implementing linearized soft-core potential
    """
    # NOTE: Main expression must use ONLY built-in variables (r, q1, q2, sigma1, sigma2, epsilon1, epsilon2)
    # and global parameters. Custom intermediate variables CANNOT appear in the main expression.
    # All intermediate variables are defined AFTER the main expression for documentation only.

    if has_nbfix_terms:
        # NBFIX mode: Use tabulated functions for LJ parameters
        # Formula: E_LJ = C12/r^12 - C6/r^6
        # where C6 = bcoef(type1, type2) and C12 = acoef(type1, type2)^2
        #
        # CRITICAL: OpenMM's backward variable definition (defining variables after using them)
        # only works for INDEPENDENT variables. Variables CANNOT reference other backward-defined
        # variables. This was discovered through extensive testing (tests/test_sonnet/charmm/nbfix_gapsys/).
        #
        # Solution: Only use backward definition for c6 and c12 (which depend only on built-in functions).
        # All other calculations must be inlined into the main expression.

        # Helper: inline calculations
        q_prod = "q1 * q2"
        soft_factor = "(1.0 - lambda_val)^0.1666667"

        # sig_eff and eps_eff in terms of c6/c12
        sig_eff = "(c12/c6)^0.1666667"
        eps_eff = "(c6*c6)/(4.0*c12)"

        # LJ softcore switching distance
        rsw_lj = f"({sig_eff} * alpha_lj * 1.1224 * {soft_factor})"
        rsw_lj_safe = f"max({rsw_lj}, 1.0e-6)"
        rsw_lj_sq = f"({rsw_lj_safe} * {rsw_lj_safe})"
        u_lj = f"({sig_eff} / {rsw_lj_safe})"

        # LJ softcore parameters at r_sw
        V_lj_sw = f"(4.0 * {eps_eff} * ({u_lj}^12 - {u_lj}^6))"
        F_lj_sw = f"((24.0 * {eps_eff} / {rsw_lj_safe}) * (2.0 * {u_lj}^12 - {u_lj}^6))"
        dF_lj_sw = f"((24.0 * {eps_eff} / {rsw_lj_sq}) * (-26.0 * {u_lj}^12 + 7.0 * {u_lj}^6))"

        # Coulomb softcore switching distance
        rsw_q = f"(alpha_coul * (1.0 + sigma_coul * abs({q_prod})) * {soft_factor})"
        rsw_q_safe = f"max({rsw_q}, 1.0e-6)"
        rsw_q_sq = f"({rsw_q_safe} * {rsw_q_safe})"
        C_pre = f"(ONE_4PI_EPS0 * {q_prod})"

        # Coulomb softcore parameters at r_sw
        V_q_sw = f"({C_pre} / {rsw_q_safe})"
        F_q_sw = f"({C_pre} / {rsw_q_sq})"
        dF_q_sw = f"(-2.0 * {C_pre} / ({rsw_q_sq} * {rsw_q_safe}))"

        # Build fully inlined energy expression
        energy_expr = f"""
lambda_val * (
    select(step(r - {rsw_lj_safe}),
           c12/(r^12) - c6/(r^6),
           {V_lj_sw} + {F_lj_sw} * ({rsw_lj} - r) - 0.5 * {dF_lj_sw} * ({rsw_lj} - r)^2)
    +
    select(step(r - {rsw_q_safe}),
           ONE_4PI_EPS0 * {q_prod} / r,
           {V_q_sw} + {F_q_sw} * ({rsw_q} - r) - 0.5 * {dF_q_sw} * ({rsw_q} - r)^2)
);

c6 = bcoef(type1, type2);
c12 = (acoef(type1, type2))^2
"""

        nb_force = mm.CustomNonbondedForce(energy_expr)

        # Global parameters
        nb_force.addGlobalParameter("lambda_val", current_lambda)
        nb_force.addGlobalParameter("alpha_lj", alpha_lj)
        nb_force.addGlobalParameter("alpha_coul", alpha_coul)
        nb_force.addGlobalParameter("sigma_coul", sigma_coul)
        nb_force.addGlobalParameter("ONE_4PI_EPS0", 138.935456)  # kJ/mol/nm/e^2

        # Per-particle parameters for NBFIX mode
        nb_force.addPerParticleParameter("type")  # Type index for NBFIX lookup
        nb_force.addPerParticleParameter("q")  # Charge

    else:
        # Standard mode: Use combination rules
        # Helper expressions for inline use
        e = "sqrt(epsilon1 * epsilon2)"
        sig = "0.5 * (sigma1 + sigma2)"
        q_prod = "q1 * q2"
        soft = "(1.0 - lambda_val)^0.1666667"

        # LJ parameters (inline)
        rsw_lj_inline = f"{sig} * alpha_lj * 1.1224 * {soft}"
        rsw_lj_safe_inline = f"max({rsw_lj_inline}, 1.0e-6)"
        rsw_lj_sq_inline = f"{rsw_lj_safe_inline}^2"
        u_lj_inline = f"{sig} / {rsw_lj_safe_inline}"
        V_lj_sw_inline = f"4.0 * {e} * ({u_lj_inline}^12 - {u_lj_inline}^6)"
        F_lj_sw_inline = (
            f"(24.0 * {e} / {rsw_lj_safe_inline}) * (2.0 * {u_lj_inline}^12 - {u_lj_inline}^6)"
        )
        dF_lj_sw_inline = f"(24.0 * {e} / {rsw_lj_sq_inline}) * (-26.0 * {u_lj_inline}^12 + 7.0 * {u_lj_inline}^6)"

        # Coulomb parameters (inline)
        rsw_q_inline = f"alpha_coul * (1.0 + sigma_coul * abs({q_prod})) * {soft}"
        rsw_q_safe_inline = f"max({rsw_q_inline}, 1.0e-6)"
        rsw_q_sq_inline = f"{rsw_q_safe_inline}^2"
        C_pre_inline = f"ONE_4PI_EPS0 * {q_prod}"
        V_q_sw_inline = f"{C_pre_inline} / {rsw_q_safe_inline}"
        F_q_sw_inline = f"{C_pre_inline} / {rsw_q_sq_inline}"
        dF_q_sw_inline = f"-2.0 * {C_pre_inline} / ({rsw_q_sq_inline} * {rsw_q_safe_inline})"

        energy_expr = f"""
lambda_val * (
    select(step(r - {rsw_lj_safe_inline}),
           4.0 * {e} * ({sig}/r)^12 - 4.0 * {e} * ({sig}/r)^6,
           {V_lj_sw_inline} + {F_lj_sw_inline} * ({rsw_lj_inline} - r) - 0.5 * {dF_lj_sw_inline} * ({rsw_lj_inline} - r)^2)
    +
    select(step(r - {rsw_q_safe_inline}),
           ONE_4PI_EPS0 * {q_prod} / r,
           {V_q_sw_inline} + {F_q_sw_inline} * ({rsw_q_inline} - r) - 0.5 * {dF_q_sw_inline} * ({rsw_q_inline} - r)^2)
);

/* Documentation: LJ parameters (NOT used in main expression above) */
e = {e};
sig = {sig};
q_prod = {q_prod};
soft_factor = {soft};

rsw_lj = alpha_lj * 1.1224 * sig * soft_factor;
rsw_lj_safe = max(rsw_lj, 1.0e-6);
rsw_lj_sq = rsw_lj_safe * rsw_lj_safe;
u_lj = sig / rsw_lj_safe;
V_lj_sw  = 4.0 * e * (u_lj^12 - u_lj^6);
F_lj_sw  = (24.0 * e / rsw_lj_safe) * (2.0 * u_lj^12 - u_lj^6);
dF_lj_sw = (24.0 * e / rsw_lj_sq) * (-26.0 * u_lj^12 + 7.0 * u_lj^6);

/* Documentation: Coulomb parameters (NOT used in main expression above) */
rsw_q = alpha_coul * (1.0 + sigma_coul * abs(q_prod)) * soft_factor;
rsw_q_safe = max(rsw_q, 1.0e-6);
rsw_q_sq = rsw_q_safe * rsw_q_safe;
C_pre = ONE_4PI_EPS0 * q_prod;
V_q_sw  = C_pre / rsw_q_safe;
F_q_sw  = C_pre / rsw_q_sq;
dF_q_sw = -2.0 * C_pre / (rsw_q_sq * rsw_q_safe);
"""

        nb_force = mm.CustomNonbondedForce(energy_expr)

        # Global parameters
        nb_force.addGlobalParameter("lambda_val", current_lambda)
        nb_force.addGlobalParameter("alpha_lj", alpha_lj)
        nb_force.addGlobalParameter("alpha_coul", alpha_coul)
        nb_force.addGlobalParameter("sigma_coul", sigma_coul)
        nb_force.addGlobalParameter("ONE_4PI_EPS0", 138.935456)  # kJ/mol/nm/e^2

        # Per-particle parameters (sigma1, epsilon1 for each particle)
        # Note: Using sigma/epsilon from combination rule 2 (sigma/epsilon directly)
        # For combination rule 1/3, we need to convert C6/C12 to sigma/epsilon
        nb_force.addPerParticleParameter("q")  # Charge
        nb_force.addPerParticleParameter("sigma")  # LJ sigma (for rule 2) or derived from C6/C12
        nb_force.addPerParticleParameter("epsilon")  # LJ epsilon

    # Set cutoff
    nb_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    nb_force.setCutoffDistance(nonbonded_cutoff.value_in_unit(unit.nanometer))

    # Use switching function only if switchDistance is provided (matching standard mode behavior)
    if switchDistance is not None:
        nb_force.setUseSwitchingFunction(True)
        if unit.is_quantity(switchDistance):
            nb_force.setSwitchingDistance(switchDistance.value_in_unit(unit.nanometer))
        else:
            nb_force.setSwitchingDistance(switchDistance)
    else:
        nb_force.setUseSwitchingFunction(False)

    # LRC is disabled for Gapsys linearized potential
    # Gapsys linearized potential doesn't support analytical long-range correction
    nb_force.setUseLongRangeCorrection(False)

    return nb_force
