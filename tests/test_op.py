"""Unit tests for vulcan.op: numerics of transport, time stepping, rates and bookkeeping.

These tests build a small (nz = 10) model state through the same set-up sequence as
`run_vulcan`, for the default CHO network, without running the integration.
"""

from __future__ import annotations

import copy
import logging
import os
import pickle

import numpy as np
import pytest
import scipy.linalg

from vulcan import Config, build_atm, op, store
from vulcan.chem_funs import ni, nr, spec_list

NZ = 10
R_ROS2 = 1.0 + 1.0 / np.sqrt(2.0)  # gamma^-1 of the 2nd-order Rosenbrock scheme


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope='module')
def base_state(tmp_path_factory):
    """Initial state of the default HD 189733b model on a 10-layer grid (no integration)."""
    tmp = tmp_path_factory.mktemp('vulcan_op')
    cfg = Config()
    cfg.nz = NZ
    cfg.ini_mix = 'const_mix'  # avoids the FastChem executable
    cfg.output_dir = str(tmp / 'out') + '/'
    cfg.plot_dir = str(tmp / 'out' / 'plot') + '/'
    cfg.plot_TP = False
    cfg.use_photo = False
    cfg.save_type = ['pickle', 'plain']

    var, atm, para = store.Variables(cfg), store.AtmData(cfg), store.Parameters(cfg)
    make_atm = build_atm.Atm(cfg)
    output = op.Output(cfg)
    atm = make_atm.f_pico(atm)
    atm = make_atm.load_TPK(atm)
    make_atm.mol_diff(atm)
    rate = op.ReadRate(cfg)
    var = rate.read_rate(var, atm)
    var = rate.rev_rate(var, atm)
    var = rate.remove_rate(var)
    ini = build_atm.InitialAbun(cfg)
    var = ini.ini_y(var, atm)
    var = ini.ele_sum(var)
    atm = make_atm.f_mu_dz(var, atm, output)
    make_atm.BC_flux(atm)
    return cfg, var, atm, para


@pytest.fixture
def state(base_state):
    """Independent copy of the base state that a test may mutate."""
    return copy.deepcopy(base_state)


@pytest.fixture
def mixed_state(state):
    """State with an asymmetric, layer-dependent composition and chemistry switched off.

    Chemistry is switched off (all k = 0) so that the right-hand side and the Jacobian
    contain transport only.
    """
    cfg, var, atm, para = state
    rng = np.random.default_rng(0)
    y = rng.uniform(0.01, 1.0, (NZ, ni)) * atm.n_0[:, None] / ni
    var.y = y
    var.ymix = y / y.sum(axis=1, keepdims=True)
    var.y_prev = y.copy()
    var.dt = 1.0e3
    var.k = {i: np.zeros(NZ) for i in range(1, nr + 1)}
    return cfg, var, atm, para


def _column_weights(atm):
    """Control-volume widths implied by the flux-form discretisation of the transport.

    Interior layers use the mean of the two adjacent interface spacings; the boundary
    layers use the single adjacent spacing. With these weights a conservative flux-form
    operator telescopes, so sum_j w_j * (dn_j/dt) = 0 for zero-flux boundaries.
    """
    w = np.empty(NZ)
    w[0], w[-1] = atm.dzi[0], atm.dzi[-1]
    w[1:-1] = 0.5 * (atm.dzi[:-1] + atm.dzi[1:])
    return w


def _relative_column_change(w, tendency):
    """|sum_j w_j dn_j/dt| relative to the gross weighted tendency, per species."""
    gross = np.sum(np.abs(w[:, None] * tendency), axis=0)
    return np.abs(w @ tendency) / np.where(gross > 0, gross, 1.0)


# ---------------------------------------------------------------------------
# safe_rm
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_safe_rm_removes_file_and_tree(tmp_path):
    """safe_rm deletes a plain file and a whole directory tree, reporting success."""
    f = tmp_path / 'a.txt'
    f.write_text('x')
    d = tmp_path / 'tree'
    (d / 'sub').mkdir(parents=True)
    (d / 'sub' / 'b.txt').write_text('y')

    assert op.safe_rm(str(f)) is True
    assert not f.exists()
    assert op.safe_rm(str(d)) is True
    assert not d.exists()


@pytest.mark.unit
def test_safe_rm_missing_and_empty_paths_are_not_removed(tmp_path):
    """A missing path and an empty path are reported as not removed (falsy)."""
    assert op.safe_rm(str(tmp_path / 'does_not_exist')) is False
    # documented as returning bool, currently returns None for ''
    assert not op.safe_rm('')


@pytest.mark.unit
def test_safe_rm_refuses_git_repository(tmp_path, caplog):
    """A directory containing a .git folder is never deleted."""
    repo = tmp_path / 'repo'
    (repo / '.git').mkdir(parents=True)
    with caplog.at_level(logging.WARNING):
        assert op.safe_rm(str(repo)) is False
    assert (repo / '.git').is_dir()
    assert 'Git repository' in caplog.text


@pytest.mark.unit
def test_safe_rm_refuses_current_directory(tmp_path, monkeypatch):
    """The current working directory is never deleted."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'keep.txt').write_text('x')
    assert op.safe_rm(str(tmp_path)) is False
    assert (tmp_path / 'keep.txt').is_file()


# ---------------------------------------------------------------------------
# Tridiagonal solver
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize('nf', [2, 7])
def test_thomas_vec_matches_dense_solve(nf):
    """The Thomas algorithm solves a tridiagonal system with several right-hand sides.

    Reference: numpy.linalg.solve on the assembled dense matrix. a[0] and c[-1] lie
    outside the matrix and are set to NaN to show that they are ignored.
    """
    rng = np.random.default_rng(1)
    a = rng.uniform(-1, 1, nf)
    c = rng.uniform(-1, 1, nf)
    b = 3.0 + rng.uniform(0, 1, nf)  # diagonally dominant, well conditioned
    a[0], c[-1] = np.nan, np.nan
    d = rng.uniform(-5, 5, (nf, 3))

    dense = np.diag(b) + np.diag(a[1:], -1) + np.diag(c[:-1], 1)
    expected = np.linalg.solve(dense, d)

    x = op.ODESolver.thomas_vec(a, b, c, d)
    # well-conditioned algebra: machine-precision tolerance
    np.testing.assert_allclose(x, expected, rtol=1e-12, atol=1e-12)


# ---------------------------------------------------------------------------
# Transport operators
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_eddy_diffusion_vanishes_for_well_mixed_column(mixed_state):
    """Eddy diffusion acts on mixing-ratio gradients, so a well-mixed column is steady.

    A column with height-independent mixing ratios but a number density falling with
    height must give zero tendency; diffusing number density instead would not.
    """
    cfg, var, atm, para = mixed_state
    solver = op.Ros2(cfg)
    x = np.linspace(1.0, 3.0, ni)
    x /= x.sum()
    y_mixed = np.outer(atm.n_0, x)

    scale = np.abs(solver.diffdf_no_mol(var.y, atm)).max()
    tendency = solver.diffdf_no_mol(y_mixed, atm)
    # round-off of O(10) terms of size `scale`
    np.testing.assert_allclose(tendency, 0.0, atol=1e-13 * scale)


@pytest.mark.unit
@pytest.mark.parametrize('with_advection', [False, True])
def test_eddy_diffusion_conserves_column(mixed_state, with_advection):
    """With zero-flux boundaries, eddy diffusion and advection conserve each species' column."""
    cfg, var, atm, para = mixed_state
    if with_advection:
        atm.vz = np.linspace(-5.0, 5.0, NZ - 1)  # cm/s, both signs to exercise upwinding
    solver = op.Ros2(cfg)
    tendency = solver.diffdf_no_mol(var.y, atm)
    # telescoping sum: machine precision
    assert _relative_column_change(_column_weights(atm), tendency).max() < 1e-13


@pytest.mark.unit
@pytest.mark.parametrize('settling', [False, True])
def test_molecular_diffusion_conserves_column_uniform_gravity(mixed_state, settling):
    """Molecular diffusion (and settling) conserve each species' column when g is uniform."""
    cfg, var, atm, para = mixed_state
    atm.g[:] = atm.g.mean()
    solver = op.Ros2(cfg)
    if settling:
        cfg.use_settling = True
        atm.vs = np.zeros((NZ - 1, ni))
        atm.vs[:, spec_list.index('H2O')] = -np.linspace(1.0, 50.0, NZ - 1)
        tendency = solver.diffdf_settling(var.y, atm)
    else:
        tendency = solver.diffdf(var.y, atm)
    assert _relative_column_change(_column_weights(atm), tendency).max() < 1e-13


@pytest.mark.unit
def test_settling_moves_species_downward(mixed_state):
    """Downward settling (vs < 0) removes a species from the top layer and adds it at the bottom.

    Only the settling species may change relative to the no-settling operator.
    """
    cfg, var, atm, para = mixed_state
    solver = op.Ros2(cfg)
    isp = spec_list.index('H2O')
    atm.vs = np.zeros((NZ - 1, ni))
    atm.vs[:, isp] = -10.0  # cm/s, downward

    no_settling = solver.diffdf(var.y, atm)
    cfg.use_settling = True
    extra = solver.diffdf_settling(var.y, atm) - no_settling

    assert extra[-1, isp] < 0  # top layer loses
    assert extra[0, isp] > 0  # bottom layer gains
    others = np.delete(extra, isp, axis=1)
    np.testing.assert_allclose(others, 0.0, atol=1e-12 * np.abs(no_settling).max())


@pytest.mark.unit
@pytest.mark.parametrize(
    'rhs_name, lhs_name, flags',
    [
        ('diffdf_no_mol', 'lhs_jac_no_mol', {'use_moldiff': False}),
        ('diffdf', 'lhs_jac_tot', {}),
        ('diffdf_settling', 'lhs_jac_settling', {'use_settling': True}),
    ],
)
def test_transport_jacobian_matches_rhs(mixed_state, rhs_name, lhs_name, flags):
    """The Rosenbrock matrix is consistent with the transport right-hand side.

    The transport terms are linear in y for fixed total density, f(y) = D y, and the solver
    matrix is lhs = I/(r dt) - J with J = D when chemistry is off. Hence
    (I/(r dt) - lhs) @ y must reproduce f(y) for every solver variant. Advection and
    settling velocities are non-zero so that those terms are exercised.
    """
    cfg, var, atm, para = mixed_state
    for key, val in flags.items():
        setattr(cfg, key, val)
    atm.vz = np.linspace(-5.0, 5.0, NZ - 1)
    atm.vs = np.zeros((NZ - 1, ni))
    atm.vs[:, spec_list.index('H2O')] = -np.linspace(1.0, 50.0, NZ - 1)
    solver = op.Ros2(cfg)

    rhs = getattr(solver, rhs_name)(var.y, atm).flatten()
    lhs = getattr(solver, lhs_name)(var, atm)
    c0 = 1.0 / (R_ROS2 * var.dt)
    jac = c0 * np.eye(NZ * ni) - lhs

    # cancellation between terms of size max|rhs|: ~1e-12 relative is round-off
    np.testing.assert_allclose(jac @ var.y.flatten(), rhs, atol=1e-10 * np.abs(rhs).max())


# ---------------------------------------------------------------------------
# Banded linear solver
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_solver_matrix_bandwidth_is_number_of_species(state):
    """The Rosenbrock matrix with chemistry and transport has bandwidth exactly ni.

    Chemistry couples species within a layer (|row - col| < ni) and transport couples each
    species only to itself in the adjacent layers (|row - col| = ni). band_lu relies on
    this: any entry further from the diagonal would be silently dropped.
    """
    cfg, var, atm, para = state
    var.dt = 1.0e3
    lhs = op.Ros2(cfg).lhs_jac_tot(var, atm)
    rows, cols = np.nonzero(lhs)
    offset = np.abs(rows - cols)
    assert offset.max() == ni


@pytest.mark.unit
def test_band_solve_matches_dense_solve(mixed_state):
    """One band LU factorisation solves the system for several right-hand sides.

    Reference: numpy.linalg.solve on the dense matrix. Chemistry is off, so the matrix is
    the well-conditioned transport operator plus I/(r dt).
    """
    cfg, var, atm, para = mixed_state
    lhs = op.Ros2(cfg).lhs_jac_tot(var, atm)
    rng = np.random.default_rng(2)
    lu = op.Ros2.band_lu(lhs, ni)
    for _ in range(2):
        b = rng.uniform(-1.0, 1.0, NZ * ni)
        expected = np.linalg.solve(lhs, b)
        # well-conditioned algebra: machine-precision tolerance
        np.testing.assert_allclose(op.Ros2.band_solve(lu, b), expected, rtol=1e-12, atol=0)


@pytest.mark.unit
def test_band_lu_rejects_singular_and_non_finite_matrices():
    """band_lu raises like scipy.linalg.solve_banded for singular or non-finite input."""
    nb, nlev = 3, 4
    a = np.eye(nb * nlev)
    a[5, 5] = 0.0  # zero row and column: exactly singular
    with pytest.raises(np.linalg.LinAlgError):
        op.Ros2.band_lu(a, nb)

    a = np.eye(nb * nlev)
    a[2, 3] = np.nan
    with pytest.raises(ValueError):
        op.Ros2.band_lu(a, nb)


# ---------------------------------------------------------------------------
# Rosenbrock time step
# ---------------------------------------------------------------------------


def _ros2_step(cfg, var, atm, para, dt):
    v = copy.deepcopy(var)
    p = copy.deepcopy(para)
    v.dt = dt
    v.y_prev = v.y.copy()
    v, p = op.Ros2(cfg).solver(v, atm, p)
    return v, p


@pytest.mark.unit
def test_ros2_is_second_order_against_matrix_exponential(mixed_state):
    """One Ros2 step converges to the exact solution with local error O(dt^3).

    With chemistry off, no advection and eddy diffusion only, the total density is
    stationary, so dy/dt = D y with a constant matrix D. The exact solution is
    expm(D dt) y0 (scipy.linalg.expm, independent of VULCAN). For a 2nd-order method,
    halving dt reduces the one-step error by 2^3 = 8.
    """
    cfg, var, atm, para = mixed_state
    cfg.use_moldiff = False
    solver = op.Ros2(cfg)
    # D is the transport Jacobian (validated against the RHS in another test)
    c0 = 1.0 / (R_ROS2 * var.dt)
    D = c0 * np.eye(NZ * ni) - solver.lhs_jac_no_mol(var, atm)

    # the asymptotic regime needs |lambda_max| dt << 1; check it rather than assume it
    lam_max = np.abs(np.linalg.eigvals(D)).max()
    dts = (1250.0, 625.0, 312.5)
    assert lam_max * dts[0] < 0.1

    errors = []
    for dt in dts:
        exact = scipy.linalg.expm(D * dt) @ var.y.flatten()
        v, _ = _ros2_step(cfg, var, atm, para, dt)
        errors.append(np.abs(v.y.flatten() - exact).max() / np.abs(exact).max())

    ratios = np.array(errors[:-1]) / np.array(errors[1:])
    # asymptotic ratio 8 (a 1st-order method gives 4, 3rd-order 16); higher-order terms
    # still lower it by a few per cent at these steps
    assert np.all((ratios > 7.0) & (ratios < 9.0)), ratios
    assert errors[0] < 1e-8


@pytest.mark.unit
def test_ros2_step_conserves_column_and_keeps_well_mixed_state(mixed_state):
    """A Ros2 step conserves every column, and a well-mixed column is a fixed point."""
    cfg, var, atm, para = mixed_state
    cfg.use_moldiff = False
    w = _column_weights(atm)

    v, p = _ros2_step(cfg, var, atm, para, 1e5)
    # each stage solves (I/(r dt) - D) k = f with w^T D = 0, so w^T k = 0 to round-off
    np.testing.assert_allclose(w @ v.y, w @ var.y, rtol=1e-12)
    np.testing.assert_allclose(v.ymix.sum(axis=1), 1.0, rtol=1e-12)
    assert p.delta >= 0

    x = np.linspace(1.0, 3.0, ni)
    x /= x.sum()
    var.y = np.outer(atm.n_0, x)
    var.ymix = var.y / var.y.sum(axis=1, keepdims=True)
    v, _ = _ros2_step(cfg, var, atm, para, 1e5)
    np.testing.assert_allclose(v.y, var.y, rtol=1e-12)


@pytest.mark.unit
@pytest.mark.parametrize(
    'delta_over_rtol, expected_factor',
    [
        (0.25, 0.9 * 2.0),  # 0.9 * sqrt(4) = 1.8, inside the [0.5, 2] clamp
        (100.0, 0.5),  # 0.9 * 0.1 = 0.09, clamped to dt_var_min
        (0.0, 2.0),  # delta = 0 treated as 0.01 rtol -> 9, clamped to dt_var_max
    ],
)
def test_step_size_follows_error_controller(state, delta_over_rtol, expected_factor):
    """The step grows/shrinks by 0.9 sqrt(rtol/delta), clamped to [dt_var_min, dt_var_max]."""
    cfg, var, atm, para = state
    var.dt = 1.0e3
    para.delta = delta_over_rtol * cfg.rtol
    out = op.Ros2(cfg).step_size(
        var,
        para,
        dt_var_min=cfg.dt_var_min,
        dt_var_max=cfg.dt_var_max,
        dt_min=cfg.dt_min,
        dt_max=cfg.dt_max,
    )
    assert out.dt == pytest.approx(1.0e3 * expected_factor, rel=1e-12)


@pytest.mark.unit
def test_step_size_respects_dt_bounds(state):
    """The new step never leaves [dt_min, dt_max]."""
    cfg, var, atm, para = state
    solver = op.Ros2(cfg)
    para.delta = 1e-6 * cfg.rtol
    var.dt = 0.9e3
    assert solver.step_size(var, para, 0.5, 2.0, 1e-8, 1e3).dt == pytest.approx(1e3)
    para.delta = 1e6 * cfg.rtol
    var.dt = 1.5e-8
    assert solver.step_size(var, para, 0.5, 2.0, 1e-8, 1e3).dt == pytest.approx(1e-8)


# ---------------------------------------------------------------------------
# Clipping, element conservation and step acceptance
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_loss_counts_atoms_from_composition(state):
    """Element totals follow the chemical formulae: H = 2 H2O + 4 CH4 + 2 H2, etc."""
    cfg, var, atm, para = state
    n_h2o, n_ch4, n_co2, n_h2 = 3.0e12, 5.0e11, 2.0e10, 7.0e15
    y = np.zeros((NZ, ni))
    y[:, spec_list.index('H2O')] = n_h2o
    y[:, spec_list.index('CH4')] = n_ch4
    y[:, spec_list.index('CO2')] = n_co2
    y[:, spec_list.index('H2')] = n_h2
    var.y = y
    var.atom_ini = {'H': 1.0e17, 'O': 1.0e13, 'C': 1.0e13}

    var = op.Ros2(cfg).loss(var)

    expected = {
        'H': NZ * (2 * n_h2o + 4 * n_ch4 + 2 * n_h2),
        'O': NZ * (n_h2o + 2 * n_co2),
        'C': NZ * (n_ch4 + n_co2),
    }
    for atom, total in expected.items():
        # integer stoichiometry times exact inputs: machine precision
        assert var.atom_sum[atom] == pytest.approx(total, rel=1e-12)
        assert var.atom_loss[atom] == pytest.approx(
            (total - var.atom_ini[atom]) / var.atom_ini[atom], rel=1e-12
        )


@pytest.mark.unit
def test_clip_zeroes_negatives_and_renormalises(state):
    """Negative densities are zeroed and the mixing ratios renormalised.

    Values in [nega_cut, pos_cut) are zeroed and their magnitude is counted in nega_y.
    Larger negatives are not counted, but are still zeroed because their (negative)
    mixing ratio is below mtol.
    """
    cfg, var, atm, para = state
    y = np.abs(var.y) + 1.0
    i_co = spec_list.index('CO')
    y[2, i_co] = -0.5  # inside [nega_cut=-1, 0): clipped
    y[4, i_co] = -5.0  # below nega_cut: not counted, zeroed by the mtol rule
    var.y = y
    var.ymix = y / y.sum(axis=1, keepdims=True)
    para.nega_y = 0.0

    var, para = op.Ros2(cfg).clip(var, para, atm, pos_cut=0, nega_cut=-1)

    assert var.y[2, i_co] == 0.0
    assert var.y[4, i_co] == 0.0
    assert np.all(var.y >= 0)
    assert para.nega_y == pytest.approx(0.5, rel=1e-12)
    np.testing.assert_allclose(var.ymix.sum(axis=1), 1.0, rtol=1e-12)


def _ready_for_acceptance(var, para, cfg):
    var.y_prev = var.y.copy()
    var.atom_loss_prev = dict.fromkeys(cfg.atom_list, 0.0)
    var.atom_loss = dict.fromkeys(cfg.atom_list, 0.0)
    para.delta = 0.1 * cfg.rtol
    var.dt = 1.0e3


@pytest.mark.unit
def test_step_ok_rejects_each_failure_mode(state):
    """A step is accepted only if y >= 0, the element budget is stable and delta <= rtol."""
    cfg, var, atm, para = state
    solver = op.Ros2(cfg)
    _ready_for_acceptance(var, para, cfg)
    assert solver.step_ok(var, para, loss_eps=cfg.loss_eps, rtol=cfg.rtol)

    para.delta = 2.0 * cfg.rtol
    assert not solver.step_ok(var, para, loss_eps=cfg.loss_eps, rtol=cfg.rtol)

    _ready_for_acceptance(var, para, cfg)
    var.atom_loss['C'] = 2.0 * cfg.loss_eps
    assert not solver.step_ok(var, para, loss_eps=cfg.loss_eps, rtol=cfg.rtol)

    _ready_for_acceptance(var, para, cfg)
    var.y = var.y.copy()
    var.y[0, 0] = -1.0
    assert not solver.step_ok(var, para, loss_eps=cfg.loss_eps, rtol=cfg.rtol)


@pytest.mark.unit
def test_step_reject_truncation_error_resets_and_shrinks_step(state):
    """A step with delta > rtol is undone and retried with dt * dt_var_min."""
    cfg, var, atm, para = state
    solver = op.Ros2(cfg)
    _ready_for_acceptance(var, para, cfg)
    y_prev = var.y_prev.copy()
    var.y = var.y * 3.0
    para.delta = 2.0 * cfg.rtol

    stop = solver.step_reject(var, para, loss_eps=cfg.loss_eps, rtol=cfg.rtol)

    assert stop is False
    assert para.delta_count == 1
    np.testing.assert_array_equal(var.y, y_prev)
    assert var.dt == pytest.approx(1.0e3 * cfg.dt_var_min, rel=1e-12)


@pytest.mark.unit
def test_step_reject_negative_reports_the_right_species(state, caplog):
    """A negative solution is counted and logged with its layer and species name."""
    cfg, var, atm, para = state
    cfg.use_print_prog = True
    solver = op.Ros2(cfg)
    _ready_for_acceptance(var, para, cfg)
    var.y = var.y.copy()
    var.y[3, spec_list.index('C2H2')] = -2.0

    with caplog.at_level(logging.WARNING):
        solver.step_reject(var, para, loss_eps=cfg.loss_eps, rtol=cfg.rtol)

    assert para.nega_count == 1
    assert "['C2H2']" in caplog.text
    assert 'from levels: [3]' in caplog.text


@pytest.mark.unit
def test_step_reject_element_loss_branch(state, caplog):
    """A step rejected for element non-conservation is counted as a loss rejection."""
    cfg, var, atm, para = state
    cfg.use_print_prog = True
    solver = op.Ros2(cfg)
    _ready_for_acceptance(var, para, cfg)
    var.atom_loss['H'] = 1.0

    with caplog.at_level(logging.WARNING):
        solver.step_reject(var, para, loss_eps=cfg.loss_eps, rtol=cfg.rtol)

    assert para.loss_count == 1
    assert 'Element conservation' in caplog.text


@pytest.mark.unit
def test_step_reject_gives_up_at_minimum_step(state):
    """Once dt would drop below dt_min, negatives are clipped and the step is forced through."""
    cfg, var, atm, para = state
    solver = op.Ros2(cfg)
    _ready_for_acceptance(var, para, cfg)
    var.y_prev = var.y.copy()
    var.y_prev[1, 1] = -3.0  # the state restored by the rejection is negative
    var.dt = cfg.dt_min  # halving takes it below dt_min
    para.delta = 2.0 * cfg.rtol

    assert solver.step_reject(var, para, loss_eps=cfg.loss_eps, rtol=cfg.rtol) is True
    assert var.dt == cfg.dt_min
    assert np.all(var.y >= 0)


# ---------------------------------------------------------------------------
# Rate coefficients
# ---------------------------------------------------------------------------


def _reaction_index(var, equation):
    for i, eq in var.Rf.items():
        if eq == equation:
            return i
    pytest.skip(f"'{equation}' is not in the network")


@pytest.mark.unit
def test_low_temperature_limits_apply_only_below_thresholds(state):
    """Low-T caps replace k only at or below each reaction's temperature threshold.

    Caps and thresholds are those cited in the code from Moses et al. (2005); they are
    not independently verified here. CH3 + H uses the Lindemann form k0/(1 + k0 M/kinf)
    with k0 = 6e-29 cm6/s and kinf = 2.06e-10 T^-0.4 cm3/s.
    """
    cfg, var, atm, para = state
    # temperatures straddle the thresholds 200, 277.5 and 300 K
    atm.Tco = np.array([150.0, 199.0, 201.0, 250.0, 277.0, 278.0, 299.0, 301.0, 400.0, 900.0])
    atm.M = np.logspace(19, 13, NZ)

    i_ch4 = _reaction_index(var, 'H + CH3 + M -> CH4 + M')
    i_c2h5 = _reaction_index(var, 'H + C2H4 + M -> C2H5 + M')
    i_c2h6 = _reaction_index(var, 'H + C2H5 + M -> C2H6 + M')
    sentinel = 1.0
    for i in (i_ch4, i_c2h5, i_c2h6):
        var.k[i] = np.full(NZ, sentinel)
    untouched = {i: var.k[i].copy() for i in (1, 3, i_ch4 + 1)}

    var = op.ReadRate(cfg).lim_lowT_rates(var, atm)

    T, M = atm.Tco, atm.M
    k0, kinf = 6e-29, 2.06e-10 * T**-0.4
    exp_ch4 = np.where(T <= 277.5, k0 / (1.0 + k0 * M / kinf), sentinel)
    exp_c2h5 = np.where(T <= 300.0, 3.7e-30, sentinel)
    exp_c2h6 = np.where(T <= 200.0, 2.49e-27, sentinel)

    np.testing.assert_allclose(var.k[i_ch4], exp_ch4, rtol=1e-12)
    np.testing.assert_allclose(var.k[i_c2h5], exp_c2h5, rtol=1e-12)
    np.testing.assert_allclose(var.k[i_c2h6], exp_c2h6, rtol=1e-12)
    for i, k in untouched.items():
        np.testing.assert_array_equal(var.k[i], k)


@pytest.mark.unit
def test_remove_rate_zeroes_listed_reactions(state):
    """Reactions in remove_list have zero rate coefficients, now and when re-evaluated."""
    cfg, var, atm, para = state
    cfg.remove_list = [1, 2]
    assert np.any(var.k[1] > 0)
    var = op.ReadRate(cfg).remove_rate(var)
    for i in (1, 2):
        np.testing.assert_array_equal(var.k[i], 0.0)
        np.testing.assert_array_equal(var.k_fun[i](atm.Tco, atm.M), 0.0)
    assert np.any(var.k[3] > 0)


# ---------------------------------------------------------------------------
# Stopping criteria
# ---------------------------------------------------------------------------


def _integration(cfg):
    return op.Integration(op.Ros2(cfg), op.Output(cfg), cfg)


def _history(var, para, n, last_change=0.0):
    """Fill the stored history with n steps of a constant state, optionally changing the last."""
    para.count = n
    var.t_time = list(np.logspace(0, 8, n))
    var.t = var.t_time[-1]
    var.y_time = [var.y.copy() for _ in range(n)]
    var.y_time[-1] = var.y * (1.0 + last_change)
    var.ymix = var.y / var.y.sum(axis=1, keepdims=True)
    var.aflux_change = 0.0


@pytest.mark.unit
def test_stop_reports_convergence(state):
    """A run whose abundances no longer change, after count_min steps, ends as converged."""
    cfg, var, atm, para = state
    _history(var, para, cfg.count_min + 10)
    assert _integration(cfg).stop(var, para, atm)
    assert para.end_case == 1
    assert var.longdy == 0.0


@pytest.mark.unit
def test_stop_does_not_converge_while_abundances_change(state):
    """A 50% change over the convergence window is not steady state; the run continues."""
    cfg, var, atm, para = state
    _history(var, para, cfg.count_min + 10, last_change=0.5)
    assert not _integration(cfg).stop(var, para, atm)
    assert para.end_case == 0
    # longdy = max |dy| / n_0 / ymix = 0.5 * ysum / n_0, maximised over layers
    expected = 0.5 * np.max(var.y.sum(axis=1) / atm.n_0)
    assert var.longdy == pytest.approx(expected, rel=1e-12)


@pytest.mark.unit
def test_stop_at_runtime_and_step_limits(state):
    """Exceeding runtime gives end_case 2; exceeding count_max gives end_case 3."""
    cfg, var, atm, para = state
    integ = _integration(cfg)
    para.start_time = 0.0

    var.t, para.count = 2.0 * cfg.runtime, 5  # too few steps to test convergence
    assert integ.stop(var, para, atm)
    assert para.end_case == 2

    para.end_case = 0
    var.t, para.count = 1.0, cfg.count_max + 1
    assert integ.stop(var, para, atm)
    assert para.end_case == 3


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_save_out_writes_pickle_and_plain(state):
    """The pickle round-trips the final abundances; the plain files are labelled text."""
    cfg, var, atm, para = state
    cfg.use_plot_evo = False
    cfg.use_plot_end = False
    cfg.save_type = ['pickle', 'plain']
    cfg.out_name = 'unit.pkl'  # a known extension is replaced, not doubled
    out = op.Output(cfg)
    var.y_time, var.t_time = [], []

    out.save_out(var, atm, para)

    pkl = os.path.join(cfg.output_dir, 'unit.pkl')
    with open(pkl, 'rb') as hdl:
        data = pickle.load(hdl)
    assert data['variable']['species'] == spec_list
    np.testing.assert_array_equal(data['variable']['ymix'], var.ymix)
    np.testing.assert_array_equal(data['atm']['Tco'], atm.Tco)

    for ext in ('.var', '.atm', '.par'):
        assert os.path.isfile(os.path.join(cfg.output_dir, 'unit' + ext))
    with open(os.path.join(cfg.output_dir, 'unit.var')) as hdl:
        assert hdl.readline().strip() == '# species'
    assert not os.path.exists(os.path.join(cfg.output_dir, 'unit.nc'))


@pytest.mark.unit
def test_output_rejects_unknown_save_type(state):
    """An unknown output format is a configuration error, raised before the run starts."""
    cfg, var, atm, para = state
    cfg.save_type = 'hdf5'
    with pytest.raises(ValueError, match='Invalid save_type'):
        op.Output(cfg)
