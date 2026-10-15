"""Quick end-to-end tests of `run_vulcan` on coarse grids.

Each run goes through the whole pipeline (grid, initial abundances, rates, photochemistry,
Ros2 integration to convergence, output) on 10-20 layers, taking a few seconds. All runs
start from well-mixed abundances (ini_mix = 'const_mix'), so FastChem is not needed. The
assertions are physical properties of the converged state, checked against references that
do not come from VULCAN's kinetics:

- detailed balance of every reversible reaction at thermochemical equilibrium;
- geometric dilution of the stellar spectrum, (R_star / a)^2, at the top of the atmosphere.

The default CHO network is used throughout, since `chem_funs.py` is generated once per
session (tests/conftest.py) and other modules bind to it at import.
"""

from __future__ import annotations

import pickle

import numpy as np
import pytest

from vulcan import Config, chem_funs, run_vulcan
from vulcan.chem_funs import spec_list
from vulcan.phy_const import au, r_sun


def _config(out_dir, **overrides):
    """Default config writing only a pickle to out_dir, with plotting off."""
    cfg = Config()
    cfg.output_dir = str(out_dir) + '/'
    cfg.plot_dir = str(out_dir / 'plot') + '/'
    cfg.plot_TP = False
    cfg.use_live_plot = False
    cfg.use_plot_end = False
    cfg.use_plot_evo = False
    cfg.use_print_prog = False
    cfg.save_type = ['pickle']
    for key, val in overrides.items():
        setattr(cfg, key, val)
    return cfg


def _run(cfg):
    run_vulcan(cfg, make_network=False)
    with open(cfg.output_dir + cfg.out_name + '.pkl', 'rb') as hdl:
        return pickle.load(hdl)


# ---------------------------------------------------------------------------
# Thermochemical equilibrium: hot isothermal column, no photochemistry, no transport
# ---------------------------------------------------------------------------

# 1500 K between 100 and 0.1 bar, starting from CH4 and H2O with no CO: the kinetics must
# convert much of the carbon to CO at low pressure, a large change in composition.
EQ_SETTINGS = dict(
    nz=10,
    atm_type='isothermal',
    Tiso=1500.0,
    P_b=1e8,
    P_t=1e5,
    use_photo=False,
    use_moldiff=False,
    Kzz_prof='const',
    const_Kzz=0.0,  # no transport: each layer relaxes to its own local equilibrium
    diff_esc=[],
    ini_mix='const_mix',
    const_mix={'CH4': 2.5e-4, 'H2O': 5.0e-4, 'He': 0.16, 'H2': 1.0 - 0.16 - 7.5e-4},
)


@pytest.fixture(scope='module')
def equilibrium_run(tmp_path_factory):
    """Run from the well-mixed start to convergence."""
    return _run(_config(tmp_path_factory.mktemp('vulcan_eq'), **EQ_SETTINGS))


@pytest.mark.integration
def test_equilibrium_satisfies_detailed_balance(equilibrium_run):
    """At thermochemical equilibrium each reversible reaction is individually balanced.

    A steady state only requires the net tendency of each species to vanish; equilibrium
    requires forward rate = reverse rate for every reaction pair. Rates are recomputed here
    from the stored k and number densities as k * prod(reactant densities), with M the total
    density, independently of the generated chemdf.
    """
    run = equilibrium_run
    assert run['parameter']['end_case'] == 1  # converged, not stopped by a limit
    var, atm = run['variable'], run['atm']
    y, k, x = var['y'], var['k'], var['ymix']

    def density_product(i, side):
        """Product of the reactant (side 0) or product (side 1) densities of reaction i."""
        prod = np.ones_like(atm['M'])
        for sp in chem_funs.re_wM_dict[i][side]:
            prod = prod * (atm['M'] if sp == 'M' else y[:, spec_list.index(sp)])
        return prod

    n_checked = 0
    for i in range(1, chem_funs.nr + 1, 2):  # forward reactions are odd, reverse is i + 1
        if not (np.any(k[i]) and np.any(k[i + 1])):
            continue  # irreversible or switched off
        species = [sp for sp in sum(chem_funs.re_wM_dict[i], []) if sp != 'M']
        # layers where every participant is above VULCAN's convergence floor (mtol_conv)
        layers = np.all([x[:, spec_list.index(sp)] > 1e-20 for sp in species], axis=0)
        if not layers.any():
            continue
        # the reverse reaction i + 1 consumes the products of reaction i
        forward = (k[i] * density_product(i, 0))[layers]
        reverse = (k[i + 1] * density_product(i, 1))[layers]
        # tolerance: converged to the 1e-2 steady-state criterion on abundances of up to four
        # species per rate; 2e-5 is seen in practice
        np.testing.assert_allclose(forward, reverse, rtol=1e-3, err_msg=f'reaction {i}')
        n_checked += 1
    assert n_checked >= 100


@pytest.mark.integration
def test_repeated_run_in_one_process_is_bitwise_identical(equilibrium_run, tmp_path):
    """A second run with an identical, fresh config reproduces the first exactly.

    PROTEUS calls run_vulcan repeatedly in one process, so no state may leak between runs.
    """
    first = equilibrium_run
    second = _run(_config(tmp_path, **EQ_SETTINGS))
    assert second['parameter']['count'] == first['parameter']['count']
    np.testing.assert_array_equal(second['variable']['y'], first['variable']['y'])


# ---------------------------------------------------------------------------
# Photochemistry: default HD 189733b column on a coarse grid
# ---------------------------------------------------------------------------


@pytest.fixture(scope='module')
def photochem_run(tmp_path_factory):
    """Default model on 20 layers, started from the default well-mixed abundances."""
    tmp = tmp_path_factory.mktemp('vulcan_photo')
    cfg = _config(tmp, nz=20, ini_mix='const_mix')
    return cfg, _run(cfg)


@pytest.mark.integration
def test_photochemical_run_converges_to_a_valid_state(photochem_run):
    """The run converges, mixing ratios sum to one, densities are non-negative and the
    column is in hydrostatic balance (n_i = n_0 * x_i in every layer)."""
    _, run = photochem_run
    var, atm = run['variable'], run['atm']
    assert run['parameter']['end_case'] == 1
    np.testing.assert_allclose(var['ymix'].sum(axis=1), 1.0, rtol=1e-12)  # round-off
    assert np.all(var['y'] >= 0)
    np.testing.assert_allclose(var['y'], atm['n_0'][:, None] * var['ymix'], rtol=1e-12)


@pytest.mark.integration
def test_top_of_atmosphere_flux_is_diluted_stellar_spectrum(photochem_run):
    """The direct beam entering the top layer is the stellar surface flux times
    (R_star / a)^2, interpolated onto the wavelength bins (geometric dilution)."""
    cfg, run = photochem_run
    var = run['variable']
    raw = np.genfromtxt(cfg.sflux_file, skip_header=1, names=['lambda', 'flux'])
    dilution = (cfg.r_star * r_sun / (au * cfg.orbit_radius)) ** 2
    expected = np.interp(var['bins'], raw['lambda'], raw['flux'], left=0, right=0) * dilution
    assert expected.max() > 0
    # linear interpolation of the same data: machine precision
    np.testing.assert_allclose(var['sflux'][-1], expected, rtol=1e-12, atol=0)


@pytest.mark.integration
def test_direct_beam_is_attenuated_with_depth(photochem_run):
    """Optical depth grows downward, so the direct beam never increases with depth, and
    photolysis rates are non-negative everywhere and positive at the top."""
    _, run = photochem_run
    var = run['variable']
    tau, sflux = var['tau'], var['sflux']  # index 0 is the bottom
    assert np.all(np.diff(tau, axis=0) <= 0)
    assert np.all(np.diff(sflux, axis=0) >= 0)
    assert sflux[0].sum() < sflux[-1].sum()  # some light is absorbed in the column
    j_rates = np.array(list(var['J_sp'].values()))  # (n_branch, nz)
    assert np.all(j_rates >= 0)
    assert np.any(j_rates[:, -1] > 0)


# ---------------------------------------------------------------------------
# Invalid input
# ---------------------------------------------------------------------------

BAD_RUN = dict(nz=8, ini_mix='const_mix', use_photo=False, count_max=20)


@pytest.mark.integration
def test_short_isothermal_run_completes(tmp_path):
    """Control for the invalid-input tests: the same short set-up with T = 1000 K runs and
    stops at the step limit, so the failures below are caused by the bad value alone."""
    cfg = _config(tmp_path, atm_type='isothermal', Tiso=1000.0, **BAD_RUN)
    run = _run(cfg)
    assert run['parameter']['end_case'] == 3  # maximal number of steps reached
    assert np.all(np.isfinite(run['variable']['y']))


@pytest.mark.integration
@pytest.mark.parametrize('tiso', [0.0, -100.0])
def test_non_positive_temperature_raises(tmp_path, tiso):
    """A temperature <= 0 K must stop the run instead of producing output.

    Currently caught only when the resulting NaN rate coefficients reach the linear solver.
    """
    cfg = _config(tmp_path, atm_type='isothermal', Tiso=tiso, **BAD_RUN)
    with pytest.raises(ValueError):
        run_vulcan(cfg, make_network=False)
    assert not (tmp_path / 'example.pkl').exists()


@pytest.mark.integration
@pytest.mark.parametrize(
    'setting, error',
    [({'ini_mix': 'bogus'}, OSError), ({'atm_type': 'bogus'}, ValueError)],
)
def test_unknown_setup_option_raises(tmp_path, setting, error):
    """An unrecognised initial-composition or atmosphere option is rejected."""
    cfg = _config(tmp_path, **{**BAD_RUN, **setting})
    with pytest.raises(error):
        run_vulcan(cfg, make_network=False)
