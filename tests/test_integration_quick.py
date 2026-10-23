"""Quick end-to-end tests of `run_vulcan` on coarse grids.

Each run goes through the whole pipeline (grid, initial abundances, rates, photochemistry,
Ros2 integration to convergence, output) on 10-20 layers, taking a few seconds. The
assertions are physical properties of the converged state, checked against references that
do not come from VULCAN's kinetics:

- thermochemical equilibrium from FastChem, an independent Gibbs-minimisation code with its
  own thermodynamic data (VULCAN's reverse rates use NASA9 polynomials);
- detailed balance of every reversible reaction at equilibrium;
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

ELEMENTS = ('H', 'O', 'C', 'He')


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


def _element_totals(y):
    """Number of atoms of each element in the column (sum over layers), from composition data."""
    from vulcan.build_atm import compo, compo_row

    return {
        el: sum(
            compo[compo_row.index(sp)][el] * y[:, i] for i, sp in enumerate(spec_list)
        ).sum()
        for el in ELEMENTS
    }


# ---------------------------------------------------------------------------
# Thermochemical equilibrium: hot isothermal column, no photochemistry, no transport
# ---------------------------------------------------------------------------

# 1500 K between 100 and 0.1 bar: CH4 dominates carbon at the bottom and CO at the top, so
# the equilibrium varies strongly with pressure. FastChem's other elements (N, S, ...) are
# scaled to negligible abundance so the element budget matches the CHO network.
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
    use_solar=False,
    C_H=2.5e-4,
    O_H=5.0e-4,
    fastchem_met_scale=1e-10,
)


def _fastchem_equilibrium(out_dir):
    """FastChem equilibrium mixing ratios on the EQ_SETTINGS grid, shape (nz, ni)."""
    from vulcan import build_atm, store

    cfg = _config(out_dir, ini_mix='eq', **EQ_SETTINGS)
    var, atm = store.Variables(cfg), store.AtmData(cfg)
    make_atm = build_atm.Atm(cfg)
    atm = make_atm.f_pico(atm)
    atm = make_atm.load_TPK(atm)
    var = build_atm.InitialAbun(cfg).ini_y(var, atm)
    return var.y / atm.n_0[:, None]


def _well_mixed_start(x_eq):
    """H2, He, H2O, CH4 mixture with the same element ratios as the equilibrium x_eq.

    Built per H atom: CH4 holds all carbon, H2O all oxygen, and H2 the remaining hydrogen.
    CO, the main carbon carrier at the top of the equilibrium column, starts at zero.
    """
    el = _element_totals(x_eq[:1])
    c_h, o_h, he_h = el['C'] / el['H'], el['O'] / el['H'], el['He'] / el['H']
    h2 = (1.0 - 4.0 * c_h - 2.0 * o_h) / 2.0
    total = c_h + o_h + he_h + h2
    return {'CH4': c_h / total, 'H2O': o_h / total, 'He': he_h / total, 'H2': h2 / total}


@pytest.fixture(scope='module')
def equilibrium_runs(tmp_path_factory):
    """FastChem reference plus two runs, one started at it and one from a well-mixed state."""
    tmp = tmp_path_factory.mktemp('vulcan_eq')
    x_eq = _fastchem_equilibrium(tmp / 'ref')
    const_mix = _well_mixed_start(x_eq)
    from_eq = _run(_config(tmp / 'from_eq', ini_mix='eq', **EQ_SETTINGS))
    from_mix = _run(
        _config(tmp / 'from_mix', ini_mix='const_mix', const_mix=const_mix, **EQ_SETTINGS)
    )
    return x_eq, const_mix, from_eq, from_mix


def _assert_matches_equilibrium(run, x_eq):
    assert run['parameter']['end_case'] == 1  # converged, not stopped by a limit
    x = run['variable']['ymix']
    # species and layers with a well-defined equilibrium abundance; trace species far below
    # this are not converged (VULCAN's own check ignores mixing ratios < mtol_conv = 1e-20)
    mask = x_eq > 1e-10
    assert len(set(np.nonzero(mask)[1])) >= 10  # non-vacuous: many species compared
    # rtol: the steady-state criterion allows ~1% residual drift (yconv_cri = 0.01); the two
    # codes' thermodynamic data differ by much less (~4e-4 seen at 1500 K)
    np.testing.assert_allclose(x[mask], x_eq[mask], rtol=1e-2)


@pytest.mark.integration
def test_fastchem_equilibrium_is_a_fixed_point_of_the_kinetics(equilibrium_runs):
    """Started at FastChem's equilibrium, the forward and reverse kinetics leave it unchanged.

    Reverse rates are forward rates divided by equilibrium constants from VULCAN's NASA9
    Gibbs energies; a wrong equilibrium constant would drive the column away from FastChem.
    """
    x_eq, _, from_eq, _ = equilibrium_runs
    _assert_matches_equilibrium(from_eq, x_eq)


@pytest.mark.integration
def test_kinetics_relax_to_fastchem_equilibrium(equilibrium_runs):
    """From a well-mixed CH4/H2O start with no CO, the kinetics reach FastChem's equilibrium.

    This exercises the full Ros2 integration through a large change in composition: CO goes
    from zero to the dominant carbon carrier at 0.1 bar.
    """
    x_eq, const_mix, _, from_mix = equilibrium_runs
    i_co = spec_list.index('CO')
    assert 'CO' not in const_mix
    assert x_eq[-1, i_co] > 0.5 * _element_totals(x_eq[-1:])['C']  # CO dominant at the top
    _assert_matches_equilibrium(from_mix, x_eq)


@pytest.mark.integration
def test_equilibrium_satisfies_detailed_balance(equilibrium_runs):
    """At thermochemical equilibrium each reversible reaction is individually balanced.

    A steady state only requires the net tendency of each species to vanish; equilibrium
    requires forward rate = reverse rate for every reaction pair. Rates are recomputed here
    from the stored k and number densities as k * prod(reactant densities), with M the total
    density, independently of the generated chemdf.
    """
    _, _, _, run = equilibrium_runs
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
def test_repeated_run_in_one_process_is_bitwise_identical(equilibrium_runs, tmp_path):
    """A second run with an identical, fresh config reproduces the first exactly.

    PROTEUS calls run_vulcan repeatedly in one process, so no state may leak between runs.
    """
    _, const_mix, _, first = equilibrium_runs
    second = _run(_config(tmp_path, ini_mix='const_mix', const_mix=const_mix, **EQ_SETTINGS))
    assert second['parameter']['count'] == first['parameter']['count']
    np.testing.assert_array_equal(second['variable']['y'], first['variable']['y'])


# ---------------------------------------------------------------------------
# Photochemistry: default HD 189733b column on a coarse grid
# ---------------------------------------------------------------------------


@pytest.fixture(scope='module')
def photochem_run(tmp_path_factory):
    """Default model on 20 layers; r_star and orbit_radius are set explicitly because the
    dilution pin of the top-of-atmosphere test depends on them."""
    tmp = tmp_path_factory.mktemp('vulcan_photo')
    cfg = _config(tmp, nz=20, r_star=0.805, orbit_radius=0.03142)
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
@pytest.mark.physics_invariant
@pytest.mark.reference_pinned
def test_top_of_atmosphere_flux_is_diluted_stellar_spectrum(photochem_run, validation_plot):
    """The direct beam at the top interface is the stellar surface flux times
    (R_star / a)^2, interpolated onto the wavelength bins (geometric dilution).

    Reference: the stellar input file read with numpy and interpolated with np.interp,
    independently of the scipy interpolation in build_atm.read_sflux; an input-to-output
    identity of the top boundary on a real run, not a radiative-transfer benchmark.
    """
    cfg, run = photochem_run
    var = run['variable']
    raw = np.genfromtxt(cfg.sflux_file, skip_header=1, names=['lambda', 'flux'])
    dilution = (cfg.r_star * r_sun / (au * cfg.orbit_radius)) ** 2
    expected = np.interp(var['bins'], raw['lambda'], raw['flux'], left=0, right=0) * dilution
    assert expected.max() > 0
    # pinned from the IAU 2015 nominal R_sun (6.957e8 m) and the IAU 2012 au, not phy_const
    assert dilution == pytest.approx(0.0141962, rel=1e-5)
    # linear interpolation of the same data: machine precision
    rtol = 1e-12
    np.testing.assert_allclose(var['sflux'][-1], expected, rtol=rtol, atol=0)
    # discrimination: R_star / a without the square is 8.4 times larger in every bin
    wrong = expected / np.sqrt(dilution)
    assert np.all(np.abs(var['sflux'][-1] - wrong) > rtol * wrong)

    def figure():
        import matplotlib.pyplot as plt

        rel = np.abs(var['sflux'][-1] - expected) / expected
        fig, (top, bottom) = plt.subplots(2, 1, sharex=True, figsize=(6.4, 5.6))
        top.plot(var['bins'][::25], expected[::25], label='input file')
        top.plot(var['bins'][::100], var['sflux'][-1][::100], 'o', ms=3, label='VULCAN')
        top.set(yscale='log', ylabel='Flux\n[erg/s/cm2/nm]')
        top.legend()
        bottom.plot(var['bins'][::50], np.maximum(rel[::50], 1e-18), '.', ms=2)
        bottom.axhline(rtol, color=plt.rcParams['text.color'], ls='--')
        bottom.set(
            yscale='log', ylim=(1e-18, 1e-10), yticks=[1e-16, 1e-12], xlabel='Wavelength [nm]'
        )
        bottom.set_ylabel('Rel. difference')
        for ax in (top, bottom):
            ax.minorticks_off()
        fig.tight_layout()
        return fig

    validation_plot(figure, 'toa_flux_diluted_stellar_spectrum')


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
    """A temperature <= 0 K is rejected when the atmosphere is built, before any rate
    coefficient is evaluated (so no NaN arithmetic, and no output)."""
    cfg = _config(tmp_path, atm_type='isothermal', Tiso=tiso, **BAD_RUN)
    with pytest.raises(ValueError, match='Temperature must be finite and positive'):
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
