"""Unit tests for vulcan.plot: composition parsing, colours, labels and the plot routines."""

from __future__ import annotations

import logging
import re
from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
import pytest

from vulcan import plot

SPECIES = ['H2', 'H', 'H2O', 'CH4', 'CO', 'He']
NZ = 6


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close('all')


@pytest.fixture
def keep_figures(monkeypatch):
    """Stop the plot routines closing their figures, so the plotted data can be inspected."""
    monkeypatch.setattr(plot.plt, 'close', lambda *args: None)


@pytest.fixture
def model(tmp_path):
    """Minimal cfg/var/atm/para for a 6-layer column with an asymmetric composition."""
    pco = np.logspace(8, 2, NZ)  # dyn cm-2, bottom to top
    pico = np.logspace(8.2, 1.8, NZ + 1)
    n_0 = pco / (1.380649e-16 * np.linspace(1500.0, 900.0, NZ))
    rng = np.random.default_rng(3)
    ymix = rng.uniform(0.05, 1.0, (NZ, len(SPECIES)))
    ymix /= ymix.sum(axis=1, keepdims=True)
    zco = np.linspace(0.0, 6.0e7, NZ + 1)  # cm

    cfg = SimpleNamespace(
        plot_spec=['H2O', 'CH4', 'CO'],
        plot_height=False,
        use_condense=False,
        condense_sp=[],
        P_b=1e8,
        P_t=1e2,
        plot_dir=str(tmp_path) + '/',
        plot_fmt='png',
        plot_dpi=50,
    )
    var = SimpleNamespace(
        ymix=ymix,
        y=ymix * n_0[:, None],
        t=1.5e7,
        longdydt=3e-5,
        t_time=np.logspace(0, 7, 5),
        y_time=np.stack([ymix * n_0[:, None] * (1 + 0.1 * i) for i in range(5)]),
    )
    atm = SimpleNamespace(
        pco=pco,
        pico=pico,
        n_0=n_0,
        Tco=np.linspace(1500.0, 900.0, NZ),
        zco=zco,
        zmco=0.5 * (zco[1:] + zco[:-1]),
        Kzz=np.logspace(8, 10, NZ - 1),
        sat_mix={'H2O': np.full(NZ, 1e-4)},
    )
    para = SimpleNamespace(count=42, pic_count=0)
    return cfg, var, atm, para


# ---------------------------------------------------------------------------
# Composition, labels and colours
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize(
    'mol, expected',
    [
        ('H2O', {'H': 2, 'O': 1}),
        ('CH3COOH', {'C': 2, 'H': 4, 'O': 2}),  # repeated elements are summed
        ('C10H22', {'C': 10, 'H': 22}),  # multi-digit counts
        ('NaCl', {'Na': 1, 'Cl': 1}),  # two-letter elements
        ('He', {'He': 1}),
        ('CH2_1', {'C': 1, 'H': 2}),  # VULCAN excited-state suffix is ignored
        ('S8_l_s', {'S': 8}),  # VULCAN condensate suffix is ignored
    ],
)
def test_mol_to_ele_counts_atoms(mol, expected):
    """Element counts follow the chemical formula."""
    assert plot.mol_to_ele(mol) == expected

@pytest.mark.unit
@pytest.mark.parametrize(
    'gas, expected',
    [('H2O', 'H$_2$O'), ('C2H6', 'C$_2$H$_6$'), ('He', 'He')],
)
def test_latexify_subscripts_digits(gas, expected):
    """Digits in a formula become LaTeX subscripts; letters are unchanged."""
    assert plot.latexify(gas) == expected


@pytest.mark.unit
def test_get_colour_uses_preset_table():
    """Gases with a preset colour use it."""
    assert plot.get_colour('H2O') == '#027FB1'
    assert plot.get_colour('CO2') == '#D24901'


@pytest.mark.unit
@pytest.mark.parametrize('a, b', [('O2', 'O3'), ('NO2', 'N2O4'), ('C2H6', 'CH3')])
def test_generated_colour_depends_on_element_ratios_only(a, b):
    """Generated colours mix element colours by composition, normalised to the brightest
    channel, so formulae with the same element ratios get the same colour."""
    assert plot.get_colour(a) == plot.get_colour(b)


@pytest.mark.unit
@pytest.mark.parametrize('gas', ['C2H6', 'HCN', 'NH2', 'OCS', 'H2SO4', 'HO2', 'C6H6'])
def test_generated_colour_is_valid_and_not_near_white(gas):
    """Generated colours are valid hex strings, and dimmed if close to white (sum <= 3*220)."""
    colour = plot.get_colour(gas)
    assert re.fullmatch(r'#[0-9a-f]{6}', colour)
    rgb = [int(colour[i : i + 2], 16) for i in (1, 3, 5)]
    assert sum(rgb) <= 3 * 220
    assert max(rgb) > 0


@pytest.mark.unit
@pytest.mark.parametrize('gas', ['ZrO2', 'h2o'])
def test_unknown_composition_falls_back(gas, caplog):
    """An unknown element or an unparsable name gets the fallback colour, with a warning."""
    with caplog.at_level(logging.WARNING):
        assert plot.get_colour(gas) == '#ff00ff'
    assert 'fallback' in caplog.text


# ---------------------------------------------------------------------------
# Plot routines
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_plot_update_plots_mixing_ratio_against_pressure_in_bar(model, tmp_path, keep_figures):
    """Profiles are mixing ratio versus pressure in bar (1 bar = 1e6 dyn cm-2).

    Without an explicit path the plot is saved as _recent and copied to a numbered frame.
    """
    cfg, var, atm, para = model
    plot.plot_update(cfg, SPECIES, var, atm, para)

    ax = plt.gcf().axes[0]
    lines = {ln.get_label(): ln for ln in ax.get_lines()}
    for sp in cfg.plot_spec:
        ln = lines[plot.latexify(sp)]
        np.testing.assert_allclose(ln.get_xdata(), var.ymix[:, SPECIES.index(sp)], rtol=0)
        np.testing.assert_allclose(ln.get_ydata(), atm.pco / 1e6, rtol=1e-15)
    assert ax.get_yscale() == 'log'
    assert ax.yaxis_inverted()  # pressure increases downward

    assert (tmp_path / '_recent.png').stat().st_size > 0
    assert (tmp_path / '00000.png').stat().st_size > 0
    assert para.pic_count == 1


@pytest.mark.unit
def test_plot_update_height_axis_and_saturation_curve(model, tmp_path, keep_figures):
    """With plot_height, profiles use layer-centre heights in km and condensates show a
    dashed saturation curve."""
    cfg, var, atm, para = model
    cfg.plot_height = True
    cfg.use_condense = True
    cfg.condense_sp = ['H2O']
    fpath = tmp_path / 'mix.png'
    plot.plot_update(cfg, SPECIES, var, atm, para, fpath=str(fpath))

    ax = plt.gcf().axes[0]
    lines = {ln.get_label(): ln for ln in ax.get_lines()}
    np.testing.assert_allclose(lines['CH$_4$'].get_ydata(), atm.zmco / 1e5, rtol=1e-15)
    sat = lines['H$_2$O sat']
    assert sat.get_linestyle() == '--'
    np.testing.assert_allclose(sat.get_xdata(), atm.sat_mix['H2O'])
    assert fpath.stat().st_size > 0
    assert not (tmp_path / '00000.png').exists()  # no frame copy for an explicit path


@pytest.mark.unit
def test_plot_evo_plots_top_layer_mixing_ratio_history(model, tmp_path, keep_figures):
    """The evolution plot shows y/n_0 in the selected layer (default: top) versus time."""
    cfg, var, atm, para = model
    plot.plot_evo(cfg, SPECIES, var, atm)

    ax = plt.gcf().axes[0]
    ln = {ln.get_label(): ln for ln in ax.get_lines()}['CO']
    isp = SPECIES.index('CO')
    np.testing.assert_allclose(ln.get_xdata(), var.t_time)
    np.testing.assert_allclose(ln.get_ydata(), var.y_time[:, -1, isp] / atm.n_0[-1], rtol=1e-14)
    assert (tmp_path / 'evo.png').stat().st_size > 0


@pytest.mark.unit
def test_plot_evo_inter_writes_figure(model, tmp_path):
    """The interrupted-run evolution plot is written to plot_dir."""
    cfg, var, atm, para = model
    plot.plot_evo_inter(cfg, SPECIES, var, atm)
    assert (tmp_path / 'evo.png').stat().st_size > 0


@pytest.mark.unit
def test_plot_TP_places_Kzz_on_interfaces(model, tmp_path, keep_figures):
    """Temperature is plotted at layer centres and Kzz at the interior interfaces."""
    cfg, var, atm, para = model
    plot.plot_TP(cfg, atm)

    ax1, ax2 = plt.gcf().axes
    np.testing.assert_allclose(ax1.get_lines()[0].get_ydata(), atm.pco / 1e6)
    np.testing.assert_allclose(ax2.get_lines()[0].get_xdata(), atm.Kzz)
    np.testing.assert_allclose(ax2.get_lines()[0].get_ydata(), atm.pico[1:-1] / 1e6)
    assert (tmp_path / 'TPK_initial.png').stat().st_size > 0


@pytest.mark.unit
@pytest.mark.parametrize('routine', ['plot_update', 'plot_evo', 'plot_evo_inter', 'plot_TP'])
def test_plot_routines_leave_no_open_figures(model, routine):
    """Each plot routine closes its figures, so repeated live plotting does not leak memory."""
    cfg, var, atm, para = model
    plt.close('all')
    if routine == 'plot_TP':
        plot.plot_TP(cfg, atm)
    elif routine == 'plot_update':
        plot.plot_update(cfg, SPECIES, var, atm, para)
    else:
        getattr(plot, routine)(cfg, SPECIES, var, atm)
    assert plt.get_fignums() == []
