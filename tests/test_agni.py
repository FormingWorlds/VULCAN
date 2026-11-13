"""Tests for vulcan.agni: unit conversion and composition hand-off, and a smoke run with AGNI.

The unit tests use a stand-in for the AGNI atmosphere with the fields the functions read, so
they run without Julia. The integration test sets up a real AGNI atmosphere and is skipped
when juliacall or the AGNI folder is missing.
"""

from __future__ import annotations

import importlib.util
import os
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from vulcan import Config, agni, build_atm, store
from vulcan.chem_funs import spec_list
from vulcan.climate import regrid_kzz, regrid_temperature


def _fake_atmos():
    """AGNI-like atmosphere: 4 layers from 1e4 to 1e1 Pa, Kzz on the 5 level edges (SI)."""
    p = np.array([1e1, 1e2, 1e3, 1e4])  # Pa, top down as in AGNI
    pl = np.array([3e0, 3e1, 3e2, 3e3, 3e4])
    return SimpleNamespace(
        p=p,
        pl=pl,
        tmp=np.array([200.0, 300.0, 400.0, 500.0]),
        Kzz=np.array([1.0, 2.0, 3.0, 4.0, 5.0]),  # m2 s-1
        gas_names=['H2O', 'H2'],
        gas_vmr={'H2O': np.full(4, 0.5), 'H2': np.full(4, 0.5)},
        gas_ovmr={'H2O': np.full(4, 0.5), 'H2': np.full(4, 0.5)},
    )


@pytest.mark.unit
def test_climate_update_converts_agni_si_units_to_cgs():
    """Pressures go from Pa to dyn cm-2 (x10) and Kzz from m2 s-1 to cm2 s-1 (x1e4).

    The inverse factors (x0.1, x1e-4) would put the profiles 2 and 8 orders of magnitude off.
    """
    atmos = _fake_atmos()
    update = agni.climate_update(atmos)
    np.testing.assert_allclose(update.p, atmos.p * 10.0, rtol=1e-15)
    np.testing.assert_allclose(update.p_kzz, atmos.pl * 10.0, rtol=1e-15)
    np.testing.assert_allclose(update.kzz, atmos.Kzz * 1e4, rtol=1e-15)
    np.testing.assert_array_equal(update.tmp, atmos.tmp)
    assert update.kzz[0] / (atmos.Kzz[0] * 1e-4) == pytest.approx(1e8)


@pytest.mark.unit
def test_set_composition_maps_vulcan_profiles_onto_agni_layers():
    """AGNI layers at 1e3 and 1e1 Pa take VULCAN's mixing ratios at 1e4 and 1e2 dyn cm-2.

    VULCAN's pressures are in dyn cm-2 (10 per Pa); reading them as Pa would take the
    composition from a level ten times deeper, where H2O is 0.8 instead of 0.2.
    """
    atmos = _fake_atmos()
    pco = np.array([1e5, 1e4, 1e3, 1e2])  # dyn cm-2, bottom up as in VULCAN
    species = ['H2O', 'H2']
    ymix = np.array([[0.8, 0.2], [0.6, 0.4], [0.4, 0.6], [0.2, 0.8]])

    agni.set_composition(atmos, pco, species, ymix)

    # AGNI layer at 1e1 Pa = 1e2 dyn cm-2 is VULCAN's top layer; at 1e3 Pa = 1e4 dyn cm-2
    assert atmos.gas_vmr['H2O'][0] == pytest.approx(0.2, rel=1e-12)
    assert atmos.gas_vmr['H2O'][2] == pytest.approx(0.6, rel=1e-12)
    np.testing.assert_allclose(atmos.gas_vmr['H2O'] + atmos.gas_vmr['H2'], 1.0, rtol=1e-14)
    np.testing.assert_array_equal(atmos.gas_ovmr['H2O'], atmos.gas_vmr['H2O'])


@pytest.mark.unit
def test_agni_dir_defaults_to_the_vulcan_folder():
    """An empty `agni_dir` means VULCAN/AGNI; a set one is used as given."""
    cfg = Config()
    assert agni.agni_dir(cfg) == agni.paths.AGNI_DIR
    cfg.agni_dir = '/somewhere/AGNI'
    assert agni.agni_dir(cfg) == '/somewhere/AGNI'


@pytest.mark.unit
def test_init_refuses_rayleigh_with_greygas():
    """AGNI cannot combine grey-gas radiation with Rayleigh scattering; refuse it up front."""
    cfg = Config()
    cfg.spectral_file, cfg.use_rayleigh = 'greygas', True
    with pytest.raises(ValueError, match='Rayleigh'):
        agni.init_agni_atmos(cfg, None, None)
    assert cfg.use_rayleigh  # the config is left as it was


def _agni_available() -> bool:
    """juliacall is installed and AGNI is in VULCAN/AGNI; checked without starting Julia."""
    if importlib.util.find_spec('juliacall') is None:
        return False
    return os.path.isfile(os.path.join(agni.paths.AGNI_DIR, 'agni.jl'))


def _julia_starts() -> bool:
    """Julia starts through juliacall, tried in a subprocess.

    A Julia runtime that fails to load ends the whole Python process, which no skip in
    this process could catch.
    """
    cmd = [sys.executable, '-c', 'from juliacall import Main; Main.seval("1")']
    try:
        return subprocess.run(cmd, capture_output=True, timeout=600).returncode == 0
    except subprocess.TimeoutExpired:
        return False


@pytest.mark.integration
@pytest.mark.skipif(not _agni_available(), reason='needs juliacall and AGNI in VULCAN/AGNI')
def test_agni_climate_keeps_the_temperature_and_returns_a_positive_kzz(tmp_path):
    """With `solve_rce = False` AGNI returns VULCAN's temperature and a finite, positive Kzz.

    This sets up a grey-gas AGNI atmosphere on a 10-layer VULCAN grid, so it checks that the
    calls to AGNI match its present interface. The temperature is linear in ln p, which the
    round trip through AGNI's grid keeps exactly inside the grid; reading VULCAN's dyn cm-2
    as Pa would shift it by 40 K x ln 10 = 92 K.
    """
    if not _julia_starts():
        pytest.skip('Julia does not start through juliacall in this environment')
    cfg = Config()
    cfg.nz = 10
    cfg.ini_mix = 'const_mix'
    cfg.output_dir = str(tmp_path) + '/'
    cfg.plot_dir = str(tmp_path / 'plot') + '/'
    cfg.plot_TP = False
    cfg.use_live_plot = False
    cfg.use_photo = False
    cfg.agni_nlev = 20
    cfg.save_type = ['pickle']
    os.makedirs(cfg.output_dir, exist_ok=True)

    var, atm = store.Variables(cfg), store.AtmData(cfg)
    make_atm = build_atm.Atm(cfg)
    atm = make_atm.f_pico(atm)
    atm = make_atm.load_TPK(atm)
    atm.Tco = 1500.0 + 40.0 * np.log(atm.pco / atm.pco[0])  # 1500 K at the bottom, cooler above
    var = build_atm.InitialAbun(cfg).ini_y(var, atm)
    var.ymix = var.y / var.y.sum(axis=1, keepdims=True)

    climate = agni.AgniClimate(cfg, atm, var)
    try:
        update = climate(atm.pco.copy(), atm.Tco.copy(), list(spec_list), var.ymix.copy())
        # the energy solver accepts every keyword VULCAN passes (no MethodError) and leaves a
        # physical profile, whether or not it converges from this start
        solved = agni._solve_energy(climate.atmos, cfg)
        t_solved = np.asarray(climate.atmos.tmp, dtype=float)
    finally:
        climate.close()
    assert solved in (True, False)
    assert np.all(np.isfinite(t_solved))
    assert np.all(t_solved > 0)

    # inside AGNI's layer centres the profile comes back exactly; the end layers fall between
    # an AGNI layer and AGNI's held end value
    t_back = regrid_temperature(update, atm.pco)
    np.testing.assert_allclose(t_back[1:-1], atm.Tco[1:-1], rtol=1e-9)
    kzz = regrid_kzz(update, atm.pico[1:-1])
    assert np.all(np.isfinite(kzz))
    assert np.all(kzz > 0)
