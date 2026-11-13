"""Unit tests for vulcan.climate: the climate update, its regridding and the composition map."""

from __future__ import annotations

import numpy as np
import pytest

from vulcan.climate import (
    ClimateUpdate,
    interp_logp,
    map_composition,
    max_relative_change,
    regrid_kzz,
    regrid_temperature,
)

pytestmark = pytest.mark.unit

# A climate grid from 1e8 to 1e2 dyn cm-2 (100 bar to 0.1 mbar), ordered top down as AGNI is
P_CLIM = np.logspace(2, 8, 25)
# A coarse grid, two decades between points, on which the wrong interpolation errs clearly
P_COARSE = np.logspace(2, 8, 4)


def _update(tmp=None, kzz=None, p_kzz=None):
    tmp = np.linspace(300.0, 1500.0, P_CLIM.size) if tmp is None else tmp
    return ClimateUpdate(p=P_CLIM, tmp=tmp, p_kzz=p_kzz, kzz=kzz)


def test_regrid_temperature_is_exact_for_a_profile_linear_in_log_pressure():
    """T = a + b ln p is reproduced exactly between the grid points of the climate model.

    On a grid two decades apart, a linear interpolation in p instead of ln p misses this
    profile by about 75 K midway (in ln p) between two points.
    """
    a, b = 2000.0, -40.0
    update = ClimateUpdate(p=P_COARSE, tmp=a + b * np.log(P_COARSE))
    pco = np.logspace(2.5, 7.5, 11)  # VULCAN layers between the climate grid points
    exact = a + b * np.log(pco)
    np.testing.assert_allclose(regrid_temperature(update, pco), exact, rtol=1e-12)

    linear_in_p = np.interp(pco, P_COARSE, a + b * np.log(P_COARSE))
    assert np.max(np.abs(linear_in_p - exact)) > 50.0


def test_regrid_temperature_holds_end_values_outside_the_climate_grid():
    """Above and below the climate grid the temperature of its top and bottom layer holds."""
    tmp = np.linspace(200.0, 2000.0, P_CLIM.size)  # rises with pressure
    update = _update(tmp=tmp)
    out = regrid_temperature(update, np.array([1.0, 1e10]))
    assert out[0] == pytest.approx(200.0, rel=1e-12)  # above the top: top value
    assert out[1] == pytest.approx(2000.0, rel=1e-12)  # below the bottom: bottom value


def test_regrid_kzz_is_exact_for_a_power_law():
    """Kzz = K0 (p / p0)^-0.5, as above a convective zone, is reproduced exactly.

    On a grid two decades apart, interpolating Kzz itself instead of its logarithm gives
    0.55 K1 instead of 0.32 K1 midway between two points, about 70% too high.
    """
    kzz = 1e5 * (P_COARSE / 1e6) ** -0.5
    update = _update(kzz=kzz, p_kzz=P_COARSE)
    pico_inner = np.logspace(2.5, 7.5, 11)
    exact = 1e5 * (pico_inner / 1e6) ** -0.5
    np.testing.assert_allclose(regrid_kzz(update, pico_inner), exact, rtol=1e-10)

    linear = interp_logp(P_COARSE, kzz, pico_inner)
    assert np.max(np.abs(linear / exact - 1.0)) > 0.5


def test_regrid_kzz_with_zeros_interpolates_linearly_and_keeps_zeros():
    """A Kzz profile with zeros (no mixing) is interpolated linearly and stays at zero there."""
    kzz = np.where(P_CLIM > 1e5, 0.0, 1e6)
    update = _update(kzz=kzz, p_kzz=P_CLIM)
    out = regrid_kzz(update, np.array([1e3, 1e7]))
    assert out[0] == pytest.approx(1e6, rel=1e-12)
    assert out[1] == 0.0
    assert np.all(np.isfinite(out))


def test_regrid_kzz_refuses_an_update_without_kzz():
    """`Kzz_prof = 'climate'` needs the hook to return kzz; a missing profile is an error."""
    with pytest.raises(ValueError, match='no kzz'):
        regrid_kzz(_update(), np.array([1e3, 1e4]))


@pytest.mark.parametrize(
    'kwargs, match',
    [
        (dict(tmp=np.full(P_CLIM.size, -10.0)), 'positive'),
        (dict(tmp=np.full(P_CLIM.size - 1, 300.0)), 'at least 2 values'),
        (dict(kzz=np.ones(P_CLIM.size)), 'together'),
        (dict(kzz=-np.ones(P_CLIM.size), p_kzz=P_CLIM), 'negative'),
        (dict(kzz=np.full(P_CLIM.size, np.nan), p_kzz=P_CLIM), 'finite'),
    ],
    ids=[
        'negative-temperature',
        'length-mismatch',
        'kzz-without-pressure',
        'negative-kzz',
        'nan-kzz',
    ],
)
def test_climate_update_rejects_unphysical_profiles(kwargs, match):
    """A temperature <= 0 K, a negative or non-finite Kzz, or mismatched arrays are refused."""
    with pytest.raises(ValueError, match=match):
        _update(**kwargs)
    assert _update().tmp.min() == pytest.approx(300.0)  # the valid default is accepted


def test_map_composition_renormalises_and_keeps_gases_vulcan_lacks():
    """Shared gases take VULCAN's profile; others keep theirs; each level sums to 1.

    VULCAN has H2O, H2 and the radical OH (absent from the climate model); the climate
    model has H2O, H2 and Ar (absent from VULCAN). On each level the mapped mixing ratios
    are (x_H2O, x_H2, x_Ar) / (x_H2O + x_H2 + x_Ar), which drops OH from the sum.
    """
    p_src = np.array([1e6, 1e4, 1e2])
    species = ['H2O', 'H2', 'OH']
    ymix = np.array([[0.6, 0.39, 0.01], [0.3, 0.6, 0.1], [0.1, 0.6, 0.3]])
    p_tgt = np.array([1e6, 1e2])
    current = {
        'H2O': np.array([0.5, 0.5]),
        'H2': np.array([0.45, 0.45]),
        'Ar': np.array([0.05, 0.05]),
    }

    out = map_composition(p_src, species, ymix, p_tgt, current)

    assert set(out) == {'H2O', 'H2', 'Ar'}
    total = np.sum(list(out.values()), axis=0)
    np.testing.assert_allclose(total, 1.0, rtol=1e-14)
    # bottom level: (0.6, 0.39, 0.05) / 1.04
    np.testing.assert_allclose(
        [out['H2O'][0], out['H2'][0], out['Ar'][0]],
        np.array([0.6, 0.39, 0.05]) / 1.04,
        rtol=1e-12,
    )
    # top level: (0.1, 0.6, 0.05) / 0.75; without the division H2O would be 0.1
    assert out['H2O'][1] == pytest.approx(0.1 / 0.75, rel=1e-12)


def test_map_composition_refuses_an_empty_level():
    """A level with no gas of the climate model left is an error, not a division by zero."""
    with pytest.raises(ValueError, match='empty'):
        map_composition(
            np.array([1e6, 1e2]),
            ['OH'],
            np.array([[1.0], [1.0]]),
            np.array([1e6]),
            {'H2O': np.array([0.0])},
        )


def test_max_relative_change_is_the_largest_layer_change():
    """max |dT| / T over the layers; a 30 K change at 300 K (10%) beats 50 K at 1000 K (5%)."""
    t_old = np.array([300.0, 1000.0])
    t_new = np.array([330.0, 1050.0])
    assert max_relative_change(t_old, t_new) == pytest.approx(0.1, rel=1e-12)
    assert max_relative_change(t_old, t_old) == 0.0


def test_climate_update_copies_the_hook_arrays():
    """Changing the hook's arrays after it returned does not change the update.

    AGNI's arrays are Julia memory that its next solve overwrites; a view would follow it.
    """
    p, tmp = P_CLIM.copy(), np.full(P_CLIM.size, 500.0)
    update = ClimateUpdate(p=p, tmp=tmp)
    tmp[:] = 100.0
    p[:] = 1.0
    np.testing.assert_array_equal(update.tmp, 500.0)
    np.testing.assert_array_equal(update.p, P_CLIM)
