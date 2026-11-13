"""Unit tests for vulcan.vulcan: the checks on the climate coupling settings."""

from __future__ import annotations

import pytest

from vulcan import Config
from vulcan.climate import ClimateUpdate
from vulcan.vulcan import check_climate_config

pytestmark = pytest.mark.unit


def _hook(pco, tco, species, ymix):
    return ClimateUpdate(p=pco, tmp=tco)


def test_check_climate_config_accepts_the_defaults_and_a_coupled_setup():
    """The default config (no coupling) and Kzz from a climate hook are both valid."""
    cfg = Config()
    check_climate_config(cfg)
    cfg.Kzz_prof = 'climate'
    check_climate_config(cfg, climate=_hook)
    cfg.agni_call_frq = 50
    check_climate_config(cfg)
    assert cfg.T_cross_sp == []


@pytest.mark.parametrize(
    'settings, with_hook, match',
    [
        (dict(Kzz_prof='climate'), False, 'needs a climate hook'),
        (dict(T_cross_sp=['H2O']), True, 'T_cross_sp'),
        (dict(fix_species=['H2O_l_s']), True, 'fix_species'),
        (dict(agni_T_rtol=0.0), False, 'agni_T_rtol'),
        (dict(agni_max_calls=0), False, 'agni_max_calls'),
    ],
    ids=[
        'kzz-without-climate',
        'temperature-dependent-cross-sections',
        'fixed-species',
        'zero-rtol',
        'no-calls',
    ],
)
def test_check_climate_config_refuses_settings_that_cannot_work(settings, with_hook, match):
    """Each unusable setting is refused with a message that names it."""
    cfg = Config()
    for key, val in settings.items():
        setattr(cfg, key, val)
    with pytest.raises(ValueError, match=match):
        check_climate_config(cfg, climate=_hook if with_hook else None)
    assert cfg.agni_call_frq == 0  # the check changes no setting
