"""
---

Interface between VULCAN and a climate model that updates the temperature and the eddy
diffusion profile during an integration.

A climate hook is any callable `hook(pco, Tco, species, ymix)` that returns a
`ClimateUpdate`. VULCAN passes its layer-centre pressures `pco` [dyn cm-2] from the bottom of
the column up, its temperature `Tco` [K] on those layers, the list of species names and
their volume mixing ratios `ymix` [layer, species]. The hook returns the new temperature
profile, and optionally a new eddy diffusion profile, on its own pressure grid; VULCAN maps
them onto its grid with `regrid_temperature` and `regrid_kzz`. `vulcan.agni.AgniClimate` is
the hook for AGNI; a coupled framework can pass its own hook to `vulcan.run_vulcan`.

All quantities are in CGS units, as everywhere in VULCAN.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class ClimateUpdate:
    """Profiles returned by a climate hook, each on the hook's own pressure grid.

    The arrays are copied, so they do not change with the hook's own arrays.

    Attributes
    ----------
    p : np.ndarray
        Pressures of `tmp` [dyn cm-2], in any order.
    tmp : np.ndarray
        Temperature [K].
    p_kzz : np.ndarray or None
        Pressures of `kzz` [dyn cm-2], in any order. None when `kzz` is None.
    kzz : np.ndarray or None
        Eddy diffusion coefficient [cm2 s-1], or None to leave VULCAN's profile unchanged.
    """

    p: np.ndarray
    tmp: np.ndarray
    p_kzz: np.ndarray | None = None
    kzz: np.ndarray | None = None

    def __post_init__(self):
        self.p, self.tmp = _checked_pair(self.p, self.tmp, 'tmp')
        if not np.all(self.tmp > 0):
            raise ValueError('Climate temperature must be positive in every layer')
        if (self.kzz is None) != (self.p_kzz is None):
            raise ValueError('Climate kzz and p_kzz must be given together')
        if self.kzz is not None:
            self.p_kzz, self.kzz = _checked_pair(self.p_kzz, self.kzz, 'kzz')
            if not np.all(self.kzz >= 0):
                raise ValueError('Climate kzz must not be negative')


def _checked_pair(p, val, name: str) -> tuple[np.ndarray, np.ndarray]:
    """Return copies of a pressure and a value profile as float arrays, after checking them.

    The copies keep the update intact when the hook's own arrays change or are freed later,
    as the arrays of a Julia model are when it solves again or is deallocated.
    """
    p = np.array(p, dtype=float).ravel()
    val = np.array(val, dtype=float).ravel()
    if p.size != val.size or p.size < 2:
        raise ValueError(
            f'Climate {name} needs at least 2 values on as many pressures, '
            f'got {val.size} values on {p.size} pressures'
        )
    if not (np.all(np.isfinite(p)) and np.all(p > 0) and np.all(np.isfinite(val))):
        raise ValueError(f'Climate {name} profile must be finite, on positive pressures')
    return p, val


def interp_logp(p_src, val_src, p_tgt) -> np.ndarray:
    """Interpolate linearly in log pressure; outside the source range the end values hold.

    `p_src` and `p_tgt` share any one pressure unit; `p_src` may be in any order.
    """
    p_src = np.asarray(p_src, dtype=float)
    val_src = np.asarray(val_src, dtype=float)
    order = np.argsort(p_src)
    return np.interp(
        np.log(np.asarray(p_tgt, dtype=float)), np.log(p_src[order]), val_src[order]
    )


def regrid_temperature(update: ClimateUpdate, pco: np.ndarray) -> np.ndarray:
    """Temperature of a climate update on VULCAN's layer-centre pressures `pco` [K].

    Linear in log pressure. Above or below the climate grid the temperature of its top or
    bottom layer holds, so no value is extrapolated.
    """
    return interp_logp(update.p, update.tmp, np.asarray(pco, dtype=float))


def regrid_kzz(update: ClimateUpdate, pico_inner: np.ndarray) -> np.ndarray:
    """Eddy diffusion of a climate update on VULCAN's inner interfaces `pico[1:-1]` [cm2 s-1].

    Linear in log pressure of log Kzz when Kzz is positive everywhere, which reproduces a
    power law in pressure exactly; linear in Kzz otherwise. The end values hold outside the
    climate grid.
    """
    if update.kzz is None:
        raise ValueError('The climate update carries no kzz profile')
    p_tgt = np.asarray(pico_inner, dtype=float)
    if np.all(update.kzz > 0):
        return np.exp(interp_logp(update.p_kzz, np.log(update.kzz), p_tgt))
    return interp_logp(update.p_kzz, update.kzz, p_tgt)


def map_composition(
    p_src: np.ndarray,
    species: list[str],
    ymix: np.ndarray,
    p_tgt: np.ndarray,
    current: dict[str, np.ndarray],
) -> dict[str, np.ndarray]:
    """Map VULCAN mixing ratios onto the gases and pressures of a climate model.

    A gas of the climate model that VULCAN also has takes VULCAN's profile, interpolated
    linearly in log pressure (end values hold outside VULCAN's grid). A gas that VULCAN does
    not have keeps its profile from `current`. The mixing ratios are then divided by their
    sum on each level, so that they sum to 1 again: VULCAN species absent from the climate
    model, such as radicals and condensates, are left out of that sum.

    Parameters
    ----------
    p_src : np.ndarray
        VULCAN pressures of `ymix` [any unit, the same as `p_tgt`].
    species : list of str
        VULCAN species names, the columns of `ymix`.
    ymix : np.ndarray
        VULCAN volume mixing ratios [layer, species].
    p_tgt : np.ndarray
        Pressures of the climate model's layers, in the same unit as `p_src`.
    current : dict
        The climate model's present mixing ratio profile of each of its gases, on `p_tgt`.

    Returns
    -------
    dict
        The new mixing ratio profile of each gas of the climate model, on `p_tgt`.
    """
    p_src = np.asarray(p_src, dtype=float)
    p_tgt = np.asarray(p_tgt, dtype=float)
    ymix = np.asarray(ymix, dtype=float)
    out = {}
    for gas, prof in current.items():
        if gas in species:
            col = np.clip(ymix[:, species.index(gas)], 0.0, None)
            out[gas] = interp_logp(p_src, col, p_tgt)
        else:
            out[gas] = np.asarray(prof, dtype=float).copy()
    total = np.sum(list(out.values()), axis=0)
    if not np.all(total > 0):
        raise ValueError('The mapped composition is empty on at least one level')
    return {gas: prof / total for gas, prof in out.items()}


def max_relative_change(t_old: np.ndarray, t_new: np.ndarray) -> float:
    """Largest relative temperature change over the layers, max |T_new - T_old| / T_old."""
    t_old = np.asarray(t_old, dtype=float)
    return float(np.max(np.abs(np.asarray(t_new, dtype=float) - t_old) / t_old))
