"""
---

Climate hook that couples VULCAN to the AGNI radiative-convective model.

`AgniClimate` sets up an AGNI atmosphere on VULCAN's pressure range and, at each climate
update, hands it VULCAN's mixing ratios, solves for the temperature profile (or, with
`solve_rce = False`, keeps VULCAN's temperature) and returns the temperature and AGNI's eddy
diffusion profile as a `vulcan.climate.ClimateUpdate`. `set_composition` and
`climate_update` work on any AGNI atmosphere, so a coupled framework that builds its own
AGNI atmosphere can use them in its own hook.

AGNI works in SI units and VULCAN in CGS: pressures are converted with `PA_TO_DYN` and
diffusion coefficients with `M2_TO_CM2` at this boundary.

Imports
---

- [`vulcan.climate`](climate.md): the climate-hook interface.
- [`vulcan.config`](config.md): `Config` class for handling configuration settings.
- [`vulcan.paths`](paths.md): Paths to various files and directories used in VULCAN.
"""

from __future__ import annotations

import logging
import os

import numpy as np
from scipy.integrate import trapezoid

from . import paths
from .climate import ClimateUpdate, interp_logp, map_composition
from .config import Config
from .phy_const import au, r_sun, spacer

log = logging.getLogger('fwl.' + __name__)

AGNI_LOGFILE_NAME = 'agni_recent.log'
PA_TO_DYN = 10.0  # 1 Pa = 10 dyn cm-2
M2_TO_CM2 = 1.0e4  # 1 m2 s-1 = 1e4 cm2 s-1


def _julia():
    """The Julia `Main` module, imported only when AGNI is used."""
    from juliacall import Main as jl

    return jl


def agni_dir(vulcan_cfg: Config) -> str:
    """AGNI folder: `agni_dir` from the config, or `VULCAN/AGNI` when that is empty."""
    return vulcan_cfg.agni_dir or paths.AGNI_DIR


def activate_julia(vulcan_cfg: Config):
    """
    Load AGNI into Julia, unless a caller (such as PROTEUS) has loaded it already.

    When AGNI is already loaded, its Julia environment and logger are left as they are.
    Otherwise the Julia environment in `agni_dir` is activated, headless plotting is set
    up, AGNI is imported and its log goes to `agni_recent.log` in `Config.output_dir`.

    Parameters
    ----------
    vulcan_cfg : Config
        VULCAN configuration; `agni_dir` and `output_dir` are used.

    Returns
    -------
    juliacall module
        The Julia `Main` module.

    Raises
    ------
    FileNotFoundError
        If AGNI is not loaded yet and `agni.jl` is not in `agni_dir`.
    """
    jl = _julia()
    if bool(jl.seval('isdefined(Main, :AGNI)')):
        log.info('AGNI is already loaded; using it as it is')
        return jl

    root = agni_dir(vulcan_cfg)
    if not os.path.isfile(os.path.join(root, 'agni.jl')):
        raise FileNotFoundError(os.path.join(root, 'agni.jl'))

    log.info(f'Activating Julia environment {root}')
    jl.seval('using Pkg')
    jl.Pkg.activate(root)
    jl.seval('ENV["GKSwstype"] = "100"')
    jl.seval('using Plots')
    jl.seval('default(label=nothing, dpi=250)')
    jl.seval('import AGNI')

    # AGNI keeps this log file open for the whole run
    logpath = os.path.join(vulcan_cfg.output_dir, AGNI_LOGFILE_NAME)
    jl.AGNI.setup_logging(logpath, 2)
    log.debug("AGNI will log to '%s'" % logpath)
    return jl


def _instellation(vulcan_cfg: Config) -> float:
    """Bolometric instellation at the planet [W m-2], from the spectrum at the stellar surface."""
    sflux = np.loadtxt(vulcan_cfg.sflux_file).T  # nm, erg s-1 cm-2 nm-1
    flux = trapezoid(sflux[1], x=sflux[0])  # erg s-1 cm-2 at the stellar surface
    flux *= (vulcan_cfg.r_star * r_sun) ** 2 / (vulcan_cfg.orbit_radius * au) ** 2
    return float(flux * 1.0e-3)  # erg s-1 cm-2 to W m-2


def init_agni_atmos(vulcan_cfg: Config, atm, var):
    """
    Set up and allocate an AGNI atmosphere on VULCAN's pressure range.

    The composition starts from VULCAN's column-median mixing ratios, leaving out species
    with '_' in their name (condensates and excited states), and the temperature from
    VULCAN's profile (`set_temperature`). The instellation comes from `_instellation`
    [W m-2]. The boundary pressures `atm.pico[0]` and `atm.pico[-1]` are converted from
    dyn cm-2 to bar, gravity and planet radius from CGS to SI, and `Config.sl_angle` from
    radians to degrees.

    A `runtime.sf` spectral file in `Config.output_dir` is used unchanged when it exists;
    otherwise AGNI makes it from `Config.spectral_file` in `agni_dir` and the stellar
    spectrum. FastChem files go to `agni_fastchem` in `Config.output_dir`.

    Parameters
    ----------
    vulcan_cfg : Config
        VULCAN configuration.
    atm : store.AtmData
        Atmosphere data; reads `pico` [dyn cm-2], `pco` [dyn cm-2] and `Tco` [K].
    var : store.Variables
        Model variables; reads `ymix`.

    Returns
    -------
    AGNI.atmosphere.Atmos_t
        The allocated AGNI atmosphere.

    Raises
    ------
    ValueError
        If `spectral_file = 'greygas'` is combined with `use_rayleigh`.
    RuntimeError
        If AGNI fails to set up or allocate the atmosphere.
    """
    if vulcan_cfg.spectral_file == 'greygas' and vulcan_cfg.use_rayleigh:
        raise ValueError('AGNI does not support Rayleigh scattering with greygas')
    from .chem_funs import spec_list

    jl = _julia()
    log.debug('New AGNI atmosphere')
    atmos = jl.AGNI.atmosphere.Atmos_t()

    sflux_path = vulcan_cfg.sflux_file
    try_spfile = os.path.join(vulcan_cfg.output_dir, 'runtime.sf')
    if os.path.exists(try_spfile):
        input_sf, input_star = try_spfile, ''  # AGNI made it in an earlier run
    elif vulcan_cfg.spectral_file == 'greygas':
        input_sf, input_star = 'greygas', sflux_path
    else:
        input_sf = os.path.join(agni_dir(vulcan_cfg), vulcan_cfg.spectral_file)
        input_star = sflux_path

    # condensates and excited states carry an underscore
    vol_dict = {
        g: float(np.median(var.ymix[:, i])) for i, g in enumerate(spec_list) if '_' not in g
    }

    fastchem_work = os.path.join(vulcan_cfg.output_dir, 'agni_fastchem')
    os.makedirs(fastchem_work, exist_ok=True)

    succ = jl.AGNI.atmosphere.setup_b(
        atmos,
        agni_dir(vulcan_cfg),
        vulcan_cfg.output_dir,
        input_sf,
        _instellation(vulcan_cfg),
        float(vulcan_cfg.f_diurnal),
        0.0,
        float(vulcan_cfg.sl_angle * 180.0 / np.pi),
        float(vulcan_cfg.Tsurf_guess),
        float(vulcan_cfg.gs * 0.01),  # cm s-2 to m s-2
        float(vulcan_cfg.Rp * 0.01),  # cm to m
        int(vulcan_cfg.agni_nlev),
        float(atm.pico[0] * 1e-6),  # dyn cm-2 to bar
        float(atm.pico[-1] * 1e-6),
        vol_dict,
        '',
        flag_rayleigh=bool(vulcan_cfg.use_rayleigh),
        flag_cloud=False,
        overlap_method='ee',
        albedo_s=float(vulcan_cfg.surf_albedo),
        surface_material='greybody',
        condensates=[],
        use_all_gases=False,
        fastchem_work=fastchem_work,
        real_gas=False,
        skin_d=0.01,
        skin_k=2.0,
        tmp_magma=float(vulcan_cfg.Tsurf_guess),
        tmp_floor=50.0,
        flux_int=float(vulcan_cfg.agni_flux_int),
    )
    if not succ:
        raise RuntimeError('Could not initialise AGNI atmos object')
    if not jl.AGNI.atmosphere.allocate_b(atmos, input_star):
        raise RuntimeError('Could not allocate AGNI atmos object')

    set_temperature(atmos, atm.pco, atm.Tco)
    return atmos


def deallocate_atmos(atmos):
    """
    Deallocate an AGNI atmosphere.

    Parameters
    ----------
    atmos : AGNI.atmosphere.Atmos_t
        The atmosphere to deallocate.
    """
    _julia().AGNI.atmosphere.deallocate_b(atmos)


def set_temperature(atmos, pco, tco):
    """
    Set the AGNI temperature profile from VULCAN's.

    The temperature is interpolated in log pressure onto AGNI's layers and levels (see
    `vulcan.climate.interp_logp`); the surface temperature is that of the last AGNI
    level. AGNI's layer properties are then recomputed.

    Parameters
    ----------
    atmos : AGNI.atmosphere.Atmos_t
        AGNI atmosphere, updated in place.
    pco : np.ndarray
        VULCAN layer-centre pressures [dyn cm-2].
    tco : np.ndarray
        VULCAN temperatures on `pco` [K].
    """
    p_layers = np.asarray(atmos.p, dtype=float) * PA_TO_DYN
    p_levels = np.asarray(atmos.pl, dtype=float) * PA_TO_DYN
    atmos.tmp[:] = interp_logp(pco, tco, p_layers)
    atmos.tmpl[:] = interp_logp(pco, tco, p_levels)
    atmos.tmp_surf = float(atmos.tmpl[-1])
    _julia().AGNI.atmosphere.calc_layer_props_b(atmos)


def set_composition(atmos, pco, species: list[str], ymix):
    """
    Give an AGNI atmosphere VULCAN's mixing ratio profiles (see `map_composition`).

    A gas that VULCAN lacks keeps its profile from before condensation (`gas_ovmr`).

    Parameters
    ----------
    atmos : AGNI.atmosphere.Atmos_t
        AGNI atmosphere; `gas_vmr` and `gas_ovmr` are updated in place.
    pco : np.ndarray
        VULCAN layer-centre pressures [dyn cm-2].
    species : list of str
        VULCAN species names, the columns of `ymix`.
    ymix : np.ndarray
        VULCAN volume mixing ratios, shape (nz, ni) [1].
    """
    # gas_ovmr, not the rained-out gas_vmr, so a gas VULCAN lacks keeps its full inventory
    current = {str(g): np.asarray(atmos.gas_ovmr[g], dtype=float) for g in atmos.gas_names}
    p_layers = np.asarray(atmos.p, dtype=float) * PA_TO_DYN
    for gas, prof in map_composition(pco, species, ymix, p_layers, current).items():
        atmos.gas_vmr[gas][:] = prof
        atmos.gas_ovmr[gas][:] = prof


def climate_update(atmos) -> ClimateUpdate:
    """
    The temperature and eddy diffusion of an AGNI atmosphere, in VULCAN's CGS units.

    Temperature is on AGNI's layer centres and Kzz on its level edges. Pressures are
    converted from Pa to dyn cm-2 and Kzz from m2 s-1 to cm2 s-1.

    Parameters
    ----------
    atmos : AGNI.atmosphere.Atmos_t
        AGNI atmosphere.

    Returns
    -------
    ClimateUpdate
        `p` [dyn cm-2] and `tmp` [K] on the layers, `p_kzz` [dyn cm-2] and `kzz`
        [cm2 s-1] on the levels.
    """
    return ClimateUpdate(
        p=np.asarray(atmos.p, dtype=float) * PA_TO_DYN,
        tmp=np.asarray(atmos.tmp, dtype=float),
        p_kzz=np.asarray(atmos.pl, dtype=float) * PA_TO_DYN,
        kzz=np.asarray(atmos.Kzz, dtype=float) * M2_TO_CM2,
    )


def _solve_energy(atmos, vulcan_cfg: Config) -> bool:
    """Solve for the energy-conserving temperature profile, with a second, wider attempt."""
    jl = _julia()
    modplot = 1 if vulcan_cfg.use_live_plot else 0
    for attempt in (1, 2):
        log.info('Attempt %d' % attempt)
        wide = attempt == 2
        success = bool(
            jl.AGNI.solver.solve_energy_b(
                atmos,
                sol_type=int(vulcan_cfg.agni_sol_type),
                method=1,
                chem=False,
                conduct=False,
                convect=True,
                sens_heat=True,
                latent=False,
                rainout=True,
                max_steps=70,
                max_runtime=900.0,
                conv_atol=float(vulcan_cfg.agni_atol),
                conv_rtol=float(vulcan_cfg.agni_rtol),
                ls_method=1 if wide else 2,
                ls_increase=1.1 if wide else 1.04,
                dx_max=600.0 if wide else 300.0,
                easy_start=not bool(atmos.is_solved),
                perturb_all=True,
                save_frames=False,
                modplot=modplot,
            )
        )
        if success:
            log.info('Attempt %d succeeded' % attempt)
            return True
        log.warning('Attempt %d failed' % attempt)
    log.error('Maximum attempts when executing AGNI')
    return False


class AgniClimate:
    """Climate hook that solves the temperature profile with AGNI (see `vulcan.climate`).

    With `solve_rce = True` AGNI solves for radiative-convective equilibrium with the
    surface condition `agni_sol_type`. With `solve_rce = False` VULCAN's temperature is
    returned unchanged, on VULCAN's own grid, and AGNI only computes the convective fluxes
    and the eddy diffusion profile from it.
    """

    def __init__(self, vulcan_cfg: Config, atm, var):
        """
        Parameters
        ----------
        vulcan_cfg : Config
            VULCAN configuration.
        atm : store.AtmData
            Atmosphere data, passed to `init_agni_atmos`.
        var : store.Variables
            Model variables, passed to `init_agni_atmos`.
        """
        self.cfg = vulcan_cfg
        activate_julia(vulcan_cfg)
        self.atmos = init_agni_atmos(vulcan_cfg, atm, var)

    def __call__(self, pco, tco, species, ymix) -> ClimateUpdate:
        """
        Run AGNI with VULCAN's present composition and return the new profiles.

        Parameters
        ----------
        pco : np.ndarray
            VULCAN layer-centre pressures [dyn cm-2].
        tco : np.ndarray
            VULCAN temperatures on `pco` [K].
        species : list of str
            VULCAN species names, the columns of `ymix`.
        ymix : np.ndarray
            VULCAN volume mixing ratios, shape (nz, ni) [1].

        Returns
        -------
        ClimateUpdate
            With `solve_rce`, AGNI's temperature and Kzz (see `climate_update`); otherwise
            `pco` and `tco` unchanged, with AGNI's Kzz [cm2 s-1].

        Raises
        ------
        RuntimeError
            If AGNI does not converge with `solve_rce`.
        """
        log.info('Running climate calculation...')
        set_composition(self.atmos, pco, species, ymix)
        if self.cfg.solve_rce:
            log.info('Using nonlinear solver to conserve fluxes')
            if not _solve_energy(self.atmos, self.cfg):
                raise RuntimeError('AGNI did not converge at this climate update')
            log.info(spacer)
            return climate_update(self.atmos)

        # Keep the temperature on VULCAN's own grid: a round trip through AGNI's grid would
        # smooth it at every update. AGNI only fills Kzz from it.
        log.info('Using the present temperature profile')
        set_temperature(self.atmos, pco, tco)
        _julia().AGNI.energy.calc_fluxes_b(self.atmos, convective=True)
        update = climate_update(self.atmos)
        log.info(spacer)
        return ClimateUpdate(p=pco, tmp=tco, p_kzz=update.p_kzz, kzz=update.kzz)

    def close(self):
        """Deallocate the AGNI atmosphere."""
        deallocate_atmos(self.atmos)
