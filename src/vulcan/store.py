"""
----

Module that contains classes that store all the variables used in VULCAN.

Copyright (C) 2016 Shang-Min Tsai (Shami)

There are three classes: `Variables`, `AtmData`, and `Parameters`, which store the main physical variables, atmospheric variables, and numerical parameters, respectively.

Imports
---
- [`vulcan.config`](config.md): `Config` class for handling configuration settings.
- `vulcan.chem_funs`: `ni`, `spec_list`: number of species and the species list of the network.

Array conventions shared by all three classes: index 0 of the vertical axis is the bottom
of the atmosphere; `nz` is the number of layers (`Config.nz`) and `ni` the number of species.
Quantities with `nz - 1` values sit on the interior interfaces between layers.

"""

from __future__ import annotations

import numpy as np

from vulcan.config import Config

from .chem_funs import ni, spec_list  # number of species and reactions in the network


class Variables(object):
    """
    Composition, rate coefficients and time-stepping state of a run.

    One instance is passed through the whole run (usually as `var` or `data_var`) and is
    updated in place by the solver. Photochemistry attributes such as `bins`, `sflux_top`,
    `tau`, `aflux` and `cross` are added later by `op.ReadRate.make_bins_read_cross` and
    `op.ODESolver.compute_tau` / `compute_flux`.

    Attributes
    ----------
    y : np.ndarray
        Number densities, shape (nz, ni) [cm-3].
    y_prev : np.ndarray
        Number densities at the previous accepted step, shape (nz, ni) [cm-3].
    y_ini : np.ndarray
        Initial number densities, shape (nz, ni) [cm-3].
    ymix : np.ndarray
        Mixing ratios, shape (nz, ni) [dimensionless]. With condensation, normalised by the
        sum over the gas species only (`AtmData.gas_indx`).
    k : dict
        Rate coefficients keyed by reaction number, each of shape (nz,): [s-1] for
        unimolecular and photochemical reactions, [cm3 s-1] for two-body reactions and
        [cm6 s-1] for three-body reactions.
    t : float
        Integration time [s].
    dt : float
        Current time step [s]; starts at `Config.dttry`.
    dy : float
        Largest relative change of y over the last step [dimensionless].
    dy_prev : float
        `dy` of the step before [dimensionless].
    dydt : float
        `dy / dt` [s-1].
    longdy : float
        Largest change of the mixing ratios, relative to themselves, over the convergence
        look-back window (dn in eq. 10 of Tsai et al. 2017) [dimensionless].
    longdydt : float
        `longdy` divided by the length of that window (dn/dt in eq. 11 of Tsai et al. 2017)
        [s-1].
    aflux_change : float
        Largest relative change of the actinic flux at the last photochemistry update
        [dimensionless]; stays 0 without photochemistry.
    y_time, t_time : list
        `y` [cm-3] and `t` [s] after each accepted step.
    atom_loss_time : list
        `atom_loss` values after each accepted step [dimensionless].
    atom_ini : dict
        Per element, the initial sum over all species and layers of the number of atoms
        times the number density [cm-3]. This is a sum over layers, not a column density.
    atom_sum : dict
        The same sum for the current `y` [cm-3].
    atom_loss : dict
        Relative change `(atom_sum - atom_ini) / atom_ini` per element [dimensionless].
    atom_loss_prev : dict
        `atom_loss` at the previous step [dimensionless].
    atom_conden : dict
        Per element; set to 0 by `build_atm.InitialAbun.ele_sum` and not updated elsewhere.
    Rf : dict
        Reaction equation string keyed by reaction number.
    Rindx : dict
        Reaction number as written in the network file, keyed by the internal number.
    a, n, E : dict
        Arrhenius parameters of k = a T^n exp(-E/T), keyed by reaction number: a in the
        units of k, n [dimensionless], E [K].
    a_inf, n_inf, E_inf : dict
        The same parameters for the high-pressure limit of three-body reactions.
    k_fun : dict
        Callables `k_fun[i](T, M)` returning k[i] from temperature [K] and total number
        density [cm-3].
    kinf_fun : dict
        Callables returning the high-pressure-limit rate coefficient [cm3 s-1].
    k_fun_new : dict
        Pressure-dependent `k_fun` of three-body reactions with a high-pressure limit;
        merged into `k_fun` by `op.ReadRate.read_rate`.
    k_inf : dict
        Not filled by the current code.
    photo_sp : set
        Species with photodissociation reactions in the network.
    pho_rate_index : dict
        Reaction number keyed by (species, branch index) for photodissociation.
    n_branch : dict
        Number of photodissociation branches per species.
    ion_sp : set
        Species with photoionisation reactions in the network.
    ion_rate_index, ion_branch : dict
        As `pho_rate_index` and `n_branch`, for photoionisation.
    charge_list : list
        Species with non-zero charge, excluding the electron 'e'.
    def_bin_min, def_bin_max : float
        Wavelength range of the stellar spectrum clipped to 2-700 nm [nm]; narrowed further
        by `op.ReadRate.make_bins_read_cross`.
    threshold : dict
        Photodissociation threshold wavelength per species [nm], set by
        `op.ReadRate.make_bins_read_cross`.
    cross_T_sp_list : dict
        Temperatures with tabulated cross sections, per species [K].
    conden_re_list : list
        Reaction numbers of the condensation reactions.
    var_save : list of str
        Names of the attributes written to the output file.
    var_evol_save : list of str
        Names of the time series written to the output when `Config.save_evolution` is set.
    """

    def __init__(self, vulcan_cfg: Config):
        """
        Parameters
        ----------
        vulcan_cfg : Config
            Run configuration; reads `nz`, `dttry`, `sflux_file`, `use_photo`,
            `T_cross_sp` and `use_ion`.
        """

        nz = vulcan_cfg.nz

        self.k = {}  # rate coefficients: k[1] is the rate constant of R1 reaction at every level (same shape as Tco and pco)
        self.y = np.zeros(
            (nz, ni)
        )  # current number density in the shape of (number of vertical levels, number of species)
        self.y_prev = np.zeros((nz, ni))  # number density at the previous step
        self.ymix = np.zeros((nz, ni))  # current mixing ratios
        self.y_ini = np.zeros((nz, ni))  # the initial number density
        self.t = 0  # integration time
        self.dt = vulcan_cfg.dttry  # time step
        self.dy = 1.0  # the max change of y from the previous step to current step
        self.dy_prev = 1.0  # the max change of y from i-2 step to i-1 step
        self.dydt = 1.0  # the max change of dy/dt
        self.longdy = 1.0  # the max change of y over a period of time
        # for checking the steady state (dn in (10) in Tsai et al 2017)
        self.longdydt = 1.0  # the max change of dy/dt over a period of time
        # (dn/dt in (11) in Tsai et al 2017)
        # self.slope_min = 0 #

        self.dy_time = []  # storing dy at each step
        self.dydt_time = []
        self.atim_loss_time = []
        self.ymix_time = []  # storing the mixing ratio at each step
        self.y_time = []  # storing the number density at each step
        self.t_time = []  # storing the time  at each step
        self.dt_time = []  # storing the time step  at each step
        self.atom_loss_time = []  # storing the loss of atoms at each step

        self.atom_ini = {}
        self.atom_sum = {}
        self.atom_loss = {}
        self.atom_loss_prev = {}
        self.atom_conden = {}

        self.Rf = {}  # reaction equation
        self.Rindx = {}
        (
            self.a,
            self.n,
            self.E,
            self.a_inf,
            self.n_inf,
            self.E_inf,
        ) = [{} for i in range(6)]
        self.k_fun, self.k_inf = [{} for i in range(2)]
        self.photo_sp = set()
        self.pho_rate_index, self.n_branch, self.wavelen = {}, {}, {}
        self.ion_rate_index, self.ion_branch, self.ion_wavelen, self.ion_br_ratio = (
            {},
            {},
            {},
            {},
        )
        self.charge_list, self.ion_sp = (
            [],
            set(),
        )  # charge_list: list of species with non-zero charge; ion_sp: species subjected to photoionisation

        self.kinf_fun = {}
        self.k_fun_new = {}

        # the max change of the actinic flux (for convergence)
        # if photochemistry is off, the value remaines 0 for checking convergence
        self.aflux_change = 0.0

        # The temporary wavelegth range (nm) given by the stellar flux
        # It will later be adjusted in make_bins_read_cross in op.py considering all molecules the absorbe photons, taking the smaller range of the two
        sflux_data = np.genfromtxt(
            vulcan_cfg.sflux_file, dtype=float, skip_header=1, names=['lambda', 'flux']
        )

        # setting the spectral bins based on the stellar spectrum, bun not shorter than 2nm ann not longer than 700 nm. This will further be ajusted in op.py while reading cross sections
        self.def_bin_min = max(sflux_data['lambda'][0], 2.0)
        self.def_bin_max = min(sflux_data['lambda'][-1], 700.0)

        # Define what variables to save in the output file!
        self.var_save = [
            'k',
            'y',
            'ymix',
            'y_ini',
            't',
            'dt',
            'longdy',
            'longdydt',
            'atom_ini',
            'atom_sum',
            'atom_loss',
            'atom_conden',
            'aflux_change',
            'Rf',
        ]
        if vulcan_cfg.use_photo:
            self.var_save.extend(
                [
                    'nbin',
                    'bins',
                    'dbin1',
                    'dbin2',
                    'tau',
                    'sflux',
                    'aflux',
                    'cross',
                    'cross_scat',
                    'cross_J',
                    'J_sp',
                    'n_branch',
                ]
            )
            if vulcan_cfg.T_cross_sp:
                self.var_save.extend(['cross_J', 'cross_T'])
            if vulcan_cfg.use_ion:
                self.var_save.extend(
                    [
                        'charge_list',
                        'ion_sp',
                        'cross_Jion',
                        'Jion_sp',
                        'ion_wavelen',
                        'ion_branch',
                        'ion_br_ratio',
                    ]
                )
        # 'ion_list' stores all the non-neutral species in build.atm whereas 'ion_sp' is for the species that actually have ionisation reactions in the network
        self.var_evol_save = ['y_time', 't_time']
        self.conden_re_list = []

        # new for rading ratios
        self.threshold = {}
        # list of avaliable temperatures of cross sections
        self.cross_T_sp_list = {}

        # TEST
        self.v_ratio = np.ones(nz)


class AtmData(object):
    """
    Atmospheric structure: pressure, temperature, heights, mixing and boundary conditions.

    One instance (usually `atm` or `data_atm`) is filled by the methods of `build_atm.Atm`
    and read by the solver. Layer i lies between the interfaces i and i + 1.

    Attributes
    ----------
    pco : np.ndarray
        Pressure at the layer centres, log-spaced from `P_b` to `P_t`, shape (nz,)
        [dyn cm-2].
    pico : np.ndarray
        Pressure at the interfaces, shape (nz + 1,) [dyn cm-2]; set by `Atm.f_pico`.
    Tco : np.ndarray
        Temperature at the layer centres, shape (nz,) [K]; set by `Atm.load_TPK`.
    M, n_0 : np.ndarray
        Total number density p / (k_B T), shape (nz,) [cm-3]. `M` is used as the third body.
    Kzz : np.ndarray
        Eddy diffusion coefficient at the interior interfaces, shape (nz - 1,) [cm2 s-1].
    vz : np.ndarray
        Vertical advection velocity at the interior interfaces, shape (nz - 1,)
        [cm s-1]; positive is upward.
    dz : np.ndarray
        Layer thickness, shape (nz,) [cm].
    dzi : np.ndarray
        Distance between the centres of adjacent layers, shape (nz - 1,) [cm].
    zco : np.ndarray
        Height of the interfaces, shape (nz + 1,) [cm]. Zero at the bottom, or for a gas
        giant with `P_b` >= 1 bar at the interface closest to 1 bar.
    zmco : np.ndarray
        Height of the layer centres, shape (nz,) [cm].
    Hp : np.ndarray
        Pressure scale height, shape (nz,) [cm].
    mu : np.ndarray
        Mean molar mass, shape (nz,) [g mol-1].
    ms : np.ndarray
        Molar mass of each species, shape (ni,) [g mol-1].
    gs : float
        Gravity at the surface, or at 1 bar for a gas giant [cm s-2].
    g : np.ndarray
        Gravity in each layer, shape (nz,) [cm s-2].
    Dzz : np.ndarray
        Molecular diffusion coefficient at the interior interfaces, shape (nz - 1, ni)
        [cm2 s-1]; set by `Atm.mol_diff`.
    Dzz_cen : np.ndarray
        Shape (nz, ni); not used by the current code.
    vm : np.ndarray
        Molecular diffusion velocity read by `op.ODESolver.diffdf_settling_vm`, shape
        (nz, ni) [cm s-1]. Nothing in the current code sets it, so it stays zero.
    vs : np.ndarray
        Settling velocity of the non-gaseous species at the interior interfaces, shape
        (nz - 1, ni) [cm s-1]; negative is downward.
    alpha : np.ndarray
        Thermal diffusion factor of each species, shape (ni,) [dimensionless]; set by
        `Atm.mol_diff` for the chosen `atm_base`.
    top_flux : np.ndarray
        Flux at the top boundary, shape (ni,) [cm-2 s-1]; positive into the atmosphere,
        negative for escape.
    bot_flux : np.ndarray
        Flux at the bottom boundary, shape (ni,) [cm-2 s-1]; positive into the atmosphere.
    bot_vdep : np.ndarray
        Deposition velocity at the bottom boundary, shape (ni,) [cm s-1]; the surface loss
        is `bot_vdep * y[0]`.
    bot_fix_sp : np.ndarray
        Fixed bottom mixing ratios read by `Atm.BC_flux`, shape (ni,) [dimensionless]. The
        solver takes the fixed values from `Config.use_fix_sp_bot` instead.
    sat_p : dict
        Saturation vapour pressure per condensing species, shape (nz,) [dyn cm-2].
    sat_mix : dict
        Saturation mixing ratio `sat_p / pco` per condensing species, capped at 1 and
        scaled by `Config.humidity` for H2O, shape (nz,) [dimensionless].
    gas_indx : list of int
        Species indices that are not in `Config.non_gas_sp`.
    fix_sp_indx : dict
        Per species in `Config.fix_species`, the indices of that species in the flattened
        (nz * ni) array of y.
    fix_e_indx : np.ndarray
        The same indices for the electron; only with `use_ion`.
    r_p : dict
        Particle radius per condensate [cm].
    rho_p : dict
        Particle density per condensate [g cm-3].
    """

    def __init__(self, vulcan_cfg: Config):
        """
        Parameters
        ----------
        vulcan_cfg : Config
            Run configuration; reads `nz`, `P_b`, `P_t`, `gs`, `non_gas_sp`, `fix_species`,
            `use_ion`, `use_condense`, `r_p` and `rho_p`.
        """

        nz = vulcan_cfg.nz

        self.pco = np.logspace(
            np.log10(vulcan_cfg.P_b), np.log10(vulcan_cfg.P_t), nz
        )  # pressure grids
        self.pico = np.empty(nz + 1)  # pressure grids at the interface
        self.dz = np.zeros(nz)  # height grids
        self.dzi = np.zeros(nz - 1)  # height grids at the interface
        self.zco = np.zeros(nz + 1)  # not used in calculation
        self.zmco = np.empty(nz)
        self.Tco = np.empty(nz)  # temperature grids
        self.Kzz = np.zeros(nz - 1)
        self.vz = np.zeros(nz - 1)  # vertical velocity
        self.M = np.empty(nz)  # the number density of the third body
        self.n_0 = np.empty(nz)  # the total number density
        self.Hp = np.empty(nz)  # the scale height
        self.mu = np.empty(nz)  # mean molecular weight
        self.ms = np.empty(ni)  # molecular weight for every species
        self.Dzz = np.zeros((nz - 1, ni))  # molecular diffusion (nz,ni)
        self.Dzz_cen = np.zeros((nz, ni))  # molecular diffusion (nz,ni)
        self.vm = np.zeros((nz, ni))  # molecular diffusion (nz,ni)
        self.vs = np.zeros((nz - 1, ni))  # the settling velocity
        self.alpha = np.zeros(
            ni
        )  # thermal diffusion factor = -0.25 for every species except for H and H2 (defined in mol_diff() in build_atm.py)
        self.gs = vulcan_cfg.gs  # the gravitational acceleration at the surface or at 1 bar
        self.g = np.zeros(nz)  # g(z)

        self.top_flux = np.zeros(ni)  # the assigned flux at the top boundary
        self.bot_flux = np.zeros(ni)  # the assigned flux at the bottom boundary
        self.bot_vdep = np.zeros(
            ni
        )  # the deposition velocity at the bottom boundary (surface sink)
        self.bot_fix_sp = np.zeros(
            ni
        )  # the fixed mixing ratios for certain species at the bottom boundary

        self.sat_p = {}
        self.sat_mix = {}
        # excluding non-gaseous species while computing ymix from y
        self.gas_indx = [_ for _ in range(ni) if spec_list[_] not in vulcan_cfg.non_gas_sp]

        self.fix_sp_indx = {}
        if hasattr(vulcan_cfg, 'fix_species'):  # if fix_species is defined in vulcan_cfg
            for sp in vulcan_cfg.fix_species:
                self.fix_sp_indx[sp] = np.arange(
                    spec_list.index(sp), spec_list.index(sp) + ni * nz, ni
                )

        # turning of df(e) while using charge conservation. it'll be very slow if not doing this
        if vulcan_cfg.use_ion:
            self.fix_e_indx = np.arange(
                spec_list.index('e'), spec_list.index('e') + ni * nz, ni
            )

        # TEST condensation excluding non-gaseous species
        self.r_p, self.rho_p = {}, {}
        if vulcan_cfg.use_condense:
            for sp in vulcan_cfg.r_p.keys():
                self.r_p[sp] = vulcan_cfg.r_p[sp]
            for sp in vulcan_cfg.rho_p.keys():
                self.rho_p[sp] = vulcan_cfg.rho_p[sp]


class Parameters(object):
    """
    Counters and switches of the numerical integration.

    Attributes
    ----------
    count : int
        Number of steps taken.
    delta : float
        Relative truncation error of the last step [dimensionless]; a step with
        `delta > Config.rtol` is rejected.
    nega_count : int
        Steps rejected because of negative number densities.
    loss_count : int
        Steps rejected because the element conservation changed too much.
    delta_count : int
        Steps rejected because of the truncation error.
    small_y : float
        Running sum of the small positive number densities clipped to 0 [cm-3].
    nega_y : float
        Running sum of the magnitudes of the negative number densities clipped to 0
        [cm-3].
    end_case : int
        Why the integration stopped: 0 not set, 1 steady state, 2 runtime limit,
        3 step limit (`count_max`).
    solver_str : str
        Name of the solver method that `op.Ros2.one_step` calls.
    switch_final_photo_frq : bool
        True once the photochemistry update interval has switched to
        `Config.final_update_photo_frq` near convergence.
    where_varies_most : np.ndarray
        `longdy / ymix` per layer and species at the last convergence check, shape
        (nz, ni) [dimensionless].
    pic_count : int
        Number of live-plot frames written.
    fix_species_start : bool
        True once the condensed species in `Config.fix_species` are held fixed.
    tableau20 : list of tuple
        The 20 "Tableau 20" colours as RGB in [0, 1].
    """

    def __init__(self, vulcan_cfg: Config):
        """
        Parameters
        ----------
        vulcan_cfg : Config
            Run configuration; reads `nz`.
        """

        nz = vulcan_cfg.nz

        self.nega_y = 0
        self.small_y = 0
        self.delta = 0
        self.count = 0  # number of steps taken
        self.nega_count = 0
        self.loss_count = 0
        self.delta_count = 0
        self.end_case = 0
        self.solver_str = ''  # for assigning the name of solver
        self.switch_final_photo_frq = False
        self.where_varies_most = np.zeros(
            (nz, ni)
        )  # recording from where and what species flucating from convergence
        self.pic_count = 0  # for live plotting

        # TEST
        self.fix_species_start = False
        # self.conden_water_start = False

        # These are the "Tableau 20" colors as RGB.
        self.tableau20 = [
            (31, 119, 180),
            (255, 127, 14),
            (44, 160, 44),
            (214, 39, 40),
            (148, 103, 189),
            (140, 86, 75),
            (227, 119, 194),
            (127, 127, 127),
            (188, 189, 34),
            (23, 190, 207),
            (174, 199, 232),
            (255, 187, 120),
            (152, 223, 138),
            (255, 152, 150),
            (197, 176, 213),
            (196, 156, 148),
            (247, 182, 210),
            (199, 199, 199),
            (219, 219, 141),
            (158, 218, 229),
        ]

        # Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
        for i in range(len(self.tableau20)):
            r, g, b = self.tableau20[i]
            self.tableau20[i] = (r / 255.0, g / 255.0, b / 255.0)
