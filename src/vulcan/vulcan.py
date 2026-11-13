#!/usr/bin/env python3

"""

----

Main module for running VULCAN.
Contains the `run_vulcan` function, which executes the main VULCAN workflow, and the `run_cli` function, which serves as the entry point when running the script directly.

Imports
----
- [`vulcan.config`](config.md): `Config` class for handling configuration settings.
- [`vulcan.paths`](paths.md): Paths to various files and directories used in VULCAN.
- [`vulcan.logs`](logs.md): Logger setup for VULCAN.
- [`vulcan.make_chem_funs`](make_chem_funs.md): `make_all`, generates all chemical functions from the chemical network.

"""

from __future__ import annotations

# Import system modules
import logging
import os
import sys
import time

from vulcan import paths
from vulcan.config import Config
from vulcan.logs import setup_logger
from vulcan.make_chem_funs import make_all

# Get number of threads from OMP_NUM_THREADS variable, if set
NUM_THREADS = 2
NUM_THREADS = max(1, int(os.environ.get('OMP_NUM_THREADS', NUM_THREADS)))
for k in ('MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[k] = str(NUM_THREADS)

# disable matplotlib debug logging
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)

# initialise logger for vulcan
log = logging.getLogger('fwl.' + __name__)


def check_climate_config(vulcan_cfg: Config, climate=None) -> None:
    """Refuse climate coupling settings that cannot work, before any work is done.

    Raises
    ------
    ValueError
        If `Kzz_prof = 'climate'` has no climate to come from; if temperature-dependent
        cross sections (read once, at the start) or fixed species (held at densities of the
        old temperature) are requested with a climate; or if `agni_T_rtol` or
        `agni_max_calls` is not positive.
    """
    coupled = climate is not None or vulcan_cfg.agni_call_frq > 0
    if vulcan_cfg.Kzz_prof == 'climate' and not coupled:
        raise ValueError("Kzz_prof = 'climate' needs a climate hook or agni_call_frq > 0")
    if coupled and vulcan_cfg.T_cross_sp:
        raise ValueError('T_cross_sp must be empty when the climate updates the temperature')
    if coupled and vulcan_cfg.fix_species:
        raise ValueError('fix_species must be empty when the climate updates the temperature')
    if not vulcan_cfg.agni_T_rtol > 0:
        raise ValueError(f'agni_T_rtol must be positive, got {vulcan_cfg.agni_T_rtol}')
    if not vulcan_cfg.agni_max_calls >= 1:
        raise ValueError(f'agni_max_calls must be at least 1, got {vulcan_cfg.agni_max_calls}')


def run_vulcan(vulcan_cfg: Config, make_network: bool = True, climate=None):
    """Run VULCAN with the given configuration.

    Parameters
    ----------
    vulcan_cfg : Config
        The configuration object containing all settings for the VULCAN run.
    make_network : bool, optional
        Whether to remake the chemical network (chem_funs.py). Default is True.
    climate : callable, optional
        Climate hook (see `vulcan.climate`) that updates the temperature, and with
        `Kzz_prof = 'climate'` the eddy diffusion, before the integration, every
        `agni_call_frq` steps and at each steady state. Without one, `agni_call_frq > 0`
        couples to AGNI through `vulcan.agni.AgniClimate`.
    """

    check_climate_config(vulcan_cfg, climate)

    # make network?
    if make_network:
        log.debug('Making chem_funs.py ...')
        make_all(vulcan_cfg)
    else:
        log.debug('Skip making chem_funs.py')

    log.info('Running VULCAN')

    # Import modules for running VULCAN
    from . import build_atm, op, store

    ### read in the basic chemistry data
    with open(paths.COM_FILE, 'r') as f:
        columns = f.readline()  # reading in the first line
        num_ele = (
            len(columns.split()) - 2
        )  # number of elements (-2 for removing "species" and "mass")
    type_list = ['int' for i in range(num_ele)]
    type_list.insert(0, 'U20')
    type_list.append('float')

    ### create the instances for storing the variables and parameters
    data_var = store.Variables(vulcan_cfg)
    data_atm = store.AtmData(vulcan_cfg)
    data_para = store.Parameters(vulcan_cfg)

    # record starting CPU time
    data_para.start_time = time.time()

    # create atmosphere object
    make_atm = build_atm.Atm(vulcan_cfg)

    # for plotting and printing
    output = op.Output(vulcan_cfg)

    # write config to disk
    vulcan_cfg.write_file()

    # construct pico
    data_atm = make_atm.f_pico(data_atm)

    # construct Tco and Kzz
    data_atm = make_atm.load_TPK(data_atm)

    # construct Dzz (molecular diffusion)
    # Only setting up ms (the species molecular weight) if vulcan_cfg.use_moldiff == False
    make_atm.mol_diff(data_atm)

    # calculating the saturation pressure
    if vulcan_cfg.use_condense:
        make_atm.sp_sat(data_atm)

    # for reading rates
    rate = op.ReadRate(vulcan_cfg)

    # read-in network and calculating forward rates
    data_var = rate.read_rate(data_var, data_atm)

    # for low-T rates e.g. Jupiter
    if vulcan_cfg.use_lowT_limit_rates:
        data_var = rate.lim_lowT_rates(data_var, data_atm)

    # reversing rates
    data_var = rate.rev_rate(data_var, data_atm)

    # removing rates
    data_var = rate.remove_rate(data_var)

    # initialing y and ymix (the number density and the mixing ratio of every species)
    ini_abun = build_atm.InitialAbun(vulcan_cfg)
    data_var = ini_abun.ini_y(data_var, data_atm)

    # storing the initial total number of atmos
    data_var = ini_abun.ele_sum(data_var)

    # calculating mean molecular weight, dz, and dzi and plotting TP
    data_atm = make_atm.f_mu_dz(data_var, data_atm, output)

    # specify the BC
    make_atm.BC_flux(data_atm)

    # ============== Execute VULCAN  ==============
    # time-steping in the while loop until conv() returns True or count > count_max

    # setting the numerical solver to the desinated one in vulcan_cfg
    solver = op.Ros2(vulcan_cfg)

    # Setting up for photo chemistry
    if vulcan_cfg.use_photo:
        rate.make_bins_read_cross(data_var, data_atm)
        # rate.read_cross(data_var)
        make_atm.read_sflux(data_var, data_atm)

        # computing the optical depth (tau), flux, and the photolisys rates (J) for the first time
        solver.compute_tau(data_var, data_atm)
        solver.compute_flux(data_var, data_atm)
        solver.compute_J(data_var, data_atm)
        # they will be updated in op.Integration by the assigned frequence

        # removing rates
        data_var = rate.remove_rate(data_var)

    # Assgining the specific solver corresponding to different B.C.s
    integ = op.Integration(solver, output, vulcan_cfg)
    solver.naming_solver(data_para)

    # Running the integration loop
    log.info(f'Starting VULCAN integration with {os.environ["OMP_NUM_THREADS"]} threads...')
    # climate hook: the caller's, or AGNI's when agni_call_frq > 0
    own_climate = climate is None and vulcan_cfg.agni_call_frq > 0
    try:
        if own_climate:
            from .agni import AgniClimate

            climate = AgniClimate(vulcan_cfg, data_atm, data_var)
        if climate is not None:
            integ.update_climate(climate, data_var, data_atm, data_para, make_atm)
        integ(data_var, data_atm, data_para, make_atm, climate=climate)
    finally:
        if own_climate and climate is not None:
            climate.close()

    # Save result to disk
    output.save_out(data_var, data_atm, data_para)


def run_cli():
    """Entry point for the script when run directly."""

    print('Starting VULCAN from command line')

    # Make config
    vulcan_cfg = Config()

    # Make output folder
    if not os.path.exists(vulcan_cfg.output_dir):
        os.makedirs(vulcan_cfg.output_dir)

    # Setup basic logging
    logpath = os.path.join(vulcan_cfg.output_dir, 'vulcan.log')
    setup_logger(logpath=logpath, logterm=True, level=vulcan_cfg.log_level)
    log = logging.getLogger('fwl.' + __name__)

    # Remake chem_funs by default. Disabled when passing -n flag.
    make_network = bool('-n' not in sys.argv)

    # Run model
    run_vulcan(vulcan_cfg, make_network)


if __name__ == '__main__':
    run_cli()
