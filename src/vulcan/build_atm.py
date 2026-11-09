"""
---

Module for building the atmospheric structure and calculating the initial mixing ratios for VULCAN.

Imports
---
- [`vulcan.config`](config.md): `Config` class for handling configuration settings.
- [`vulcan.paths`](paths.md): Paths to various files and directories used in VULCAN.
- [`vulcan.phy_const`](phy_const.md): Physical constants used in VULCAN.
- `chem_funs`: `ni`, `spec_list`: number of species and reactions in the network, and the list of species.

"""

from __future__ import annotations

import logging
import os
import pickle
import subprocess
from shutil import copyfile

import numpy as np
import scipy
import scipy.optimize as sop
from scipy import interpolate

from vulcan import plot
from vulcan.config import Config
from vulcan.paths import COM_FILE, FASTCHEM_DIR
from vulcan.phy_const import Navo, au, kb, r_sun

from .chem_funs import ni  # number of species and reactions in the network
from .chem_funs import spec_list as species

log = logging.getLogger('fwl.' + __name__)

### read in the basic chemistry data
with open(COM_FILE, 'r') as f:
    columns = f.readline()  # reading in the first line
    num_ele = (
        len(columns.split()) - 2
    )  # number of elements (-2 for removing "species" and "mass")
type_list = ['int' for i in range(num_ele)]
type_list.insert(0, 'U20')
type_list.append('float')
compo = np.genfromtxt(COM_FILE, names=True, dtype=type_list)
compo_row = list(compo['species'])


class InitialAbun(object):
    """
    Initial composition of the atmosphere, chosen by `Config.ini_mix`.

    The options are 'eq' (chemical equilibrium from FastChem), 'vulcan_ini' (the result of
    a previous run), 'table' (mixing ratios from a text file), 'const_mix' (constant mixing
    ratios) and 'const_lowt' (H2, H2O, CH4, NH3 and He solved from the elemental ratios).
    """

    def __init__(self, vulcan_cfg: Config):
        """
        Parameters
        ----------
        vulcan_cfg : Config
            Run configuration. `ini_mix` is converted to lower case in place.
        """
        self.cfg = vulcan_cfg
        self.ini_m = [0.9, 0.1, 0.0, 0.0, 0]  # initial guess
        self.atom_list = vulcan_cfg.atom_list
        self.cfg.ini_mix = self.cfg.ini_mix.lower()  # normalise string to lowercase

    def abun_lowT(self, x):
        """
        Residuals of the low-temperature composition, for `scipy.optimize.fsolve`.

        Parameters
        ----------
        x : sequence of float
            Mixing ratios of H2, H2O, CH4, He and NH3 [dimensionless].

        Returns
        -------
        tuple of float
            Residuals of x1 + x2 + x3 + x4 = 1 and of the O/H, C/H, He/H and N/H number
            ratios (`Config.O_H`, `C_H`, `He_H`, `N_H`).
        """
        O_H, C_H, He_H, N_H = self.cfg.O_H, self.cfg.C_H, self.cfg.He_H, self.cfg.N_H
        x1, x2, x3, x4, x5 = x
        f1 = x1 + x2 + x3 + x4 - 1.0
        f2 = x2 - (2 * x1 + 2 * x2 + 4 * x3 + 3 * x5) * O_H
        f3 = x3 - (2 * x1 + 2 * x2 + 4 * x3 + 3 * x5) * C_H
        f4 = x4 - (2 * x1 + 2 * x2 + 4 * x3 + 3 * x5) * He_H
        f5 = x5 - (2 * x1 + 2 * x2 + 4 * x3 + 3 * x5) * N_H
        return f1, f2, f3, f4, f5

    def abun_highT(self, x):
        """
        Residuals of the high-temperature composition, for `scipy.optimize.fsolve`.

        Not called by the current code; `ini_mol` uses `abun_lowT` only.

        Parameters
        ----------
        x : sequence of float
            Mixing ratios of H2, H2O, CO, He and N2 [dimensionless].

        Returns
        -------
        tuple of float
            Residuals of x1 + x2 + x3 + x4 = 1 and of the O/H, C/H, He/H and N/H number
            ratios.
        """
        O_H, C_H, He_H, N_H = self.cfg.O_H, self.cfg.C_H, self.cfg.He_H, self.cfg.N_H
        x1, x2, x3, x4, x5 = x
        f1 = x1 + x2 + x3 + x4 - 1.0
        f2 = x2 + x3 - (2 * x1 + 2 * x2) * O_H
        f3 = x3 - (2 * x1 + 2 * x2) * C_H
        f4 = x4 - (2 * x1 + 2 * x2) * He_H
        f5 = x5 * 2 - (2 * x1 + 2 * x2) * N_H
        return f1, f2, f3, f4, f5

    def ini_mol(self):
        """
        Solve `abun_lowT` for `ini_mix = 'const_lowt'`.

        Returns
        -------
        np.ndarray or None
            Mixing ratios of H2, H2O, CH4, He and NH3 [dimensionless], or None for any
            other `ini_mix`.
        """
        if self.cfg.ini_mix == 'const_lowt':
            return np.array(sop.fsolve(self.abun_lowT, self.ini_m))

    def ini_fc(self, data_var, data_atm):
        """
        Run FastChem for the equilibrium composition on the current T-P profile.

        Writes the FastChem inputs inside `paths.FASTCHEM_DIR`: it copies the parameter file
        with or without ions over `input/parameters.dat`, writes the elemental abundances to
        `input/element_abundances_vulcan.dat` and the profile to `input/vulcan_TP/
        vulcan_TP.dat` (pressure in bar, T in K). With `Config.use_solar = False` the
        elements in `atom_list` take the ratios `<X>_H` from the configuration and the other
        FastChem elements are scaled by `fastchem_met_scale`. FastChem writes its result to
        `output/vulcan_EQ.dat`, which `ini_y` reads.

        Parameters
        ----------
        data_var : store.Variables
            Not used.
        data_atm : store.AtmData
            Reads `pco` [dyn cm-2] and `Tco` [K].

        Raises
        ------
        RuntimeError
            If the FastChem executable is missing.
        """

        # check that fastchem is compiled in the correct place
        if not os.path.isfile(os.path.join(FASTCHEM_DIR, 'fastchem')):
            raise RuntimeError(
                f'FastChem cannot be found. Try compiling it by running `make` inside {FASTCHEM_DIR}'
            )

        # reading-in the default elemental abundances from Lodders 2009
        # depending on including ion or not (whether there is e- in the fastchem elemental abundance dat)
        solar_ele = os.path.join(FASTCHEM_DIR, 'input', 'solar_element_abundances.dat')
        if self.cfg.use_ion:
            copyfile(
                os.path.join(FASTCHEM_DIR, 'input', 'parameters_ion.dat'),
                os.path.join(FASTCHEM_DIR, 'input', 'parameters.dat'),
            )
        else:
            copyfile(
                os.path.join(FASTCHEM_DIR, 'input', 'parameters_wo_ion.dat'),
                os.path.join(FASTCHEM_DIR, 'input', 'parameters.dat'),
            )

        with open(solar_ele, 'r') as f:
            new_str = ''
            ele_list = list(self.cfg.atom_list)
            ele_list.remove('H')

            fc_list = [
                'C',
                'N',
                'O',
                'S',
                'P',
                'Si',
                'Ti',
                'V',
                'Cl',
                'K',
                'Na',
                'Mg',
                'F',
                'Ca',
                'Fe',
            ]

            if self.cfg.use_solar:
                new_str = f.read()  # read in as a string
                log.info('Initializing with the default solar abundance.')

            else:  # using costomized elemental abundances
                log.info('Initializing with the customized elemental abundance:')
                log.info('{:4}'.format('H') + str('1.'))
                for line in f.readlines():
                    li = line.split()
                    sp = li[0].strip()

                    if sp in ele_list:
                        # read-in vulcan_cfg.sp_H
                        sp_abun = getattr(self.cfg, sp + '_H')
                        fc_abun = 12.0 + np.log10(sp_abun)
                        line = sp + '\t' + '{0:.4f}'.format(fc_abun) + '\n'
                        log.info('{:4}'.format(sp) + '{0:.4E}'.format(sp_abun))

                    elif sp in fc_list:  # other elements included in fastchem but not in VULCAN
                        sol_ratio = li[1].strip()
                        if hasattr(
                            self.cfg, 'fastchem_met_scale'
                        ):  # vulcan_cfg.fastchem_met_scale exists
                            met_scale = self.cfg.fastchem_met_scale
                        else:
                            met_scale = 1.0
                            log.warning(
                                'fastchem_met_scale not specified in vulcan_cfg. Using solar metallicity for other elements not included in vulcan.'
                            )

                        new_ratio = float(sol_ratio) + np.log10(met_scale)
                        line = sp + '\t' + '{0:.4f}'.format(new_ratio) + '\n'

                    new_str += line

            # make the new elemental abundance file
            with open(
                os.path.join(FASTCHEM_DIR, 'input', 'element_abundances_vulcan.dat'), 'w'
            ) as f:
                f.write(new_str)

        # write a T-P text file for fast_chem
        vulcan_TP_dir = os.path.join(FASTCHEM_DIR, 'input', 'vulcan_TP')
        os.makedirs(vulcan_TP_dir, exist_ok=True)
        with open(os.path.join(vulcan_TP_dir, 'vulcan_TP.dat'), 'w') as f:
            ost = '#p (bar)    T (K)\n'
            for n, p in enumerate(data_atm.pco):  # p in bar in fast_chem
                ost += (
                    '{:.3e}'.format(p / 1.0e6) + '\t' + '{:.1f}'.format(data_atm.Tco[n]) + '\n'
                )
            ost = ost[:-1]
            f.write(ost)

        # run fastchem
        try:
            fc_cmd = [
                os.path.join(FASTCHEM_DIR, 'fastchem'),
                os.path.join(FASTCHEM_DIR, 'input', 'config.input'),
            ]
            log.debug(f'Fastchem command: {fc_cmd}')
            process = subprocess.Popen(
                fc_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=FASTCHEM_DIR,
            )

            def check_io():
                """Log each line FastChem has written to stdout, until none is waiting."""
                while True:
                    out_std = process.stdout.readline().decode().strip()
                    if out_std:
                        log.log(logging.INFO, out_std)
                    else:
                        break

            # keep checking stdout/stderr until the child exits
            while process.poll() is None:
                check_io()

        except subprocess.CalledProcessError:
            raise RuntimeError(
                f'FastChem cannot run properly. Try compiling it by running `make` inside {FASTCHEM_DIR}'
            )

    def ini_y(self, data_var, data_atm):
        """
        Set the initial number densities and mixing ratios for `Config.ini_mix`.

        With `use_condense`, also sets the saturation mixing ratios and, with
        `use_ini_cold_trap`, caps each condensing species at saturation and holds it above
        the cold trap at the saturation mixing ratio of the cold-trap level.

        Parameters
        ----------
        data_var : store.Variables
            Filled in place: `y` [cm-3], `y_ini` [cm-3], `ymix` and, with `use_ion`,
            `charge_list`.
        data_atm : store.AtmData
            Reads `pco` [dyn cm-2], `Tco` [K], `M` and `n_0` [cm-3], and with condensation
            `sat_p` [dyn cm-2]; sets `sat_mix` and `conden_min_lev`.

        Returns
        -------
        store.Variables
            `data_var`.

        Raises
        ------
        ValueError
            If a 'table' file has a different number of layers, or `use_ion` is set and the
            network has no charged species.
        IOError
            If `ini_mix` is not recognised.
        """

        nz = self.cfg.nz
        ini_mol = self.ini_mol()
        ini = np.zeros(ni)
        y_ini = data_var.y
        gas_tot = data_atm.M
        charge_list = []  # list of charged species excluding echarge_list

        if self.cfg.ini_mix == 'eq':
            self.ini_fc(data_var, data_atm)
            fc = np.genfromtxt(
                os.path.join(FASTCHEM_DIR, 'output', 'vulcan_EQ.dat'),
                names=True,
                dtype=None,
                skip_header=0,
            )
            for sp in species:
                # Atomic P hack (genfromtxt gets the Pressure as index otherwise)
                if sp == 'P':
                    y_ini[:, species.index(sp)] = (
                        fc['P_1'] * gas_tot
                    )  # this also changes data_var.y because the address of y array has passed to y_ini
                    # print(sp,y_ini[:,species.index(sp)]/gas_tot)
                    continue
                if sp in fc.dtype.names:
                    y_ini[:, species.index(sp)] = (
                        fc[sp] * gas_tot
                    )  # this also changes data_var.y because the address of y array has passed to y_ini

                else:
                    log.warning(sp + ' not included in fastchem.')

                if self.cfg.use_ion:
                    if compo[compo_row.index(sp)]['e'] != 0:
                        charge_list.append(sp)

            # remove the fc output
            subprocess.call(['rm', 'vulcan_EQ.dat'], cwd=os.path.join(FASTCHEM_DIR, 'output'))

        elif self.cfg.ini_mix == 'vulcan_ini':
            log.info('Initializing with compositions from pickle file ' + self.cfg.vul_ini)
            with open(self.cfg.vul_ini, 'rb') as handle:
                vul_data = pickle.load(handle)

            # y_ini = np.copy(vul_data['variable']['y'])
            # data_var.y = np.copy(y_ini)
            for sp in species:
                if sp in vul_data['variable']['species']:
                    y_ini[:, species.index(sp)] = vul_data['variable']['y'][
                        :, vul_data['variable']['species'].index(sp)
                    ]
                else:
                    log.warning(sp + ' not included in the input mixing ratios file.')

            # if vulcan_cfg.use_ion : charge_list = vul_data['variable']['charge_list']
            if self.cfg.use_ion:
                for sp in species:
                    if compo[compo_row.index(sp)]['e'] != 0:
                        charge_list.append(sp)

        elif self.cfg.ini_mix == 'table':
            log.info('Initializing with compositions from ASCII file ' + self.cfg.vul_ini)
            table = np.genfromtxt(self.cfg.vul_ini, names=True, dtype=None, skip_header=1)
            if not len(data_atm.pco) == len(table['Pressure']):
                raise ValueError(
                    'The input profile has different layers than the current setting'
                )
            for sp in species:
                try:
                    arr = data_atm.n_0 * table[sp]
                except (ValueError, KeyError):
                    arr = np.zeros(len(data_atm.pco))
                data_var.y[:, species.index(sp)] = arr

        elif self.cfg.ini_mix == 'const_mix':
            log.info(
                'Initializing with constant (well-mixed) abundances: ' + str(self.cfg.const_mix)
            )
            for sp in self.cfg.const_mix.keys():
                y_ini[:, species.index(sp)] = (
                    gas_tot * self.cfg.const_mix[sp]
                )  # this also changes data_var.y
            if self.cfg.use_ion:
                for sp in species:
                    if compo[compo_row.index(sp)]['e'] != 0:
                        charge_list.append(sp)

        else:
            for i in range(nz):
                if self.cfg.ini_mix == 'const_lowt':
                    y_ini[i, :] = ini
                    y_ini[i, species.index('H2')] = ini_mol[0] * gas_tot[i]
                    y_ini[i, species.index('H2O')] = ini_mol[1] * gas_tot[i]
                    y_ini[i, species.index('CH4')] = ini_mol[2] * gas_tot[i]
                    y_ini[i, species.index('NH3')] = ini_mol[4] * gas_tot[i]
                    # assign rest of the particles to He
                    y_ini[i, species.index('He')] = gas_tot[i] - np.sum(y_ini[i, :])

                else:
                    raise IOError(
                        'Initial mixing ratios unknown. Check the setting in vulcan_cfg.py.'
                    )

        if self.cfg.use_condense:
            for sp in self.cfg.condense_sp:
                data_atm.sat_mix[sp] = data_atm.sat_p[sp] / data_atm.pco

                # fixed 2022
                data_atm.sat_mix[sp] = np.minimum(1.0, data_atm.sat_mix[sp])

                if sp == 'H2O':
                    data_atm.sat_mix[sp] *= self.cfg.humidity

                if self.cfg.use_ini_cold_trap:
                    if self.cfg.ini_mix != 'table' and self.cfg.ini_mix != 'vul_ini':
                        # the level where condensation starts
                        conden_bot = np.argmax(
                            data_atm.n_0 * data_atm.sat_mix[sp]
                            <= data_var.y[:, species.index(sp)]
                        )
                        # conden_status: ture if the partial p >= the saturation p
                        sat_rho = data_atm.n_0 * data_atm.sat_mix[sp]
                        conden_status = data_var.y[:, species.index(sp)] >= sat_rho

                        # take the min between the mixing ratio and the saturation mixing ratio
                        data_var.y[:, species.index(sp)] = np.minimum(
                            data_atm.n_0 * data_atm.sat_mix[sp],
                            data_var.y[:, species.index(sp)],
                        )

                        if list(
                            data_var.y[conden_status, species.index(sp)]
                        ):  # if it condenses
                            min_sat = np.amin(
                                data_atm.sat_mix[sp][conden_status]
                            )  # the mininum value of the saturation p within the saturation region
                            conden_min_lev = np.where(data_atm.sat_mix[sp] == min_sat)[0][0]

                            data_atm.conden_min_lev = conden_min_lev

                            log.debug(
                                sp
                                + ' condensed from nz = '
                                + str(conden_bot)
                                + ' to the minimum level nz = '
                                + str(conden_min_lev)
                                + ' (cold trap)'
                            )
                            # data_var.y[conden_min_lev:,species.index(sp)] = (y_ini[conden_min_lev,species.index(sp)]/data_atm.n_0[conden_min_lev]) *data_atm.n_0[conden_min_lev:]
                            data_var.y[conden_min_lev:, species.index(sp)] = (
                                data_atm.sat_mix[sp][conden_min_lev]
                                * data_atm.n_0[conden_min_lev:]
                            )

        # re-normalisation
        # TEST
        # Excluding the non-gaseous species
        if self.cfg.use_condense:
            exc_conden = [_ for _ in range(ni) if species[_] not in self.cfg.non_gas_sp]
            ysum = np.sum(y_ini[:, exc_conden], axis=1).reshape((-1, 1))
        else:
            ysum = np.sum(y_ini, axis=1).reshape((-1, 1))

        data_var.y_ini = np.copy(y_ini)
        data_var.ymix = y_ini / ysum

        if self.cfg.use_ion:
            # if the charge_list is empty (no species with nonzero charges include)
            if not charge_list:
                raise ValueError(
                    'vulcan_cfg.use_ion = True but the network with ions is not supplied.'
                )
            else:
                if 'e' in charge_list:
                    charge_list.remove('e')
                data_var.charge_list = charge_list

        return data_var

    def ele_sum(self, data_var):
        """
        Record the initial element totals for the conservation check.

        Sets `atom_ini[atom]`, the sum over species and layers of the number of atoms times
        the number density [cm-3], and resets `atom_loss` and `atom_conden` to 0.

        Parameters
        ----------
        data_var : store.Variables
            Reads `y`; filled in place.

        Returns
        -------
        store.Variables
            `data_var`.
        """

        for atom in self.atom_list:
            data_var.atom_ini[atom] = np.sum(
                [compo[compo_row.index(species[i])][atom] * data_var.y[:, i] for i in range(ni)]
            )
            data_var.atom_loss[atom] = 0.0
            data_var.atom_conden[atom] = 0.0

        return data_var


class Atm(object):
    """
    Build the atmospheric structure stored in `store.AtmData`.

    The methods are called in this order by `vulcan.run_vulcan`: `f_pico`, `load_TPK`,
    `mol_diff`, `sp_sat` (with condensation), then after the initial abundances `f_mu_dz`
    and `BC_flux`, and `read_sflux` with photochemistry.
    """

    def __init__(self, vulcan_cfg: Config):
        """
        Parameters
        ----------
        vulcan_cfg : Config
            Run configuration; the profile settings are copied to attributes.
        """
        self.cfg = vulcan_cfg
        self.gs = vulcan_cfg.gs  # gravity
        self.P_b = vulcan_cfg.P_b
        self.P_t = vulcan_cfg.P_t
        self.type = vulcan_cfg.atm_type
        self.use_Kzz = vulcan_cfg.use_Kzz
        self.Kzz_prof = vulcan_cfg.Kzz_prof
        self.const_Kzz = vulcan_cfg.const_Kzz
        self.use_vz = vulcan_cfg.use_vz
        self.vz_prof = vulcan_cfg.vz_prof
        self.const_vz = vulcan_cfg.const_vz
        self.use_settling = vulcan_cfg.use_settling
        self.non_gas_sp = vulcan_cfg.non_gas_sp

    def f_pico(self, data_atm):
        """
        Set the interface pressures from the layer-centre pressures.

        Each interior interface takes the geometric mean of the two adjacent centres; the
        bottom and top interfaces are extrapolated by the same ratio in log p.

        Parameters
        ----------
        data_atm : store.AtmData
            Reads `pco` [dyn cm-2]; sets `pico`, shape (nz + 1,) [dyn cm-2].

        Returns
        -------
        store.AtmData
            `data_atm`.
        """

        pco = data_atm.pco

        # construct pico
        pco_up1 = np.roll(pco, 1)
        pi = (pco * pco_up1) ** 0.5
        pi[0] = pco[0] ** 1.5 * pco[1] ** (-0.5)
        pi = np.append(pi, pco[-1] ** 1.5 * pco[-2] ** (-0.5))

        # store pico
        data_atm.pico = pi
        # data_atm.pco = pco

        return data_atm

    def load_TPK(self, data_atm):
        """
        Set the temperature, eddy diffusion and vertical velocity profiles.

        `Config.atm_type` sets T: 'isothermal' (`Tiso`), 'analytical' (`TP_H14` with
        `para_anaTP`), 'file' (interpolated in p from `atm_file`, columns Pressure
        [dyn cm-2], Temp [K] and, with `Kzz_prof = 'file'`, Kzz [cm2 s-1]; constant beyond
        the file range), 'vulcan_ini' (from a previous run's pickle) or 'table' (p and T from
        `vul_ini`). `Kzz_prof` sets Kzz on the interior interfaces: 'const' (`const_Kzz`),
        'file', 'JM16' (1e5 (300 / p[mbar])^0.5 cm2 s-1, at least `K_deep`) or 'Pfunc'
        (`K_max` (`K_p_lev`[bar] / p[bar])^0.4, at least `K_max`). `vz_prof` sets vz:
        'const' (`const_vz`) or 'file' (column vz of `atm_file`, which needs
        `atm_type = 'file'`). Kzz and vz are zero when `use_Kzz` and `use_vz` are off.

        Parameters
        ----------
        data_atm : store.AtmData
            Reads `pco` and `pico` [dyn cm-2]; sets `Tco` [K], `Kzz` [cm2 s-1], `vz`
            [cm s-1], and `M` and `n_0` = p / (k_B T) [cm-3]. 'table' also replaces `pco`.

        Returns
        -------
        store.AtmData
            `data_atm`.

        Raises
        ------
        ValueError
            If `atm_type` is not recognised, a 'table' file has a different number of
            layers, or T is not finite and positive in every layer.
        IOError
            If `Kzz_prof` or `vz_prof` is not recognised.
        """

        PTK_fun = {}
        nz = self.cfg.nz

        # IF switches for TP types
        if self.type == 'isothermal':
            data_atm.Tco = np.repeat(self.cfg.Tiso, nz)
            # data_atm.Kzz = np.repeat(self.const_Kzz,nz-1)
            # data_atm.vz = np.repeat(self.const_vz,nz-1)

        elif self.type == 'analytical':
            # plotting T-P on the fly
            para_atm = self.cfg.para_anaTP

            # return the P-T function
            PTK_fun['pT'] = lambda pressure: self.TP_H14(pressure, *para_atm)
            data_atm.Tco = PTK_fun['pT'](data_atm.pco)
            # data_atm.Kzz = np.repeat(self.const_Kzz,nz-1)
            # data_atm.vz = np.repeat(self.const_vz,nz-1)

        # for atm_type = 'file' and also Kzz_prof = 'file
        elif self.type == 'file':
            if self.Kzz_prof == 'file':
                atm_table = np.genfromtxt(
                    self.cfg.atm_file, names=True, dtype=None, skip_header=1
                )
                p_file, T_file, Kzz_file = (
                    atm_table['Pressure'],
                    atm_table['Temp'],
                    atm_table['Kzz'],
                )

            else:
                atm_table = np.genfromtxt(
                    self.cfg.atm_file, names=True, dtype=None, skip_header=1
                )
                p_file, T_file = atm_table['Pressure'], atm_table['Temp']

            if max(p_file) < data_atm.pco[0] or min(p_file) > data_atm.pco[-1]:
                log.warning(
                    'P_b and P_t assgined in vulcan.cfg are out of range of the input. Constant extension is used.'
                )

            PTK_fun['pT'] = interpolate.interp1d(
                p_file,
                T_file,
                assume_sorted=False,
                bounds_error=False,
                fill_value=(T_file[np.argmin(p_file)], T_file[np.argmax(p_file)]),
            )
            # store Tco in data_atm
            try:
                data_atm.Tco = PTK_fun['pT'](data_atm.pco)

            # for SciPy earlier than v0.18.0
            except ValueError:
                PTK_fun['pT'] = interpolate.interp1d(
                    p_file,
                    T_file,
                    assume_sorted=False,
                    bounds_error=False,
                    fill_value=T_file[np.argmin(p_file)],
                )
                data_atm.Tco = PTK_fun['pT'](data_atm.pco)

            if self.use_Kzz and self.Kzz_prof == 'file':
                PTK_fun['pK'] = interpolate.interp1d(
                    p_file,
                    Kzz_file,
                    assume_sorted=False,
                    bounds_error=False,
                    fill_value=(Kzz_file[np.argmin(p_file)], Kzz_file[np.argmax(p_file)]),
                )
                # store Kzz in data_atm
                try:
                    data_atm.Kzz = PTK_fun['pK'](data_atm.pico[1:-1])
                # for SciPy earlier than v0.18.0
                except ValueError:
                    PTK_fun['pK'] = interpolate.interp1d(
                        p_file,
                        Kzz_file,
                        assume_sorted=False,
                        bounds_error=False,
                        fill_value=Kzz_file[np.argmin(p_file)],
                    )
                    data_atm.Kzz = PTK_fun['pK'](data_atm.pico[1:-1])

            elif self.Kzz_prof == 'const':
                data_atm.Kzz = np.repeat(self.const_Kzz, nz - 1)

        elif self.type == 'vulcan_ini':
            log.info('Initializing PT from the prvious run ' + self.cfg.vul_ini)
            with open(self.cfg.vul_ini, 'rb') as handle:
                vul_data = pickle.load(handle)

            data_atm.Tco = vul_data['atm']['Tco']

        elif self.type == 'table':
            log.info('Initializing PT from the prvious run ' + self.cfg.vul_ini)
            table = np.genfromtxt(self.cfg.vul_ini, names=True, dtype=None, skip_header=1)
            if not len(data_atm.pco) == len(table['Pressure']):
                raise ValueError(
                    'The input profile has different layers than the current setting...'
                )
            data_atm.pco = table['Pressure']
            data_atm.Tco = table['Temp']

        else:
            raise ValueError(f'atm_type {self.type} not recongized.')

        # IF switches for Kzz types
        if self.Kzz_prof == 'const':
            data_atm.Kzz = np.repeat(self.const_Kzz, nz - 1)
        elif self.Kzz_prof == 'JM16':  # Kzz profiles assumed in Moses et al.2016
            data_atm.Kzz = 1e5 * (300.0 / (data_atm.pico[1:-1] * 1e-3)) ** 0.5
            data_atm.Kzz = np.maximum(self.cfg.K_deep, data_atm.Kzz)
        elif self.Kzz_prof == 'Pfunc':  # Kzz profiles assumed in Tsai 2020
            data_atm.Kzz = (
                self.cfg.K_max * (self.cfg.K_p_lev * 1e6 / (data_atm.pico[1:-1])) ** 0.4
            )
            data_atm.Kzz = np.maximum(self.cfg.K_max, data_atm.Kzz)

        elif self.Kzz_prof == 'file':
            pass  # already defined within atm_type = 'file
        else:
            raise IOError(
                '"Kzz_prof" (the type of Kzz profile) cannot be recongized.\nPlease assign it as "file" or "const" or "JM16" in vulcan_cfg.'
            )

        # IF switches for Vz types
        if self.vz_prof == 'const':
            data_atm.vz = np.repeat(self.const_vz, nz - 1)
        elif self.vz_prof == 'file':
            inter_vz = interpolate.interp1d(
                atm_table['Pressure'],
                atm_table['vz'],
                assume_sorted=False,
                bounds_error=False,
                fill_value=0,
            )
            data_atm.vz = inter_vz(data_atm.pico[1:-1])
        else:
            raise IOError(
                '"vz_prof" cannot be recongized.\nPlease assign it as "file" or "const" in vulcan_cfg.'
            )

        if not self.use_Kzz:
            # store Kzz in data_atm
            data_atm.Kzz = np.zeros(nz - 1)
        if not self.use_vz:
            data_atm.vz = np.zeros(nz - 1)

        # every rate coefficient and the number density depend on T: reject unphysical
        # profiles here rather than letting NaNs spread through the integration
        if not np.all(np.isfinite(data_atm.Tco) & (data_atm.Tco > 0)):
            raise ValueError(
                'Temperature must be finite and positive in every layer; got T from '
                f'{np.nanmin(data_atm.Tco)} to {np.nanmax(data_atm.Tco)} K '
                f'(atm_type = {self.type})'
            )

        # calculating and storing M(the third body)
        data_atm.M = data_atm.pco / (kb * data_atm.Tco)
        data_atm.n_0 = data_atm.M.copy()

        return data_atm

    # T(P) profile in Heng et al. 2014 (126)
    def TP_H14(self, pco, *args_analytical):
        """
        Analytical T-P profile of Heng et al. (2014), their eq. 126.

        The collision-induced absorption term is switched off. The shortwave scattering
        parameter sets the albedo, (1 - beta_S) / (1 + beta_S), which reduces T_irr by
        (1 - albedo)^(1/4). Gravity is `Config.gs` throughout.

        Parameters
        ----------
        pco : np.ndarray
            Pressure [dyn cm-2].
        *args_analytical : float
            T_int [K], T_irr [K], kappa_L and kappa_S (longwave and shortwave opacities)
            [cm2 g-1], beta_S and beta_L [dimensionless], as in `Config.para_anaTP`.

        Returns
        -------
        np.ndarray
            Temperature [K], same shape as `pco`.
        """

        # convert args_analytical tuple to a list so we can modify it
        T_int, T_irr, ka_0, ka_s, beta_s, beta_l = list(args_analytical)

        g = self.cfg.gs
        P_b = self.cfg.P_b

        # albedo(beta_s) also affects T_irr
        albedo = (1.0 - beta_s) / (1.0 + beta_s)
        T_irr *= (1 - albedo) ** 0.25
        eps_L = 3.0 / 8
        eps_L3 = 1.0 / 3
        ka_CIA = 0
        m = pco / g
        m_0 = P_b / g
        ka_l = ka_0 + ka_CIA * m / m_0
        term1 = (
            T_int**4
            / 4
            * (1 / eps_L + m / (eps_L3 * beta_l**2) * (ka_0 + ka_CIA * m / (2 * m_0)))
        )
        term2 = 1 / (2 * eps_L) + scipy.special.expn(2, ka_s * m / beta_s) * (
            ka_s / (ka_l * beta_s) - (ka_CIA) * m * beta_s / (eps_L3 * ka_s * m_0 * beta_l**2)
        )
        term3 = (
            ka_0
            * beta_s
            / (eps_L3 * ka_s * beta_l**2)
            * (1.0 / 3 - scipy.special.expn(4, ka_s * m / beta_s))
        )
        term4 = 0.0  # related to CIA
        T = (term1 + T_irr**4 / 8 * (term2 + term3 + term4)) ** 0.25

        return T

    def mol_mass(self, sp):
        """
        Molar mass of a species from the composition file `paths.COM_FILE` [g mol-1].
        """
        return compo['mass'][compo_row.index(sp)]

    def mean_mass(self, var, atm, ni):
        """
        Set the mean molar mass `atm.mu` [g mol-1] from the mixing ratios `var.ymix`.

        Parameters
        ----------
        var : store.Variables
            Reads `ymix`.
        atm : store.AtmData
            Sets `mu`, shape (nz,) [g mol-1].
        ni : int
            Number of species.

        Returns
        -------
        store.AtmData
            `atm`.
        """
        nz = self.cfg.nz
        mu = np.zeros(nz)
        for i in range(ni):
            mu += self.mol_mass(species[i]) * var.ymix[:, i]
        atm.mu = mu
        return atm

    def f_mu_dz(self, data_var, data_atm, output):
        """
        Set the mean molar mass, gravity, scale height and the height grid.

        Integrates hydrostatic balance layer by layer, dz = H_p ln(p_i / p_{i+1}) with
        H_p = k_B T N_A / (mu g) and g = g_s (R_p / (R_p + z))^2. For a gas giant with
        `P_b` >= 1 bar, z = 0 and g = g_s at the interface closest to 1 bar; otherwise at
        the bottom. With `use_settling`, also sets the Stokes settling velocity of the
        non-gaseous species, -2/9 rho_p r_p^2 g / eta, with the viscosity eta of the
        `atm_base` gas from Cloutman (Table 1); CO2 uses the N2 values. Plots the T-P
        profile when `plot_TP` is set and warns when T is outside 200-6000 K, the range of
        the thermodynamic data.

        Parameters
        ----------
        data_var : store.Variables
            Reads `ymix`.
        data_atm : store.AtmData
            Reads `Tco` [K] and `pico` [dyn cm-2]; sets `mu` [g mol-1], `g` and `gs`
            [cm s-2], `Hp`, `dz`, `dzi`, `zco`, `zmco` [cm], `pref_indx`, with
            `use_moldiff` `Ti` [K] and `Hpi` [cm] on the interfaces, and with
            `use_settling` `vs` [cm s-1].
        output : op.Output
            Not used.

        Returns
        -------
        store.AtmData
            `data_atm`.

        Raises
        ------
        ValueError
            If a non-gaseous species has no particle radius or density.
        """

        nz = self.cfg.nz
        dz, zco = data_atm.dz, data_atm.zco  # pressure defined at interfaces
        Tco, pico = data_atm.Tco.copy(), data_atm.pico.copy()
        Hp = data_atm.Hp

        if (
            not self.cfg.rocky and self.P_b >= 1e6
        ):  # if the lower BC greater than 1bar for gas giants
            # Find the index of pico closest to 1bar
            pref_indx = min(range(nz + 1), key=lambda i: abs(np.log10(pico[i]) - 6.0))
        else:
            pref_indx = 0

        # updating and storing mu
        data_atm = self.mean_mass(data_var, data_atm, ni)
        mu = data_atm.mu
        gs = self.gs
        gz = data_atm.g
        # updating and storing mu
        data_atm = self.mean_mass(data_var, data_atm, ni)

        for i in range(pref_indx, nz):
            if i == pref_indx:
                gz[i] = gs
                Hp[i] = kb * Tco[i] / (mu[i] / Navo * gs)
            else:
                gz[i] = gs * (self.cfg.Rp / (self.cfg.Rp + zco[i])) ** 2
                Hp[i] = kb * Tco[i] / (mu[i] / Navo * gz[i])
            dz[i] = Hp[i] * np.log(
                pico[i] / pico[i + 1]
            )  # pico[i+1] has a lower P than pico[i] (higer height)
            zco[i + 1] = zco[i] + dz[i]  # zco is set zero at 1bar for gas giants

        # for pref_indx != zero
        if not pref_indx == 0:
            for i in range(pref_indx - 1, -1, -1):
                gz[i] = gs * (self.cfg.Rp / (self.cfg.Rp + zco[i + 1])) ** 2
                Hp[i] = kb * Tco[i] / (mu[i] / Navo * gz[i])
                dz[i] = Hp[i] * np.log(pico[i] / pico[i + 1])
                zco[i] = zco[i + 1] - dz[i]  # from i+1 propogating down to i

        zmco = 0.5 * (zco + np.roll(zco, -1))
        zmco = zmco[:-1]
        dzi = 0.5 * (dz + np.roll(dz, 1))
        dzi = dzi[1:]
        # for the j grid, dzi[j] from the grid above and dz[j-1] from the grid below

        # for the molecular diffsuion
        if self.cfg.use_moldiff:
            Ti = 0.5 * (Tco + np.roll(Tco, -1))
            data_atm.Ti = Ti[:-1]
            Hpi = 0.5 * (Hp + np.roll(Hp, -1))
            data_atm.Hpi = Hpi[:-1]

        # updating and storing dz and dzi
        # data_atm.dz = dz
        data_atm.dzi = dzi
        # data_atm.zco = zco
        data_atm.zmco = zmco
        data_atm.g, data_atm.gs = gz, gs
        data_atm.pref_indx = pref_indx

        if self.use_settling:
            # TESTing settling velocity
            # based on L. D. Cloutman: A Database of Selected Transport Coefficients for Combustion Studies (Table 1.)
            if self.cfg.atm_base == 'N2':
                na = 1.52
                a = 1.186e-5
                b = 86.54
            elif self.cfg.atm_base == 'H2':
                na = 1.67
                a = 1.936e-6
                b = 2.187
            elif self.cfg.atm_base == 'CO2':
                log.warning('NO CO2 viscosity yet! (using N2 instead)')
                na = 1.52
                a = 1.186e-5
                b = 86.54
            elif self.cfg.atm_base == 'H2O':
                na = 1.5
                a = 1.6e-5
                b = 0
            elif self.cfg.atm_base == 'O2':
                na = 1.46
                a = 2.294e-5
                b = 164.4

            dmu = a * data_atm.Tco**na / (b + data_atm.Tco)  # g cm-1 s-1 dynamic viscosity

            for sp in self.cfg.non_gas_sp:
                try:
                    rho_p = data_atm.rho_p[sp]
                    r_p = data_atm.r_p[sp]

                except (ValueError, KeyError):
                    raise ValueError(sp + ' has not been prescribed size and density!')

                # Calculating the setteling (terminal) velocity
                gi = 0.5 * (data_atm.g + np.roll(data_atm.g, -1))
                gi = gi[:-1]
                data_atm.vs[:, species.index(sp)] = -1.0 * (
                    2.0 / 9 * rho_p * r_p**2 * gi / dmu[1:]
                )

        # plot T-P profile
        if self.cfg.plot_TP:
            plot.plot_TP(self.cfg, data_atm)

        # print warning when T exceeds the valid range of Gibbs free energy (NASA polynomials)
        if np.any(np.logical_or(data_atm.Tco < 200, data_atm.Tco > 6000)):
            log.warning('Temperatures exceed the valid range of Gibbs free energy.\n')

        return data_atm

    def read_sflux(self, var, atm):
        """
        Read the stellar spectrum and put the flux at the top of the atmosphere on the bins.

        The flux at the stellar surface from `Config.sflux_file` (wavelength [nm], flux
        [erg cm-2 s-1 nm-1]) is scaled by (R_star / a)^2 and interpolated linearly onto
        `var.bins`, with zero outside the file range. Logs the fraction of the energy in
        the bin range that the binned flux keeps, both integrated with the trapezoidal rule.

        Parameters
        ----------
        var : store.Variables
            Reads `bins` [nm]; sets `sflux_top` [erg cm-2 s-1 nm-1] and
            `sflux_din12_indx`, the first bin at `Config.dbin_12trans` (-1 if none).
        atm : store.AtmData
            Sets `sflux_raw`, the spectrum as read.
        """
        atm.sflux_raw = np.genfromtxt(
            self.cfg.sflux_file, dtype=float, skip_header=1, names=['lambda', 'flux']
        )

        # for values outside the boundary => fill_value = 0
        bins = var.bins

        dbin1 = self.cfg.dbin1
        dbin2 = self.cfg.dbin2

        inter_sflux = interpolate.interp1d(
            atm.sflux_raw['lambda'],
            atm.sflux_raw['flux']
            * (self.cfg.r_star * r_sun / (au * self.cfg.orbit_radius)) ** 2,
            bounds_error=False,
            fill_value=0,
        )
        for n, ld in enumerate(var.bins):
            var.sflux_top[n] = inter_sflux(ld)
            if ld == self.cfg.dbin_12trans:
                var.sflux_din12_indx = n
            # not converting to actinic flux yet *1/(hc/ld)

        # Check for energy conservation
        # finding the index for the left & right pts that match var.bins in the raw data
        raw_flux = (
            atm.sflux_raw['flux']
            * (self.cfg.r_star * r_sun / (au * self.cfg.orbit_radius)) ** 2
        )
        raw_left_indx = np.searchsorted(atm.sflux_raw['lambda'], bins[0], side='right')
        raw_right_indx = np.searchsorted(atm.sflux_raw['lambda'], bins[-1], side='right') - 1

        sum_orgin = 0
        # for checking the trapezoidal error in energy conservation
        for n in range(raw_left_indx, raw_right_indx):
            sum_orgin += (
                0.5
                * (raw_flux[n] + raw_flux[n + 1])
                * (atm.sflux_raw['lambda'][n + 1] - atm.sflux_raw['lambda'][n])
            )
        sum_orgin += (
            0.5
            * (inter_sflux(bins[0]) + raw_flux[raw_left_indx])
            * (atm.sflux_raw['lambda'][raw_left_indx] - bins[0])
        )
        sum_orgin += (
            0.5
            * (inter_sflux(bins[-1]) + raw_flux[raw_right_indx])
            * (bins[-1] - atm.sflux_raw['lambda'][raw_right_indx])
        )

        # dbin_12trans outside the bins
        if 'sflux_din12_indx' not in vars(var).keys():
            var.sflux_din12_indx = -1

        sum_bin = dbin1 * np.sum(var.sflux_top[: var.sflux_din12_indx])
        sum_bin -= dbin1 * 0.5 * (var.sflux_top[0] + var.sflux_top[var.sflux_din12_indx - 1])
        sum_bin += dbin2 * np.sum(var.sflux_top[var.sflux_din12_indx :])
        sum_bin -= dbin2 * 0.5 * (var.sflux_top[var.sflux_din12_indx] + var.sflux_top[-1])

        log.debug(
            'The stellar flux is interpolated onto uniform grid of '
            + str(self.cfg.dbin1)
            + ' (<'
            + str(self.cfg.dbin_12trans)
            + ' nm) and '
            + str(self.cfg.dbin2)
            + ' (>='
            + str(self.cfg.dbin_12trans)
            + ' nm)'
            + ' and conserving '
            + '{:.2f}'.format(100 * sum_bin / sum_orgin)
            + ' %'
            + ' energy.'
        )

    def mol_diff(self, atm):
        """
        Set the molar masses and the molecular diffusion coefficients.

        Without `use_moldiff`, only `ms` is set. Otherwise the coefficient of each species
        in the `atm_base` gas is a reference binary coefficient A T^s / n (H2: CH4-H2 from
        Moses 2000a; N2 and O2: CH4-N2 and CH4-O2 from Banks & Kockarts 1973; CO2: H2-CO2),
        scaled by the square root of the ratio of reduced masses, and evaluated with T and n
        averaged onto the interior interfaces. Non-gaseous species get zero. Also sets the
        thermal diffusion factors.

        Parameters
        ----------
        atm : store.AtmData
            Reads `Tco` [K] and `n_0` [cm-3]; sets `ms` [g mol-1], `Dzz` [cm2 s-1] and
            `alpha` [dimensionless].

        Raises
        ------
        IOError
            If `atm_base` is not H2, N2, O2 or CO2.
        """
        Tco = atm.Tco
        n_0 = atm.n_0

        # using the value defined on the interface
        Tco_i = np.delete((Tco + np.roll(Tco, 1)) * 0.5, 0)
        n0_i = np.delete((n_0 + np.roll(n_0, 1)) * 0.5, 0)

        if not self.cfg.use_moldiff:
            for i in range(len(species)):
                # this is required even without molecular weight
                atm.ms[i] = compo[compo_row.index(species[i])][-1]
            return

        if self.cfg.atm_base == 'H2':
            Dzz_gen = lambda T, n_tot, mi: (
                2.2965e17 * T**0.765 / n_tot * (16.04 / mi * (mi + 2.016) / 18.059) ** 0.5
            )  # from Moses 2000a

            # scaling with (15.27) in [Aeronomy part B by Banks & Kockarts(1973)]
            # *( m_ref/mi*(mi+ m_base)/m_tot )**0.5 (m_ref is the molecular mass of the ref-base e.g. CH4 in CH4-H2 in Moses 2000a)

            # # thermal diffusion factor (>0 means (heavier) components diffuse toward colder rigions)
            if 'H' in species:
                atm.alpha[species.index('H')] = -0.1  # simplified from Moses 2000a
            if 'He' in species:
                atm.alpha[species.index('He')] = 0.145
            for sp in species:
                if self.mol_mass(sp) > 4.0:
                    atm.alpha[species.index(sp)] = 0.25

        elif (
            self.cfg.atm_base == 'N2'
        ):  # use CH4-N2 in Aeronomy [Banks ] as a reference to scale by the molecular mass
            Dzz_gen = lambda T, n_tot, mi: (
                7.34e16 * T**0.75 / n_tot * (16.04 / mi * (mi + 28.014) / 44.054) ** 0.5
            )

            # # thermal diffusion factor (>0 means (heavier) components diffuse toward colder rigions)
            if 'H' in species:
                atm.alpha[species.index('H')] = -0.25
            if 'H2' in species:
                atm.alpha[species.index('H2')] = -0.25
            if 'He' in species:
                atm.alpha[species.index('He')] = -0.25
            if 'Ar' in species:
                atm.alpha[species.index('Ar')] = 0.17

        elif (
            self.cfg.atm_base == 'O2'
        ):  # use CH4-O2 in Aeronomy [Banks ] as a reference to scale by the molecular mass
            Dzz_gen = lambda T, n_tot, mi: (
                7.51e16 * T**0.759 / n_tot * (16.04 / mi * (mi + 32) / 48.04) ** 0.5
            )

            # # thermal diffusion factor (>0 means (heavier) components diffuse toward colder rigions)
            if 'H' in species:
                atm.alpha[species.index('H')] = -0.25
            if 'H2' in species:
                atm.alpha[species.index('H2')] = -0.25
            if 'He' in species:
                atm.alpha[species.index('He')] = -0.25
            if 'Ar' in species:
                atm.alpha[species.index('Ar')] = 0.17

        elif (
            self.cfg.atm_base == 'CO2'
        ):  # use H2-CO2 in Hu seager as a reference to scale by the molecular mass
            Dzz_gen = lambda T, n_tot, mi: (
                2.15e17 * T**0.750 / n_tot * (2.016 / mi * (mi + 44.001) / 46.017) ** 0.5
            )

            # # thermal diffusion factor (>0 means (heavier) components diffuse toward colder rigions)
            if 'H' in species:
                atm.alpha[species.index('H')] = -0.25
            if 'H2' in species:
                atm.alpha[species.index('H2')] = -0.25
            if 'He' in species:
                atm.alpha[species.index('He')] = -0.25
            if 'Ar' in species:
                atm.alpha[species.index('Ar')] = 0.17

        else:
            raise IOError('\n Unknow atm_base!')

        for i in range(len(species)):
            # input should be float or in the form of nz-long 1D array
            atm.Dzz[:, i] = Dzz_gen(Tco_i, n0_i, self.mol_mass(species[i]))
            # constructing the molecular weight for every species
            # this is required even without molecular weight
            atm.ms[i] = compo[compo_row.index(species[i])][-1]

        # setting the molecuar diffusion of the non-gaseous species to zero
        for sp in [_ for _ in self.cfg.non_gas_sp if _ in species]:
            atm.Dzz[:, species.index(sp)] = 0

    def BC_flux(self, atm):
        """
        Read the boundary conditions from the files named in the configuration.

        With `use_topflux`, `top_BC_flux_file` gives a species and its flux per line. With
        `use_botflux`, `bot_BC_flux_file` gives a species, its flux and its deposition
        velocity; with `use_fix_sp_bot`, the fourth column of that file is the fixed
        mixing ratio.

        Parameters
        ----------
        atm : store.AtmData
            Sets `top_flux` and `bot_flux` [cm-2 s-1] (positive into the atmosphere),
            `bot_vdep` [cm s-1] and `bot_fix_sp` [dimensionless].
        """
        # read in the const top BC
        if self.cfg.use_topflux:
            log.info('Using the prescribed constant top flux.')
            with open(self.cfg.top_BC_flux_file) as f:
                for line in f.readlines():
                    if not line.startswith('#') and line.strip():
                        li = line.split()
                        atm.top_flux[species.index(li[0])] = li[1]

        # read in the const bottom BC
        if self.cfg.use_botflux:
            log.info('Using the prescribed constant bottom flux.')
            with open(self.cfg.bot_BC_flux_file) as f:
                for line in f.readlines():
                    if not line.startswith('#') and line.strip():
                        li = line.split()
                        atm.bot_flux[species.index(li[0])] = li[1]
                        atm.bot_vdep[species.index(li[0])] = li[2]

        # using fixed-mixing-ratio BC
        if self.cfg.use_fix_sp_bot:
            log.info('Using the prescribed fixed bottom mixing ratios.')
            with open(self.cfg.bot_BC_flux_file) as f:
                for line in f.readlines():
                    if not line.startswith('#') and line.strip():
                        li = line.split()
                        atm.bot_fix_sp[species.index(li[0])] = li[3]

    # TEST condensation
    def sp_sat(self, atm):
        """
        Set the saturation vapour pressure of each species in `Config.condense_sp`.

        Fits by species: H2O over ice and water (Ackerman & Marley 2001), NH3 (Weast 1971),
        H2SO4 (Ayers), S2 and S8 (Zahnle 2017, refitted from Lyons 2008), C (NIST) and
        H2S over ice and liquid (Giauque & Blue 1936).

        Parameters
        ----------
        atm : store.AtmData
            Reads `Tco` [K]; sets `sat_p[sp]`, shape (nz,) [dyn cm-2].

        Raises
        ------
        IOError
            If a species has no saturation data.
        """
        # the list that the data has been coded
        sat_sp_list = ['H2O', 'NH3', 'H2SO4', 'S2', 'S8', 'C', 'H2S']

        nz = self.cfg.nz

        for sp in self.cfg.condense_sp:
            if sp not in sat_sp_list:
                raise IOError(
                    'No saturation vapor data for '
                    + sp
                    + '. Check the sp_sat function in build_atm.py'
                )

            T = np.copy(atm.Tco)

            if sp == 'H2O':
                # T is in C
                T -= 273.0
                # from Seinfeld & Pandis 2006, P in mbar in the book
                # a_water = (
                #     6.107799961,
                #     4.436518521e-1,
                #     1.428945805e-2,
                #     2.650648471e-4,
                #     3.031240396e-6,
                #     2.034080948e-8,
                #     6.136820929e-11,
                # )
                # a_ice = (
                #     6.109177956,
                #     5.034698970e-1,
                #     1.886013408e-2,
                #     4.176223716e-4,
                #     5.824720280e-6,
                #     4.838803174e-8,
                #     1.838826904e-10,
                # )

                # saturate_p_1 = (T<0)*( a_ice[0] + a_ice[1]*T + a_ice[2]*T**2 + a_ice[3]*T**3 + a_ice[4]*T**4 + a_ice[5]*T**5 + a_ice[6]*T**6 ) +\
                #  (T>0)*(a_water[0] + a_water[1]*T + a_water[2]*T**2 + a_water[3]*T**3 + a_water[4]*T**4 + a_water[5]*T**5 + a_water[6]*T**6)

                # ice from Ackerman&Marley (2001)
                c0 = 6111.5
                c1 = 23.036
                c2 = -333.7
                c3 = 279.82
                # water from Ackerman&Marley (2001)
                w0 = 6112.1
                w1 = 18.729
                w2 = -227.3
                w3 = 257.87

                saturate_p = (T < 0) * (c0 * np.exp((c1 * T + T**2 / c2) / (T + c3))) + (
                    T > 0
                ) * (w0 * np.exp((w1 * T + T**2 / w2) / (T + w3)))

                atm.sat_p[sp] = saturate_p
                # atm.sat_p[sp] = c0 * np.exp( (c1*T + T**2/c2)/(T + c3) )
                # atm.sat_p[sp] = moses_ice

            elif sp == 'NH3':
                # from Weast (1971) in bar
                c0 = 10.53
                c1 = -2161.0
                c2 = -86596.0
                saturate_p = np.exp(c0 + c1 / T + c2 / T**2)
                atm.sat_p[sp] = saturate_p * 1.0e6

            elif sp == 'H2SO4':
                # change to Kulmala later
                p_ayers = np.e ** (-10156.0 / T + 16.259)  # in atm
                atm.sat_p[sp] = p_ayers * 1.01325 * 1e6  # arm to cgs

            elif sp == 'S2':
                atm.sat_p[sp] = np.zeros(nz)
                # from Zahnle 2017 (refitted from Lyons 2008)
                atm.sat_p[sp][T < 413] = (
                    np.exp(27.0 - 18500.0 / T[T < 413]) * 1e6
                )  # in bar => cgs
                atm.sat_p[sp][T >= 413] = np.exp(16.1 - 14000.0 / T[T >= 413]) * 1e6

            elif sp == 'S4':
                atm.sat_p[sp] = np.zeros(nz)
                # from Lyons 2008
                atm.sat_p[sp] = 10 ** (6.0028 - 6047.5 / T) * 1.01325e6  # atm to vgs

            elif sp == 'S8':
                atm.sat_p[sp] = np.zeros(nz)
                # from Zahnle 2017 (refitted from Lyons 2008)
                atm.sat_p[sp][T < 413] = (
                    np.exp(20.0 - 11800.0 / T[T < 413]) * 1e6
                )  # in bar => cgs
                atm.sat_p[sp][T >= 413] = np.exp(9.6 - 7510.0 / T[T >= 413]) * 1e6

            elif sp == 'C':
                atm.sat_p[sp] = np.zeros(nz)
                # from NIST (dyna cm^-2)
                a = 3.27860e01
                b = -8.65139e04
                c = 4.80395e-01
                atm.sat_p[sp] = np.exp(a + b / (atm.Tco + c))

            elif sp == 'H2S':
                # from Giauque and Blue(1936) in cmHg (adoped in Atreya's book)
                h2s_ice_log10 = -1329.0 / T + 9.28588 - 0.0051263 * T  # 164.9 <= T <= 187.6
                h2s_l_log10 = -1145.0 / T + 7.94746 - 0.00322 * T  # 187.6 < T <= 213.2

                saturate_p = 10 ** ((T <= 187.6) * h2s_ice_log10 + (T > 187.6) * h2s_l_log10)
                atm.sat_p[sp] = saturate_p * 0.001333 * 1.0e6
