# Variables and functions to help with plotting functions
from __future__ import annotations

import logging
import os
import re
import shutil

import matplotlib.pyplot as plt
import numpy as np

log = logging.getLogger('fwl.' + __name__)


def mol_to_ele(mol: str):
    """
    Count the atoms of each element in a molecule name.

    Anything other than element symbols and counts, such as the '_l_s' of a condensate or
    a charge sign, is ignored.

    Parameters
    ----------
    mol : str
        Molecule name, for example 'H2O' or 'CH3OH'.

    Returns
    -------
    dict of str to int
        Number of atoms of each element, for example {'H': 2, 'O': 1}.

    Raises
    ------
    ValueError
        If the name does not start with an upper-case letter or contains no element.
    """

    # Validate
    if not str(mol[0]).isupper():
        raise ValueError(f"Molecule name '{mol}' is invalid")

    # Get atoms
    decomp = re.findall(r'([A-Z][a-z]?)(\d*)', mol)
    elems = {}
    for ev in decomp:
        if ev[1] == '':
            val = 1
        else:
            val = int(ev[1])
        elems[str(ev[0])] = elems.get(ev[0], 0) + val

    # Check that what we got is reasonable
    if not elems:
        raise ValueError(f"Could not decompose molecule '{mol}'")

    return elems


# Standard plotting colours
_preset_colours = {
    # Default colour
    '_fallback': '#ff00ff',
    # Volatile gases
    'H2O': '#027FB1',
    'CO2': '#D24901',
    'H2': '#008C01',
    'CH4': '#C720DD',
    'CO': '#D1AC02',
    'N2': '#870036',
    'S2': '#FF8FA1',
    'SO2': '#00008B',
    'H2S': '#2eff90',
    'NH3': '#675200',
    'SiO': '#b2df8a',
    'MgO': '#a6cee3',
    'SiO2': '#1b9e77',
    'FeO': '#7f0000',
    'TiO': '#9C45BC',
    'TiO2': '#421e39',
    # Volatile elements
    'H': '#0000cc',
    'C': '#ff0000',
    'O': '#88bb22',
    'N': '#ffaa00',
    'S': '#ff22ff',
    'P': '#33ccff',
    # refractory elements
    'Fe': '#888888',
    'Si': '#aa2277',
    'Mg': '#996633',
    'Na': '#ccff00',
    'Ti': '#0d2959',
    'V': '#ff6600',
    'K': '#2e8b57',
    'Al': '#5e2172',
    'Ca': '#167A9B',
    # Noble gases (each is both a gas species and an element)
    'He': '#7f7f7f',
    'Ne': '#e6550d',
    'Ar': '#3182bd',
    'Kr': '#31a354',
    'Xe': '#756bb1',
    # GLobal energy fluxes
    'OLR': '#dc143c',
    'ASF': '#4169e1',
    'sct': '#2e8b57',
    'tidal': '#daa520',
    'radio': '#C720DD',
    'star': '#FF8FA1',
    # Atmosphere energy fluxes (copied from AGNI)
    'flux_r': '#c0c0c0',
    'flux_n': '#000000',
    'flux_c': '#6495ed',
    'flux_t': '#ff4400',
    'flux_o': '#66CD00',
    'flux_p': '#ecb000',
    # Model components
    'atm': '#444444',
    'int': '#ff7f0e',
    'core': '#aa3077',
    'atm_bkg': '#f2faff',
    'int_bkg': '#fffaf2',
    'cor_bkg': '#efefef',
    # Clouds and aerosols
    'cloud': '#027FB1',
    'frsoot': '#876e44',
    'soot': '#63401a',
    'agsoot': '#251607',
    'nitrate': '#ffaa00',
    'ash': '#970505',
    'sulph': '#ff22ff',
    'biogenic': '#33ccff',
    'bioms1': '#33ccaa',
    'delta': '#6a9064',
    'dustdiv1': '#fcfcd7',
    'dustdiv2': '#fbfca1',
    'dustdiv3': '#f6f873',
    'dustdiv4': '#fcff5d',
    'dustdiv5': '#fcff33',
    'dustdiv6': '#fbff21',
    'naclflm': '#fca1ff',
    'nacljet': '#ee29f5',
}


def _generate_colour(gas: str):
    """
    Make a colour for a gas by mixing the preset colours of its elements.

    The element colours are weighted by the number of atoms, then brightened and capped so
    that the colour stays visible on white. The method follows AGNI's phys.jl.

    Parameters
    ----------
    gas : str
        Molecule name.

    Returns
    -------
    str
        Colour as '#rrggbb'; the fallback magenta if the name cannot be split into
        elements or an element has no preset colour.
    """

    # Break into atoms
    try:
        atoms = mol_to_ele(gas)
    except ValueError:
        log.warning(f"Using fallback colour for '{gas}'")
        return _preset_colours['_fallback']

    # Red, green, blue components
    red = 0.0
    gre = 0.0
    blu = 0.0

    # For each atom, add the contriution of its rgb components
    for e in atoms.keys():
        if e not in _preset_colours.keys():
            log.warning(f"Using fallback colour for '{gas}'")
            return _preset_colours['_fallback']
        red += int(_preset_colours[e][1:2], base=16) * atoms[e]
        gre += int(_preset_colours[e][3:4], base=16) * atoms[e]
        blu += int(_preset_colours[e][5:6], base=16) * atoms[e]

    # Normalisation constant
    norm = max((red, gre, blu))

    # Normalise colours to 0-255
    red = int(255 * red / norm)
    gre = int(255 * gre / norm)
    blu = int(255 * blu / norm)

    # Prevent the colour getting too close to white, which is hard to see on a plot
    c_lim = 220
    if red + gre + blu > c_lim * 3:
        scale = c_lim / 255.0
        red = int(scale * red)
        gre = int(scale * gre)
        blu = int(scale * blu)

    # Convert to hex string
    colour = f'#{red:02x}{gre:02x}{blu:02x}'

    return colour


def get_colour(thing: str):
    """
    Colour for a gas, element, flux or model component in the plots.

    Parameters
    ----------
    thing : str
        Name to look up in the preset colours, for example 'H2O', 'OLR' or 'atm'.

    Returns
    -------
    str
        Colour as '#rrggbb': the preset colour if there is one, otherwise one mixed from
        the elements of the name (`_generate_colour`).
    """

    # Try getting it from dictionary
    try:
        colour = _preset_colours[thing]

    # Otherwise, generate it systematically
    except KeyError:
        colour = _generate_colour(thing)

    return colour


def latexify(gas: str):
    """
    Write a gas name with its digits as LaTeX subscripts, for plot labels.

    Each digit becomes its own subscript, so 'H2O' gives 'H$_2$O' and 'C10' gives
    'C$_1$$_0$'.

    Parameters
    ----------
    gas : str
        Gas name.

    Returns
    -------
    str
        The label for matplotlib.
    """

    out = ''
    for c in gas:
        c = str(c)
        if c.isnumeric():
            out += r'$_%s$' % c
        else:
            out += c
    return out


def plot_update(cfg, species, var, atm, para, fpath=None):
    """
    Plot the mixing ratios of `Config.plot_spec` and the temperature at the current step.

    Mixing ratios are drawn against pressure [bar] or, with `Config.plot_height`, height
    [km], with the saturation mixing ratios dashed and the temperature [K] dotted on a
    second x axis.

    Parameters
    ----------
    cfg : Config
        Run configuration.
    species : list of str
        Species names, the columns of `var.ymix`.
    var : store.Variables
        Reads `ymix`, `t` [s] and `longdydt` [s-1].
    atm : store.AtmData
        Reads `pco`, `zco`, `zmco`, `Tco` and, with condensation, `sat_mix`.
    para : store.Parameters
        Reads `count`; increments `pic_count`.
    fpath : str, optional
        File to save to. If None, the plot goes to `_recent.<plot_fmt>` in `plot_dir` and
        is copied to `<pic_count>.<plot_fmt>` (five digits) for the live-plot frames.
    """

    log.debug('Plotting mixing ratios')

    colors = [
        'b',
        'r',
        'c',
        'm',
        'y',
        'k',
        'orange',
        'pink',
        'grey',
        'darkred',
        'darkblue',
        'salmon',
        'chocolate',
        'mediumspringgreen',
        'steelblue',
        'plum',
        'hotpink',
    ]

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    color_index = 0
    for sp in cfg.plot_spec:
        sp_lab = latexify(sp)

        if color_index == len(colors):  # when running out of colors
            colors.append(tuple(np.random.rand(3)))

        color = get_colour(sp)

        if not cfg.plot_height:
            (line,) = ax.plot(
                var.ymix[:, species.index(sp)], atm.pco / 1.0e6, color=color, label=sp_lab
            )
            if cfg.use_condense and sp in cfg.condense_sp:
                ax.plot(
                    atm.sat_mix[sp],
                    atm.pco / 1.0e6,
                    color=color,
                    label=sp_lab + ' sat',
                    ls='--',
                )

            ax.set_yscale('log')
            ax.invert_yaxis()
            ax.set_ylabel('Pressure [bar]')
            ax.set_ylim((cfg.P_b / 1.0e6, cfg.P_t / 1.0e6))

        else:  # plotting with height
            (line,) = ax.plot(
                var.ymix[:, species.index(sp)], atm.zmco / 1.0e5, color=color, label=sp_lab
            )
            if cfg.use_condense and sp in cfg.condense_sp:
                ax.plot(
                    atm.sat_mix[sp],
                    atm.zco[1:] / 1.0e5,
                    color=color,
                    label=sp_lab + ' sat',
                    ls='--',
                )

            ax.set_ylim((atm.zco[0] / 1e5, atm.zco[-1] / 1e5))
            ax.set_ylabel('Height [km]')

    title = 'i = %5d steps     t = %.2e seconds' % (para.count, var.t)
    title += '\n dy/dt = %.2e' % (var.longdydt)
    ax.set_title(title)
    ax.set_xscale('log')
    ax.set_xlabel('Volume mixing ratio')
    ax.set_xlim(1e-16, 1.2)
    ax.legend(fontsize=10, labelspacing=0.2, loc='upper left', bbox_to_anchor=(1.0, 1.0))

    # temperature profile
    axt = ax.twiny()
    axt.set_xlabel('Temperature [K]')
    if cfg.plot_height:
        axt_yarr = atm.zmco / 1.0e5
    else:
        axt_yarr = atm.pco / 1.0e6
    axt.plot(atm.Tco, axt_yarr, color='k', lw=0.5, ls='dotted')

    if fpath is None:
        save_fpath = os.path.join(cfg.plot_dir, f'_recent.{cfg.plot_fmt}')
        copy_fpath = os.path.join(cfg.plot_dir, f'{para.pic_count:05d}.{cfg.plot_fmt}')
    else:
        save_fpath = fpath
        copy_fpath = None

    log.debug(f'Plotting to {save_fpath}')
    fig.savefig(save_fpath, dpi=cfg.plot_dpi, bbox_inches='tight')
    plt.close('all')

    if copy_fpath:
        shutil.copyfile(save_fpath, copy_fpath)

    para.pic_count += 1


def plot_evo(cfg, species, var, atm, plot_j=-1, plot_ymin=1e-16, dn=1):
    """
    Plot the mixing ratios of `Config.plot_spec` in one layer against time.

    The mixing ratios are the saved number densities divided by the present `atm.n_0`.
    The figure is saved as `evo.<plot_fmt>`; the name is appended to `plot_dir`
    directly, so `plot_dir` must end with a path separator.

    Parameters
    ----------
    cfg : Config
        Run configuration; reads `plot_spec`, `plot_dir`, `plot_fmt` and `plot_dpi`.
    species : list of str
        Species names, the last axis of `var.y_time`.
    var : store.Variables
        Reads `y_time` [cm-3], an array of shape (number of saved steps, nz, ni), and
        `t_time` [s].
    atm : store.AtmData
        Reads `n_0` [cm-3].
    plot_j : int, optional
        Layer to plot; the default -1 is the top layer.
    plot_ymin : float, optional
        Lower limit of the mixing-ratio axis [dimensionless].
    dn : int, optional
        Plot every `dn`-th saved step.
    """

    plot_spec = cfg.plot_spec
    plot_dir = cfg.plot_dir
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    ymix_time = np.array(var.y_time / atm.n_0[:, np.newaxis])

    for i, sp in enumerate(cfg.plot_spec):
        col = get_colour(sp)
        lbl = latexify(sp)
        ax.plot(var.t_time[::dn], ymix_time[::dn, plot_j, species.index(sp)], c=col, label=lbl)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Volume mixing ratios')
    ax.set_ylim((plot_ymin, 1.2))
    ax.legend(fontsize=10, labelspacing=0.2, loc='upper left', bbox_to_anchor=(1.0, 1.0))
    fig.savefig(plot_dir + f'evo.{cfg.plot_fmt}', dpi=cfg.plot_dpi)
    plt.close('all')


def plot_evo_inter(cfg, species, var, atm, plot_j=-1, plot_ymin=1e-16, dn=1):
    """
    Plot the mixing ratios of `Config.plot_spec` in one layer against time, mid-run.

    Not called by the current code.

    Like `plot_evo`, but for a run that was interrupted; it converts `var.t_time` to an
    array in place and overwrites `evo.<plot_fmt>` in `plot_dir`.

    Parameters
    ----------
    cfg : Config
        Run configuration; reads `plot_spec`, `plot_dir`, `plot_fmt` and `plot_dpi`.
    species : list of str
        Species names, the last axis of `var.y_time`.
    var : store.Variables
        Reads `y_time` [cm-3] and `t_time` [s].
    atm : store.AtmData
        Reads `n_0` [cm-3].
    plot_j : int, optional
        Layer to plot; the default -1 is the top layer.
    plot_ymin : float, optional
        Lower limit of the mixing-ratio axis [dimensionless].
    dn : int, optional
        Plot every `dn`-th saved step.
    """
    var.t_time = np.array(var.t_time)
    ymix_time = np.array(var.y_time / atm.n_0[:, np.newaxis])

    plot_spec = cfg.plot_spec
    plot_dir = cfg.plot_dir
    plt.figure('evolution')

    for i, sp in enumerate(cfg.plot_spec):
        plt.plot(
            var.t_time[::dn],
            ymix_time[::dn, plot_j, species.index(sp)],
            c=plt.cm.rainbow(float(i) / len(plot_spec)),
            label=sp,
        )

    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.xlabel('time')
    plt.ylabel('mixing ratios')
    plt.ylim((plot_ymin, 1.0))
    plt.legend(frameon=0, prop={'size': 14}, loc='best')
    plt.savefig(plot_dir + f'evo.{cfg.plot_fmt}', dpi=cfg.plot_dpi)
    plt.close('all')


def plot_TP(cfg, atm):
    """
    Plot the initial temperature and Kzz profiles.

    Temperature [K] is drawn as a solid line and Kzz [cm2 s-1] as a dashed line on a
    second x axis, against pressure [bar] or, with `Config.plot_height`, height [km].
    The figure is saved as `TPK_initial.<plot_fmt>` in `Config.plot_dir`; the name is
    appended to `plot_dir` directly, so `plot_dir` must end with a path separator.

    Parameters
    ----------
    cfg : Config
        Run configuration.
    atm : store.AtmData
        Atmosphere data; reads `Tco`, `Kzz`, `pco`, `pico` and `zmco`.
    """
    plot_dir = cfg.plot_dir
    # plt.figure('TPK')
    fig, ax1 = plt.subplots()
    ax2 = ax1.twiny()  # ax1 and ax2 share y-axis

    if not cfg.plot_height:
        ax1.semilogy(atm.Tco, atm.pco / 1.0e6, c='black')
        ax2.loglog(atm.Kzz, atm.pico[1:-1] / 1.0e6, c='k', ls='--')
        plt.gca().invert_yaxis()
        plt.ylim((cfg.P_b / 1.0e6, cfg.P_t / 1.0e6))
        ax1.set_ylabel('Pressure (bar)')

    else:  # plotting with height
        ax1.plot(atm.Tco, atm.zmco / 1.0e5, c='black')
        ax2.semilogx(atm.Kzz, atm.zmco[1:] / 1.0e5, c='k', ls='--')
        ax1.set_ylabel('Height (km)')

    # plt.xlabel("Temperature (K)")
    ax1.set_xlabel('Temperature (K)')
    ax2.set_xlabel(r'K$_{zz}$ (cm$^2$s$^{-1}$)')

    fig.savefig(plot_dir + f'TPK_initial.{cfg.plot_fmt}', dpi=cfg.plot_dpi)
    plt.close('all')
