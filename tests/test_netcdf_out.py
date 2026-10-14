from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

nc = pytest.importorskip('netCDF4')

from vulcan.netcdf_out import FILL, compute_reaction_rates, write_netcdf  # noqa: E402

# Synthetic network: 2A + M <-> B + M, and photolysis B -> C (reverse switched off)
SPECIES = ['A', 'B', 'C']
NR = 4
RE_WM_DICT = {
    1: [['A', 'A', 'M'], ['B', 'M']],
    2: [['B', 'M'], ['A', 'A', 'M']],
    3: [['B'], ['C']],
    4: [['C'], ['B']],
}
NZ = 4
NBIN = 5


def _state(use_photo=True, save_evolution=False):
    """Build minimal var, atm, para and cfg containers for a 4-layer atmosphere."""
    rng = np.random.default_rng(42)
    pco = np.logspace(6, 2, NZ)
    Tco = np.linspace(1500.0, 800.0, NZ)
    M = pco / (1.380649e-16 * Tco)
    y = rng.uniform(0.1, 1.0, (NZ, len(SPECIES))) * M[:, None]

    k = {
        1: np.full(NZ, 1e-32),
        2: np.full(NZ, 1e-12),
        3: np.linspace(1e-6, 1e-4, NZ),
        4: np.zeros(NZ),
    }

    var = SimpleNamespace(
        y=y,
        ymix=y / y.sum(axis=1, keepdims=True),
        y_ini=y * 0.9,
        k=k,
        t=1e10,
        dt=1e8,
        longdy=1e-3,
        longdydt=1e-15,
        atom_ini={'H': 1.0, 'O': 2.0},
        atom_sum={'H': 1.0, 'O': 2.0},
        atom_loss={'H': 0.0, 'O': 0.0},
        atom_conden={'H': 0.0, 'O': 0.0},
        pho_rate_index={('B', 1): 3},
        ion_rate_index={},
        conden_re_list=[],
        bins=np.linspace(10.0, 300.0, NBIN),
        aflux=np.ones((NZ, NBIN)),
        sflux_top=np.ones(NBIN),
        sflux=np.ones((NZ + 1, NBIN)),
        dflux_u=np.zeros((NZ + 1, NBIN)),
        dflux_d=np.zeros((NZ + 1, NBIN)),
        tau=np.zeros((NZ + 1, NBIN)),
        t_time=np.array([0.0, 1.0, 2.0]),
        y_time=np.stack([y, y, y]),
    )
    atm = SimpleNamespace(
        pco=pco,
        pico=np.logspace(6.2, 1.8, NZ + 1),
        Tco=Tco,
        zmco=np.linspace(0.0, 3e7, NZ),
        zco=np.linspace(-5e6, 3.5e7, NZ + 1),
        Kzz=np.full(NZ - 1, 1e9),
        vz=np.zeros(NZ - 1),
        mu=np.full(NZ, 2.3),
        g=np.full(NZ, 2140.0),
        Hp=np.full(NZ, 1e7),
        n_0=M.copy(),
        M=M,
        ms=np.array([1.0, 2.0, 2.0]),
    )
    para = SimpleNamespace(end_case=1, count=1234, start_time=None)
    cfg = SimpleNamespace(
        nz=NZ,
        atom_list=['H', 'O'],
        network='synthetic',
        use_photo=use_photo,
        save_evolution=save_evolution,
        save_evo_frq=10,
        Rp=7e9,
        gs=2140.0,
        rocky=False,
        agni_call_frq=0,
    )
    return var, atm, para, cfg


@pytest.mark.unit
def test_reaction_rates_match_hand_computed_mass_action():
    var, atm, _, _ = _state()
    rates = compute_reaction_rates(var.k, var.y, atm.M, RE_WM_DICT, SPECIES, NR)

    A, B, C = var.y[:, 0], var.y[:, 1], var.y[:, 2]
    assert rates.shape == (NZ, NR)
    np.testing.assert_allclose(rates[:, 0], var.k[1] * A**2 * atm.M, rtol=1e-14)
    np.testing.assert_allclose(rates[:, 1], var.k[2] * B * atm.M, rtol=1e-14)
    np.testing.assert_allclose(rates[:, 2], var.k[3] * B, rtol=1e-14)
    np.testing.assert_array_equal(rates[:, 3], 0.0)
    # cm-3 s-1 must not depend on the k dict being mutated
    assert var.k[1][0] == 1e-32


@pytest.mark.unit
def test_reaction_rates_reproduce_net_chemical_tendency():
    """Sum of stoichiometry x rate must equal an independently coded dn/dt."""
    var, atm, _, _ = _state()
    rates = compute_reaction_rates(var.k, var.y, atm.M, RE_WM_DICT, SPECIES, NR)

    nu = np.zeros((NR, len(SPECIES)))
    for i in range(1, NR + 1):
        for sp in RE_WM_DICT[i][0]:
            if sp != 'M':
                nu[i - 1, SPECIES.index(sp)] -= 1
        for sp in RE_WM_DICT[i][1]:
            if sp != 'M':
                nu[i - 1, SPECIES.index(sp)] += 1
    dndt = rates @ nu

    A, B, C = var.y[:, 0], var.y[:, 1], var.y[:, 2]
    net1 = var.k[1] * A**2 * atm.M - var.k[2] * B * atm.M
    net3 = var.k[3] * B - var.k[4] * C
    expected = np.stack([-2 * net1, net1 - net3, net3], axis=1)
    np.testing.assert_allclose(dndt, expected, rtol=1e-12)
    # element (mass) conservation: 2 m_A = m_B = m_C in this toy network
    np.testing.assert_allclose(dndt @ atm.ms, 0.0, atol=1e-12 * np.abs(expected).max())


@pytest.mark.unit
def test_write_netcdf_roundtrip(tmp_path):
    var, atm, para, cfg = _state(use_photo=True, save_evolution=True)
    fpath = tmp_path / 'out.nc'
    write_netcdf(str(fpath), var, atm, para, cfg, SPECIES, NR, RE_WM_DICT)

    with nc.Dataset(fpath) as ds:
        assert ds.dimensions['layer'].size == NZ
        assert ds.dimensions['level'].size == NZ + 1
        assert ds.dimensions['reaction'].size == NR
        assert ds.dimensions['wavelength'].size == NBIN
        assert ds.dimensions['time'].size == 3

        assert list(ds['species'][:]) == SPECIES
        assert ds['reaction_equation'][0] == 'A + A + M -> B + M'
        assert ds['reaction_equation'][1] == 'B + M -> A + A + M'
        assert list(ds['reaction_type'][:]) == [
            'thermal',
            'thermal',
            'photolysis',
            'photolysis',
        ]
        np.testing.assert_array_equal(ds['reaction_is_reverse'][:], [0, 1, 0, 1])
        np.testing.assert_array_equal(ds['reaction_order'][:], [3, 2, 1, 1])

        np.testing.assert_allclose(ds['number_density'][:], var.y)
        np.testing.assert_allclose(ds['mixing_ratio'][:].sum(axis=1), 1.0, rtol=1e-12)
        np.testing.assert_allclose(ds['rate_coefficient'][:, 2], var.k[3])
        np.testing.assert_allclose(
            ds['reaction_rate'][:],
            compute_reaction_rates(var.k, var.y, atm.M, RE_WM_DICT, SPECIES, NR),
        )
        assert ds['pressure'].units == 'dyn cm-2'
        assert ds['actinic_flux'].units == 'photons cm-2 s-1 nm-1'

        # Kzz is defined only on interior interfaces
        kzz = ds['Kzz'][:]
        assert np.ma.is_masked(kzz[0]) and np.ma.is_masked(kzz[-1])
        np.testing.assert_allclose(kzz[1:-1], atm.Kzz)

        assert ds.termination == 'converged to steady state'
        assert ds.steps == 1234
        assert ds.rocky == 0
        assert 'wall_time_s' not in ds.ncattrs()
        assert '"network": "synthetic"' in ds.config


@pytest.mark.unit
def test_write_netcdf_without_photo_or_evolution(tmp_path):
    var, atm, para, cfg = _state(use_photo=False, save_evolution=False)
    fpath = tmp_path / 'out.nc'
    write_netcdf(str(fpath), var, atm, para, cfg, SPECIES, NR, RE_WM_DICT)

    with nc.Dataset(fpath) as ds:
        assert 'wavelength' not in ds.dimensions
        assert 'actinic_flux' not in ds.variables
        assert 'time' not in ds.dimensions
        assert ds['Kzz']._FillValue == FILL


def _output_stub(save_type, out_name='example', output_humanread=None):
    op = pytest.importorskip('vulcan.op', reason='needs a generated chem_funs.py')
    out = object.__new__(op.Output)
    out.cfg = SimpleNamespace(save_type=save_type, out_name=out_name, output_dir='/x/out')
    if output_humanread is not None:
        out.cfg.output_humanread = output_humanread
    return out


@pytest.mark.unit
def test_save_types_accepts_str_or_list():
    assert _output_stub('NetCDF')._save_types() == ['netcdf']
    assert _output_stub(['pickle', 'netcdf', 'pickle'])._save_types() == ['pickle', 'netcdf']
    with pytest.raises(ValueError, match='Invalid save_type'):
        _output_stub('hdf5')._save_types()


@pytest.mark.unit
def test_save_types_maps_deprecated_output_humanread():
    assert _output_stub('netcdf', output_humanread=True)._save_types() == ['netcdf', 'plain']
    assert _output_stub('pickle', output_humanread=False)._save_types() == ['pickle']


@pytest.mark.unit
@pytest.mark.parametrize('out_name', ['vulcan', 'vulcan.pkl', 'vulcan.nc'])
def test_out_path_replaces_known_extension(out_name):
    out = _output_stub('netcdf', out_name=out_name)
    assert out._out_path('.nc') == '/x/out/vulcan.nc'
    assert out._out_path('.pkl') == '/x/out/vulcan.pkl'
