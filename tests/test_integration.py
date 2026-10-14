from __future__ import annotations

import os

os.environ['OMP_NUM_THREADS'] = '1'  # noqa
import pickle
from pathlib import Path

import numpy as np
import pytest

from vulcan import Config, run_vulcan


@pytest.mark.integration
@pytest.mark.timeout(1800)
def test_default_config_run_creates_default_outputs():
    repo_root = Path(__file__).resolve().parents[1]
    fastchem_bin = repo_root / 'fastchem_vulcan' / 'fastchem'

    assert fastchem_bin.is_file(), (
        f'Missing FastChem binary at {fastchem_bin}; compile fastchem_vulcan first.'
    )

    # expected output files and config dump file
    output_file = repo_root / 'output' / 'example.pkl'
    nc_file = repo_root / 'output' / 'example.nc'
    cfg_dump = repo_root / 'output' / 'vulcan_cfg.txt'
    output_file.unlink(missing_ok=True)
    nc_file.unlink(missing_ok=True)
    cfg_dump.unlink(missing_ok=True)

    # Set up the config
    vulcan_cfg = Config()
    vulcan_cfg.save_type = ['netcdf', 'pickle']

    # Run VULCAN
    run_vulcan(vulcan_cfg, make_network=True)

    # Check output files exist
    assert output_file.is_file()
    assert nc_file.is_file()

    # Check config exists
    assert cfg_dump.is_file()

    # Check plots folder exists
    assert (repo_root / 'output' / 'plot').is_dir()

    # Check mix.png and evo.png exist
    assert (repo_root / 'output' / 'plot' / 'mix.png').is_file()
    assert (repo_root / 'output' / 'plot' / 'evo.png').is_file()

    check_netcdf_consistency(nc_file, output_file)


def check_netcdf_consistency(nc_file, pkl_file):
    """The NetCDF output must agree with the pickle and with the model's own RHS."""
    import netCDF4 as nc

    from vulcan import chem_funs

    with open(pkl_file, 'rb') as hdl:
        pkl = pickle.load(hdl)

    with nc.Dataset(nc_file) as ds:
        species = list(ds['species'][:])
        y = np.asarray(ds['number_density'][:])
        ymix = np.asarray(ds['mixing_ratio'][:])
        M = np.asarray(ds['third_body_density'][:])
        k_arr = np.asarray(ds['rate_coefficient'][:])
        rates = np.asarray(ds['reaction_rate'][:])
        Tco = np.asarray(ds['temperature'][:])

    assert species == pkl['variable']['species']
    np.testing.assert_array_equal(ymix, pkl['variable']['ymix'])
    np.testing.assert_array_equal(Tco, pkl['atm']['Tco'])
    np.testing.assert_allclose(ymix.sum(axis=1), 1.0, rtol=1e-6)

    # Net tendency from the stored reaction rates equals the generated chemdf
    nr = chem_funs.nr
    nu = np.zeros((nr, len(species)))
    for i in range(1, nr + 1):
        for sp in chem_funs.re_dict[i][0]:
            nu[i - 1, species.index(sp)] -= 1
        for sp in chem_funs.re_dict[i][1]:
            nu[i - 1, species.index(sp)] += 1
    k = {i: k_arr[:, i - 1] for i in range(1, nr + 1)}
    expected = chem_funs.chemdf(y, M, k)
    scale = np.abs(rates) @ np.abs(nu)  # gross production + loss per species
    np.testing.assert_array_less(np.abs(rates @ nu - expected), 1e-10 * scale + 1e-300)
