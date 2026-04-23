from __future__ import annotations

import os

os.environ['OMP_NUM_THREADS'] = '1'  # noqa
from pathlib import Path

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

    # expected output file and config dump file
    output_file = repo_root / 'output' / 'example.pkl'
    cfg_dump = repo_root / 'output' / 'vulcan_cfg.txt'
    output_file.unlink(missing_ok=True)
    cfg_dump.unlink(missing_ok=True)

    # Set environment variables for the subprocess
    env = os.environ.copy()
    env['MPLBACKEND'] = 'Agg'

    # Set up the config
    vulcan_cfg = Config()

    # Run VULCAN
    run_vulcan(vulcan_cfg, make_network=True)

    # Check output file exists
    assert output_file.is_file()

    # Check config exists
    assert cfg_dump.is_file()

    # Check plots folder exists
    assert (repo_root / 'output' / 'plot').is_dir()

    # Check mix.png and evo.png exist
    assert (repo_root / 'output' / 'plot' / 'mix.png').is_file()
    assert (repo_root / 'output' / 'plot' / 'evo.png').is_file()
