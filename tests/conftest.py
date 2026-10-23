from __future__ import annotations

import os

import pytest

os.environ.setdefault('OMP_NUM_THREADS', '1')

import matplotlib  # noqa: E402

matplotlib.use('Agg')  # tests never open GUI windows


def pytest_configure(config):
    """Generate chem_funs.py for the default network before any test module is collected.

    `vulcan.op`, `vulcan.build_atm` and `vulcan.store` import the generated module at
    import time, so it must exist before test modules that import them are collected (a
    fresh checkout, as in CI, does not have it).

    It is regenerated even if present: the argument order of the generated rate functions
    depends on the per-process string hash seed, so a file written by another process can
    disagree with one regenerated later in this process (e.g. by the integration test).
    Code already bound to the old `chemdf` would then call the new rate functions with
    arguments in the wrong order. Generating it here keeps every regeneration in this
    session identical. Generation takes a few seconds.
    """
    from vulcan import Config
    from vulcan.make_chem_funs import make_all

    make_all(Config())


@pytest.hookimpl(trylast=True)
def pytest_sessionstart(session):
    """Print the data summary of every rendered validation plot, for the platform comparison."""
    import json
    import warnings

    from proteus_mpl import validation

    render = validation.render

    def logged(*args, **kwargs):
        files, summary = render(*args, **kwargs)
        warnings.warn(f'PILOT-SUMMARY {args[1]} {json.dumps(summary)}', stacklevel=2)
        return files, summary

    validation.render = logged
