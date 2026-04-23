from __future__ import annotations


def test_package_imports():
    import vulcan

    assert isinstance(vulcan.__version__, str)
    assert vulcan.__version__


def test_cli_module_imports():
    from vulcan import run_vulcan

    assert callable(run_vulcan)
