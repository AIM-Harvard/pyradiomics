from __future__ import annotations

import importlib.metadata

import radiomics as m


def test_version():
    assert importlib.metadata.version("pyradiomics") == m.__version__
