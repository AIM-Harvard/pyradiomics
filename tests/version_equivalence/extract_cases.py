"""Shared extraction routine for the equivalence test suite.

`extract_all(...)` runs the BROAD param file (Original+LoG+Wavelet x all
feature classes) on the named case(s) and returns
``{case_id: {feature_name: float}}`` plus an env-provenance block. Used both by
``generate_baseline.py`` (to freeze a reference) and ``test_version_equivalence.py``
(to extract under the current env for comparison).
"""
from __future__ import annotations

import sys
import tempfile
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

HERE = Path(__file__).resolve().parent
PARAMS = HERE / "params_equivalence.yaml"


def _extractor():
    from radiomics import featureextractor

    return featureextractor.RadiomicsFeatureExtractor(str(PARAMS))


def _features(result):
    out = {}
    for k, v in result.items():
        if k.startswith("diagnostics_"):
            continue
        try:
            out[k] = float(v)
        except (TypeError, ValueError):
            continue
    return out


def _phantom_case():
    from phantom import make_phantom  # local import; same dir

    img, msk = make_phantom()
    return ("phantom", img, msk)


def _brain1_case():
    import radiomics

    img, msk = radiomics.getTestCase("brain1", tempfile.gettempdir())
    if img is None or msk is None:
        raise RuntimeError("brain1 test case unavailable (no network/cache)")
    return ("brain1", img, msk)


def list_cases(include_brain1: bool):
    cases = [_phantom_case()]
    if include_brain1:
        try:
            cases.append(_brain1_case())
        except Exception as exc:  # noqa: BLE001
            print(f"[warn] skipping brain1: {exc}", file=sys.stderr)
    return cases


def env_block():
    import numpy as np
    import radiomics
    import SimpleITK as sitk

    return {
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "pyradiomics": radiomics.__version__,
        "SimpleITK": sitk.Version.VersionString(),
    }


def extract_all(include_brain1: bool = True):
    sys.path.insert(0, str(HERE))
    ext = _extractor()
    cases = {}
    for cid, img, msk in list_cases(include_brain1):
        cases[cid] = _features(ext.execute(img, msk, label=1))
    return {"env": env_block(), "cases": cases}
