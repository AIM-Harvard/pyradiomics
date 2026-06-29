"""Cross-version feature-equivalence test.

Extracts the broad feature panel (Original+LoG+Wavelet x all classes) in the
CURRENT environment and asserts every feature matches a frozen reference
baseline within tolerance. This is the missing scientific-trust guarantee for
the modernized (scikit-build-core + numpy 2 + cp3.10-3.13) build: it proves the
new stack reproduces the established pyradiomics feature values.

Run:
  pytest test_version_equivalence.py -v
Select a specific baseline:
  PYRADIOMICS_BASELINE=baseline/baseline_3.0.1.json pytest -v

Tolerances:
  RTOL  relative, default 1e-6 (well above the ~1e-14 float-reassociation noise
        observed between numpy 1.x and 2.x, far below any biologically/
        statistically meaningful change).
  ATOL  absolute floor for near-zero reference values, default 1e-8.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from extract_cases import extract_all

HERE = Path(__file__).resolve().parent
RTOL = float(os.environ.get("PYRADIOMICS_RTOL", "1e-6"))
ATOL = float(os.environ.get("PYRADIOMICS_ATOL", "1e-8"))


def _baseline_path() -> Path:
    env = os.environ.get("PYRADIOMICS_BASELINE")
    if env:
        p = Path(env)
        return p if p.is_absolute() else HERE / p
    # default: the single *.json under baseline/ (or baseline.json)
    bdir = HERE / "baseline"
    preferred = bdir / "baseline.json"
    if preferred.exists():
        return preferred
    cands = sorted(bdir.glob("*.json"))
    if not cands:
        pytest.skip("no baseline JSON found in baseline/ — run generate_baseline.py first")
    return cands[0]


@pytest.fixture(scope="module")
def baseline():
    return json.loads(_baseline_path().read_text())


@pytest.fixture(scope="module")
def current(baseline):
    # match the baseline's case set: only attempt brain1 if the baseline has it
    include_brain1 = "brain1" in baseline["cases"]
    return extract_all(include_brain1=include_brain1)


def test_envs_differ_is_ok(baseline, current):
    # informational: the whole point is comparing across different stacks
    print("\nbaseline env:", baseline["env"])
    print("current  env:", current["env"])


def _case_ids(baseline, current):
    return sorted(set(baseline["cases"]) & set(current["cases"]))


def test_case_coverage(baseline, current):
    common = _case_ids(baseline, current)
    assert common, "no shared cases between baseline and current extraction"
    for cid in common:
        assert baseline["cases"][cid], f"baseline case {cid} has no features"


def test_feature_coverage(baseline, current):
    for cid in _case_ids(baseline, current):
        b = set(baseline["cases"][cid])
        c = set(current["cases"][cid])
        missing = b - c
        assert not missing, f"{cid}: {len(missing)} baseline features missing in current: {sorted(missing)[:5]}"


def test_feature_values_equivalent(baseline, current):
    failures = []
    n_compared = 0
    for cid in _case_ids(baseline, current):
        bcase, ccase = baseline["cases"][cid], current["cases"][cid]
        for feat, bval in bcase.items():
            if feat not in ccase:
                continue
            cval = ccase[feat]
            n_compared += 1
            denom = max(abs(bval), abs(cval), ATOL)
            rel = abs(bval - cval) / denom
            if rel > RTOL:
                failures.append((cid, feat, bval, cval, rel))
    print(f"\ncompared {n_compared} feature values @ rtol={RTOL:g} across "
          f"{len(_case_ids(baseline, current))} case(s)")
    if failures:
        msg = "\n".join(
            f"  {cid}:{feat}  base={b:.6g} cur={c:.6g} rel={r:.2e}"
            for cid, feat, b, c, r in sorted(failures, key=lambda x: -x[4])[:20]
        )
        pytest.fail(f"{len(failures)} feature(s) exceed rtol={RTOL:g}:\n{msg}")
