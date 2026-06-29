# pyradiomics cross-version feature-equivalence test suite

**Purpose.** Prove that a modernized pyradiomics build (scikit-build-core +
**numpy 2** + CPython 3.10–3.13 wheels) reproduces the established feature values
of the last upstream release (**3.0.1 / numpy 1.x**) within floating-point
tolerance. This is the *scientific-trust* guarantee that the modernization does
not silently shift IBSI feature values — the verification piece that accompanies
[@steubk](https://github.com/steubk)'s build-system work
(`steubk/pyradiomics`, branch `fix/python312-compatibility`).

> **Credit.** The numpy-2 / scikit-build-core / cibuildwheel modernization is
> @steubk's. This suite adds the feature-equivalence verification on top of it;
> it is intended to accompany an upstream PR with @steubk as co-author, not to
> fork-and-claim.

## What it checks

A deliberately **broad** panel so the whole numpy-2 migration surface is
exercised:

- **Image types**: Original + LoG (σ=1.0) + Wavelet → C-extension filtering,
  scipy LoG, and PyWavelets paths.
- **Feature classes**: all 7 (firstorder, glcm, glrlm, glszm, gldm, ngtdm,
  shape) → every `_cmatrices` / `_cshape` C kernel plus the pure-numpy
  first-order path.
- **Cases** (both fully public / reproducible — no private data):
  - a deterministic synthetic phantom (`phantom.py`, seed-stable, no download), and
  - the canonical `brain1` test case (auto-fetched via `radiomics.getTestCase`;
    skipped gracefully if offline).
- = **944 features × 2 cases = 1888 comparisons.**

Each feature is compared to a frozen reference baseline at relative tolerance
`RTOL` (default `1e-6`). The cross-stack difference observed in practice is
~1e-14 (numpy 1↔2 float reassociation on reductions such as Skewness / Kurtosis
/ Variance / GLCM ClusterShade) — orders of magnitude inside tolerance.

## Files

| file | role |
|---|---|
| `params_equivalence.yaml` | broad extraction params (all classes × 3 image types) |
| `phantom.py` | deterministic public synthetic image+mask |
| `extract_cases.py` | shared extraction routine (env-agnostic) |
| `generate_baseline.py` | freeze a reference baseline (run under the reference env) |
| `test_version_equivalence.py` | pytest: current env vs frozen baseline |
| `baseline/baseline_3.0.1.json` | frozen reference (pyradiomics 3.0.1 / numpy 1.26.4 / py3.10.19) |

## Usage

```bash
# 1. (one time, under the reference env) freeze the baseline
python generate_baseline.py --out baseline/baseline_3.0.1.json

# 2. (under any new build) verify equivalence
pip install pytest
PYRADIOMICS_BASELINE=baseline/baseline_3.0.1.json pytest test_version_equivalence.py -v -s

# offline-only (skip brain1): generate with --no-brain1; the test auto-matches the case set
```

Env knobs: `PYRADIOMICS_BASELINE`, `PYRADIOMICS_RTOL` (default 1e-6),
`PYRADIOMICS_ATOL` (default 1e-8).

## Result on the steubk build

| baseline | current | comparisons | verdict |
|---|---|---|---|
| 3.0.1 / numpy 1.26.4 / py3.10.19 | 3.1.0+py312fix / numpy 2.5.0 / py3.13.11 | 1888 | ✅ all within rtol 1e-6 |
