#!/usr/bin/env python3
"""Freeze a reference feature baseline for the equivalence test suite.

Run this ONCE under a trusted reference environment (e.g. the upstream-released
pyradiomics 3.0.1 + numpy 1.x) to produce a baseline JSON that
``test_version_equivalence.py`` then checks any new build against.

Usage:
  python generate_baseline.py [--no-brain1] [--out baseline/baseline.json]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from extract_cases import extract_all

HERE = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-brain1", action="store_true", help="phantom only")
    ap.add_argument("--out", default=str(HERE / "baseline" / "baseline.json"))
    args = ap.parse_args()

    payload = extract_all(include_brain1=not args.no_brain1)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, sort_keys=True))

    print("env:", payload["env"])
    for cid, feats in payload["cases"].items():
        print(f"  {cid}: {len(feats)} features")
    print("wrote", out)


if __name__ == "__main__":
    main()
