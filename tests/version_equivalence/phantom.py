"""Deterministic synthetic CT-like phantom for offline equivalence testing.

Builds a fixed, seed-stable 3-D textured volume + spherical ROI mask as
SimpleITK images. No network, no bundled data — fully reproducible across
machines and pyradiomics versions, so any feature-value difference is
attributable to the library/numpy stack alone.
"""
from __future__ import annotations

import numpy as np
import SimpleITK as sitk

SEED = 20260628
SHAPE = (40, 64, 64)          # (z, y, x)
SPACING = (1.5, 1.0, 1.0)     # anisotropic, like a thin-slice CT


def make_phantom():
    rng = np.random.default_rng(SEED)
    z, y, x = SHAPE
    zz, yy, xx = np.mgrid[0:z, 0:y, 0:x]

    # smooth anatomical gradient + a couple of blobs + reproducible texture
    base = (
        40.0
        + 0.8 * xx
        + 0.5 * yy
        + 25.0 * np.exp(-(((xx - 28) ** 2 + (yy - 30) ** 2 + (zz - 18) ** 2) / 120.0))
        - 15.0 * np.exp(-(((xx - 40) ** 2 + (yy - 36) ** 2 + (zz - 22) ** 2) / 80.0))
    )
    texture = rng.normal(0.0, 6.0, size=SHAPE)
    img_arr = (base + texture).astype(np.float32)

    # spherical ROI fully inside the FOV
    cz, cy, cx, r = 20, 32, 32, 14
    mask_arr = (
        ((zz - cz) ** 2 + (yy - cy) ** 2 + (xx - cx) ** 2) <= r ** 2
    ).astype(np.uint8)

    img = sitk.GetImageFromArray(img_arr)
    msk = sitk.GetImageFromArray(mask_arr)
    img.SetSpacing(SPACING)
    msk.SetSpacing(SPACING)
    return img, msk
