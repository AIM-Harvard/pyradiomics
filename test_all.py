# test_all.py
from __future__ import annotations

import numpy as np
import torch
import SimpleITK as sitk

# 1️⃣ 官方 pyradiomics
from radiomics import featureextractor as pyrad_featureextractor

# 2️⃣ 你的 tensor 版 radiomics_torch
from radiomics_torch import featureextractor as torch_featureextractor


# =========================================================
# 工具函数：numpy <-> SimpleITK
# =========================================================

def numpy_to_sitk(image_np: np.ndarray,
                  spacing_zyx=(1.0, 1.0, 1.0),
                  origin=None,
                  direction=None) -> sitk.Image:
    """
    image_np: shape (Z, Y, X)
    spacing_zyx: (dz, dy, dx)
    SimpleITK spacing 顺序是 (sx, sy, sz) = (dx, dy, dz)
    """
    img = sitk.GetImageFromArray(image_np.astype(np.float32))

    if spacing_zyx is not None:
        dz, dy, dx = spacing_zyx
        img.SetSpacing((float(dx), float(dy), float(dz)))

    if origin is not None:
        img.SetOrigin(tuple(origin))
    if direction is not None:
        img.SetDirection(tuple(direction))

    return img


def numpy_mask_to_sitk(mask_np: np.ndarray,
                       spacing_zyx=(1.0, 1.0, 1.0),
                       origin=None,
                       direction=None) -> sitk.Image:
    mask = sitk.GetImageFromArray(mask_np.astype(np.uint8))

    if spacing_zyx is not None:
        dz, dy, dx = spacing_zyx
        mask.SetSpacing((float(dx), float(dy), float(dz)))

    if origin is not None:
        mask.SetOrigin(tuple(origin))
    if direction is not None:
        mask.SetDirection(tuple(direction))

    return mask


# =========================================================
# radiomics_torch：开启所有特征（不含 shape2D）
# =========================================================

def radiomics_torch_all(image_t: torch.Tensor,
                        mask_t: torch.Tensor,
                        spacing_zyx=(1.0, 1.0, 1.0),
                        settings=None,
                        label: int = 1) -> dict:
    """
    用 radiomics_torch 跑所有 feature class
    （firstorder / shape / glcm / glrlm / glszm / gldm / ngtdm），不启用 shape2D。
    """
    image_t = image_t.to(torch.float32)

    # mask 中 label 区域
    mask_label_t = torch.zeros_like(mask_t, dtype=torch.int16)
    mask_label_t[mask_t > 0] = label

    if settings is None:
        settings = {
            "resampledPixelSpacing": None,  # 不重采样
            "normalize": False,
            "normalizeScale": 1,
            "label": label,
            "binWidth": 25,
            "force2D": False,
            "distances": [1],
        }

    extractor = torch_featureextractor.RadiomicsFeatureExtractor(**settings)

    # 不直接 enableAllFeatures，避免 shape2D
    extractor.disableAllFeatures()
    for cls in ["firstorder", "shape", "glcm", "gldm", "glrlm", "glszm", "ngtdm"]:
        extractor.enableFeatureClassByName(cls)

    extractor.disableAllImageTypes()
    extractor.enableImageTypeByName("Original")

    print("\n[radiomics_torch] Enabled features:", extractor.enabledFeatures)

    result = extractor.execute(image_t, mask_label_t)

    # 关键修正点：任何能 float(v) 的都收下，支持 torch.Tensor / numpy / python 标量
    out: dict[str, float] = {}
    for k, v in result.items():
        try:
            out[k] = float(v)
        except Exception:
            # 例如有些 value 是 list/dict/对象，直接跳过
            continue

    return out


# =========================================================
# 官方 pyradiomics：开启所有特征（同样不含 shape2D）
# =========================================================

def radiomics_pyradiomics_all(image_t: torch.Tensor,
                              mask_t: torch.Tensor,
                              spacing_zyx=(1.0, 1.0, 1.0),
                              settings=None,
                              label: int = 1) -> dict:
    """
    用官方 pyradiomics + SimpleITK 跑所有 feature class（不启用 shape2D）。
    """
    image_np = image_t.detach().cpu().numpy().astype(np.float32)
    mask_np = (mask_t.detach().cpu().numpy() > 0).astype(np.uint8) * label

    image_sitk = numpy_to_sitk(image_np, spacing_zyx)
    mask_sitk = numpy_mask_to_sitk(mask_np, spacing_zyx)

    if settings is None:
        settings = {
            "resampledPixelSpacing": [1, 1, 1],
            "interpolator": sitk.sitkBSpline,
            "normalize": False,
            "normalizeScale": 1,
            "label": label,
            "binWidth": 25,
            "force2D": False,
            "distances": [1],
        }

    extractor = pyrad_featureextractor.RadiomicsFeatureExtractor(**settings)

    extractor.disableAllFeatures()
    for cls in ["firstorder", "shape", "glcm", "gldm", "glrlm", "glszm", "ngtdm"]:
        extractor.enableFeatureClassByName(cls)

    extractor.disableAllImageTypes()
    extractor.enableImageTypeByName("Original")

    print("\n[pyradiomics_official] Enabled features:", extractor.enabledFeatures)

    result = extractor.execute(image_sitk, mask_sitk)

    out: dict[str, float] = {}
    for k, v in result.items():
        try:
            out[k] = float(v)
        except Exception:
            continue

    return out


# =========================================================
# 对比所有特征
# =========================================================

def compare_all_features(image_t: torch.Tensor,
                         mask_t: torch.Tensor,
                         spacing_zyx=(1.0, 1.0, 1.0),
                         pyrad_torch_settings=None,
                         pyrad_sitk_settings=None,
                         print_table: bool = True):
    """
    对比：
      - radiomics_torch (all feature classes except shape2D)
      - 官方 pyradiomics (同样的 feature classes)
    """
    res_torch = radiomics_torch_all(
        image_t, mask_t, spacing_zyx=spacing_zyx, settings=pyrad_torch_settings
    )
    res_pyrad = radiomics_pyradiomics_all(
        image_t, mask_t, spacing_zyx=spacing_zyx, settings=pyrad_sitk_settings
    )

    if print_table:
        all_keys = sorted(set(res_torch.keys()) | set(res_pyrad.keys()))

        print("\n=================== All Features Comparison (no shape2D) ===================")
        print("{:<60s} {:>15s} {:>15s} {:>15s}".format(
            "Feature",
            "Torch",
            "Official",
            "Torch - Off",
        ))
        print("-" * 115)

        def fmt(x):
            return "None" if x is None else f"{x:.6g}"

        for k in all_keys:
            vt = res_torch.get(k, None)
            vo = res_pyrad.get(k, None)
            diff = None if vt is None or vo is None else vt - vo

            print("{:<60s} {:>15s} {:>15s} {:>15s}".format(
                k,
                fmt(vt),
                fmt(vo),
                fmt(diff),
            ))

    return res_torch, res_pyrad


# =========================================================
# demo
# =========================================================

if __name__ == "__main__":
    # 构造一个简单 16^3 体数据
    Z = Y = X = 16
    img = torch.zeros((Z, Y, X), dtype=torch.float32)
    msk = torch.zeros((Z, Y, X), dtype=torch.uint8)

    # ROI: 4:12 的立方体 => 8x8x8 = 512 个 voxel
    # 随机离散灰度，让所有矩阵类特征有点结构
    roi = torch.randint(low=1, high=20, size=(8, 8, 8), dtype=torch.int32).to(torch.float32)
    img[4:12, 4:12, 4:12] = roi
    msk[4:12, 4:12, 4:12] = 1

    spacing = (1.0, 1.0, 1.0)

    res_torch, res_pyrad = compare_all_features(
        img, msk, spacing_zyx=spacing,
        pyrad_torch_settings=None,
        pyrad_sitk_settings=None,
        print_table=True,
    )

    print("\nTorch feature count:", len(res_torch))
    print("Official feature count:", len(res_pyrad))
