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
                  spacing_zyx=None,
                  origin=None,
                  direction=None) -> sitk.Image:
    """
    image_np: shape (Z, Y, X)
    spacing_zyx: (dz, dy, dx)
    SITK spacing 顺序是 (sx, sy, sz) = (dx, dy, dz)
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
                       spacing_zyx=None,
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
# radiomics_torch 版 shape
# =========================================================

def shape_radiomics_torch(image_t: torch.Tensor,
                          mask_t: torch.Tensor,
                          spacing_zyx=(1.0, 1.0, 1.0),
                          settings=None,
                          label: int = 1) -> dict:
    """
    用 radiomics_torch 跑 shape 特征。
    """
    image_t = image_t.to(torch.float32)

    # mask 中 label 区域
    mask_label_t = torch.zeros_like(mask_t, dtype=torch.int16)
    mask_label_t[mask_t > 0] = label

    if settings is None:
        settings = {
            "resampledPixelSpacing": None,  # 不重采样，直接用原 spacing
            "normalize": False,
            "normalizeScale": 1,
            "label": label,
        }

    extractor = torch_featureextractor.RadiomicsFeatureExtractor(**settings)
    extractor.disableAllImageTypes()
    extractor.enableImageTypeByName("Original")
    extractor.disableAllFeatures()
    extractor.enableFeatureClassByName("shape")

    print("radiomics_torch enabled shape features:", extractor.enabledFeatures)

    # radiomics_torch: 直接 tensor
    result = extractor.execute(image_t, mask_label_t)

    # 只保留 shape 相关特征
    return {
        k: float(v) for k, v in result.items()
        if "shape" in k
    }


# =========================================================
# 官方 pyradiomics 版 shape
# =========================================================

def shape_pyradiomics_sitk(image_t: torch.Tensor,
                           mask_t: torch.Tensor,
                           spacing_zyx=(1.0, 1.0, 1.0),
                           settings=None,
                           label: int = 1) -> dict:
    """
    用官方 pyradiomics + SimpleITK 跑 shape。
    """
    image_np = image_t.detach().cpu().numpy().astype(np.float32)
    mask_np = (mask_t.detach().cpu().numpy() > 0).astype(np.uint8) * label

    image_sitk = numpy_to_sitk(image_np, spacing_zyx)
    mask_sitk = numpy_mask_to_sitk(mask_np, spacing_zyx)

    if settings is None:
        settings = {
            "resampledPixelSpacing": [1, 1, 1],  # 一般 shape 用原始 spacing 也可以
            "interpolator": sitk.sitkBSpline,
            "normalize": False,
            "normalizeScale": 1,
            "label": label,
        }

    extractor = pyrad_featureextractor.RadiomicsFeatureExtractor(**settings)
    extractor.disableAllImageTypes()
    extractor.enableImageTypeByName("Original")
    extractor.disableAllFeatures()
    extractor.enableFeatureClassByName("shape")

    print("pyradiomics (official) enabled shape features:", extractor.enabledFeatures)

    result = extractor.execute(image_sitk, mask_sitk)

    return {
        k: float(v) for k, v in result.items()
        if "shape" in k
    }


# =========================================================
# shape 对比
# =========================================================

def compare_shape_all(image_t: torch.Tensor,
                      mask_t: torch.Tensor,
                      spacing_zyx=(1.0, 1.0, 1.0),
                      pyrad_torch_settings=None,
                      pyrad_sitk_settings=None,
                      print_table: bool = True):
    """
    对比：
      - radiomics_torch shape
      - 官方 pyradiomics shape
    """
    sh_rtorch = shape_radiomics_torch(
        image_t, mask_t, spacing_zyx=spacing_zyx, settings=pyrad_torch_settings
    )
    sh_pyrad = shape_pyradiomics_sitk(
        image_t, mask_t, spacing_zyx=spacing_zyx, settings=pyrad_sitk_settings
    )

    if print_table:
        all_keys = sorted(set(sh_rtorch.keys()) | set(sh_pyrad.keys()))

        print("\n===== Shape comparison =====")
        print("{:<45s} {:>15s} {:>15s} {:>15s}".format(
            "Feature",
            "RadTorch",
            "RadOfficial",
            "RT - RO",
        ))
        print("-" * 100)

        def fmt(x):
            return "None" if x is None else f"{x:.6g}"

        for k in all_keys:
            vrt = sh_rtorch.get(k, None)
            vro = sh_pyrad.get(k, None)
            diff_rt_ro = None if vrt is None or vro is None else vrt - vro

            print("{:<45s} {:>15s} {:>15s} {:>15s}".format(
                k,
                fmt(vrt),
                fmt(vro),
                fmt(diff_rt_ro),
            ))

    return sh_rtorch, sh_pyrad


# =========================================================
# demo
# =========================================================

if __name__ == "__main__":
    # 构造一个简单 16^3 体数据：
    # 只要 mask 一致，shape 只和几何有关，灰度值无所谓
    Z = Y = X = 16
    img = torch.zeros((Z, Y, X), dtype=torch.float32)
    msk = torch.zeros((Z, Y, X), dtype=torch.uint8)

    # ROI: 4:12 的立方体 => 8x8x8 = 512 个 voxel
    img[4:12, 4:12, 4:12] = 100.0  # 给 ROI 一个非零值，防止某些实现做特殊处理
    msk[4:12, 4:12, 4:12] = 1

    spacing = (1.0, 1.0, 1.0)

    sh_rtorch, sh_pyrad = compare_shape_all(
        img, msk, spacing_zyx=spacing,
        pyrad_torch_settings=None,
        pyrad_sitk_settings=None,
        print_table=True,
    )

    print("\nShape radiomics_torch keys:", sh_rtorch.keys())
    print("Shape pyradiomics (official) keys:", sh_pyrad.keys())