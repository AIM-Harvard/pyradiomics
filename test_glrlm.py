# test_glrlm.py
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
# （和 test_glcm.py / test_gldm.py 一样）
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
# radiomics_torch 版 GLRLM
# =========================================================

def glrlm_radiomics_torch(image_t: torch.Tensor,
                          mask_t: torch.Tensor,
                          spacing_zyx=(1.0, 1.0, 1.0),
                          settings=None,
                          label: int = 1) -> dict:
    """
    用 radiomics_torch 跑 GLRLM 特征。
    """
    image_t = image_t.to(torch.float32)

    # mask 中 label 区域
    mask_label_t = torch.zeros_like(mask_t, dtype=torch.int16)
    mask_label_t[mask_t > 0] = label

    if settings is None:
        # 尽量跟官方 pyradiomics 对齐（可以根据项目需要调整）
        settings = {
            "resampledPixelSpacing": None,  # 不重采样
            "normalize": False,
            "normalizeScale": 1,
            "label": label,
            "binWidth": 25,
            "force2D": False,
            # GLRLM 相关参数（如果你的 radiomics_torch 有自定义，可以在这里补）
            # 例如：'GLRLM_bins': Nr 之类的；否则用默认
        }

    extractor = torch_featureextractor.RadiomicsFeatureExtractor(**settings)
    extractor.disableAllImageTypes()
    extractor.enableImageTypeByName("Original")
    extractor.disableAllFeatures()
    extractor.enableFeatureClassByName("glrlm")

    print("radiomics_torch enabled glrlm features:", extractor.enabledFeatures)

    result = extractor.execute(image_t, mask_label_t)

    # 只保留 GLRLM 相关特征
    return {
        k: float(v) for k, v in result.items()
        if "glrlm" in k
    }


# =========================================================
# 官方 pyradiomics 版 GLRLM
# =========================================================

def glrlm_pyradiomics_sitk(image_t: torch.Tensor,
                           mask_t: torch.Tensor,
                           spacing_zyx=(1.0, 1.0, 1.0),
                           settings=None,
                           label: int = 1) -> dict:
    """
    用官方 pyradiomics + SimpleITK 跑 GLRLM。
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
            "symmetricalGLCM": True,  # 这个跟 GLCM 一致即可，GLRLM 本身不会用到
        }

    extractor = pyrad_featureextractor.RadiomicsFeatureExtractor(**settings)
    extractor.disableAllImageTypes()
    extractor.enableImageTypeByName("Original")
    extractor.disableAllFeatures()
    extractor.enableFeatureClassByName("glrlm")

    print("pyradiomics (official) enabled glrlm features:", extractor.enabledFeatures)

    result = extractor.execute(image_sitk, mask_sitk)

    return {
        k: float(v) for k, v in result.items()
        if "glrlm" in k
    }


# =========================================================
# GLRLM 对比
# =========================================================

def compare_glrlm_all(image_t: torch.Tensor,
                      mask_t: torch.Tensor,
                      spacing_zyx=(1.0, 1.0, 1.0),
                      pyrad_torch_settings=None,
                      pyrad_sitk_settings=None,
                      print_table: bool = True):
    """
    对比：
      - radiomics_torch GLRLM
      - 官方 pyradiomics GLRLM
    """
    glrlm_rtorch = glrlm_radiomics_torch(
        image_t, mask_t, spacing_zyx=spacing_zyx, settings=pyrad_torch_settings
    )
    glrlm_pyrad = glrlm_pyradiomics_sitk(
        image_t, mask_t, spacing_zyx=spacing_zyx, settings=pyrad_sitk_settings
    )

    if print_table:
        all_keys = sorted(set(glrlm_rtorch.keys()) | set(glrlm_pyrad.keys()))

        print("\n===== GLRLM comparison =====")
        print("{:<55s} {:>15s} {:>15s} {:>15s}".format(
            "Feature",
            "RadTorch",
            "RadOfficial",
            "RT - RO",
        ))
        print("-" * 110)

        def fmt(x):
            return "None" if x is None else f"{x:.6g}"

        for k in all_keys:
            vrt = glrlm_rtorch.get(k, None)
            vro = glrlm_pyrad.get(k, None)
            diff_rt_ro = None if vrt is None or vro is None else vrt - vro

            print("{:<55s} {:>15s} {:>15s} {:>15s}".format(
                k,
                fmt(vrt),
                fmt(vro),
                fmt(diff_rt_ro),
            ))

    return glrlm_rtorch, glrlm_pyrad


# =========================================================
# demo
# =========================================================

if __name__ == "__main__":
    # 构造一个简单 16^3 体数据，与 test_glcm/test_gldm 保持一致
    Z = Y = X = 16
    img = torch.zeros((Z, Y, X), dtype=torch.float32)
    msk = torch.zeros((Z, Y, X), dtype=torch.uint8)

    # ROI: 4:12 的立方体 => 8x8x8 = 512 个 voxel
    # 为了让 GLRLM 有点结构，不是常数，这里给个线性梯度
    roi = torch.linspace(0, 255, steps=8 * 8 * 8, dtype=torch.float32).reshape(8, 8, 8)
    img[4:12, 4:12, 4:12] = roi
    msk[4:12, 4:12, 4:12] = 1

    spacing = (1.0, 1.0, 1.0)

    glrlm_rtorch, glrlm_pyrad = compare_glrlm_all(
        img, msk, spacing_zyx=spacing,
        pyrad_torch_settings=None,
        pyrad_sitk_settings=None,
        print_table=True,
    )

    print("\nGLRLM radiomics_torch keys:", glrlm_rtorch.keys())
    print("GLRLM pyradiomics (official) keys:", glrlm_pyrad.keys())
