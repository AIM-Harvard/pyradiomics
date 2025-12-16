# test_glszm.py
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
# （和 test_glcm.py / test_gldm.py / test_glrlm.py 保持一致）
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
# radiomics_torch 版 GLSZM
# =========================================================

def glszm_radiomics_torch(image_t: torch.Tensor,
                          mask_t: torch.Tensor,
                          spacing_zyx=(1.0, 1.0, 1.0),
                          settings=None,
                          label: int = 1) -> dict:
    """
    用 radiomics_torch 跑 GLSZM 特征。
    """
    image_t = image_t.to(torch.float32)

    # mask 中 label 区域
    mask_label_t = torch.zeros_like(mask_t, dtype=torch.int16)
    mask_label_t[mask_t > 0] = label

    if settings is None:
        # 尽量和官方 pyradiomics 设置保持一致（可以按项目需要调整）
        settings = {
            "resampledPixelSpacing": None,  # 不重采样
            "normalize": False,
            "normalizeScale": 1,
            "label": label,
            "binWidth": 25,
            "force2D": False,
        }

    extractor = torch_featureextractor.RadiomicsFeatureExtractor(**settings)
    extractor.disableAllImageTypes()
    extractor.enableImageTypeByName("Original")
    extractor.disableAllFeatures()
    extractor.enableFeatureClassByName("glszm")

    print("radiomics_torch enabled glszm features:", extractor.enabledFeatures)

    # radiomics_torch: 直接 tensor
    result = extractor.execute(image_t, mask_label_t)

    # 只保留 GLSZM 相关特征
    return {
        k: float(v) for k, v in result.items()
        if "glszm" in k
    }


# =========================================================
# 官方 pyradiomics 版 GLSZM
# =========================================================

def glszm_pyradiomics_sitk(image_t: torch.Tensor,
                           mask_t: torch.Tensor,
                           spacing_zyx=(1.0, 1.0, 1.0),
                           settings=None,
                           label: int = 1) -> dict:
    """
    用官方 pyradiomics + SimpleITK 跑 GLSZM。
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
            # GLSZM 跟 GLCM 的 symmetricalGLCM 没关系，但留不留都无所谓
            "symmetricalGLCM": True,
        }

    extractor = pyrad_featureextractor.RadiomicsFeatureExtractor(**settings)
    extractor.disableAllImageTypes()
    extractor.enableImageTypeByName("Original")
    extractor.disableAllFeatures()
    extractor.enableFeatureClassByName("glszm")

    print("pyradiomics (official) enabled glszm features:", extractor.enabledFeatures)

    result = extractor.execute(image_sitk, mask_sitk)

    return {
        k: float(v) for k, v in result.items()
        if "glszm" in k
    }


# =========================================================
# GLSZM 对比
# =========================================================

def compare_glszm_all(image_t: torch.Tensor,
                      mask_t: torch.Tensor,
                      spacing_zyx=(1.0, 1.0, 1.0),
                      pyrad_torch_settings=None,
                      pyrad_sitk_settings=None,
                      print_table: bool = True):
    """
    对比：
      - radiomics_torch GLSZM
      - 官方 pyradiomics GLSZM
    """
    glszm_rtorch = glszm_radiomics_torch(
        image_t, mask_t, spacing_zyx=spacing_zyx, settings=pyrad_torch_settings
    )
    glszm_pyrad = glszm_pyradiomics_sitk(
        image_t, mask_t, spacing_zyx=spacing_zyx, settings=pyrad_sitk_settings
    )

    if print_table:
        all_keys = sorted(set(glszm_rtorch.keys()) | set(glszm_pyrad.keys()))

        print("\n===== GLSZM comparison =====")
        print("{:<60s} {:>15s} {:>15s} {:>15s}".format(
            "Feature",
            "RadTorch",
            "RadOfficial",
            "RT - RO",
        ))
        print("-" * 115)

        def fmt(x):
            return "None" if x is None else f"{x:.6g}"

        for k in all_keys:
            vrt = glszm_rtorch.get(k, None)
            vro = glszm_pyrad.get(k, None)
            diff_rt_ro = None if vrt is None or vro is None else vrt - vro

            print("{:<60s} {:>15s} {:>15s} {:>15s}".format(
                k,
                fmt(vrt),
                fmt(vro),
                fmt(diff_rt_ro),
            ))

    return glszm_rtorch, glszm_pyrad


# =========================================================
# demo
# =========================================================

if __name__ == "__main__":
    # 构造一个简单 16^3 体数据
    Z = Y = X = 16
    img = torch.zeros((Z, Y, X), dtype=torch.float32)
    msk = torch.zeros((Z, Y, X), dtype=torch.uint8)

    # ROI: 4:12 的立方体 => 8x8x8 = 512 个 voxel
    # 为了让 GLSZM 有些结构，不是常数，这里给个线性梯度
    roi = torch.linspace(0, 255, steps=8 * 8 * 8, dtype=torch.float32).reshape(8, 8, 8)
    img[4:12, 4:12, 4:12] = roi
    msk[4:12, 4:12, 4:12] = 1

    spacing = (1.0, 1.0, 1.0)

    glszm_rtorch, glszm_pyrad = compare_glszm_all(
        img, msk, spacing_zyx=spacing,
        pyrad_torch_settings=None,
        pyrad_sitk_settings=None,
        print_table=True,
    )

    print("\nGLSZM radiomics_torch keys:", glszm_rtorch.keys())
    print("GLSZM pyradiomics (official) keys:", glszm_pyrad.keys())
