# test_ct_real.py
from __future__ import annotations

import numpy as np
import torch
import SimpleITK as sitk

from monai.transforms import (
    Compose, LoadImaged, Spacingd, EnsureChannelFirstd, Transposed
)

# 1️⃣ 官方 pyradiomics
from radiomics import featureextractor as pyrad_featureextractor

# 2️⃣ 你的 tensor 版 radiomics_torch
from radiomics_torch import featureextractor as torch_featureextractor


# =========================================================
# 工具函数：numpy <-> SimpleITK（如果后面想用 SITK 形式也可以用到）
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
# MONAI: 从路径读取 CT & mask -> spacing(1,1,1) -> (Z,Y,X) tensor
# =========================================================

def load_ct_with_monai_as_torch(
    image_path: str,
    mask_path: str,
    target_spacing_xyz=(1.0, 1.0, 1.0),
) -> tuple[torch.Tensor, torch.Tensor, tuple[float, float, float]]:
    """
    用 MONAI 从文件读取 CT 和 mask，并转换到 radiomics_torch 需要的格式。

    Pipeline:
      - LoadImaged(keys=["image", "mask"])
      - EnsureChannelFirstd(keys=["image", "mask"])  # (C, X, Y, Z)
      - Spacingd(pixdim=(1,1,1), mode=("bilinear","nearest"))
      - Transposed(keys=["image", "mask"], indices=(0, 3, 2, 1))
        => (C, Z, Y, X)
      - 去掉 C 维，得到 (Z, Y, X)

    返回:
      image_t: torch.float32, shape (Z, Y, X)
      mask_t: torch.int16,   shape (Z, Y, X)
      spacing_zyx: (dz, dy, dx) = (1,1,1)
    """

    transforms = Compose([
        LoadImaged(keys=["image", "mask"]),
        EnsureChannelFirstd(keys=["image", "mask"]),  # -> (1, X, Y, Z)
        Spacingd(
            keys=["image", "mask"],
            pixdim=target_spacing_xyz,   # (sx, sy, sz) = (dx, dy, dz)
            mode=("bilinear", "nearest")
        ),
        # (C, X, Y, Z) -> (C, Z, Y, X)
        Transposed(keys=["image", "mask"], indices=(0, 3, 2, 1)),
    ])

    data = transforms({"image": image_path, "mask": mask_path})

    image_np = np.asarray(data["image"])  # (1, Z, Y, X)
    mask_np = np.asarray(data["mask"])    # (1, Z, Y, X)

    # 去掉 channel 维
    image_np = image_np[0]  # (Z, Y, X)
    mask_np = mask_np[0]    # (Z, Y, X)

    image_t = torch.from_numpy(image_np.copy()).to(torch.float32)
    mask_t = torch.from_numpy(mask_np.copy()).to(torch.int16)

    # Spacingd 已经把 spacing 变成 (1,1,1)
    spacing_zyx = (1.0, 1.0, 1.0)

    return image_t, mask_t, spacing_zyx


# =========================================================
# radiomics_torch：从 path 跑所有特征（不含 shape2D）
# =========================================================

def radiomics_torch_all_from_paths(
    image_path: str,
    mask_path: str,
    label: int = 1,
) -> dict:
    """
    从 CT path + mask path 出发，用 MONAI -> torch -> radiomics_torch。
    """
    image_t, mask_t, spacing_zyx = load_ct_with_monai_as_torch(
        image_path, mask_path, target_spacing_xyz=(1.0, 1.0, 1.0)
    )

    # radiomics_torch: mask 中 label 区域
    mask_label_t = torch.zeros_like(mask_t, dtype=torch.int16)
    mask_label_t[mask_t > 0] = label

    settings = {
        # 已经在 MONAI 里 spacing 到 (1,1,1)，这里不再重采样
        "resampledPixelSpacing": None,
        "normalize": False,
        "normalizeScale": 1,
        "label": label,
        "binWidth": 25,
        "force2D": False,
        "distances": [1],
    }

    extractor = torch_featureextractor.RadiomicsFeatureExtractor(**settings)

    extractor.disableAllFeatures()
    for cls in ["firstorder", "shape", "glcm", "gldm", "glrlm", "glszm", "ngtdm"]:
        extractor.enableFeatureClassByName(cls)

    extractor.disableAllImageTypes()
    extractor.enableImageTypeByName("Original")

    print("\n[radiomics_torch] Enabled features:", extractor.enabledFeatures)

    result = extractor.execute(image_t, mask_label_t)

    out: dict[str, float] = {}
    for k, v in result.items():
        try:
            out[k] = float(v)
        except Exception:
            continue

    return out


# =========================================================
# 官方 pyradiomics：直接用 path（仍然不启用 shape2D）
# =========================================================

def radiomics_pyradiomics_all_from_paths(
    image_path: str,
    mask_path: str,
    label: int = 1,
) -> dict:
    """
    官方 pyradiomics：直接用文件路径。
    设置 resampledPixelSpacing=[1,1,1]，让它内部重采样到同样的 spacing。
    """
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

    result = extractor.execute(image_path, mask_path)

    out: dict[str, float] = {}
    for k, v in result.items():
        try:
            out[k] = float(v)
        except Exception:
            continue

    return out


# =========================================================
# 对比输出
# =========================================================

def compare_all_features_from_paths(
    image_path: str,
    mask_path: str,
    label: int = 1,
    print_table: bool = True,
):
    torch_res = radiomics_torch_all_from_paths(image_path, mask_path, label=label)
    pyrad_res = radiomics_pyradiomics_all_from_paths(image_path, mask_path, label=label)

    if print_table:
        all_keys = sorted(set(torch_res.keys()) | set(pyrad_res.keys()))

        print("\n=================== All Features Comparison on REAL CT (no shape2D) ===================")
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
            vt = torch_res.get(k, None)
            vo = pyrad_res.get(k, None)
            diff = None if vt is None or vo is None else vt - vo

            print("{:<60s} {:>15s} {:>15s} {:>15s}".format(
                k,
                fmt(vt),
                fmt(vo),
                fmt(diff),
            ))

    print("\nTorch feature count:", len(torch_res))
    print("Official feature count:", len(pyrad_res))

    return torch_res, pyrad_res


# =========================================================
# Demo: 在真实 CT 上跑（把路径改成你自己的）
# =========================================================

if __name__ == "__main__":
    # TODO: 换成你自己的 CT / mask 路径
    # 例如 NIfTI:  "/data/cases/case001_ct.nii.gz"
    #             "/data/cases/case001_mask.nii.gz"
    image_path = "/path/to/ct_image.nii.gz"
    mask_path = "/path/to/ct_mask.nii.gz"

    compare_all_features_from_paths(image_path, mask_path, label=1, print_table=True)
