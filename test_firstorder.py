import numpy as np
import torch
import SimpleITK as sitk

# 1️⃣ 原始 pyradiomics（pip 安装的那个）
from radiomics import featureextractor as pyrad_featureextractor

# 2️⃣ 你改的 tensor 版，在 /Users/donglihe/projects/pyradiomics/radiomics_torch
from radiomics_torch import featureextractor as torch_featureextractor


# =========================================================
# 1. 你的「纯 torch」版本 first-order 计算
# =========================================================

def firstorder_torch_manual(image_t: torch.Tensor,
                            mask_t: torch.Tensor,
                            spacing_zyx=None) -> dict:
    """
    你自己写的 / 将来要自己写的 pure-torch first-order 实现。
    现在放一个示例版（Mean / Var / Std / Min / Max / Median / Range / Energy / TotalEnergy）

    image_t: (Z, Y, X) float tensor
    mask_t:  (Z, Y, X) bool / 0/1 tensor
    spacing_zyx: (dz, dy, dx)   用于 TotalEnergy
    """
    if not torch.is_floating_point(image_t):
        image_t = image_t.to(torch.float32)

    mask_bool = mask_t > 0
    vals = image_t[mask_bool]
    vals = vals[~torch.isnan(vals)]

    if vals.numel() == 0:
        return {
            "firstorder_Mean_manual": float("nan"),
            "firstorder_Variance_manual": float("nan"),
        }

    mean = vals.mean()
    var = vals.var(unbiased=False)
    std = vals.std(unbiased=False)
    minimum = vals.min()
    maximum = vals.max()
    median = vals.median()
    rng = maximum - minimum

    energy = torch.sum(vals ** 2)
    if spacing_zyx is not None:
        dz, dy, dx = spacing_zyx
        voxel_volume = float(dz * dy * dx)
        total_energy = energy * voxel_volume
    else:
        total_energy = energy

    return {
        "firstorder_Mean_manual": float(mean),
        "firstorder_Variance_manual": float(var),
        "firstorder_Std_manual": float(std),
        "firstorder_Minimum_manual": float(minimum),
        "firstorder_Maximum_manual": float(maximum),
        "firstorder_Median_manual": float(median),
        "firstorder_Range_manual": float(rng),
        "firstorder_Energy_manual": float(energy),
        "firstorder_TotalEnergy_manual": float(total_energy),
    }


# =========================================================
# 2. radiomics_torch 版本（你的 tensor 版 radiomics）
# =========================================================

def firstorder_radiomics_torch(image_t: torch.Tensor,
                               mask_t: torch.Tensor,
                               spacing_zyx=(1.0, 1.0, 1.0),
                               settings: dict | None = None,
                               label: int = 1) -> dict:
    """
    用你自己改的 radiomics_torch 跑 firstorder。

    这个包的 featureextractor.execute 里是用 torch 操作的：
      boundingBox = torch.where(maskArray == label)
    所以这里直接传 tensor 进去。
    """
    image_t = image_t.to(torch.float32)

    mask_label_t = torch.zeros_like(mask_t, dtype=torch.int16)
    mask_label_t[mask_t > 0] = label

    if settings is None:
        settings = {
            "resampledPixelSpacing": None,
            "normalize": False,
            "normalizeScale": 1,
            "label": label,
            "binWidth": 25,
        }

    extractor = torch_featureextractor.RadiomicsFeatureExtractor(**settings)
    extractor.disableAllImageTypes()
    extractor.enableImageTypeByName("Original")
    extractor.disableAllFeatures()
    extractor.enableFeatureClassByName("firstorder")

    print("radiomics_torch enabled features:", extractor.enabledFeatures)

    # 🔑 radiomics_torch 版本：直接 tensor
    result = extractor.execute(image_t, mask_label_t)

    return {
        k: float(v) for k, v in result.items()
        if "firstorder" in k
    }


# =========================================================
# 3. 原版 pyradiomics（SimpleITK + C 扩展）
# =========================================================

def numpy_to_sitk(image_np: np.ndarray,
                  spacing_zyx=None,
                  origin=None,
                  direction=None) -> sitk.Image:
    """
    image_np: shape (Z, Y, X)
    SimpleITK: axis0 认为是 Z
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


def firstorder_pyradiomics_sitk(image_t: torch.Tensor,
                                mask_t: torch.Tensor,
                                spacing_zyx=(1.0, 1.0, 1.0),
                                settings: dict | None = None,
                                label: int = 1) -> dict:
    """
    用原版 pyradiomics（radiomics）+ SimpleITK 跑 firstorder。
    这里需要确保 import 到的是 pip 安装的那个 radiomics，不是你的 radiomics_torch。
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
        }

    extractor = pyrad_featureextractor.RadiomicsFeatureExtractor(**settings)
    extractor.disableAllImageTypes()
    extractor.enableImageTypeByName("Original")
    extractor.disableAllFeatures()
    extractor.enableFeatureClassByName("firstorder")

    print("pyradiomics (official) enabled features:", extractor.enabledFeatures)

    # 🔑 原版：传 SITK image
    result = extractor.execute(image_sitk, mask_sitk)

    return {
        k: float(v) for k, v in result.items()
        if "firstorder" in k
    }


# =========================================================
# 4. 三方对比函数
# =========================================================

def compare_firstorder_all(image_t: torch.Tensor,
                           mask_t: torch.Tensor,
                           spacing_zyx=(1.0, 1.0, 1.0),
                           torch_settings=None,
                           pyrad_torch_settings=None,
                           pyrad_sitk_settings=None,
                           print_table: bool = True):
    """
    同一套 image/mask 上，对比：
      - 你自己写的 firstorder_torch_manual
      - radiomics_torch (tensor 版)
      - 原版 pyradiomics (sitk + C 扩展)
    """
    # 1) 自己写的纯 torch
    fo_manual = firstorder_torch_manual(image_t, mask_t, spacing_zyx=spacing_zyx)

    # 2) radiomics_torch
    fo_rtorch = firstorder_radiomics_torch(
        image_t, mask_t, spacing_zyx=spacing_zyx, settings=pyrad_torch_settings
    )

    # 3) 原版 pyradiomics
    fo_pyrad = firstorder_pyradiomics_sitk(
        image_t, mask_t, spacing_zyx=spacing_zyx, settings=pyrad_sitk_settings
    )

    if print_table:
        # 1) 建立“核心特征名” -> 对应 key 的映射
        # manual:  firstorder_Mean_manual   -> core: "Mean"
        # rtorch:  original_firstorder_Mean -> core: "Mean"
        # pyrad:   original_firstorder_Mean -> core: "Mean"

        def core_from_manual(k: str) -> str | None:
            if k.startswith("firstorder_") and k.endswith("_manual"):
                return k[len("firstorder_"):-len("_manual")]
            return None

        def core_from_radiomics(k: str) -> str | None:
            if k.startswith("original_firstorder_"):
                return k[len("original_firstorder_"):]
            return None

        manual_by_core: dict[str, str] = {}
        for k in fo_manual.keys():
            core = core_from_manual(k)
            if core is not None:
                manual_by_core[core] = k

        rtorch_by_core: dict[str, str] = {}
        for k in fo_rtorch.keys():
            core = core_from_radiomics(k)
            if core is not None:
                rtorch_by_core[core] = k

        pyrad_by_core: dict[str, str] = {}
        for k in fo_pyrad.keys():
            core = core_from_radiomics(k)
            if core is not None:
                pyrad_by_core[core] = k

        # 2) 所有出现过的核心特征名
        all_cores = sorted(set(manual_by_core) |
                           set(rtorch_by_core) |
                           set(pyrad_by_core))

        print("\n===== First-order comparison =====")
        print("{:<40s} {:>15s} {:>15s} {:>15s} {:>15s} {:>15s}".format(
            "Feature",
            "ManualTorch",
            "RadTorch",
            "RadOfficial",
            "M - RT",
            "RT - RO",
        ))
        print("-" * 120)

        def fmt(x):
            return "None" if x is None else f"{x:.6g}"

        for core in all_cores:
            mk = manual_by_core.get(core)
            rk = rtorch_by_core.get(core)
            ok = pyrad_by_core.get(core)

            vm = fo_manual.get(mk) if mk is not None else None
            vrt = fo_rtorch.get(rk) if rk is not None else None
            vro = fo_pyrad.get(ok) if ok is not None else None

            diff_m_rt = None if vm is None or vrt is None else vm - vrt
            diff_rt_ro = None if vrt is None or vro is None else vrt - vro

            # 行名：统一显示成 original_firstorder_<Core>
            feature_name = f"original_firstorder_{core}"

            print("{:<40s} {:>15s} {:>15s} {:>15s} {:>15s} {:>15s}".format(
                feature_name,
                fmt(vm),
                fmt(vrt),
                fmt(vro),
                fmt(diff_m_rt),
                fmt(diff_rt_ro),
            ))

    return fo_manual, fo_rtorch, fo_pyrad


# =========================================================
# 5. demo
# =========================================================

if __name__ == "__main__":
    # 造一个简单 16^3 toy volume：中间 8^3 是 100，其他为 0
    Z = Y = X = 16
    img = torch.zeros((Z, Y, X), dtype=torch.float32)
    msk = torch.zeros((Z, Y, X), dtype=torch.uint8)

    img[4:8, 4:12, 4:12] = 100.0
    img[8:10, 8:10, 4:12] = 150.0
    img[10:12, 10:12, 4:12] = 214.0
    msk[4:12, 4:12, 4:12] = 1

    spacing = (1.0, 1.0, 1.0)

    fo_manual, fo_rtorch, fo_pyrad = compare_firstorder_all(
        img, msk, spacing_zyx=spacing,
        torch_settings=None,
        pyrad_torch_settings=None,
        pyrad_sitk_settings=None,
        print_table=True,
    )

    print("\nManual torch keys:", fo_manual.keys())
    print("radiomics_torch keys:", fo_rtorch.keys())
    print("pyradiomics (official) keys:", fo_pyrad.keys())
