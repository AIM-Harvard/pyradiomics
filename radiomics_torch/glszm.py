from __future__ import annotations

import torch

from radiomics_torch import base
from .cmatrices import calculate_glszm_torch, calculate_glszm_torch_vectorized
from .utils import delete_torch


class RadiomicsGLSZM(base.RadiomicsFeaturesBase):
    r"""
    （docstring 原样保留，这里省略解释，逻辑同 numpy 版）
    """

    def __init__(self, inputImage, inputMask, **kwargs):
        super().__init__(inputImage, inputMask, **kwargs)

        self.P_glszm = None
        self.imageArray = self._applyBinning(self.imageArray)

    def _initCalculation(self, voxelCoordinates=None):
        self.P_glszm = self._calculateMatrix(voxelCoordinates)

        self._calculateCoefficients()

        self.logger.debug(
            "GLSZM feature class initialized, calculated GLSZM with shape %s",
            self.P_glszm.shape,
        )

    def _calculateMatrix(self, voxelCoordinates=None):
        """
        Number of times a region with a
        gray level and voxel count occurs in an image. P_glszm[level, voxel_count] = # occurrences

        For 3D-images this concerns a 26-connected region, for 2D an 8-connected region
        """
        self.logger.debug("Calculating GLSZM matrix in Torch")
        Ng = self.coefficients["Ng"]

        # Ns = number of voxels in ROI
        # numpy 版: Ns = np.sum(self.maskArray)
        Ns = int(self.maskArray.sum().item())

        matrix_args = [
            self.imageArray,
            self.maskArray,
            Ng,
            Ns,
            self.settings.get("force2D", False),
            self.settings.get("force2Ddimension", 0),
        ]
        if self.voxelBased:
            matrix_args += [self.settings.get("kernelRadius", 1), voxelCoordinates]

        # 约定: calculate_glszm_torch 返回 Tensor, 形状 (Nvox, Ng, Ns)
        # P_glszm = calculate_glszm_torch(*matrix_args)  # (Nvox, Ng, Ns)
        P_glszm = calculate_glszm_torch_vectorized(
            self.imageArray,
            self.maskArray,
            Ng,
            Ns,
            self.settings.get("force2D", False),
            self.settings.get("force2Ddimension", 0),
            kernelRadius=0,
            voxels=None,
        )
        P_glszm = P_glszm.to(device=self.device, dtype=torch.float64)

        # ---- 删除 ROI 中不存在的灰度行（和 numpy 版保持完全一致） ----
        NgVector = list(range(1, Ng + 1))  # All possible gray values

        gray_levels_t = self.coefficients["grayLevels"]
        if isinstance(gray_levels_t, torch.Tensor):
            GrayLevels = gray_levels_t.to(torch.int64).tolist()
        else:
            GrayLevels = [int(g) for g in gray_levels_t]

        emptyGrayLevels_list = sorted(set(NgVector) - set(GrayLevels))  # 不在 ROI 中的灰度
        if len(emptyGrayLevels_list) > 0:
            emptyGrayLevels = torch.tensor(
                emptyGrayLevels_list, dtype=torch.int64, device=self.device
            )
            # 注意灰度是从 1 开始，索引从 0 开始 → 减 1
            P_glszm = delete_torch(P_glszm, emptyGrayLevels - 1, dim=1)

        return P_glszm

    def _calculateCoefficients(self):
        self.logger.debug("Calculating GLSZM coefficients")

        # P_glszm: (Nvox, Ng, Ns)
        ps = torch.sum(self.P_glszm, dim=1)  # (Nvox, Ns)
        pg = torch.sum(self.P_glszm, dim=2)  # (Nvox, Ng)

        ivector = self.coefficients["grayLevels"].to(torch.float64)  # (Ng,)
        jvector = torch.arange(
            1, self.P_glszm.shape[2] + 1, dtype=torch.float64, device=self.device
        )  # (Ns,)

        # Nz: number of zones
        Nz = torch.sum(self.P_glszm, dim=(1, 2))  # (Nvox,)
        Nz[Nz == 0] = 1

        # Np: number of voxels represented by GLSZM
        Np = torch.sum(ps * jvector[None, :], dim=1)  # (Nvox,)
        Np[Np == 0] = 1

        # 删除 zone size 不存在的列（size 上没出现）
        emptyZoneSizes_mask = torch.sum(ps, dim=0) == 0  # (Ns,)
        if emptyZoneSizes_mask.any():
            indices = torch.where(emptyZoneSizes_mask)[0]
            self.P_glszm = delete_torch(self.P_glszm, indices, dim=2)
            jvector = delete_torch(jvector, indices, dim=0)
            ps = delete_torch(ps, indices, dim=1)

        self.coefficients["Np"] = Np
        self.coefficients["Nz"] = Nz
        self.coefficients["ps"] = ps
        self.coefficients["pg"] = pg
        self.coefficients["ivector"] = ivector
        self.coefficients["jvector"] = jvector

    # =========================
    # Feature functions (Torch)
    # =========================

    def getSmallAreaEmphasisFeatureValue(self):
        ps = self.coefficients["ps"]
        jvector = self.coefficients["jvector"]
        Nz = self.coefficients["Nz"]

        return torch.sum(ps / (jvector[None, :] ** 2), dim=1) / Nz

    def getLargeAreaEmphasisFeatureValue(self):
        ps = self.coefficients["ps"]
        jvector = self.coefficients["jvector"]
        Nz = self.coefficients["Nz"]

        return torch.sum(ps * (jvector[None, :] ** 2), dim=1) / Nz

    def getGrayLevelNonUniformityFeatureValue(self):
        pg = self.coefficients["pg"]
        Nz = self.coefficients["Nz"]

        return torch.sum(pg**2, dim=1) / Nz

    def getGrayLevelNonUniformityNormalizedFeatureValue(self):
        pg = self.coefficients["pg"]
        Nz = self.coefficients["Nz"]

        return torch.sum(pg**2, dim=1) / (Nz**2)

    def getSizeZoneNonUniformityFeatureValue(self):
        ps = self.coefficients["ps"]
        Nz = self.coefficients["Nz"]

        return torch.sum(ps**2, dim=1) / Nz

    def getSizeZoneNonUniformityNormalizedFeatureValue(self):
        ps = self.coefficients["ps"]
        Nz = self.coefficients["Nz"]

        return torch.sum(ps**2, dim=1) / (Nz**2)

    def getZonePercentageFeatureValue(self):
        Nz = self.coefficients["Nz"]
        Np = self.coefficients["Np"]

        return Nz / Np

    def getGrayLevelVarianceFeatureValue(self):
        ivector = self.coefficients["ivector"]
        Nz = self.coefficients["Nz"]
        pg = self.coefficients["pg"] / Nz[:, None]  # 归一化

        u_i = torch.sum(pg * ivector[None, :], dim=1, keepdim=True)
        return torch.sum(pg * (ivector[None, :] - u_i) ** 2, dim=1)

    def getZoneVarianceFeatureValue(self):
        jvector = self.coefficients["jvector"]
        Nz = self.coefficients["Nz"]
        ps = self.coefficients["ps"] / Nz[:, None]  # 归一化

        u_j = torch.sum(ps * jvector[None, :], dim=1, keepdim=True)
        return torch.sum(ps * (jvector[None, :] - u_j) ** 2, dim=1)

    def getZoneEntropyFeatureValue(self):
        eps = torch.finfo(self.P_glszm.dtype).eps
        Nz = self.coefficients["Nz"]
        p_glszm = self.P_glszm / Nz[:, None, None]  # 归一化

        return -torch.sum(p_glszm * torch.log2(p_glszm + eps), dim=(1, 2))

    def getLowGrayLevelZoneEmphasisFeatureValue(self):
        pg = self.coefficients["pg"]
        ivector = self.coefficients["ivector"]
        Nz = self.coefficients["Nz"]

        return torch.sum(pg / (ivector[None, :] ** 2), dim=1) / Nz

    def getHighGrayLevelZoneEmphasisFeatureValue(self):
        pg = self.coefficients["pg"]
        ivector = self.coefficients["ivector"]
        Nz = self.coefficients["Nz"]

        return torch.sum(pg * (ivector[None, :] ** 2), dim=1) / Nz

    def getSmallAreaLowGrayLevelEmphasisFeatureValue(self):
        ivector = self.coefficients["ivector"]
        jvector = self.coefficients["jvector"]
        Nz = self.coefficients["Nz"]

        return (
            torch.sum(
                self.P_glszm
                / ((ivector[None, :, None] ** 2) * (jvector[None, None, :] ** 2)),
                dim=(1, 2),
            )
            / Nz
        )

    def getSmallAreaHighGrayLevelEmphasisFeatureValue(self):
        ivector = self.coefficients["ivector"]
        jvector = self.coefficients["jvector"]
        Nz = self.coefficients["Nz"]

        return (
            torch.sum(
                self.P_glszm
                * (ivector[None, :, None] ** 2)
                / (jvector[None, None, :] ** 2),
                dim=(1, 2),
            )
            / Nz
        )

    def getLargeAreaLowGrayLevelEmphasisFeatureValue(self):
        ivector = self.coefficients["ivector"]
        jvector = self.coefficients["jvector"]
        Nz = self.coefficients["Nz"]

        return (
            torch.sum(
                self.P_glszm
                * (jvector[None, None, :] ** 2)
                / (ivector[None, :, None] ** 2),
                dim=(1, 2),
            )
            / Nz
        )

    def getLargeAreaHighGrayLevelEmphasisFeatureValue(self):
        ivector = self.coefficients["ivector"]
        jvector = self.coefficients["jvector"]
        Nz = self.coefficients["Nz"]

        return (
            torch.sum(
                self.P_glszm
                * (ivector[None, :, None] ** 2)
                * (jvector[None, None, :] ** 2),
                dim=(1, 2),
            )
            / Nz
        )
