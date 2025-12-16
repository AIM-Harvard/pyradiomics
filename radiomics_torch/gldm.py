from __future__ import annotations

import torch

from radiomics_torch import base, deprecated

from .cmatrices import calculate_gldm_torch
from .utils import delete_torch


class RadiomicsGLDM(base.RadiomicsFeaturesBase):

    def __init__(self, inputImage, inputMask, **kwargs):
        super().__init__(inputImage, inputMask, **kwargs)

        self.gldm_a = kwargs.get("gldm_a", 0)

        self.P_gldm = None
        self.imageArray = self._applyBinning(self.imageArray)

    def _initCalculation(self, voxelCoordinates=None):
        self.P_gldm = self._calculateMatrix(voxelCoordinates)

        self.logger.debug(
            "Feature class initialized, calculated GLDM with shape %s",
            self.P_gldm.shape,
        )

    def _calculateMatrix(self, voxelCoordinates=None):
        Ng = self.coefficients["Ng"]

        matrix_args = [
            self.imageArray,
            self.maskArray,
            torch.tensor(self.settings.get("distances", [1]), device=self.device),
            Ng,
            self.gldm_a,
            self.settings.get("force2D", False),
            self.settings.get("force2Ddimension", 0),
        ]
        if self.voxelBased:
            matrix_args += [self.settings.get("kernelRadius", 1), voxelCoordinates]

        P_gldm = calculate_gldm_torch(*matrix_args)  # shape (Nv, Ng, Nd)

        # ---- 关键：用 Python int 做集合运算，而不是直接对 torch.Tensor 调用 set() ----
        NgVector = list(range(1, Ng + 1))  # All possible gray values

        gray_levels_t = self.coefficients["grayLevels"]
        # 确保拿到的是一串 Python int
        if isinstance(gray_levels_t, torch.Tensor):
            GrayLevels = gray_levels_t.to(torch.int64).tolist()
        else:
            # 兼容 numpy / list 的情况
            GrayLevels = [int(g) for g in gray_levels_t]

        emptyGrayLevels_list = sorted(set(NgVector) - set(GrayLevels))  # 真的“没出现过”的灰度
        if len(emptyGrayLevels_list) > 0:
            emptyGrayLevels = torch.tensor(
                emptyGrayLevels_list, dtype=torch.int64, device=self.device
            )
            # 删除这些灰度对应的行（注意 -1 转成 0-based index）
            P_gldm = delete_torch(P_gldm, emptyGrayLevels - 1, dim=1)

        # 下面保持不变（只是顺手加点注释）
        jvector = torch.arange(1, P_gldm.shape[2] + 1, dtype=torch.float64, device=self.device)

        # shape (Nv, Nd)
        pd = torch.sum(P_gldm, dim=1)
        # shape (Nv, Ng')
        pg = torch.sum(P_gldm, dim=2)

        # Delete columns that dependence sizes not present in the ROI
        empty_sizes = torch.sum(pd, dim=0)   # (Nd,)
        indices = torch.where(empty_sizes == 0)[0]

        if indices.numel() > 0:
            P_gldm = delete_torch(P_gldm, indices, dim=2)
            jvector = delete_torch(jvector, indices, dim=0)
            pd = delete_torch(pd, indices, dim=1)

        Nz = torch.sum(pd, dim=1)  # Nz per kernel, shape (Nv,)
        Nz[Nz == 0] = 1

        self.coefficients["Nz"] = Nz
        self.coefficients["pd"] = pd
        self.coefficients["pg"] = pg
        self.coefficients["ivector"] = self.coefficients["grayLevels"].to(torch.float64)
        self.coefficients["jvector"] = jvector

        return P_gldm


    def getSmallDependenceEmphasisFeatureValue(self):
        pd = self.coefficients["pd"]
        jvector = self.coefficients["jvector"]
        Nz = self.coefficients["Nz"]  # Nz = Np, see class docstring

        return torch.sum(pd / (jvector[None, :] ** 2), dim=1) / Nz

    def getLargeDependenceEmphasisFeatureValue(self):
        pd = self.coefficients["pd"]
        jvector = self.coefficients["jvector"]
        Nz = self.coefficients["Nz"]

        return torch.sum(pd * (jvector[None, :] ** 2), dim=1) / Nz

    def getGrayLevelNonUniformityFeatureValue(self):
        pg = self.coefficients["pg"]
        Nz = self.coefficients["Nz"]

        return torch.sum(pg**2, dim=1) / Nz

    @deprecated
    def getGrayLevelNonUniformityNormalizedFeatureValue(self):
        msg = (
            "GLDM - Gray Level Non-Uniformity Normalized is mathematically equal to First Order - "
            "Uniformity, see http://pyradiomics.readthedocs.io/en/latest/removedfeatures.html for more"
            "details"
        )
        raise DeprecationWarning(msg)

    def getDependenceNonUniformityFeatureValue(self):
        pd = self.coefficients["pd"]
        Nz = self.coefficients["Nz"]

        return torch.sum(pd**2, dim=1) / Nz

    def getDependenceNonUniformityNormalizedFeatureValue(self):
        pd = self.coefficients["pd"]
        Nz = self.coefficients["Nz"]

        return torch.sum(pd**2, dim=1) / Nz**2

    def getGrayLevelVarianceFeatureValue(self):
        ivector = self.coefficients["ivector"]
        Nz = self.coefficients["Nz"]
        pg = (
            self.coefficients["pg"] / Nz[:, None]
        )  # divide by Nz to get the normalized matrix

        u_i = torch.sum(pg * ivector[None, :], dim=1, keepdims=True)
        return torch.sum(pg * (ivector[None, :] - u_i) ** 2, dim=1)

    def getDependenceVarianceFeatureValue(self):
        jvector = self.coefficients["jvector"]
        Nz = self.coefficients["Nz"]
        pd = (
            self.coefficients["pd"] / Nz[:, None]
        )  # divide by Nz to get the normalized matrix

        u_j = torch.sum(pd * jvector[None, :], dim=1, keepdims=True)
        return torch.sum(pd * (jvector[None, :] - u_j) ** 2, dim=1)

    def getDependenceEntropyFeatureValue(self):
        eps = torch.finfo(self.P_gldm.dtype).eps
        Nz = self.coefficients["Nz"]
        p_gldm = (
            self.P_gldm / Nz[:, None, None]
        )  # divide by Nz to get the normalized matrix

        return -torch.sum(p_gldm * torch.log2(p_gldm + eps), dim=(1, 2))

    @deprecated
    def getDependencePercentageFeatureValue(self):
        msg = (
            "GLDM - Dependence Percentage always computes 1, "
            "see http://pyradiomics.readthedocs.io/en/latest/removedfeatures.html for more details"
        )
        raise DeprecationWarning(msg)

    def getLowGrayLevelEmphasisFeatureValue(self):
        pg = self.coefficients["pg"]
        ivector = self.coefficients["ivector"]
        Nz = self.coefficients["Nz"]

        return torch.sum(pg / (ivector[None, :] ** 2), dim=1) / Nz

    def getHighGrayLevelEmphasisFeatureValue(self):
        pg = self.coefficients["pg"]
        ivector = self.coefficients["ivector"]
        Nz = self.coefficients["Nz"]

        return torch.sum(pg * (ivector[None, :] ** 2), dim=1) / Nz

    def getSmallDependenceLowGrayLevelEmphasisFeatureValue(self):
        ivector = self.coefficients["ivector"]
        jvector = self.coefficients["jvector"]
        Nz = self.coefficients["Nz"]

        return (
            torch.sum(
                self.P_gldm
                / ((ivector[None, :, None] ** 2) * (jvector[None, None, :] ** 2)),
                dim=(1, 2),
            )
            / Nz
        )

    def getSmallDependenceHighGrayLevelEmphasisFeatureValue(self):
        ivector = self.coefficients["ivector"]
        jvector = self.coefficients["jvector"]
        Nz = self.coefficients["Nz"]

        return (
            torch.sum(
                self.P_gldm
                * (ivector[None, :, None] ** 2)
                / (jvector[None, None, :] ** 2),
                dim=(1, 2),
            )
            / Nz
        )

    def getLargeDependenceLowGrayLevelEmphasisFeatureValue(self):
        ivector = self.coefficients["ivector"]
        jvector = self.coefficients["jvector"]
        Nz = self.coefficients["Nz"]

        return (
            torch.sum(
                self.P_gldm
                * (jvector[None, None, :] ** 2)
                / (ivector[None, :, None] ** 2),
                dim=(1, 2),
            )
            / Nz
        )

    def getLargeDependenceHighGrayLevelEmphasisFeatureValue(self):
        ivector = self.coefficients["ivector"]
        jvector = self.coefficients["jvector"]
        Nz = self.coefficients["Nz"]

        return (
            torch.sum(
                self.P_gldm
                * ((jvector[None, None, :] ** 2) * (ivector[None, :, None] ** 2)),
                dim=(1, 2),
            )
            / Nz
        )
