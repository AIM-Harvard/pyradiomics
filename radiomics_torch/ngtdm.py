from __future__ import annotations

import torch
from radiomics_torch import base
from .cmatrices import calculate_ngtdm_torch


class RadiomicsNGTDM(base.RadiomicsFeaturesBase):

    def __init__(self, inputImage, inputMask, **kwargs):
        super().__init__(inputImage, inputMask, **kwargs)

        self.P_ngtdm = None

        # 你保证 imageArray/maskArray 已是 torch.Tensor
        self.imageArray = self._applyBinning(self.imageArray)
        self.maskArray = self.maskArray.bool()

    def _initCalculation(self, voxelCoordinates=None):
        self.P_ngtdm = self._calculateMatrix(voxelCoordinates)
        self._calculateCoefficients()

    # ------------------------------------------------
    # NGTDM Matrix 计算（Torch）
    # ------------------------------------------------
    def _calculateMatrix(self, voxelCoordinates=None):
        distances = torch.as_tensor(
            self.settings.get("distances", [1]),
            dtype=torch.int32,
            device=self.imageArray.device,
        )

        matrix_args = [
            self.imageArray,
            self.maskArray,
            distances,
            int(self.coefficients["Ng"]),
            self.settings.get("force2D", False),
            self.settings.get("force2Ddimension", 0),
        ]

        if self.voxelBased:
            kernel_radius = int(self.settings.get("kernelRadius", 1))
            matrix_args += [kernel_radius, voxelCoordinates]

        P_ngtdm = calculate_ngtdm_torch(*matrix_args)  # (Nvox, Ng, 3)

        # 删除空灰度级
        valid_gray = torch.sum(P_ngtdm[:, :, 0], dim=0) != 0
        P_ngtdm = P_ngtdm[:, valid_gray, :]

        return P_ngtdm

    # ------------------------------------------------
    # 系数计算（Torch）
    # ------------------------------------------------
    def _calculateCoefficients(self):
        n_i = self.P_ngtdm[:, :, 0]
        s_i = self.P_ngtdm[:, :, 1]
        i_vec = self.P_ngtdm[:, :, 2]

        Nvp = torch.sum(n_i, dim=1)      # (Nvox,)
        self.coefficients["Nvp"] = Nvp

        self.coefficients["p_i"] = n_i / Nvp.unsqueeze(1)
        self.coefficients["s_i"] = s_i
        self.coefficients["ivector"] = i_vec

        self.coefficients["Ngp"] = torch.sum(n_i > 0, dim=1)

        # p_i = 0 的 index，用于 mask 屏蔽
        self.coefficients["p_zero"] = torch.where(self.coefficients["p_i"] == 0)

    # ------------------------------------------------
    # Feature: Coarseness
    # ------------------------------------------------
    def getCoarsenessFeatureValue(self):
        p_i = self.coefficients["p_i"]
        s_i = self.coefficients["s_i"]

        x = torch.sum(p_i * s_i, dim=1)
        out = torch.empty_like(x)

        non_zero = x != 0
        zero = ~non_zero

        out[non_zero] = 1.0 / x[non_zero]
        out[zero] = 1e6
        return out

    # ------------------------------------------------
    # Feature: Contrast
    # ------------------------------------------------
    def getContrastFeatureValue(self):
        Ngp = self.coefficients["Ngp"]
        Nvp = self.coefficients["Nvp"]
        p_i = self.coefficients["p_i"]
        s_i = self.coefficients["s_i"]
        i = self.coefficients["ivector"]

        div = Ngp * (Ngp - 1)

        p_i_i = p_i.unsqueeze(2)
        p_i_j = p_i.unsqueeze(1)
        diff_sq = (i.unsqueeze(2) - i.unsqueeze(1)) ** 2

        A = torch.sum(p_i_i * p_i_j * diff_sq, dim=(1, 2))
        B = torch.sum(s_i, dim=1) / Nvp

        contrast = A * B

        non_zero = div != 0
        zero = ~non_zero

        out = torch.zeros_like(contrast)
        out[non_zero] = contrast[non_zero] / div[non_zero]
        out[zero] = 0.0

        return out

    # ------------------------------------------------
    # Feature: Busyness
    # ------------------------------------------------
    def getBusynessFeatureValue(self):
        p_i = self.coefficients["p_i"]
        s_i = self.coefficients["s_i"]
        i = self.coefficients["ivector"]
        p_zero = self.coefficients["p_zero"]

        i_pi = i * p_i

        absdiff = torch.abs(i_pi.unsqueeze(2) - i_pi.unsqueeze(1))

        # 屏蔽 p=0 的行列
        absdiff[p_zero[0], :, p_zero[1]] = 0
        absdiff[p_zero[0], p_zero[1], :] = 0

        denom = torch.sum(absdiff, dim=(1, 2))
        numer = torch.sum(p_i * s_i, dim=1)

        out = torch.zeros_like(numer)
        non_zero = denom != 0
        out[non_zero] = numer[non_zero] / denom[non_zero]
        return out

    # ------------------------------------------------
    # Feature: Complexity
    # ------------------------------------------------
    def getComplexityFeatureValue(self):
        Nvp = self.coefficients["Nvp"]
        p_i = self.coefficients["p_i"]
        s_i = self.coefficients["s_i"]
        i = self.coefficients["ivector"]
        p_zero = self.coefficients["p_zero"]

        pi_si = p_i * s_i
        numerator = pi_si.unsqueeze(2) + pi_si.unsqueeze(1)

        # 屏蔽 p=0 的项
        numerator[p_zero[0], :, p_zero[1]] = 0
        numerator[p_zero[0], p_zero[1], :] = 0

        divisor = p_i.unsqueeze(2) + p_i.unsqueeze(1)
        divisor_mask = divisor == 0
        divisor = divisor.clone()
        divisor[divisor_mask] = 1  # numerator 也为 0，所以合理

        diff = torch.abs(i.unsqueeze(2) - i.unsqueeze(1))

        comp = torch.sum(diff * numerator / divisor, dim=(1, 2))
        return comp / Nvp

    # ------------------------------------------------
    # Feature: Strength
    # ------------------------------------------------
    def getStrengthFeatureValue(self):
        p_i = self.coefficients["p_i"]
        s_i = self.coefficients["s_i"]
        i = self.coefficients["ivector"]
        sum_s = torch.sum(s_i, dim=1)
        p_zero = self.coefficients["p_zero"]

        pi_pj = p_i.unsqueeze(2) + p_i.unsqueeze(1)
        diff_sq = (i.unsqueeze(2) - i.unsqueeze(1)) ** 2

        strength = pi_pj * diff_sq

        strength[p_zero[0], :, p_zero[1]] = 0
        strength[p_zero[0], p_zero[1], :] = 0

        numer = torch.sum(strength, dim=(1, 2))

        out = torch.zeros_like(numer)
        non_zero = sum_s != 0
        out[non_zero] = numer[non_zero] / sum_s[non_zero]
        return out
