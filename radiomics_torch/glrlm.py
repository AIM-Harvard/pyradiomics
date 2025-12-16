from __future__ import annotations

import torch

from radiomics_torch import base

from .cmatrices import calculate_glrlm_torch, glrlm_vectorized_3d
from .utils import torch_delete


class RadiomicsGLRLM(base.RadiomicsFeaturesBase):
    def __init__(self, inputImage, inputMask, **kwargs):
        super().__init__(inputImage, inputMask, **kwargs)

        self.weightingNorm = kwargs.get(
            "weightingNorm"
        )  # manhattan, euclidean, infinity

        self.P_glrlm = None
        self.imageArray = self._applyBinning(self.imageArray)

    def _initCalculation(self, voxelCoordinates=None):
        self.P_glrlm = self._calculateMatrix(voxelCoordinates)

        self._calculateCoefficients()

        self.logger.debug(
            "GLRLM feature class initialized, calculated GLRLM with shape %s",
            self.P_glrlm.shape,
        )

    def _calculateMatrix(self, voxelCoordinates=None):
        self.logger.debug("Calculating GLRLM matrix in C")

        Ng = self.coefficients["Ng"]
        Nr = torch.tensor(max(self.imageArray.shape), device=self.device)

        matrix_args = [
            self.imageArray,
            self.maskArray,
            Ng,
            Nr,
            self.settings.get("force2D", False),
            self.settings.get("force2Ddimension", 0),
        ]
        if self.voxelBased:
            matrix_args += [self.settings.get("kernelRadius", 1), voxelCoordinates]

        P_glrlm, angles = calculate_glrlm_torch(
            *matrix_args
        )  # shape (Nvox, Ng, Nr, Na)

        self.logger.debug("Process calculated matrix")

        # Delete rows that specify gray levels not present in the ROI
        NgVector = range(1, Ng + 1)  # All possible gray values
        GrayLevels = self.coefficients["grayLevels"].tolist()  # Gray values present in ROI
        emptyGrayLevels = torch.tensor(
            list(set(NgVector) - set(GrayLevels)), dtype=torch.int64, device=self.device
        )  # Gray values NOT present in ROI

        P_glrlm = torch_delete(P_glrlm, emptyGrayLevels - 1, axis=1)

        # Optionally apply a weighting factor
        if self.weightingNorm is not None:
            self.logger.debug("Applying weighting (%s)", self.weightingNorm)

            pixelSpacing = self.spacing
            weights = torch.empty(len(angles), device=self.device)
            for a_idx, a in enumerate(angles):
                if self.weightingNorm == "infinity":
                    weights[a_idx] = max(torch.abs(a) * pixelSpacing)
                elif self.weightingNorm == "euclidean":
                    weights[a_idx] = torch.sqrt(torch.sum((torch.abs(a) * pixelSpacing) ** 2))
                elif self.weightingNorm == "manhattan":
                    weights[a_idx] = torch.sum(torch.abs(a) * pixelSpacing)
                elif self.weightingNorm == "no_weighting":
                    weights[a_idx] = 1
                else:
                    self.logger.warning(
                        'weigthing norm "%s" is unknown, weighting factor is set to 1',
                        self.weightingNorm,
                    )
                    weights[a_idx] = 1

            P_glrlm = torch.sum(P_glrlm * weights[None, None, None, :], dim=3, keepdim=True)

        Nr = torch.sum(P_glrlm, dim=(1, 2))

        # Delete empty angles if no weighting is applied
        if P_glrlm.shape[3] > 1:
            emptyAngles = torch.where(torch.sum(Nr, 0) == 0)
            if len(emptyAngles[0]) > 0:  # One or more angles are 'empty'
                self.logger.debug(
                    "Deleting %d empty angles:\n%s",
                    len(emptyAngles[0]),
                    angles[emptyAngles],
                )
                P_glrlm = torch_delete(P_glrlm, emptyAngles, axis=3)
                Nr = torch_delete(Nr, emptyAngles, axis=1)
            else:
                self.logger.debug("No empty angles")

        Nr[Nr == 0] = torch.nan  # set sum to numpy.spacing(1) if sum is 0?
        self.coefficients["Nr"] = Nr

        return P_glrlm

    def _calculateCoefficients(self):
        self.logger.debug("Calculating GLRLM coefficients")

        pr = torch.sum(self.P_glrlm, dim=1)  # shape (Nvox, Nr, Na)
        pg = torch.sum(self.P_glrlm, dim=2)  # shape (Nvox, Ng, Na)

        ivector = self.coefficients["grayLevels"].to(dtype=torch.float64, device=self.device)  # shape (Ng,)
        jvector = torch.arange(
            1, self.P_glrlm.shape[2] + 1, dtype=torch.float64, device=self.device
        )  # shape (Nr,)

        # Delete columns that run lengths not present in the ROI
        emptyRunLenghts = torch.where(torch.sum(pr, dim=(0, 2)) == 0)[0]
        self.P_glrlm = torch_delete(self.P_glrlm, emptyRunLenghts, 2)
        jvector = torch_delete(jvector, emptyRunLenghts, 0)
        pr = torch_delete(pr, emptyRunLenghts, 1)

        self.coefficients["pr"] = pr
        self.coefficients["pg"] = pg
        self.coefficients["ivector"] = ivector
        self.coefficients["jvector"] = jvector

    def getShortRunEmphasisFeatureValue(self):
        pr = self.coefficients["pr"]
        jvector = self.coefficients["jvector"]
        Nr = self.coefficients["Nr"]

        sre = torch.sum((pr / (jvector[None, :, None] ** 2)), dim=1) / Nr
        return torch.nanmean(sre, dim=1)

    def getLongRunEmphasisFeatureValue(self):
        pr = self.coefficients["pr"]
        jvector = self.coefficients["jvector"]
        Nr = self.coefficients["Nr"]

        lre = torch.sum((pr * (jvector[None, :, None] ** 2)), dim=1) / Nr
        return torch.nanmean(lre, dim=1)

    def getGrayLevelNonUniformityFeatureValue(self):
        pg = self.coefficients["pg"]
        Nr = self.coefficients["Nr"]

        gln = torch.sum((pg**2), dim=1) / Nr
        return torch.nanmean(gln, dim=1)

    def getGrayLevelNonUniformityNormalizedFeatureValue(self):
        pg = self.coefficients["pg"]
        Nr = self.coefficients["Nr"]

        glnn = torch.sum(pg**2, dim=1) / (Nr**2)
        return torch.nanmean(glnn, dim=1)

    def getRunLengthNonUniformityFeatureValue(self):
        pr = self.coefficients["pr"]
        Nr = self.coefficients["Nr"]

        rln = torch.sum((pr**2), dim=1) / Nr
        return torch.nanmean(rln, dim=1)

    def getRunLengthNonUniformityNormalizedFeatureValue(self):
        pr = self.coefficients["pr"]
        Nr = self.coefficients["Nr"]

        rlnn = torch.sum((pr**2), dim=1) / Nr**2
        return torch.nanmean(rlnn, dim=1)

    def getRunPercentageFeatureValue(self):
        pr = self.coefficients["pr"]
        jvector = self.coefficients["jvector"]
        Nr = self.coefficients["Nr"]

        Np = torch.sum(pr * jvector[None, :, None], dim=1)  # shape (Nvox, Na)

        rp = Nr / Np
        return torch.nanmean(rp, dim=1)

    def getGrayLevelVarianceFeatureValue(self):
        ivector = self.coefficients["ivector"]
        Nr = self.coefficients["Nr"]
        pg = (
            self.coefficients["pg"] / Nr[:, None, :]
        )  # divide by Nr to get the normalized matrix

        u_i = torch.sum(pg * ivector[None, :, None], dim=1, keepdims=True)
        glv = torch.sum(pg * (ivector[None, :, None] - u_i) ** 2, dim=1)
        return torch.nanmean(glv, dim=1)

    def getRunVarianceFeatureValue(self):
        jvector = self.coefficients["jvector"]
        Nr = self.coefficients["Nr"]
        pr = (
            self.coefficients["pr"] / Nr[:, None, :]
        )  # divide by Nr to get the normalized matrix

        u_j = torch.sum(pr * jvector[None, :, None], dim=1, keepdims=True)
        rv = torch.sum(pr * (jvector[None, :, None] - u_j) ** 2, dim=1)
        return torch.nanmean(rv, dim=1)

    def getRunEntropyFeatureValue(self):
        eps = torch.finfo(torch.float32).eps
        Nr = self.coefficients["Nr"]
        p_glrlm = (
            self.P_glrlm / Nr[:, None, None, :]
        )  # divide by Nr to get the normalized matrix

        re = -torch.sum(p_glrlm * torch.log2(p_glrlm + eps), dim=(1, 2))
        return torch.nanmean(re, dim=1)

    def getLowGrayLevelRunEmphasisFeatureValue(self):
        pg = self.coefficients["pg"]
        ivector = self.coefficients["ivector"]
        Nr = self.coefficients["Nr"]

        lglre = torch.sum((pg / (ivector[None, :, None] ** 2)), dim=1) / Nr
        return torch.nanmean(lglre, dim=1)

    def getHighGrayLevelRunEmphasisFeatureValue(self):
        pg = self.coefficients["pg"]
        ivector = self.coefficients["ivector"]
        Nr = self.coefficients["Nr"]

        hglre = torch.sum((pg * (ivector[None, :, None] ** 2)), dim=1) / Nr
        return torch.nanmean(hglre, dim=1)

    def getShortRunLowGrayLevelEmphasisFeatureValue(self):
        ivector = self.coefficients["ivector"]
        jvector = self.coefficients["jvector"]
        Nr = self.coefficients["Nr"]

        srlgle = (
            torch.sum(
                (
                    self.P_glrlm
                    / (
                        (ivector[None, :, None, None] ** 2)
                        * (jvector[None, None, :, None] ** 2)
                    )
                ),
                dim=(1, 2),
            )
            / Nr
        )
        return torch.nanmean(srlgle, dim=1)

    def getShortRunHighGrayLevelEmphasisFeatureValue(self):
        ivector = self.coefficients["ivector"]
        jvector = self.coefficients["jvector"]
        Nr = self.coefficients["Nr"]

        srhgle = (
            torch.sum(
                (
                    self.P_glrlm
                    * (ivector[None, :, None, None] ** 2)
                    / (jvector[None, None, :, None] ** 2)
                ),
                dim=(1, 2),
            )
            / Nr
        )
        return torch.nanmean(srhgle, dim=1)

    def getLongRunLowGrayLevelEmphasisFeatureValue(self):
        ivector = self.coefficients["ivector"]
        jvector = self.coefficients["jvector"]
        Nr = self.coefficients["Nr"]

        lrlgle = (
            torch.sum(
                (
                    self.P_glrlm
                    * (jvector[None, None, :, None] ** 2)
                    / (ivector[None, :, None, None] ** 2)
                ),
                (1, 2),
            )
            / Nr
        )
        return torch.nanmean(lrlgle, dim=1)

    def getLongRunHighGrayLevelEmphasisFeatureValue(self):
        ivector = self.coefficients["ivector"]
        jvector = self.coefficients["jvector"]
        Nr = self.coefficients["Nr"]

        lrhgle = (
            torch.sum(
                (
                    self.P_glrlm
                    * (
                        (jvector[None, None, :, None] ** 2)
                        * (ivector[None, :, None, None] ** 2)
                    )
                ),
                dim=(1, 2),
            )
            / Nr
        )
        return torch.nanmean(lrhgle, dim=1)
