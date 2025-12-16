from __future__ import annotations

import torch  # type: ignore[import]

from . import base, deprecated
from .cmatrices import generate_angles_torch
from .utils import nanstd, nanmin, nanmax, nanquantile


class RadiomicsFirstOrder(base.RadiomicsFeaturesBase):
    def __init__(self, imageArray, maskArray, **kwargs):
        super().__init__(imageArray, maskArray, **kwargs)

        self.pixelSpacing = torch.tensor(kwargs.get("spacing", [1, 1, 1]), device=self.device, dtype=self.dtype)
        self.voxelArrayShift = kwargs.get("voxelArrayShift", 0)
        self.discretizedImageArray = self._applyBinning(self.imageArray.detach().clone())

    def _initVoxelBasedCalculation(self):
        super()._initVoxelBasedCalculation()

        kernelRadius = torch.tensor(self.settings.get("kernelRadius", 1), device=self.device, dtype=torch.int64)

        # Get the size of the input, which depends on whether it is in masked mode or not
        if self.masked:
            size = (
                torch.max(self.labelledVoxelCoordinates, 1).values
                - torch.min(self.labelledVoxelCoordinates, 1).values
                + 1
            )
        else:
            size = torch.tensor(self.imageArray.shape, device=self.device, dtype=torch.int64)

        # Take the minimum size along each dimension from either the size of the ROI or the kernel
        boundingBoxSize = torch.minimum(size, kernelRadius * 2 + 1)

        # Calculate the offsets, which can be used to generate a list of kernel Coordinates. Shape (Nd, Nk)
        self.kernelOffsets = generate_angles_torch(
            boundingBoxSize,
            torch.arange(1, kernelRadius + 1, device=self.device, dtype=torch.int64),
            True,  # Bi-directional
            self.settings.get("force2D", False),
            self.settings.get("force2Ddimension", 0),
        )
        center_offset = torch.zeros(
            (1, self.kernelOffsets.shape[1]),
            device=self.kernelOffsets.device,
            dtype=self.kernelOffsets.dtype,
        )
        self.kernelOffsets = torch.cat((self.kernelOffsets, center_offset), dim=0)  # add center voxel
        self.kernelOffsets = self.kernelOffsets.transpose(1, 0)  # Transpose to (Nd, Nk+1)

        self.imageArray = self.imageArray.float()
        self.imageArray[~self.maskArray] = float("nan")
        pad_width = (kernelRadius, kernelRadius) * self.imageArray.dim()
        self.imageArray = torch.nn.functional.pad(
            self.imageArray,
            pad=pad_width,
            mode="constant",
            value=float("nan"),
        )
        self.maskArray = torch.nn.functional.pad(
            self.maskArray,
            pad=pad_width,
            mode="constant",
            value=False,
        )

    def _initCalculation(self, voxelCoordinates=None):

        if voxelCoordinates is None:
            self.targetVoxelArray = (
                self.imageArray[self.maskArray].float().view(1, -1)
            )

            counts = torch.bincount(self.discretizedImageArray[self.maskArray],
            minlength=self.coefficients["Ng"] + 1)
            p_i = counts[self.coefficients["grayLevels"]].view(1, -1)
        else:
            # voxelCoordinates shape (Nd, Nvox)
            voxelCoordinates = voxelCoordinates.detach().clone() + self.settings.get(
                "kernelRadius", 1
            )  # adjust for padding
            kernelCoords = (
                self.kernelOffsets[:, None, :] + voxelCoordinates[:, :, None]
            )  # Shape (Nd, Nvox, Nk)
            kernelCoords = tuple(kernelCoords)  # shape (Nd, (Nvox, Nk))

            self.targetVoxelArray = self.imageArray[kernelCoords]  # shape (Nvox, Nk)

            p_i = torch.empty(
                (voxelCoordinates.shape[1], len(self.coefficients["grayLevels"])), device=self.device
            )  # shape (Nvox, Ng)
            for gl_idx, gl in enumerate(self.coefficients["grayLevels"]):
                p_i[:, gl_idx] = torch.nansum(
                    self.discretizedImageArray[kernelCoords] == gl, 1
                )

        sumBins = torch.sum(p_i, 1, keepdims=True).float()
        sumBins[sumBins == 0] = 1  # Prevent division by 0 errors
        p_i = p_i.float() / sumBins
        self.coefficients["p_i"] = p_i

        self.logger.debug("First order feature class initialized")

    @staticmethod
    def _moment(a, moment=1):
        if moment == 1:
            return torch.tensor(0.0, device=a.device)
        mn = torch.nanmean(a, 1, keepdims=True)
        s = torch.pow((a - mn), moment)
        return torch.nanmean(s, 1)

    def getEnergyFeatureValue(self):
        shiftedParameterArray = self.targetVoxelArray + self.voxelArrayShift

        return torch.nansum(shiftedParameterArray**2, 1)

    def getTotalEnergyFeatureValue(self):
        cubicMMPerVoxel = torch.prod(self.pixelSpacing)

        return self.getEnergyFeatureValue() * cubicMMPerVoxel

    def getEntropyFeatureValue(self):
        p_i = self.coefficients["p_i"]

        eps = torch.finfo(torch.float32).eps
        return -1.0 * torch.sum(p_i * torch.log2(p_i + eps), 1)

    def getMinimumFeatureValue(self):
        return nanmin(self.targetVoxelArray, 1).values

    def get10PercentileFeatureValue(self):
        return nanquantile(self.targetVoxelArray, 0.1, dim=1)

    def get90PercentileFeatureValue(self):
        return nanquantile(self.targetVoxelArray, 0.9, dim=1)

    def getMaximumFeatureValue(self):
        return nanmax(self.targetVoxelArray, 1).values

    def getMeanFeatureValue(self):
        return torch.nanmean(self.targetVoxelArray, 1)

    def getMedianFeatureValue(self):
        return torch.nanmedian(self.targetVoxelArray, 1).values

    def getInterquartileRangeFeatureValue(self):
        return nanquantile(self.targetVoxelArray, 0.75, dim=1) - nanquantile(
            self.targetVoxelArray, 0.25, dim=1
        )

    def getRangeFeatureValue(self):
        return (
            nanmax(self.targetVoxelArray, dim=1).values
            - nanmin(self.targetVoxelArray, dim=1).values
        )

    def getMeanAbsoluteDeviationFeatureValue(self):
        u_x = torch.nanmean(self.targetVoxelArray, dim=1, keepdims=True)
        return torch.nanmean(torch.abs(self.targetVoxelArray - u_x), dim=1)

    def getRobustMeanAbsoluteDeviationFeatureValue(self):
        prcnt10 = self.get10PercentileFeatureValue()
        prcnt90 = self.get90PercentileFeatureValue()
        percentileArray = self.targetVoxelArray.clone()

        # First get a mask for all valid voxels
        valid = ~torch.isnan(percentileArray)
        # Then, update the mask to reflect all valid voxels that are outside the the closed 10-90th percentile range
        cond = ((percentileArray - prcnt10[:, None]) < 0) | (
            (percentileArray - prcnt90[:, None]) > 0
        )
        msk = valid & cond

        # Finally, exclude the invalid voxels by setting them to torch.nan.
        percentileArray[msk] = torch.nan

        return torch.nanmean(
            torch.abs(
                percentileArray - torch.nanmean(percentileArray, dim=1, keepdims=True)
            ),
            dim=1,
        )

    def getRootMeanSquaredFeatureValue(self):
        # If no voxels are segmented, prevent division by 0 and return 0
        if self.targetVoxelArray.numel() == 0:
            return 0

        shiftedParameterArray = self.targetVoxelArray + self.voxelArrayShift
        Nvox = torch.sum(~torch.isnan(self.targetVoxelArray), dim=1).to(torch.float)
        return torch.sqrt(torch.nansum(shiftedParameterArray**2, dim=1) / Nvox)

    @deprecated
    def getStandardDeviationFeatureValue(self):
        return nanstd(self.targetVoxelArray, dim=1)

    def getSkewnessFeatureValue(self):
        m2 = self._moment(self.targetVoxelArray, 2)
        m3 = self._moment(self.targetVoxelArray, 3)

        m2[m2 == 0] = 1  # Flat Region, prevent division by 0 errors
        m3[m2 == 0] = 0  # ensure Flat Regions are returned as 0

        return m3 / m2**1.5

    def getKurtosisFeatureValue(self):
        m2 = self._moment(self.targetVoxelArray, 2)
        m4 = self._moment(self.targetVoxelArray, 4)

        m2[m2 == 0] = 1  # Flat Region, prevent division by 0 errors
        m4[m2 == 0] = 0  # ensure Flat Regions are returned as 0

        return m4 / m2**2.0

    def getVarianceFeatureValue(self):
        return nanstd(self.targetVoxelArray, 1) ** 2

    def getUniformityFeatureValue(self):
        p_i = self.coefficients["p_i"]
        return torch.nansum(p_i**2, 1)