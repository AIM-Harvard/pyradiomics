from __future__ import annotations

import torch
import SimpleITK as sitk
import torch.nn.functional as F

from . import base, deprecated
from .cshapes import calculate_coefficients_torch


class RadiomicsShape(base.RadiomicsFeaturesBase):
    def __init__(self, imageArray, maskArray, **kwargs):
        super().__init__(imageArray, maskArray, **kwargs)

    def _initVoxelBasedCalculation(self):
        msg = "Shape features are not available in voxel-based mode"
        raise NotImplementedError(msg)

    def _initSegmentBasedCalculation(self):
        self.pixelSpacing = self.spacing

        # Pad inputMask to prevent index-out-of-range errors
        self.logger.debug("Padding the mask with 0s")

        pad = (1, 1, 1, 1, 1, 1)
        self.maskArray = F.pad(self.maskArray, pad, mode="constant", value=0)

        # Reassign self.maskArray using the now-padded self.inputMask
        self.maskArray = self.maskArray == self.label
        self.labelledVoxelCoordinates = torch.nonzero(self.maskArray != 0).T

        self.logger.debug("Pre-calculate Volume, Surface Area and Eigenvalues")

        # Volume, Surface Area and eigenvalues are pre-calculated
        # Compute Surface Area and volume
        self.SurfaceArea, self.Volume, self.diameters = calculate_coefficients_torch(
            self.maskArray, self.pixelSpacing
        )

        # Compute eigenvalues and -vectors
        self.Np = torch.tensor(self.labelledVoxelCoordinates.shape[1], device=self.device)
        coordinates = self.labelledVoxelCoordinates.T  # Transpose equals zip(*a)
        physicalCoordinates = coordinates * self.pixelSpacing[None, :]
        physicalCoordinates -= torch.mean(physicalCoordinates, dim=0)  # Centered at 0
        physicalCoordinates /= torch.sqrt(self.Np)
        covariance = torch.matmul(physicalCoordinates.T.clone(), physicalCoordinates)
        self.eigenValues = torch.linalg.eigvalsh(covariance)

        # Correct machine precision errors causing very small negative eigen values in case of some 2D segmentations
        machine_errors = torch.bitwise_and(self.eigenValues < 0, self.eigenValues > -1e-10)
        if torch.sum(machine_errors) > 0:
            self.logger.warning(
                "Encountered %d eigenvalues < 0 and > -1e-10, rounding to 0",
                torch.sum(machine_errors),
            )
            self.eigenValues[machine_errors] = 0

        self.eigenValues.sort()  # Sort the eigenValues from small to large

        self.logger.debug("Shape feature class initialized")

    def getMeshVolumeFeatureValue(self):
        return self.Volume

    def getVoxelVolumeFeatureValue(self):
        z, y, x = self.pixelSpacing
        return self.Np * (z * x * y)

    def getSurfaceAreaFeatureValue(self):
        return self.SurfaceArea

    def getSurfaceVolumeRatioFeatureValue(self):
        return self.SurfaceArea / self.Volume

    def getSphericityFeatureValue(self):
        return (36 * torch.pi * self.Volume**2) ** (1.0 / 3.0) / self.SurfaceArea

    @deprecated
    def getCompactness1FeatureValue(self):
        return self.Volume / (self.SurfaceArea ** (3.0 / 2.0) * torch.sqrt(torch.pi))

    @deprecated
    def getCompactness2FeatureValue(self):
        return (36.0 * torch.pi) * (self.Volume**2.0) / (self.SurfaceArea**3.0)

    @deprecated
    def getSphericalDisproportionFeatureValue(self):
        return self.SurfaceArea / (36 * torch.pi * self.Volume**2) ** (1.0 / 3.0)

    def getMaximum3DDiameterFeatureValue(self):
        return self.diameters[3]

    def getMaximum2DDiameterSliceFeatureValue(self):
        return self.diameters[0]

    def getMaximum2DDiameterColumnFeatureValue(self):
        return self.diameters[1]

    def getMaximum2DDiameterRowFeatureValue(self):
        return self.diameters[2]

    def getMajorAxisLengthFeatureValue(self):
        if self.eigenValues[2] < 0:
            self.logger.warning(
                "Major axis eigenvalue negative! (%g)", self.eigenValues[2]
            )
            return torch.nan
        return torch.sqrt(self.eigenValues[2]) * 4

    def getMinorAxisLengthFeatureValue(self):
        if self.eigenValues[1] < 0:
            self.logger.warning(
                "Minor axis eigenvalue negative! (%g)", self.eigenValues[1]
            )
            return torch.nan
        return torch.sqrt(self.eigenValues[1]) * 4

    def getLeastAxisLengthFeatureValue(self):
        if self.eigenValues[0] < 0:
            self.logger.warning(
                "Least axis eigenvalue negative! (%g)", self.eigenValues[0]
            )
            return torch.nan
        return torch.sqrt(self.eigenValues[0]) * 4

    def getElongationFeatureValue(self):
        if self.eigenValues[1] < 0 or self.eigenValues[2] < 0:
            self.logger.warning(
                "Elongation eigenvalue negative! (%g, %g)",
                self.eigenValues[1],
                self.eigenValues[2],
            )
            return torch.nan
        return torch.sqrt(self.eigenValues[1] / self.eigenValues[2])

    def getFlatnessFeatureValue(self):
        if self.eigenValues[0] < 0 or self.eigenValues[2] < 0:
            self.logger.warning(
                "Elongation eigenvalue negative! (%g, %g)",
                self.eigenValues[0],
                self.eigenValues[2],
            )
            return torch.nan
        return torch.sqrt(self.eigenValues[0] / self.eigenValues[2])
