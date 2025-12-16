from __future__ import annotations

import inspect
import logging
import math
import traceback

import torch
from numpy import ndarray
from torch import Tensor

from . import imageoperations


class RadiomicsFeaturesBase:
    def __init__(self, imageArray: Tensor | ndarray, maskArray: Tensor | ndarray, **kwargs):
        self.logger = logging.getLogger(self.__module__)
        self.logger.debug("Initializing feature class")

        if imageArray is None or maskArray is None:
            msg = "Missing input image or mask"
            raise ValueError(msg)

        self.settings = kwargs

        self.device = imageArray.device
        assert maskArray.device == self.device, "Image and mask must be on the same device"
        self.dtype = imageArray.dtype

        self.label = kwargs.get("label", 1)
        self.voxelBased = kwargs.get("voxelBased", False)

        # hardcode spacing
        self.spacing = torch.tensor([1.0, 1.0, 1.0], device=self.device)  # z, y, x 

        self.coefficients = {}

        # all features are disabled by default
        self.enabledFeatures = {}
        self.featureValues = {}

        self.featureNames = self.getFeatureNames()

        self.imageArray = torch.as_tensor(imageArray, device=self.device, dtype=self.dtype)
        self.maskArray = torch.as_tensor(maskArray, device=self.device, dtype=torch.uint8)


        if self.voxelBased:
            self._initVoxelBasedCalculation()
        else:
            self._initSegmentBasedCalculation()

    def _initSegmentBasedCalculation(self):
        self.maskArray = self.maskArray == self.label  # boolean array

    def _initVoxelBasedCalculation(self):
        self.masked = self.settings.get("maskedKernel", True)

        maskArray = self.maskArray == self.label  # boolean array
        self.labelledVoxelCoordinates = torch.nonzero(maskArray).T

        # Set up the mask array for the gray value discretization
        if self.masked:
            self.maskArray = maskArray
        else:
            # This will cause the discretization to use the entire image
            self.maskArray = torch.ones(self.imageArray.shape, device=self.device, dtype=torch.bool)

    def _initCalculation(self, voxelCoordinates=None):
        ...

    def _applyBinning(self, matrix):
        matrix, _ = imageoperations.binImage(matrix, self.maskArray, **self.settings)
        self.coefficients["grayLevels"] = torch.unique(matrix[self.maskArray])
        self.coefficients["Ng"] = self.coefficients["grayLevels"].max().item()  # max gray level in the ROI
        return matrix

    def enableFeatureByName(self, featureName, enable=True):
        if featureName not in self.featureNames:
            raise LookupError("Feature not found: " + featureName)
        if self.featureNames[featureName]:
            self.logger.warning(
                "Feature %s is deprecated, use with caution!", featureName
            )
        self.enabledFeatures[featureName] = enable

    def enableAllFeatures(self):
        for featureName, is_deprecated in self.featureNames.items():
            # only enable non-deprecated features here
            if not is_deprecated:
                self.enableFeatureByName(featureName, True)

    def disableAllFeatures(self):
        self.enabledFeatures = {}
        self.featureValues = {}

    @classmethod
    def getFeatureNames(cls):
        attributes = inspect.getmembers(cls)
        return {
            a[0][3:-12]: getattr(a[1], "_is_deprecated", False)
            for a in attributes
            if a[0].startswith("get") and a[0].endswith("FeatureValue")
        }

    def execute(self):
        if len(self.enabledFeatures) == 0:
            self.enableAllFeatures()

        if self.voxelBased:
            self._calculateVoxels()
        else:
            self._calculateSegment()

        return self.featureValues

    def _calculateVoxels(self):
        initValue = self.settings.get("initValue", 0)
        voxelBatch = self.settings.get("voxelBatch", -1)

        # Initialize the output with empty tensors
        for feature, enabled in self.enabledFeatures.items():
            if enabled:
                self.featureValues[feature] = torch.full(
                    self.imageArray.shape, initValue, device=self.device, dtype=self.dtype
                )

        # Calculate the feature values for all enabled features
        voxel_count = self.labelledVoxelCoordinates.shape[1]
        voxel_batch_idx = 0
        if voxelBatch < 0:
            voxelBatch = voxel_count
        n_batches = math.ceil(float(voxel_count) / voxelBatch)
        while voxel_batch_idx < voxel_count:
            self.logger.debug(
                "Calculating voxel batch no. %i/%i",
                voxel_batch_idx // voxelBatch + 1,
                n_batches,
            )
            voxelCoords = self.labelledVoxelCoordinates[
                :, voxel_batch_idx : voxel_batch_idx + voxelBatch
            ]
            # Calculate the feature values for the current kernel
            for success, featureName, featureValue in self._calculateFeatures(
                voxelCoords
            ):
                if success:
                    self.featureValues[featureName][
                        tuple(voxelCoords)
                    ] = featureValue

            voxel_batch_idx += voxelBatch

    def _calculateSegment(self):
        # Get the feature values using the current segment.
        for _success, featureName, featureValue in self._calculateFeatures():
            # Always store the result. In case of an error, featureValue will be NaN
            self.featureValues[featureName] = torch.squeeze(featureValue)

    def _calculateFeatures(self, voxelCoordinates=None):
        # Initialize the calculation
        # This function serves to calculate the texture matrices where applicable
        self._initCalculation(voxelCoordinates)

        self.logger.debug("Calculating features")
        for feature, enabled in self.enabledFeatures.items():
            if enabled:
                try:
                    # Use getattr to get the feature calculation methods, then use '()' to evaluate those methods
                    yield True, feature, getattr(self, f"get{feature}FeatureValue")()
                except DeprecationWarning as deprecatedFeature:
                    # Add a debug log message, as a warning is usually shown and would entail a too verbose output
                    self.logger.debug(
                        "Feature %s is deprecated: %s",
                        feature,
                        deprecatedFeature.args[0],
                    )
                except Exception:
                    self.logger.error("FAILED: %s", traceback.format_exc())
                    yield False, feature, torch.tensor(torch.nan)
