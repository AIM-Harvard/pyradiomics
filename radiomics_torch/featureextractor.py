from __future__ import annotations

import collections
import logging
from itertools import chain

import torch

from . import (
    getFeatureClasses,
    getImageTypes,
    imageoperations,
)

logger = logging.getLogger(__name__)


class RadiomicsFeatureExtractor:
    r"""
    By default, all features in all feature classes are enabled.
    By default, only `Original` input image is enabled (No filter applied).
    """

    def __init__(self, *args, **kwargs):
        self.settings = {}
        self.enabledImagetypes = {}
        self.enabledFeatures = {}

        self.featureClassNames = list(getFeatureClasses().keys())

        # Set default settings and update with and changed settings contained in kwargs
        self.settings = self._getDefaultSettings()
        logger.info("No valid config parameter, using defaults: %s", self.settings)

        self.enabledImagetypes = {"Original": {}}
        logger.info("Enabled image types: %s", self.enabledImagetypes)

        for featureClassName in self.featureClassNames:
            if featureClassName == "shape2D":  # Do not enable shape2D by default
                continue
            # skip ngtdm features for now
            self.enabledFeatures[featureClassName] = []
        print("Enabled features: %s", self.enabledFeatures)
        logger.info("Enabled features: %s", self.enabledFeatures)

        if len(kwargs) > 0:
            logger.info("Applying custom setting overrides: %s", kwargs)
            self.settings.update(kwargs)
            logger.debug("Settings: %s", self.settings)

        if self.settings.get("binCount", None) is not None:
            logger.warning(
                "Fixed bin Count enabled! However, we recommend using a fixed bin Width. See "
                "http://pyradiomics.readthedocs.io/en/latest/faq.html#radiomics-fixed-bin-width for more "
                "details"
            )

    @staticmethod
    def _getDefaultSettings():
        return {
            "minimumROIDimensions": 2,
            "minimumROISize": None,  # Skip testing the ROI size by default
            "normalize": False,
            "normalizeScale": 1,
            "removeOutliers": None,
            "resampledPixelSpacing": None,  # No resampling by default
            "interpolator": "sitkBSpline",  # Alternative: sitk.sitkBSpline
            "preCrop": False,
            "padDistance": 5,
            "distances": [1],
            "force2D": False,
            "force2Ddimension": 0,
            "resegmentRange": None,  # No resegmentation by default
            "label": 1,
            "additionalInfo": True,
        }

    def execute(
        self,
        imageArray,
        maskArray,
        label=None,
        label_channel=None,
        voxelBased=False,
    ):
        _settings = self.settings.copy()

        if label is not None:
            _settings["label"] = label
        else:
            label = _settings.get("label", 1)

        if label_channel is not None:
            _settings["label_channel"] = label_channel

        if voxelBased:
            _settings["voxelBased"] = True
            kernelRadius = _settings.get("kernelRadius", 1)
            logger.info("Starting voxel based extraction")
        else:
            kernelRadius = 0

        logger.info("Calculating features with label: %d", label)
        logger.debug("Enabled images types: %s", self.enabledImagetypes)
        logger.debug("Enabled features: %s", self.enabledFeatures)
        logger.debug("Current settings: %s", _settings)

        featureVector = collections.OrderedDict()

        coords = torch.nonzero(maskArray == label, as_tuple=False)
        z_min = int(coords[:, 0].min())
        z_max = int(coords[:, 0].max())
        y_min = int(coords[:, 1].min())
        y_max = int(coords[:, 1].max())
        x_min = int(coords[:, 2].min())
        x_max = int(coords[:, 2].max())

        boundingBox = (z_min, z_max, y_min, y_max, x_min, x_max)

        logger.debug("Image and Mask loaded and valid, starting extraction")

        if not voxelBased:
            featureVector.update(
                self.computeShape(imageArray, maskArray, boundingBox, **_settings)
            )

        logger.debug("Creating image type iterator")
        imageGenerators = []
        for imageType, customKwargs in self.enabledImagetypes.items():
            args = _settings.copy()
            args.update(customKwargs)
            msg = f'Adding image type "{imageType}" with custom settings: {customKwargs!s}'
            logger.info(msg)
            imageGenerators = chain(
                imageGenerators,
                getattr(imageoperations, f"get{imageType}Image")(imageArray, maskArray, **args),
            )

        logger.debug("Extracting features")
        # Calculate features for all (filtered) images in the generator
        for originputImage, imageTypeName, inputKwargs in imageGenerators:
            logger.info("Calculating features for %s image", imageTypeName)
            inputImage, inputMask = imageoperations.cropToTumorMask(
                originputImage, maskArray, boundingBox, padDistance=kernelRadius
            )
            featureVector.update(
                self.computeFeatures(
                    inputImage, inputMask, imageTypeName, **inputKwargs
                )
            )

        logger.debug("Features extracted")

        return featureVector

    def computeShape(self, imageArray, maskArray, boundingBox, **kwargs):
        featureVector = collections.OrderedDict()

        enabledFeatures = self.enabledFeatures

        croppedImage, croppedMask = imageoperations.cropToTumorMask(
            imageArray, maskArray, boundingBox
        )

        # Define temporary function to compute shape features
        def compute(shape_type):
            logger.info("Computing %s", shape_type)
            featureNames = enabledFeatures[shape_type]
            shapeClass = getFeatureClasses()[shape_type](
                croppedImage, croppedMask, **kwargs
            )

            if featureNames is not None:
                for feature in featureNames:
                    shapeClass.enableFeatureByName(feature)

            for featureName, featureValue in shapeClass.execute().items():
                newFeatureName = f"original_{shape_type}_{featureName}"
                featureVector[newFeatureName] = featureValue

        Nd = maskArray.dim()
        if "shape" in enabledFeatures:
            if Nd == 3:
                compute("shape")
            else:
                logger.warning(
                    "Shape features are only available 3D input (for 2D input, use shape2D). Found %iD input",
                    Nd,
                )

        return featureVector

    def computeFeatures(self, image, mask, imageTypeName, **kwargs):
        featureVector = collections.OrderedDict()
        featureClasses = getFeatureClasses()

        enabledFeatures = self.enabledFeatures

        # Calculate feature classes
        for featureClassName, featureNames in enabledFeatures.items():
            # Handle calculation of shape features separately
            if featureClassName.startswith("shape"):
                continue

            if featureClassName in featureClasses:
                logger.info("Computing %s", featureClassName)

                featureClass = featureClasses[featureClassName](image, mask, **kwargs)

                if featureNames is not None:
                    for feature in featureNames:
                        featureClass.enableFeatureByName(feature)

                for featureName, featureValue in featureClass.execute().items():
                    newFeatureName = f"{imageTypeName}_{featureClassName}_{featureName}"
                    featureVector[newFeatureName] = featureValue

        return featureVector

    def enableAllImageTypes(self):
        logger.debug("Enabling all image types")
        for imageType in getImageTypes():
            self.enabledImagetypes[imageType] = {}
        logger.debug("Enabled images types: %s", self.enabledImagetypes)

    def disableAllImageTypes(self):
        logger.debug("Disabling all image types")
        self.enabledImagetypes = {}

    def enableImageTypeByName(self, imageType, enabled=True, customArgs=None):
        if imageType not in getImageTypes():
            logger.warning("Image type %s is not recognized", imageType)
            return

        if enabled:
            if customArgs is None:
                customArgs = {}
                logger.debug(
                    "Enabling image type %s (no additional custom settings)", imageType
                )
            else:
                logger.debug(
                    "Enabling image type %s (additional custom settings: %s)",
                    imageType,
                    customArgs,
                )
            self.enabledImagetypes[imageType] = customArgs
        elif imageType in self.enabledImagetypes:
            logger.debug("Disabling image type %s", imageType)
            del self.enabledImagetypes[imageType]
        logger.debug("Enabled images types: %s", self.enabledImagetypes)

    def enableImageTypes(self, **enabledImagetypes):
        logger.debug("Updating enabled images types with %s", enabledImagetypes)
        self.enabledImagetypes.update(enabledImagetypes)
        logger.debug("Enabled images types: %s", self.enabledImagetypes)

    def enableAllFeatures(self):
        logger.debug("Enabling all features in all feature classes")
        for featureClassName in self.featureClassNames:
            self.enabledFeatures[featureClassName] = []
        logger.debug("Enabled features: %s", self.enabledFeatures)

    def disableAllFeatures(self):
        logger.debug("Disabling all feature classes")
        self.enabledFeatures = {}

    def enableFeatureClassByName(self, featureClass, enabled=True):
        if featureClass not in self.featureClassNames:
            logger.warning("Feature class %s is not recognized", featureClass)
            return

        if enabled:
            logger.debug("Enabling all features in class %s", featureClass)
            self.enabledFeatures[featureClass] = []
        elif featureClass in self.enabledFeatures:
            logger.debug("Disabling feature class %s", featureClass)
            del self.enabledFeatures[featureClass]
        logger.debug("Enabled features: %s", self.enabledFeatures)

    def enableFeaturesByName(self, **enabledFeatures):
        logger.debug("Updating enabled features with %s", enabledFeatures)
        self.enabledFeatures.update(enabledFeatures)
        logger.debug("Enabled features: %s", self.enabledFeatures)
