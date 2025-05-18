from __future__ import annotations

import pywt

from radiomics import getFeatureClasses, getImageTypes

featureClasses = getFeatureClasses()
imageTypes = getImageTypes()


def checkWavelet(value, rule_obj, path):
    if not isinstance(value, str):
        raise TypeError("Wavelet not expected type (str)")
    wavelist = pywt.wavelist()
    if value not in wavelist:
        raise ValueError(f'Wavelet "{value}" not available in pyWavelets {wavelist}')
    return True


def checkInterpolator(value, rule_obj, path):
    if value is None:
        return True
    if isinstance(value, str):
        enum = {
            "sitkNearestNeighbor",
            "sitkLinear",
            "sitkBSpline",
            "sitkGaussian",
            "sitkLabelGaussian",
            "sitkHammingWindowedSinc",
            "sitkCosineWindowedSinc",
            "sitkWelchWindowedSinc",
            "sitkLanczosWindowedSinc",
            "sitkBlackmanWindowedSinc",
        }
        if value not in enum:
            raise ValueError(
                f'Interpolator value "{value}" not valid, possible values: {enum}'
            )
    elif isinstance(value, int):
        if value < 1 or value > 10:
            raise ValueError(
                f"Intepolator value {int(value)}, must be in range of [1-10]"
            )
    else:
        raise TypeError("Interpolator not expected type (str or int)")
    return True


def checkWeighting(value, rule_obj, path):
    if value is None:
        return True
    if isinstance(value, str):
        enum = ["euclidean", "manhattan", "infinity", "no_weighting"]
        if value not in enum:
            raise ValueError(
                f'WeightingNorm value "{value}" not valid, possible values: {enum}'
            )
    else:
        raise TypeError("WeightingNorm not expected type (str or None)")
    return True


def checkFeatureClass(value, rule_obj, path):
    global featureClasses
    if value is None:
        raise TypeError("featureClass dictionary cannot be None value")
    for className, features in value.items():
        if className not in featureClasses.keys():
            msg = f"Feature Class {className} is not recognized. Available feature classes are {list(featureClasses.keys())}"
            raise ValueError(msg)
        if features is not None:
            if not isinstance(features, list):
                raise TypeError(
                    f"Value of feature class {className} not expected type (list)"
                )
            unrecognizedFeatures = set(features) - set(
                featureClasses[className].getFeatureNames()
            )
            if len(unrecognizedFeatures) > 0:
                msg = f"Feature Class {className} contains unrecognized features: {unrecognizedFeatures!s}"
                raise ValueError(msg)

    return True


def checkImageType(value, rule_obj, path):
    global imageTypes
    if value is None:
        raise TypeError("imageType dictionary cannot be None value")

    for im_type in value:
        if im_type not in imageTypes:
            msg = f"Image Type {im_type} is not recognized. Available image types are {imageTypes}"
            raise ValueError(msg)

    return True
