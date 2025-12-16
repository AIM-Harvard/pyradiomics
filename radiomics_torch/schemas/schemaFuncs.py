from __future__ import annotations

import pywt

from radiomics_torch import getFeatureClasses, getImageTypes

featureClasses = getFeatureClasses()
imageTypes = getImageTypes()


def checkWavelet(value, _rule_obj, _path):
    if not isinstance(value, str):
        msg = "Wavelet not expected type (str)"
        raise TypeError(msg)
    wavelist = pywt.wavelist()
    if value not in wavelist:
        msg = f'Wavelet "{value}" not available in pyWavelets {wavelist}'
        raise ValueError(msg)
    return True


def checkInterpolator(value, _rule_obj, _path):
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
            msg = f'Interpolator value "{value}" not valid, possible values: {enum}'
            raise ValueError(msg)
    elif isinstance(value, int):
        if value < 1 or value > 10:
            msg = f"Intepolator value {int(value)}, must be in range of [1-10]"
            raise ValueError(msg)
    else:
        msg = "Interpolator not expected type (str or int)"
        raise TypeError(msg)
    return True


def checkWeighting(value, _rule_obj, _path):
    if value is None:
        return True
    if isinstance(value, str):
        enum = ["euclidean", "manhattan", "infinity", "no_weighting"]
        if value not in enum:
            msg = f'WeightingNorm value "{value}" not valid, possible values: {enum}'
            raise ValueError(msg)
    else:
        msg = "WeightingNorm not expected type (str or None)"
        raise TypeError(msg)
    return True


def checkFeatureClass(value, _rule_obj, _path):
    if value is None:
        msg = "featureClass dictionary cannot be None value"
        raise TypeError(msg)
    for className, features in value.items():
        if className not in featureClasses:
            msg = f"Feature Class {className} is not recognized. Available feature classes are {list(featureClasses.keys())}"
            raise ValueError(msg)
        if features is not None:
            if not isinstance(features, list):
                msg = f"Value of feature class {className} not expected type (list)"
                raise TypeError(msg)
            unrecognizedFeatures = set(features) - set(
                featureClasses[className].getFeatureNames()
            )
            if len(unrecognizedFeatures) > 0:
                msg = f"Feature Class {className} contains unrecognized features: {unrecognizedFeatures!s}"
                raise ValueError(msg)

    return True


def checkImageType(value, _rule_obj, _path):
    if value is None:
        msg = "imageType dictionary cannot be None value"
        raise TypeError(msg)

    for im_type in value:
        if im_type not in imageTypes:
            msg = f"Image Type {im_type} is not recognized. Available image types are {imageTypes}"
            raise ValueError(msg)

    return True
