#!/usr/bin/env python
from __future__ import annotations

import sys

import numpy as np
import SimpleITK as sitk

from radiomics import (
    firstorder,
    getTestCase,
    glcm,
    glrlm,
    glszm,
    imageoperations,
    shape,
)

# testBinWidth = 25 this is the default bin size
# testResampledPixelSpacing = [3,3,3] no resampling for now.

# Get some test data

# Download the test case to temporary files and return it's location. If already downloaded, it is not downloaded again,
# but it's location is still returned.
imageName, maskName = getTestCase("brain1")

if (
    imageName is None or maskName is None
):  # Something went wrong, in this case PyRadiomics will also log an error
    print("Error getting testcase!")
    sys.exit()

image = sitk.ReadImage(imageName)
mask = sitk.ReadImage(maskName)

applyLog = False
applyWavelet = False

# Setting for the feature calculation.
# Currently, resampling is disabled.
# Can be enabled by setting 'resampledPixelSpacing' to a list of 3 floats (new voxel size in mm for x, y and z)
settings = {
    "binWidth": 25,
    "interpolator": sitk.sitkBSpline,
    "resampledPixelSpacing": None,
}

#
# If enabled, resample image (resampled image is automatically cropped.
#
interpolator = settings.get("interpolator")
resampledPixelSpacing = settings.get("resampledPixelSpacing")
if interpolator is not None and resampledPixelSpacing is not None:
    image, mask = imageoperations.resampleImage(image, mask, **settings)

bb, correctedMask = imageoperations.checkMask(image, mask)
if correctedMask is not None:
    mask = correctedMask
image, mask = imageoperations.cropToTumorMask(image, mask, bb)

#
# Show the first order feature calculations
#
firstOrderFeatures = firstorder.RadiomicsFirstOrder(image, mask, **settings)

firstOrderFeatures.enableFeatureByName("Mean", True)
# firstOrderFeatures.enableAllFeatures()

print("Will calculate the following first order features: ")
for f in firstOrderFeatures.enabledFeatures:
    print("  ", f)
    print(getattr(firstOrderFeatures, f"get{f}FeatureValue").__doc__)

print("Calculating first order features...")
results = firstOrderFeatures.execute()
print("done")

print("Calculated first order features: ")
for key, val in results.items():
    print("  ", key, ":", val)

#
# Show Shape features
#
shapeFeatures = shape.RadiomicsShape(image, mask, **settings)
shapeFeatures.enableAllFeatures()

print("Will calculate the following Shape features: ")
for f in shapeFeatures.enabledFeatures:
    print("  ", f)
    print(getattr(shapeFeatures, f"get{f}FeatureValue").__doc__)

print("Calculating Shape features...")
results = shapeFeatures.execute()
print("done")

print("Calculated Shape features: ")
for key, val in results.items():
    print("  ", key, ":", val)

#
# Show GLCM features
#
glcmFeatures = glcm.RadiomicsGLCM(image, mask, **settings)
glcmFeatures.enableAllFeatures()

print("Will calculate the following GLCM features: ")
for f in glcmFeatures.enabledFeatures:
    print("  ", f)
    print(getattr(glcmFeatures, f"get{f}FeatureValue").__doc__)

print("Calculating GLCM features...")
results = glcmFeatures.execute()
print("done")

print("Calculated GLCM features: ")
for key, val in results.items():
    print("  ", key, ":", val)

#
# Show GLRLM features
#
glrlmFeatures = glrlm.RadiomicsGLRLM(image, mask, **settings)
glrlmFeatures.enableAllFeatures()

print("Will calculate the following GLRLM features: ")
for f in glrlmFeatures.enabledFeatures:
    print("  ", f)
    print(getattr(glrlmFeatures, f"get{f}FeatureValue").__doc__)

print("Calculating GLRLM features...")
results = glrlmFeatures.execute()
print("done")

print("Calculated GLRLM features: ")
for key, val in results.items():
    print("  ", key, ":", val)

#
# Show GLSZM features
#
glszmFeatures = glszm.RadiomicsGLSZM(image, mask, **settings)
glszmFeatures.enableAllFeatures()

print("Will calculate the following GLSZM features: ")
for f in glszmFeatures.enabledFeatures:
    print("  ", f)
    print(getattr(glszmFeatures, f"get{f}FeatureValue").__doc__)

print("Calculating GLSZM features...")
results = glszmFeatures.execute()
print("done")

print("Calculated GLSZM features: ")
for key, val in results.items():
    print("  ", key, ":", val)

#
# Show FirstOrder features, calculated on a LoG filtered image
#
if applyLog:
    sigmaValues = np.arange(5.0, 0.0, -0.5)[::1]
    for logImage, imageTypeName, inputKwargs in imageoperations.getLoGImage(
        image, mask, sigma=sigmaValues
    ):
        logFirstorderFeatures = firstorder.RadiomicsFirstOrder(
            logImage, mask, **inputKwargs
        )
        logFirstorderFeatures.enableAllFeatures()
        results = logFirstorderFeatures.execute()
        for key, val in results.items():
            laplacianFeatureName = f"{imageTypeName}_{key}"
            print("  ", laplacianFeatureName, ":", val)
#
# Show FirstOrder features, calculated on a wavelet filtered image
#
if applyWavelet:
    for (
        decompositionImage,
        decompositionName,
        inputKwargs,
    ) in imageoperations.getWaveletImage(image, mask):
        waveletFirstOrderFeaturs = firstorder.RadiomicsFirstOrder(
            decompositionImage, mask, **inputKwargs
        )
        waveletFirstOrderFeaturs.enableAllFeatures()
        results = waveletFirstOrderFeaturs.execute()
        print("Calculated firstorder features with wavelet ", decompositionName)
        for key, val in results.items():
            waveletFeatureName = f"{decompositionName!s}_{key}"
            print("  ", waveletFeatureName, ":", val)
