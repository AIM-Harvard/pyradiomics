#!/usr/bin/env python

import sys

import SimpleITK as sitk

from radiomics import imageoperations

image = sitk.ReadImage(sys.argv[1])
mask = sitk.ReadImage(sys.argv[2])

# Resamples and crops onto bounding box defined by the label
(ii, im) = imageoperations.resampleImage(image, mask, resampledImageSpaceing=[2, 2, 2], label=1, padDistance=5)

sitk.WriteImage(ii, sys.argv[3])
sitk.WriteImage(im, sys.argv[4])
