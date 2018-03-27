#!/usr/bin/env python

from __future__ import print_function

import logging

import SimpleITK as sitk

import radiomics
from radiomics import featureextractor

# Get some test data

# Download the test case to temporary files and return it's location. If already downloaded, it is not downloaded again,
# but it's location is still returned.
imageName, maskName = radiomics.getTestCase('brain1')

if imageName is None or maskName is None:  # Something went wrong, in this case PyRadiomics will also log an error
  print('Error getting testcase!')
  exit()

# Regulate verbosity with radiomics.verbosity (default verbosity level = WARNING)
# radiomics.setVerbosity(logging.INFO)

# Get the PyRadiomics logger (default log-level = INFO)
logger = radiomics.logger
logger.setLevel(logging.DEBUG)  # set level to DEBUG to include debug log messages in log file

# Set up the handler to write out all log entries to a file
handler = logging.FileHandler(filename='testLog.txt', mode='w')
formatter = logging.Formatter("%(levelname)s:%(name)s: %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)

# Define settings for signature calculation
# These are currently set equal to the respective default values
settings = {}
settings['binWidth'] = 25
settings['resampledPixelSpacing'] = None  # [3,3,3] is an example for defining resampling (voxels with size 3x3x3mm)
settings['interpolator'] = sitk.sitkBSpline

# Initialize feature extractor
extractor = featureextractor.RadiomicsFeaturesExtractor(**settings)

# By default, only original is enabled. Optionally enable some image types:
# extractor.enableImageTypes(Original={}, LoG={}, Wavelet={})

# Disable all classes except firstorder
extractor.disableAllFeatures()

# Enable all features in firstorder
# extractor.enableFeatureClassByName('firstorder')

# Only enable mean and skewness in firstorder
extractor.enableFeaturesByName(firstorder=['Mean', 'Skewness'])

print("Calculating features")
featureVector = extractor.execute(imageName, maskName)

for featureName in featureVector.keys():
  print("Computed %s: %s" % (featureName, featureVector[featureName]))
