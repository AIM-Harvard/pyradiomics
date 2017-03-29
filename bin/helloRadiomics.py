
from __future__ import print_function

import logging
import os

import SimpleITK as sitk
import six

import radiomics
from radiomics import featureextractor


testCase = 'brain1'
dataDir = os.path.join(os.path.abspath(""), "..", "data")
imageName = os.path.join(dataDir, testCase + '_image.nrrd')
maskName = os.path.join(dataDir, testCase + '_label.nrrd')

if not os.path.exists(imageName):
  print('Error: problem finding input image', imageName)
  exit()
if not os.path.exists(maskName):
  print('Error: problem finding input labelmap', maskName)
  exit()

# Uncomment line below to output debug and info log messaging to stderr
# radiomics.debug()  # Switch on radiomics logging to stderr output from level=DEBUG (default level=WARNING)

# Alternative: regulate verbosity with radiomics.verbosity
# radiomics.setVerbosity(logging.INFO)

# Get the PyRadiomics logger (default log-level = INFO
logger = radiomics.logger
logger.setLevel(logging.DEBUG)  # set level to DEBUG to include debug log messages in log file

# Write out all log entries to a file
handler = logging.FileHandler(filename='testLog.txt', mode='w')
formatter = logging.Formatter("%(levelname)s:%(name)s: %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)

# Define settings for signature calculation
# These are currently set equal to the respective default values
kwargs = {}
kwargs['binWidth'] = 25
kwargs['resampledPixelSpacing'] = None  # [3,3,3] is an example for defining resampling (voxels with size 3x3x3mm)
kwargs['interpolator'] = sitk.sitkBSpline

# Initialize wrapperClass to generate signature
extractor = featureextractor.RadiomicsFeaturesExtractor(**kwargs)

# By default, only original is enabled. Optionally enable some filters:
# extractor.enableInputImages(Original={}, LoG={}, Wavelet={})

# Disable all classes except firstorder
extractor.disableAllFeatures()

# Enable all features in firstorder
# extractor.enableFeatureClassByName('firstorder')

# Only enable mean and skewness in firstorder
extractor.enableFeaturesByName(firstorder=['Mean', 'Skewness'])

# ############################ Progress bar section #################################################
# Uncomment this section to enable full python mode using tqdm package to print a progress bar during extraction of
# glcm and glszm. Assumes the "tqdm" package is installed.

# Comment out the following line to enable this section
"""
extractor.kwargs['enableCExtensions'] = False

import tqdm
radiomics.setVerbosity(logging.INFO)
radiomics.initProgress = tqdm.trange
# tqdm._tqdm.tqdm is the class of which an instance is returned by tqdm.trange. The update and close functions are
# called on that object
radiomics.reportProgress = tqdm._tqdm.tqdm.update
radiomics.closeProgress = tqdm._tqdm.tqdm.close

# Enable the GLCM class to show the progress bar
extractor.enableFeatureClassByName('glcm')
"""
# Comment out the preceding line to enable this section
# ############################ End of Progress bar section ###########################################

print("Active features:")
for cls, features in six.iteritems(extractor.enabledFeatures):
  if len(features) == 0:
    features = extractor.getFeatureNames(cls)
  for f in features:
    print(f)
    print(getattr(extractor.featureClasses[cls], 'get%sFeatureValue' % f).__doc__)

print("Calculating features")
featureVector = extractor.execute(imageName, maskName)

for featureName in featureVector.keys():
  print("Computed %s: %s" % (featureName, featureVector[featureName]))
