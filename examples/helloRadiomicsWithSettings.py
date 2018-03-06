#!/usr/bin/env python

from __future__ import print_function

import logging
import os

import six

import radiomics
from radiomics import featureextractor

def tqdmProgressbar():
  """
  This function will setup the progress bar exposed by the 'tqdm' package.
  Progress reporting is only used in PyRadiomics for the calculation of GLCM and GLSZM in full python mode, therefore
  enable GLCM and full-python mode to show the progress bar functionality

  N.B. This function will only work if the 'click' package is installed (not included in the PyRadiomics requirements)
  """
  global extractor
  extractor.kwargs['enableCExtensions'] = False
  # Enable the GLCM class to show the progress bar
  extractor.enableFeatureClassByName('glcm')

  radiomics.setVerbosity(logging.INFO)  # Verbosity must be at least INFO to enable progress bar

  import tqdm
  radiomics.progressReporter = tqdm.tqdm

def clickProgressbar():
  """
  This function will setup the progress bar exposed by the 'click' package.
  Progress reporting is only used in PyRadiomics for the calculation of GLCM and GLSZM in full python mode, therefore
  enable GLCM and full-python mode to show the progress bar functionality.

  Because the signature used to instantiate a click progress bar is different from what PyRadiomics expects, we need to
  write a simple wrapper class to enable use of a click progress bar. In this case we only need to change the 'desc'
  keyword argument to a 'label' keyword argument.

  N.B. This function will only work if the 'click' package is installed (not included in the PyRadiomics requirements)
  """
  global extractor

  extractor.kwargs['enableCExtensions'] = False
  # Enable the GLCM class to show the progress bar
  extractor.enableFeatureClassByName('glcm')

  radiomics.setVerbosity(logging.INFO)  # Verbosity must be at least INFO to enable progress bar

  import click

  class progressWrapper():
    def __init__(self, iterable, desc=''):
      # For a click progressbar, the description must be provided in the 'label' keyword argument.
      self.bar = click.progressbar(iterable, label=desc)

    def __iter__(self):
      return self.bar.__iter__()  # Redirect to the __iter__ function of the click progressbar

    def __enter__(self):
      return self.bar.__enter__()  # Redirect to the __enter__ function of the click progressbar

    def __exit__(self, exc_type, exc_value, tb):
      return self.bar.__exit__(exc_type, exc_value, tb)  # Redirect to the __exit__ function of the click progressbar

  radiomics.progressReporter = progressWrapper


# Get some test data

# Download the test case to temporary files and return it's location. If already downloaded, it is not downloaded again,
# but it's location is still returned.
imageName, maskName = radiomics.getTestCase('brain1')

# Get the location of the example settings file
paramsFile = os.path.abspath(os.path.join('exampleSettings', 'Params.yaml'))

if imageName is None or maskName is None:  # Something went wrong, in this case PyRadiomics will also log an error
  print('Error getting testcase!')
  exit()

# Regulate verbosity with radiomics.verbosity
# radiomics.setVerbosity(logging.INFO)

# Get the PyRadiomics logger (default log-level = INFO
logger = radiomics.logger
logger.setLevel(logging.DEBUG)  # set level to DEBUG to include debug log messages in log file

# Write out all log entries to a file
handler = logging.FileHandler(filename='testLog.txt', mode='w')
formatter = logging.Formatter("%(levelname)s:%(name)s: %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)

# Initialize feature extractor using the settings file
extractor = featureextractor.RadiomicsFeaturesExtractor(paramsFile)

# Uncomment one of these functions to show how PyRadiomics can use the 'tqdm' or 'click' package to report progress when
# running in full python mode. Assumes the respective package is installed (not included in the requirements)

# tqdmProgressbar()
# clickProgressbar()

print("Active features:")
for cls, features in six.iteritems(extractor._enabledFeatures):
  if len(features) == 0:
    features = [f for f, deprecated in six.iteritems(extractor.getFeatureNames(cls)) if not deprecated]
  for f in features:
    print(f)
    print(getattr(extractor.featureClasses[cls], 'get%sFeatureValue' % f).__doc__)

print("Calculating features")
featureVector = extractor.execute(imageName, maskName)

for featureName in featureVector.keys():
  print("Computed %s: %s" % (featureName, featureVector[featureName]))
