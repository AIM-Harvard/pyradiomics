#!/usr/bin/env python
from __future__ import annotations

import logging
import os
import sys

import radiomics
from radiomics import featureextractor, getFeatureClasses

# Get some test data

# Download the test case to temporary files and return it's location. If already downloaded, it is not downloaded again,
# but it's location is still returned.
imageName, maskName = radiomics.getTestCase("brain1")

# Get the location of the example settings file
paramsFile = os.path.abspath(os.path.join("exampleSettings", "Params.yaml"))

if (
    imageName is None or maskName is None
):  # Something went wrong, in this case PyRadiomics will also log an error
    print("Error getting testcase!")
    sys.exit()

# Regulate verbosity with radiomics.verbosity
# radiomics.setVerbosity(logging.INFO)

# Get the PyRadiomics logger (default log-level = INFO
logger = radiomics.logger
logger.setLevel(
    logging.DEBUG
)  # set level to DEBUG to include debug log messages in log file

# Write out all log entries to a file
handler = logging.FileHandler(filename="testLog.txt", mode="w")
formatter = logging.Formatter("%(levelname)s:%(name)s: %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)

# Initialize feature extractor using the settings file
extractor = featureextractor.RadiomicsFeatureExtractor(paramsFile)
featureClasses = getFeatureClasses()

print("Active features:")
for cls, orig_features in extractor.enabledFeatures.items():
    features = orig_features
    if features is None or len(features) == 0:
        features = [
            f
            for f, deprecated in featureClasses[cls].getFeatureNames().items()
            if not deprecated
        ]
    for f in features:
        print(f)
        print(getattr(featureClasses[cls], f"get{f}FeatureValue").__doc__)

print("Calculating features")
featureVector = extractor.execute(imageName, maskName)

for featureName in featureVector:
    print(f"Computed {featureName}: {featureVector[featureName]}")
