import SimpleITK as sitk
import os
import ast
import csv
import logging
import math
import numpy
from nose_parameterized import parameterized
from radiomics import imageoperations

# Get the logger. This is done outside the class, as it is needed by both the class and the custom_name_func
logger = logging.getLogger('testUtils')

def custom_name_func(testcase_func, param_num, param):
  """
  A custom test name function that will ensure that the tests are run such that they're batched with all tests for a
  given data set are run together, avoiding re-reading the data more than necessary. Tests are run in alphabetical
  order, so put the test case first. An alternate option is to right justify the test number (param_num) with zeroes
  so that the numerical and alphabetical orders are the same. Not providing this method when there are more than 10
  tests results in tests running in an order similar to:

  test_*.test_scenario_0_*

  test_*.test_scenario_10_*

  test_*.test_scenario_11_*

  ...

  test_*.test_scenario_19_*

  test_*.test_scenario_1_*

  test_*.test_scenario_20_*
  """
  global logger

  logger.debug('custom_name_func: function name = %s, param_num = {0:0>3}, param.args = %s'.format(param_num),
               testcase_func.__name__, param.args)
  return "%s_%s" % (
    testcase_func.__name__,
    parameterized.to_safe_name("_".join(str(x) for x in param.args)),
  )


class RadiomicsTestUtils:
  """
  This utility class reads in and stores the baseline files stored in 'data\baseline' (one per feature class)
  It provides utility methods to get the baseline feature value for a feature class and compare it to the result generated
  by the test.
  """
  def __init__(self):
    global logger

    self._logger = logger

    self._logger.debug('RadiomicsTestUtils')

    # the image and mask volumes
    self._image = None
    self._mask = None

    # set up file paths
    self._dataDir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
    self._baselineDir = os.path.join(self._dataDir, 'baseline')
    self._mappingDir = os.path.join(self._dataDir, 'mapping')

    self._baseline = {}
    self.readBaselineFiles()

    self._kwargs = {}
    self._featureClassName = None

    self._testCase = None

    self._results = {}
    self._diffs = {}
    for testCase in self.getTestCases():
      self._results[testCase] = {}
      self._diffs[testCase] = {}

  def setFeatureClassAndTestCase(self, className, testCase):
    """
    Set testing suite to specified testCase and feature class. Throws an assertion error if either class or test case
    are not recognized. These have to be set here together, as the settings with which the test case has to be loaded
    are defined per feature class in the baseline (extracted from provenance information).

    Only (re)loads an image/mask if the test case has changed, or the change of feature class causes a change in test
    settings.

    If feature class and test case are unchanged, nothing is reloaded and function returns False. If either feature
    class or test case is changed, function returns True.
    """
    if self._featureClassName == className and self._testCase == testCase:
      return False

    assert className in self.getFeatureClasses()
    assert testCase in self.getTestCases()

    # First set featureClass if necessary, because if settings have changed, testCase needs te be reloaded
    if self._featureClassName != className:
      self._logger.debug('Setting feature class name to %s', className)

      self._featureClassName = className

      # Check if test settings have changed
      if cmp(self._kwargs, self.getBaselineDict(className, testCase)) != 0:
        self._kwargs = self.getBaselineDict(className, testCase)
        self._testCase = None  # forces image to be reloaded (as settings have changed)

    # Next, set testCase if necessary
    if self._testCase != testCase:
      imageName = os.path.join(self._dataDir, testCase + '_image.nrrd')
      maskName = os.path.join(self._dataDir, testCase + '_label.nrrd')

      self._logger.info("Reading the image and mask for test case %s", testCase)
      self._image = sitk.ReadImage(imageName)
      self._mask = sitk.ReadImage(maskName)

      interpolator = self._kwargs.get('interpolator', sitk.sitkBSpline)
      resampledPixelSpacing = self._kwargs.get('resampledPixelSpacing', None)

      if interpolator is not None and resampledPixelSpacing is not None:
        self._image, self._mask = imageoperations.resampleImage(self._image,
                                                                self._mask,
                                                                resampledPixelSpacing,
                                                                interpolator,
                                                                self._kwargs.get('label', 1),
                                                                self._kwargs.get('padDistance', 5))
      self._image, self._mask, bb = imageoperations.cropToTumorMask(self._image, self._mask,
                                                                    self._kwargs.get('label', 1))
      self._testCase = testCase

    return True

  def getBaselineDict(self, featureClass, testCase):
    dictStr = self._baseline[featureClass][testCase].get('general_info_GeneralSettings', None)
    if dictStr is not None:
      return ast.literal_eval(str(dictStr).replace(';', ','))
    return {}

  def getTestCase(self):
    return self._testCase

  def getImage(self):
    return self._image

  def getMask(self):
    return self._mask

  def getKwargs(self):
    return self._kwargs

  def getTestCases(self):
    """
    Return all the test cases for which there are baseline information.
    """
    return self._baseline[self._baseline.keys()[0]].keys()

  def getFeatureClasses(self):
    """
    Return all the feature classes for which there are baseline information.
    """
    return self._baseline.keys()

  def readBaselineFiles(self):
    """
    Reads the 'baseline' folder contained in dataDir. All files starting with 'baseline_' are read as baseline files.
    These files should therefore be named as follows: 'baseline_<className>.csv'.
    """
    baselineFiles = [fileName for fileName in os.listdir(self._baselineDir)
                     if os.path.isfile(os.path.join(self._baselineDir, fileName)) and fileName.startswith('baseline_')]
    assert len(baselineFiles) > 0
    for baselineFile in baselineFiles:
      cls = baselineFile[9:-4]
      self._logger.debug('Reading baseline for class %s', cls)
      self._baseline[cls] = {}
      with open(os.path.join(self._baselineDir, baselineFile), 'rb') as baselineReader:
        csvReader = csv.reader(baselineReader)
        headers = csvReader.next()
        for testRow in csvReader:
          self._baseline[cls][testRow[0]] = {}
          for val_idx, val in enumerate(testRow[1:], start=1):
            self._baseline[cls][testRow[0]][headers[val_idx]] = val

  def checkResult(self, featureName, value):
    """
    Use utility methods to get and test the results against the expected baseline value for this key.
    """

    if value is None:
      self._diffs[self._testCase][featureName] = None
      self._results[self._testCase][featureName] = None
    assert (value is not None)

    if math.isnan(value):
      self._diffs[self._testCase][featureName] = numpy.nan
      self._results[self._testCase][featureName] = numpy.nan
    assert (not math.isnan(value))

    # save the result using the baseline class and feature names
    self._logger.debug('checkResults: featureName = %s', featureName)

    self._results[self._testCase][featureName] = value

    assert featureName in self._baseline[self._featureClassName][self._testCase]
    baselineValue = float(self._baseline[self._featureClassName][self._testCase][featureName])
    self._logger.debug('checkResults: for featureName %s, got baseline value = %f', featureName, baselineValue)

    if baselineValue == 0.0:
      # avoid divide by zero, the difference is either 0% if the value is also zero, or 100%
      if value - baselineValue == 0.0:
        percentDiff = 0.0
      else:
        percentDiff = 1.0
    else:
      percentDiff = abs(1.0 - (value / baselineValue))

    # save the difference
    self._diffs[self._testCase][featureName] = percentDiff

    # check for a less than three percent difference
    if (percentDiff >= 0.03):
      self._logger.error('checkResult %s, baseline value = %f, calculated = %f, diff = %f%%', featureName,
                         float(baselineValue), value, percentDiff * 100)
    assert (percentDiff < 0.03)

  def getResults(self):
    return self._results

  def getDiffs(self):
    return self._diffs

  def getDataDir(self):
    return self._dataDir

  def writeCSV(self, data, fileName):
    """
    Write out data in a csv file.
    Assumes a data structure with:

    {'id1' : {'f1':n1, 'f2':n2}, 'id2' : {'f1':n3, 'f2':n4}}
    """
    csvFile = open(fileName, 'wb')
    csvFileWriter = csv.writer(csvFile)
    # get the headers from the first row
    header = sorted(data[data.keys()[0]].keys())
    header = ['testCase'] + header
    csvFileWriter.writerow(header)
    for testCase in sorted(data.keys()):
      thisCase = data[testCase]
      thisCase['testCase'] = testCase
      row = []
      for h in header:
        row = row + [thisCase[h]]
      csvFileWriter.writerow(row)
    csvFile.close()
    self._logger.info('Wrote to file %s', fileName)
