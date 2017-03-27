
import ast
import logging
import math
import os

from nose_parameterized import parameterized
import numpy
import pandas
import SimpleITK as sitk

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
  return str("%s_%s" % (
    testcase_func.__name__,
    parameterized.to_safe_name("_".join(str(x) for x in param.args)),
  ))


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

    self._baseline = pandas.DataFrame()
    self._featureClasses = set()
    self.readBaselineFiles()

    self._kwargs = {}
    self._featureClassName = None

    self._testCase = None
    self._testedSet = set()

    self._results = pandas.DataFrame(index=self.getTestCases())
    self._diffs = pandas.DataFrame(index=self.getTestCases())

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

    # First set featureClass if necessary, because if settings have changed, testCase needs te be reloaded
    if self._featureClassName != className:
      self._logger.debug('Setting feature class name to %s', className)
      assert className in self.getFeatureClasses()

      self._featureClassName = className

      # Check if test settings have changed
      if self._kwargs != self.getBaselineSettings(className, testCase):
        self._kwargs = self.getBaselineSettings(className, testCase)
        self._testCase = None  # forces image to be reloaded (as settings have changed)

    # Next, set testCase if necessary
    if self._testCase != testCase:
      self._logger.info("Reading the image and mask for test case %s", testCase)
      assert testCase in self.getTestCases()
      self._testedSet.add(testCase)

      imageName = str(os.path.join(self._dataDir, testCase + '_image.nrrd'))
      maskName = str(os.path.join(self._dataDir, testCase + '_label.nrrd'))

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
      bb = imageoperations.checkMask(self._image, self._mask)
      self._image, self._mask = imageoperations.cropToTumorMask(self._image, self._mask, bb,
                                                                self._kwargs.get('label', 1))
      self._testCase = testCase

    return True

  def getBaselineSettings(self, featureClass, testCase):
    dictSeries = self._baseline.get('%s_general_info_GeneralSettings' % featureClass, None)
    if dictSeries is not None:
      return ast.literal_eval(dictSeries[testCase])
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
    return self._baseline.index

  def getFeatureClasses(self):
    """
    Return all the feature classes for which there are baseline information.
    """
    return self._featureClasses

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

      baseline = pandas.read_csv(os.path.join(self._baselineDir, baselineFile))
      assert 'Patient ID' in baseline
      baseline.set_index('Patient ID', inplace=True)
      baseline.columns = ['%s_%s' % (cls, feature) for feature in baseline.columns]
      self._baseline = self._baseline.join(baseline, how='outer')
      self._featureClasses.add(cls)

  def checkResult(self, featureName, value):
    """
    Use utility methods to get and test the results against the expected baseline value for this key.
    """

    longName = '%s_%s' % (self._featureClassName, featureName)
    if longName not in self._results:
      self._results = self._results.join(pandas.Series(name=longName))
      self._diffs = self._diffs.join(pandas.Series(name=longName))

    if value is None:
      self._diffs[longName][self._testCase] = None
      self._results[longName][self._testCase] = None
    assert (value is not None)

    if math.isnan(value):
      self._diffs[longName][self._testCase] = numpy.nan
      self._results[longName][self._testCase] = numpy.nan
    assert (not math.isnan(value))

    # save the result using the baseline class and feature names
    self._logger.debug('checkResults: featureName = %s', featureName)

    self._results[longName][self._testCase] = value

    assert longName in self._baseline
    baselineValue = float(self._baseline[longName][self._testCase])
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
    self._diffs[longName][self._testCase] = percentDiff

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
    Write out data in a csv file. Assumes data is a pandas DataFrame, with test cases as index (rows).
    """
    # Get the headers from the first testCase in _testedSet
    # If no tests were run, the length of _testedSet will be 0, and no files should be written
    if len(self._testedSet) > 0:
      data[data.index.isin(self._testedSet)].to_csv(fileName)
      self._logger.info('Wrote to file %s', fileName)
    else:
      self._logger.info('No test cases run, aborting file write to %s', fileName)
