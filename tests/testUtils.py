
import ast
import csv
import logging
import math
import os

from nose_parameterized import parameterized
import numpy
import SimpleITK as sitk
import six

from radiomics import getTestCase, imageoperations

# Get the logger. This is done outside the class, as it is needed by both the class and the custom_name_func
logger = logging.getLogger('radiomics.testing')


TEST_CASES = ('brain1', 'brain2', 'breast1', 'lung1', 'lung2')


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
    self._logger = logging.getLogger('radiomics.testing.utils')

    self._logger.debug('RadiomicsTestUtils')

    # the image and mask volumes
    self._image = None
    self._mask = None

    self._current_image = None
    self._current_mask = None
    self._bb = None
    self._imageType = None

    # set up file paths
    self._dataDir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
    self._baselineDir = os.path.join(self._dataDir, 'baseline')

    self._tests = set()
    self._test = None  # Test, specifies an image and mask and some configuration (settings)
    self._testCase = None  # Test image and mask to use in configured test
    self._testedSet = set()

    self._baseline = {}
    self.readBaselineFiles()

    self._current_config = {}
    self._featureClassName = None

    self._results = {}
    self._diffs = {}
    for test in self.getTests():
      self._results[test] = {}
      self._diffs[test] = {}

  def readBaselineFiles(self):
    """
    Reads the 'baseline' folder contained in dataDir. All files starting with 'baseline_' are read as baseline files.
    These files should therefore be named as follows: 'baseline_<className>.csv'.
    """
    baselineFiles = [fileName for fileName in os.listdir(self._baselineDir)
                     if os.path.isfile(os.path.join(self._baselineDir, fileName)) and fileName.startswith('baseline_')]
    assert len(baselineFiles) > 0
    for baselineFile in baselineFiles:
      newBaseline = PyRadiomicsBaseline.readBaselineFile(os.path.join(self._baselineDir, baselineFile))

      cls = newBaseline.cls
      self._logger.debug('Read baseline for class %s', cls)
      self._baseline[cls] = newBaseline
      self._tests |= newBaseline.tests

  def getTests(self):
    """
    Return all the tests for which there are baseline information.
    """
    return self._tests

  def getFeatureNames(self, className, test):
    """
    Gets all features for which a baseline value is available for the current class and test case. Returns a list
    containing the feature names (without image type and feature class specifiers, i.e. just the feature name).
    """
    if className not in self._baseline:
      return None  # No baseline available for specified class
    return self._baseline[className].getTestFeatures(test)

  def setFeatureClassAndTestCase(self, className, test):
    """
    Set testing suite to specified testCase and feature class. Throws an assertion error if either class or test case
    are not recognized. These have to be set here together, as the settings with which the test case has to be loaded
    are defined per feature class in the baseline (extracted from provenance information).

    Only (re)loads an image/mask if the test case has changed, or the change of feature class causes a change in test
    settings.

    If feature class and test case are unchanged, nothing is reloaded and function returns False. If either feature
    class or test case is changed, function returns True.
    """
    global TEST_CASES
    if self._featureClassName == className and self._test == test:
      return False

    self._test = test
    self._testedSet.add(self._test)

    # First set featureClass if necessary, because if settings have changed, testCase needs te be reloaded
    if self._featureClassName != className:
      self._logger.debug('Setting feature class name to %s', className)
      assert className in self._baseline.keys()  # Check if a baseline has been read for this class

      self._featureClassName = className

      # Check if test settings have changed
      if self._current_config != self._baseline[className].getTestConfig(test):
        self._current_config = self._baseline[className].getTestConfig(test)
        self._testCase = None  # forces image to be reloaded (as settings have changed)

    # Next, set testCase if necessary
    if self._testCase != self._current_config['TestCase']:
      self._testCase = self._current_config['TestCase']
      self._logger.info("Reading the image and mask for test case %s", self._testCase)
      assert self._current_config['TestCase'] in TEST_CASES

      imageName, maskName = getTestCase(self._testCase)

      assert imageName is not None
      assert maskName is not None

      self._image = sitk.ReadImage(imageName)
      self._mask = sitk.ReadImage(maskName)

      if 'ImageHash' in self._current_config:
        assert sitk.Hash(self._image) == self._current_config['ImageHash']
      if 'MaskHash' in self._current_config:
        assert sitk.Hash(self._mask) == self._current_config['MaskHash']

      settings = self._current_config.get('Settings', {})

      interpolator = settings.get('interpolator', sitk.sitkBSpline)
      resampledPixelSpacing = settings.get('resampledPixelSpacing', None)

      if interpolator is not None and resampledPixelSpacing is not None:
        self._image, self._mask = imageoperations.resampleImage(self._image,
                                                                self._mask,
                                                                resampledPixelSpacing,
                                                                interpolator,
                                                                settings.get('label', 1),
                                                                settings.get('padDistance', 5))
      self._bb, correctedMask = imageoperations.checkMask(self._image, self._mask, **settings)
      if correctedMask is not None:
        self._mask = correctedMask

      self._imageType = None

    return True

  def getImage(self, imageType):
    if self._imageType != imageType:
      self._applyFilter(imageType)
    return self._current_image

  def getMask(self, imageType):
    if self._imageType != imageType:
      self._applyFilter(imageType)
    return self._current_mask

  def _applyFilter(self, imageType):
    if imageType == 'original':
      self._current_image, self._current_mask = imageoperations.cropToTumorMask(self._image, self._mask, self._bb)
    else:
      raise NotImplementedError()

    self._imageType = imageType

  def getSettings(self):
    return self._current_config.get('Settings', {})

  def checkResult(self, featureName, value):
    """
    Use utility methods to get and test the results against the expected baseline value for this key.
    """

    longName = '_'.join(featureName)
    if value is None:
      self._diffs[self._test][longName] = None
      self._results[self._test][longName] = None
    assert (value is not None)

    if math.isnan(value):
      self._diffs[self._test][longName] = numpy.nan
      self._results[self._test][longName] = numpy.nan
    assert (not math.isnan(value))

    # save the result using the baseline class and feature names
    self._logger.debug('checkResults: featureName = %s', featureName)

    self._results[self._test][longName] = value

    baselineValue = self._baseline[self._featureClassName].getBaselineValue(self._test, longName)
    assert baselineValue is not None
    baselineValue = float(baselineValue)
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
    self._diffs[self._test][longName] = percentDiff

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
    # Get the headers from the first testCase in _testedSet
    # If no tests were run, the length of _testedSet will be 0, and no files should be written
    if len(self._testedSet) > 0:
      with open(fileName, 'w') as csvFile:
        csvFileWriter = csv.writer(csvFile, lineterminator='\n')
        testedCases = sorted(self._testedSet)
        header = sorted(data[testedCases[0]].keys())
        header = ['testCase'] + header
        csvFileWriter.writerow(header)
        for testCase in testedCases:
          thisCase = data[testCase]
          thisCase['testCase'] = testCase
          row = []
          for h in header:
            row = row + [thisCase.get(h, "N/A")]
          csvFileWriter.writerow(row)
        self._logger.info('Wrote to file %s', fileName)
    else:
      self._logger.info('No test cases run, aborting file write to %s', fileName)


class PyRadiomicsBaseline:

  def __init__(self, featureClassName):
    self.logger = logging.getLogger('radiomics.testing.baseline')

    self.cls = featureClassName

    self.configuration = {}
    self.baseline = {}
    self.tests = set()

  @classmethod
  def readBaselineFile(cls, baselineFile):
    featureClassName = os.path.basename(baselineFile)[9:-4]

    new_baseline = cls(featureClassName)
    new_baseline.logger.debug('Reading baseline for class %s', new_baseline.cls)

    with open(baselineFile, 'r' if six.PY3 else 'rb') as baselineReader:
      csvReader = csv.reader(baselineReader)
      tests = six.next(csvReader)[1:]

      for case in tests:
        new_baseline.configuration[case] = {}
        new_baseline.baseline[case] = {}

      for testRow in csvReader:
        for case_idx, case in enumerate(tests, start=1):
          if 'general_info' in testRow[0]:
            new_baseline.configuration[case][testRow[0]] = testRow[case_idx]
          else:
            new_baseline.baseline[case][testRow[0]] = testRow[case_idx]

      new_baseline.tests = set(tests)
    return new_baseline

  def getTestConfig(self, test):
    if test not in self.configuration:
      return {}  # This test is not present in the baseline for this class

    config = {
      'TestCase': self.configuration[test].get('general_info_TestCase', None),
      'Settings': ast.literal_eval(self.configuration[test].get('general_info_GeneralSettings', '{}')),
    }

    if 'general_info_ImageHash' in self.configuration[test]:
      config['ImageHash'] = self.configuration[test]['general_info_ImageHash']
    if 'general_info_MaskHash' in self.configuration[test]:
      config['MaskHash'] = self.configuration[test]['general_info_MaskHash']

    if config['TestCase'] is None:
      self.logger.error('Missing key "general_info_TestCase". Cannot configure!')
      return None

    return config

  def getTestFeatures(self, test):
    """
    Gets all features for which a baseline value is available for the current class and test case. Returns a list
    containing the feature names.
    """
    if test not in self.baseline:
      return None  # This test is not present in the baseline for this class
    return list(self.baseline[test].keys())

  def getBaselineValue(self, test, featureName):
    if test not in self.baseline:
      return None
    return self.baseline[test].get(featureName, None)

  def writeBaselineFile(self, baselineDir):
    baselineFile = os.path.join(baselineDir, 'baseline_%s.csv' % self.cls)
    testCases = list(self.baseline.keys())
    with open(baselineFile, 'wb') as baseline:
      csvWriter = csv.writer(baseline)
      header = ['featureName'] + testCases
      csvWriter.writerow(header)

      config = self.configuration[testCases[0]].keys()
      for c in config:
        row = [c]
        for testCase in testCases:
          row.append(str(self.configuration[testCase].get(c, '')))
        csvWriter.writerow(row)

      features = self.baseline[testCases[0]].keys()
      for f in features:
        row = [f]
        for testCase in testCases:
          row.append(str(self.baseline[testCase].get(f, '')))
        csvWriter.writerow(row)
