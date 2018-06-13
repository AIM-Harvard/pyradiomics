# to run this test, from directory above:
# setenv PYTHONPATH /path/to/pyradiomics/radiomics
# nosetests --nocapture -v tests/test_features.py

import logging
import os

from nose_parameterized import parameterized
import numpy
import SimpleITK as sitk

from radiomics import getTestCase, imageoperations

logger = logging.getLogger('radiomics.testing')
testCases = ('test_wavelet_64x64x64', 'test_wavelet_37x37x37')
baselineFile = '../data/baseline/wavelet.npy'


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
    parameterized.to_safe_name(param.args[0]),
  ))


class TestWavelet:

  def generate_scenarios():
    global logger, testCases, baselineFile

    assert os.path.isfile(baselineFile)

    baseline = numpy.load(baselineFile)
    wavelets = ('HHH', 'HHL', 'HLH', 'HLL', 'LHH', 'LHL', 'LLH', 'LLL')
    baselineDict = {wavelets[idx]: b for idx, b in enumerate(baseline)}

    for testCase in testCases:
      im_path, ma_path = getTestCase(testCase)
      assert im_path is not None
      assert ma_path is not None

      logger.debug('Loading image and mask for testCase %s', testCase)
      image = sitk.ReadImage(im_path)
      mask = sitk.ReadImage(ma_path)

      logger.debug('Checking image and mask for testCase %s', testCase)
      bb, correctedMask = imageoperations.checkMask(image, mask)
      assert bb is not None  # Check mask should pass normally

      waveletGenerator = imageoperations.getWaveletImage(image, mask)
      for wavelet_image, wavelet_name, args in waveletGenerator:
        level = wavelet_name.split('-')[1]
        yield '_'.join((testCase, wavelet_name)), image, mask, baselineDict[level]

      logger.debug('Applying preCropping with padDistance 6 to testCase %s', testCase)
      image, mask = imageoperations.cropToTumorMask(image, mask, bb, padDistance=6)

      waveletGenerator = imageoperations.getWaveletImage(image, mask)
      for wavelet_image, wavelet_name, args in waveletGenerator:
        level = wavelet_name.split('-')[1]
        yield '_'.join((testCase, 'preCropped', wavelet_name)), image, mask, baselineDict[level]

  @parameterized.expand(generate_scenarios(), testcase_func_name=custom_name_func)
  def test_scenario(self, test, image, mask, baseline):
    global logger, testUtils, featureClasses

    logger.debug('test_scenario: testCase = %s,', test)

    im_arr = sitk.GetArrayFromImage(image)
    ma_arr = sitk.GetArrayFromImage(mask) == 1  # Conver to boolean array, label = 1

    voxelArray = im_arr[ma_arr]  # 1D array of all voxels inside mask

    assert voxelArray.shape[0] == baseline.shape[0]

    diff = voxelArray - baseline

    assert numpy.max(numpy.abs(diff)) < 1e-3
