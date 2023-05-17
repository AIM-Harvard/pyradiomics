import logging
import os

import numpy
import SimpleITK as sitk

from radiomics import getTestCase, imageoperations

logger = logging.getLogger('radiomics.testing')
testCases = ('test_wavelet_64x64x64', 'test_wavelet_37x37x37')
baselineFile = os.path.join(os.path.dirname(__file__), '../data/baseline/wavelet.npy')


def pytest_generate_tests(metafunc):
  metafunc.parametrize(["testCase", "image", "mask", "baseline"], metafunc.cls.generate_scenarios())


class TestWavelet:

  @staticmethod
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

  def test_scenario(self, testCase, image, mask, baseline):
    global logger, testUtils, featureClasses

    logger.debug('test_scenario: testCase = %s,', testCase)

    im_arr = sitk.GetArrayFromImage(image)
    ma_arr = sitk.GetArrayFromImage(mask) == 1  # Conver to boolean array, label = 1

    voxelArray = im_arr[ma_arr]  # 1D array of all voxels inside mask

    assert voxelArray.shape[0] == baseline.shape[0]

    diff = voxelArray - baseline

    assert numpy.max(numpy.abs(diff)) < 1e-3
