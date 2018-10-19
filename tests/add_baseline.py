import logging
import os
import sys

import numpy
import SimpleITK as sitk
import six

from radiomics import generalinfo, getFeatureClasses, getTestCase, imageoperations
from testUtils import PyRadiomicsBaseline, RadiomicsTestUtils


class AddBaseline:

  def __init__(self):
    self.logger = logging.getLogger('radiomics.addBaseline')

    self.testUtils = RadiomicsTestUtils()

    self.testCases = sorted(self.testUtils.getTests())
    self.featureClasses = getFeatureClasses()

    dataDir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
    self.baselineDir = os.path.join(dataDir, "baseline")

  def generate_scenarios(self):
    for className, featureClass in six.iteritems(self.featureClasses):
      if not os.path.exists(os.path.join(self.baselineDir, 'baseline_%s.csv' % className)):
        self.logger.debug('generate_scenarios: featureClass = %s', className)
        for test in self.testCases:
          yield test, className

  def process_testcase(self, test, featureClassName):
    self.logger.debug('processing testCase = %s, featureClassName = %s', test, featureClassName)

    self.testUtils.setFeatureClassAndTestCase(featureClassName, test)

    testImage = self.testUtils.getImage('original')
    testMask = self.testUtils.getMask('original')

    featureClass = self.featureClasses[featureClassName](testImage, testMask, **self.testUtils.getSettings())

    featureClass.enableAllFeatures()
    featureClass.execute()

    if "_calculateMatrix" in dir(featureClass):
      cMat = getattr(featureClass, 'P_%s' % featureClassName)
      if cMat is not None:
        numpy.save(os.path.join(self.baselineDir, '%s_%s.npy' % (test, featureClassName)), cMat)

    imageTypeName = 'original'

    # Update versions to reflect which configuration generated the baseline
    generalInfo = generalinfo.GeneralInfo()

    versions = generalInfo.getGeneralInfo()
    self.new_baselines[featureClassName].configuration[test].update(versions)

    self.new_baselines[featureClassName].baseline[test] = {'%s_%s_%s' % (imageTypeName, featureClassName, key): val
                                                           for key, val in six.iteritems(featureClass.featureValues)}

  def run(self, featureClass=None):
    current_baseline = self.testUtils._baseline
    config = current_baseline[current_baseline.keys()[0]].configuration
    self.new_baselines = {}
    if featureClass is None:
      for test, newClass in self.generate_scenarios():
        if newClass not in self.new_baselines:
          self.logger.info('Adding class %s to the baseline', newClass)
          self.new_baselines[newClass] = PyRadiomicsBaseline(newClass)
          self.new_baselines[newClass].configuration = config
          # add the new baseline to test utils so it's config can be used during processing
          self.testUtils._baseline[newClass] = self.new_baselines[newClass]

        self.process_testcase(test, newClass)

      for newClass in self.new_baselines:
        self.new_baselines[newClass].writeBaselineFile(self.baselineDir)
    elif featureClass in self.featureClasses:
      if featureClass in current_baseline:
        # Re-create the baseline for specified class
        self.new_baselines[featureClass] = current_baseline[featureClass]
      else:
        # Feature class not yet present in the baseline, generate a new one
        self.new_baselines[featureClass] = PyRadiomicsBaseline(featureClass)
        self.new_baselines[featureClass].configuration = config

      for test in self.testCases:
        self.process_testcase(test, featureClass)

      self.new_baselines[featureClass].writeBaselineFile(self.baselineDir)
    else:
      self.logger.error('Feature Class %s not recognized, cannot create baseline!', featureClass)


def generateWaveletBaseline():
  logger = logging.getLogger('radiomics.addBaseline')
  baselineFile = '../data/baseline/wavelet.npy'
  testCase = 'test_wavelet_64x64x64'

  if os.path.isfile(baselineFile):  # Check if the baseline does not yet exist
    logger.info('Baseline already exists, cancelling generation')
    return

  baselineDict = {}
  wavelets = ('HHH', 'HHL', 'HLH', 'HLL', 'LHH', 'LHL', 'LLH', 'LLL')

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

    im_arr = sitk.GetArrayFromImage(image)
    ma_arr = sitk.GetArrayFromImage(mask) == 1  # Convert to boolean array, label = 1
    voxelArray = im_arr[ma_arr]  # 1D array of all voxels inside mask

    baselineDict[level] = voxelArray

  baseline = [baselineDict[wl] for wl in wavelets]
  baseline = numpy.array(baseline)
  numpy.save(baselineFile, baseline)


if __name__ == '__main__':
  generateWaveletBaseline()
  add_baseline = AddBaseline()
  if len(sys.argv) == 2:
    add_baseline.run(sys.argv[1])
  else:
    add_baseline.run()
