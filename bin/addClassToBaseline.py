import os
import csv
import radiomics
import collections
from radiomics import featureextractor


def main():
  dataDir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
  # matlabFeaturesFile = os.path.join(dataDir, "MatlabBaselineFeatures.csv")
  addBaselineFeaturesFile = os.path.join(dataDir, "AdditionalBaselineFeatures.csv")
  outputDir = os.path.join(dataDir, 'mapping')

  classMap = {}  # key is matlab class, value python class
  addBaseline = {}
  testCases = []

  with open(outputDir + os.path.sep + r'baseline2pyradiomics_featureClasses.txt') as classFile:
    csvReader = csv.reader(classFile, delimiter=':')
    for row in csvReader:
      classMap[row[0]] = row[1]
    classFile.close()

  addBaseline = readBaselineFeatures(addBaselineFeaturesFile)

  testCases = sorted(addBaseline.keys())

  kwargs = {'verbose': False,
            'binWidth': 25,
            'interpolator': None,
            'resampledPixelSpacing': None,
            'padDistance': 5,
            'voxelArrayShift': 2000,
            'symmetricalGLCM': False,
            'weightingNorm': None,
            'gldm_a': 0}

  extractor = featureextractor.RadiomicsFeaturesExtractor(**kwargs)
  featureClasses = extractor.getFeatureClassNames()
  newClasses = [cls for cls in featureClasses if cls not in classMap.values()]

  print "Adding new classes: ", newClasses

  newBaseline = {}

  # When C matrices is merged, force baseline to be build using the Python implementation of the functions.
  # radiomics.pythonMatrixCalculation(True)

  extractor.disableAllFeatures()
  extractor.disableAllInputImages()

  # Calculate new baseline on original image for all new classes
  extractor.enableInputImageByName('original')
  for cls in newClasses:
    extractor.enableFeatureClassByName(cls)

  print "Computing new baseline"
  for testCase in testCases:
    image = os.path.join(dataDir, testCase + '_image.nrrd')
    mask = os.path.join(dataDir, testCase + '_label.nrrd')
    newBaseline[testCase] = extractor.execute(image, mask)

  print "Writing new mapping"
  for cls in newClasses:
    baselineIndicesFile = os.path.join(outputDir, r'baseline_%s.txt' % (cls))
    pyradiomicsIndicesFile = os.path.join(outputDir, r'pyradiomics_%s.txt' % (cls))
    baseline2pyradiomicsIndicesFile = os.path.join(outputDir, r'baseline2pyradiomics_%s.txt' % (cls))
    featureNames = extractor.getFeaturesNames(cls)

    with open(baselineIndicesFile, 'w') as baselineIndices:
      with open(pyradiomicsIndicesFile, 'w') as pyradiomicsIndices:
        with open(baseline2pyradiomicsIndicesFile, 'w') as baseline2pyradiomicsIndices:
          for f_idx, fname in enumerate(featureNames):
            baselineIndices.write('%i:%s\n' % (f_idx, fname))
            pyradiomicsIndices.write('%i:%s\n' % (f_idx, fname))
            baseline2pyradiomicsIndices.write('%i:%i\n' % (f_idx, f_idx))
            for testCase in testCases:
              addBaseline[testCase]['%s_%s' % (cls, fname)] = newBaseline[testCase]['original_%s_%s' % (cls, fname)]
          baselineIndices.close()
          pyradiomicsIndices.close()
          baseline2pyradiomicsIndices.close()

  print "Writing new baseline"
  with open(addBaselineFeaturesFile, 'wb') as newBaselineFile:
    csvWriter = csv.writer(newBaselineFile)
    header = addBaseline[addBaseline.keys()[0]].keys()
    csvWriter.writerow(header)
    for testCase in testCases:
      row = []
      for h in header:
        row += [addBaseline[testCase][h]]
      csvWriter.writerow(row)
    newBaselineFile.close()

  print "Writing new class names"
  with open(outputDir + os.path.sep + r'baseline2pyradiomics_featureClasses.txt', 'a') as classFile:
    csvWriter = csv.writer(classFile, delimiter=':', lineterminator='\n')
    for cls in newClasses:
      csvWriter.writerow([cls, cls])
    classFile.close()


def readBaselineFeatures(baselineFeaturesFile):
  if (not os.path.exists(baselineFeaturesFile)):
    print 'Baseline features file not found %s:' % (baselineFeaturesFile)
    return {}
  baselineFeatures = collections.OrderedDict()
  csvFile = open(baselineFeaturesFile, 'rb')
  csvFileReader = csv.reader(csvFile)
  # get the column headers
  headerRow = csvFileReader.next()
  # iterate over the test cases in the file
  for testRow in csvFileReader:
    testCase = testRow[0]
    baselineFeatures[testCase] = collections.OrderedDict()
    columnIndex = 0
    for val in testRow:
      baselineFeatures[testCase][headerRow[columnIndex]] = val
      columnIndex += 1
  csvFile.close()
  return baselineFeatures


if __name__ == '__main__':
  main()
