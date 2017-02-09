from __future__ import print_function, unicode_literals, division, absolute_import
import os
import csv
from DatasetHierarchyReader import DatasetHierarchyReader


def main():
  DATA_ROOT_PATH = r''
  inputDirectory = DATA_ROOT_PATH + r''
  outputFile = DATA_ROOT_PATH + r'/FileList.csv'
  filetype = '.nrrd'

  keywordSettings = {}
  keywordSettings['image'] = ''
  keywordSettings['imageExclusion'] = 'label'
  keywordSettings['mask'] = 'label'
  keywordSettings['maskExclusion'] = ''

  print("Scanning files...")

  datasetReader = DatasetHierarchyReader(inputDirectory, filetype=filetype)
  datasetHierarchyDict = datasetReader.ReadDatasetHierarchy()

  print("Found %s patients, writing csv" % (str(len(datasetHierarchyDict.keys()))))

  try:
    with open(outputFile, 'wb') as outFile:
      cw = csv.writer(outFile, lineterminator='\n')

      for patientIndex, patientDirectory in enumerate(datasetHierarchyDict):
        patientID = os.path.basename(patientDirectory)

        for studyDirectory in datasetHierarchyDict[patientDirectory]:
          studyDate = os.path.basename(studyDirectory)

          imageFilepaths = datasetHierarchyDict[patientDirectory][studyDirectory]["reconstructions"]
          maskFilepaths = datasetHierarchyDict[patientDirectory][studyDirectory]["segmentations"]

          imageFilepath, maskFilepath = datasetReader.findImageAndLabelPair(imageFilepaths, maskFilepaths,
                                                                            keywordSettings)

          if (imageFilepath is not None) and (maskFilepath is not None):
            # ReaderName is not extracted using DatasetHierarchyReader, set it to 'N/A'
            cw.writerow([patientID, studyDate, 'N/A', imageFilepath, maskFilepath])

  except Exception as e:
    print(e)


if __name__ == '__main__':
  main()
