from __future__ import print_function

import csv
import os

from DatasetHierarchyReader import DatasetHierarchyReader


def main():
  DATA_ROOT_PATH = r'R:\TEMP'
  inputDirectory = DATA_ROOT_PATH + r'\TEST'
  outputFile = DATA_ROOT_PATH + r'/FileListNEW.csv'
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
      cw.writerow(['Patient', 'StudyDate', 'Image', 'Mask'])

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
            cw.writerow([patientID, studyDate, imageFilepath, maskFilepath])

  except Exception as exc:
    print(exc)


if __name__ == '__main__':
  main()
