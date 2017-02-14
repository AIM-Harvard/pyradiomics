# -*- coding: utf-8 -*-
from __future__ import print_function

import collections
import glob
import os

import six


class DatasetHierarchyReader(object):
  def __init__(self, inputDatasetDirectory, filetype='.nrrd'):
    self.inputDatasetDirectory = inputDatasetDirectory
    self.filetype = filetype
    self.DatabaseHierarchyDict = collections.OrderedDict()

  def setInputDatasetDirectory(self, inputDatasetDirectory):
    self.inputDatasetDirectory = inputDatasetDirectory

  def setFiletype(self, filetype):
    self.filetype = filetype

  def ReadDatasetHierarchy(self, create=False):
    patientDirectories = glob.glob(os.path.join(self.inputDatasetDirectory, '*'))

    for patientDirectory in patientDirectories:
      self.DatabaseHierarchyDict[patientDirectory] = collections.OrderedDict()
      studyDirectories = glob.glob(os.path.join(patientDirectory, '*'))

      for studyDirectory in studyDirectories:
        self.DatabaseHierarchyDict[patientDirectory][studyDirectory] = collections.OrderedDict()

        subfolders = [dirpath for dirpath in glob.glob(os.path.join(studyDirectory, '*')) if os.path.isdir(dirpath)]

        reconstructionsDirectory, images = self.readReconstructionsDirectory(studyDirectory, subfolders, create=create)
        self.DatabaseHierarchyDict[patientDirectory][studyDirectory]["reconstructions"] = images

        resourcesDirectory, resources = self.readResourcesDirectory(studyDirectory, subfolders, create=create)
        self.DatabaseHierarchyDict[patientDirectory][studyDirectory]["resources"] = resources

        segmentationsDirectory, labels = self.readSegmentationsDirectory(studyDirectory, subfolders, create=create)
        self.DatabaseHierarchyDict[patientDirectory][studyDirectory]["segmentations"] = labels

    return self.DatabaseHierarchyDict

  def readReconstructionsDirectory(self, studyDirectory, subfolders, create=False):
    images = []
    recDirectory = "NONE"
    try:
      recDirectory = [item for item in subfolders if 'reconstructions' in os.path.basename(item).lower()][0]
      images = [item for item in glob.glob(os.path.join(recDirectory, "*")) if self.filetype in os.path.basename(item)]
    except IndexError:
      if create:
        recDirectory = os.path.join(studyDirectory, "Reconstructions")
        if not os.path.exists(recDirectory):
          os.mkdir(recDirectory)
          print("\tCreated:", recDirectory)

    return recDirectory, images

  def readSegmentationsDirectory(self, studyDirectory, subfolders, create=False):
    labels = []
    segDirectory = "NONE"
    try:
      segDirectory = [item for item in subfolders if 'segmentations' in os.path.basename(item).lower()][0]
      labels = [item for item in glob.glob(os.path.join(segDirectory, "*")) if self.filetype in os.path.basename(item)]
    except IndexError:
      if create:
        segDirectory = os.path.join(studyDirectory, "Segmentations")
        if not os.path.exists(segDirectory):
          os.mkdir(segDirectory)
          print("\tCreated:", segDirectory)

    return segDirectory, labels

  def readResourcesDirectory(self, studyDirectory, subfolders, create=False):
    resources = []
    resDirectory = "NONE"
    try:
      resDirectory = [item for item in subfolders if 'resources' in os.path.basename(item).lower()][0]
      resources = [item for item in glob.glob(os.path.join(resDirectory, "*"))]
    except IndexError:
      if create:
        resDirectory = os.path.join(studyDirectory, "Resources")
        if not os.path.exists(resDirectory):
          os.mkdir(resDirectory)
          print("\tCreated:", resDirectory)

    return resDirectory, resources

  def findImageAndLabelPair(self, imageFilepaths, maskFilepaths, keywordSettings):
    """
    Accepts a list of image filepaths, a list of mask/label filepaths, and a
    dict of keyword settings in the form:

    keywordSettings['image'] = ""
    keywordSettings['imageExclusion'] = ""
    keywordSettings['mask'] = ""
    keywordSettings['maskExclusion'] = ""

    where each field is a string of words separated by commas (case and spaces do not matter).

    The output is the image filepath and mask/label filepath pair that satisfies the keyword
    conditions.
    """

    keywordSettings = {k: [str(keyword.strip()) for keyword in v.split(',')]
                       for (k, v) in six.iteritems(keywordSettings)}

    matchedImages = []
    for imageFilepath in imageFilepaths:
      imageFilename = str(os.path.basename(imageFilepath))
      if self.testString(imageFilename, keywordSettings['image'], keywordSettings['imageExclusion']):
        matchedImages.append(imageFilepath)

    matchedMasks = []
    for maskFilepath in maskFilepaths:
      maskFilename = str(os.path.basename(maskFilepath))
      if self.testString(maskFilename, keywordSettings['mask'], keywordSettings['maskExclusion']):
        matchedMasks.append(maskFilepath)

    if len(matchedImages) < 1:
      print("ERROR: No Images Matched")
    elif len(matchedImages) > 1:
      print("ERROR: Multiple Images Matched")

    if len(matchedMasks) < 1:
      print("ERROR: No Masks Matched")
    elif len(matchedMasks) > 1:
      print("ERROR: Multiple Masks Matched")

    if (len(matchedImages) == 1) and (len(matchedMasks) == 1):
      return matchedImages[0], matchedMasks[0]
    else:
      return None, None

  def testString(self, fileName, inclusionKeywords, exclusionKeywords):
    fileName = fileName.upper()
    inclusionKeywords = [keyword.upper() for keyword in inclusionKeywords if (keyword != '')]
    exclusionKeywords = [keyword.upper() for keyword in exclusionKeywords if (keyword != '')]

    result = False
    if (len(inclusionKeywords) == 0) and (len(exclusionKeywords) > 0):
      if (not any(keyword in fileName for keyword in exclusionKeywords)):
        result = True
    elif (len(inclusionKeywords) > 0) and (len(exclusionKeywords) == 0):
      if (all(keyword in fileName for keyword in inclusionKeywords)):
        result = True
    elif (len(inclusionKeywords) > 0) and (len(exclusionKeywords) > 0):
      if (all(keyword in fileName for keyword in inclusionKeywords)) and \
        (not any(keyword in fileName for keyword in exclusionKeywords)):
        result = True
    elif (len(inclusionKeywords) == 0) and (len(exclusionKeywords) == 0):
      result = True

    return result
