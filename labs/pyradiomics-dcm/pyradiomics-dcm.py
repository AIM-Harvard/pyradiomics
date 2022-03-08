import argparse
import csv
from decimal import Decimal
import distutils.spawn
import glob
import json
import logging
import os
import shutil
from subprocess import call
import sys
import tempfile

import numpy
import pandas
import pydicom

from radiomics import featureextractor

scriptlogger = logging.getLogger('radiomics.dicom')
scriptlogger.setLevel(logging.DEBUG)

def dcmImageToNRRD(inputDICOMImageDir, tempDir):
  scanNRRDFile = os.path.join(tempDir, "image.nrrd")
  if not os.path.isfile(scanNRRDFile):
    call(['plastimatch', 'convert', '--input',
          inputDICOMImageDir, '--output-img', scanNRRDFile])
  return scanNRRDFile


def dcmImageToNIfTI(inputDICOMImageDir, tempDir):
  destScanNIfTIFile = os.path.join(tempDir, "volume.nii")
  scanNIfTIFile = os.path.join(inputDICOMImageDir, "volume.nii")
  scanJSONFile = os.path.join(inputDICOMImageDir, "volume.json")
  # will save to volume.nii
  if not os.path.isfile(destScanNIfTIFile):
    cmd = ['dcm2niix', "-m", "y", "-f", "volume", inputDICOMImageDir]
    call(cmd)
    shutil.move(scanNIfTIFile, destScanNIfTIFile)
    if os.path.isfile(scanJSONFile):
      os.remove(scanJSONFile)
  return destScanNIfTIFile


# individual segments will be extracted to the destination directory into NRRD
# files, with the names assigned consecutive numbers starting from 1


def dcmSEGToNRRDs(inputSEG, tempDir):
  segmentsDir = os.path.join(tempDir, 'Segments')
  if not os.path.isdir(segmentsDir):
    os.mkdir(segmentsDir)
  call(['segimage2itkimage', '--inputDICOM',
        inputSEG, '--outputDirectory', segmentsDir])
  return glob.glob(os.path.join(segmentsDir, "*nrrd"))


def writeSR(inputSEG, inputJSON, inputDICOMImageDir, outputSR):
  cmd = [
    'tid1500writer',
    '--inputImageLibraryDirectory',
    inputDICOMImageDir,
    '--inputCompositeContextDirectory',
    os.path.split(inputSEG)[0],
    '--inputMetadata',
    inputJSON,
    '--outputDICOM',
    outputSR]
  scriptlogger.debug("Writing SR with: " + str(cmd))
  call(cmd)


def getCTSeriesUID(imageDICOMDir):
  ctFile = os.listdir(imageDICOMDir)[0]
  dcm = pydicom.read_file(os.path.join(imageDICOMDir, ctFile))
  return dcm.SeriesInstanceUID


class DICOMMetadataAccessor:
  def __init__(self, dcmFileName):
    self.dcm = pydicom.read_file(dcmFileName)

  def getInstanceUID(self):
    return self.dcm.SOPInstanceUID

  def getSeriesDescription(self):
    return self.dcm.SeriesDescription

  def getSeriesInstanceUID(self):
    return self.dcm.SeriesInstanceUID


class SEGMetadataAccessor(DICOMMetadataAccessor):
  def __init__(self, segFileName):
    DICOMMetadataAccessor.__init__(self, segFileName)

    if self.dcm.SOPClassUID != '1.2.840.10008.5.1.4.1.1.66.4':
      raise ValueError(
        "SEGMetadataAccessor: DICOM object is not Segmentation!")

  def getSegmentSegmentationTypeCode(self, segmentNumber):
    try:
      return self.dcm.SegmentSequence[segmentNumber].SegmentedPropertyTypeCodeSequence[0]
    except BaseException:
      return None

  def getTrackingIdentifier(self, segmentNumber):
    try:
      return self.dcm.SegmentSequence[segmentNumber].TrackingIdentifier
    except BaseException:
      return None

  def getTrackingUniqueIdentifier(self, segmentNumber):
    try:
      return self.dcm.SegmentSequence[segmentNumber].TrackingUID
    except BaseException:
      return None

  def getSegmentDescription(self, segmentNumber):
    try:
      return self.dcm.SegmentSequence[segmentNumber].SegmentDescription
    except BaseException:
      return None

  def getSegmentAnatomicLocationCode(self, segmentNumber):
    try:
      return self.dcm.SegmentSequence[segmentNumber].AnatomicRegionSequence[0]
    except BaseException:
      return None


class CodedValue:
  def __init__(self, value, scheme, meaning):
    self.codeValue = value
    self.codingSchemeDesignator = scheme
    self.codeMeaning = meaning

  def getDict(self):
    return {"CodeValue": self.codeValue, "CodeMeaning": self.codeMeaning,
            "CodingSchemeDesignator": self.codingSchemeDesignator}


class TID1500Metadata:

  def __init__(
    self,
    featuresDictFile,
    seriesDescription="Radiomics features"):

    self.featuresDict = self.readDictionary(featuresDictFile)

    self.m = {}
    self.m["@schema"] = "https://raw.githubusercontent.com/qiicr/dcmqi/master/doc/schemas/sr-tid1500-schema.json#"
    self.m["SeriesDescription"] = seriesDescription

    self.m["Measurements"] = []
    self.measurementGroupCount = 0

  def addMeasurementGroup(self):
    self.measurementGroupCount = self.measurementGroupCount + 1
    measurementsGroup = {}
    measurementsGroup["measurementItems"] = []
    measurementsGroup["ReferencedSegment"] = self.measurementGroupCount
    self.m["Measurements"].append(measurementsGroup)

  @staticmethod
  def readDictionary(featuresDictFile):
    return pandas.read_csv(featuresDictFile, sep='\t', low_memory=False)

  @staticmethod
  def makeHash(text, length=6):
    from base64 import b64encode
    from hashlib import sha1
    return b64encode(sha1(str.encode(text)).digest()).decode('ascii')[:length]

  def makePrivateCode(self, text):
    return CodedValue(self.makeHash(text), "99PYRADIOMICS", text).getDict()

  # returns None if prefix is not recognized, otherwise returns a tuple of
  # (measurementModifiers, derivationParameters)
  def prefix2codes(self, prefix):

    modifiers = []
    derivationParameters = []

    import re
    imageTransformationConcept = self.makePrivateCode(
      "Image transformation")

    if re.match("original", prefix):
      pass

    elif re.match("square", prefix):
      modifiers.append({"modifier": imageTransformationConcept,
                        "modifierValue": self.makePrivateCode("Square transformation")})

    elif re.match("squareroot", prefix):
      modifiers.append({"modifier": imageTransformationConcept,
                        "modifierValue": self.makePrivateCode("Square root transformation")})

    elif re.match("logarithm", prefix):
      modifiers.append({"modifier": imageTransformationConcept,
                        "modifierValue": self.makePrivateCode("Logarithm transformation")})

    elif re.match("gradient", prefix):
      modifiers.append({"modifier": imageTransformationConcept,
                        "modifierValue": self.makePrivateCode("Gradient magnitude transformation")})

    elif re.match("exponential", prefix):
      modifiers.append({"modifier": imageTransformationConcept,
                        "modifierValue": self.makePrivateCode("Exponent transformation")})

    elif re.match("exponential", prefix):
      modifiers.append({"modifier": imageTransformationConcept,
                        "modifierValue": self.makePrivateCode("Exponent transformation")})

    # parameterized processing operations
    elif re.match(r"wavelet-([HL]{2,3})", prefix):
      match = re.match(r"wavelet-([HL]{2,3})", prefix)
      modifiers.append({"modifier": imageTransformationConcept,
                        "modifierValue": self.makePrivateCode("Wavelet transformation")})
      modifiers.append({"modifier": self.makePrivateCode("Wavelet sub-band"),
                        "modifierValue": self.makePrivateCode(match.group(1))})

    elif re.match(r"log-sigma-([\d]+)-([\d]+)-([a-z]+)", prefix):
      match = re.match(r"log-sigma-([\d]+)-([\d]+)-([a-z]+)", prefix)

      units = match.group(3)
      if units == "mm":
        unitsCode = CodedValue("mm", "UCUM", "millimeters").getDict()
      elif units == "cm":
        unitsCode = CodedValue("mm", "UCUM", "centimeters").getDict()
      else:
        unitsCode = self.makePrivateCode(units)

      modifiers.append({"modifier": imageTransformationConcept,
                        "modifierValue": self.makePrivateCode("Laplacian of Gaussian")})
      derivationParameters.append({"derivationParameter": self.makePrivateCode("Kernel size"),
                                   "derivationParameterValue": str('.'.join([match.group(1), match.group(2)])),
                                   "derivationParameterUnits": unitsCode})
    else:
      # unknown prefix
      return None

    return modifiers, derivationParameters

  # adds a single measurement to the last measurement group
  def addMeasurement(
    self,
    value,
    quantityCode,
    unitsCode=CodedValue(
      "1",
      "UCUM",
      "no units"
    )):
    if self.measurementGroupCount < 1:
      scriptlogger.error(
        "Cannot add measurement - no measurement groups initialized!")
      return

    (preprocessing, featureClass, featureName) = quantityCode.split('_')

    mpTuple = self.prefix2codes(preprocessing)
    if mpTuple is None:
      return

    measurement = {}

    classSubset = self.featuresDict[self.featuresDict['pyradiomics_feature_class'] == featureClass]
    featureTuple = classSubset[classSubset['pyradiomics_feature_name'] == featureName]

    if featureTuple.empty:
      codeMeaning = featureClass + "_" + featureName
      code = self.makeHash(codeMeaning)
      measurement["quantity"] = CodedValue(
        code, "99PYRADIOMICS", codeMeaning).getDict()
      if len(code) > 16:
        scriptlogger.error("Sorry, the code value is too long!")
        sys.exit()
    else:
      measurement["quantity"] = CodedValue(
        featureTuple["IBSI_code"].values[0],
        "IBSI",
        featureTuple["IBSI_meaning"].values[0]).getDict()

    try:
      if numpy.isnan(value):
        scriptlogger.info(
          "Skipping NaN value for feature %s",
          quantityCode)
        return
    except Exception as e:
      scriptlogger.error("Exception checking for NaN: %s %s", str(e), value)
      return

    try:
      measurement["value"] = '%E' % Decimal(float(value))
    except Exception as e:
      scriptlogger.error("Exception formatting %s as Decimal: %s", value, str(e))
      scriptlogger.error("type of value: %s", type(value))

    measurement["units"] = unitsCode.getDict()

    self.m["Measurements"][-1]["measurementItems"].append(measurement)

    if len(mpTuple[0]):
      measurement["measurementModifiers"] = [m for m in mpTuple[0]]
    if len(mpTuple[1]):
      measurement["measurementDerivationParameters"] = [
        d for d in mpTuple[1]]

    return

  def saveJSONToFile(self, fileName):
    with open(fileName, 'w') as f:
      json.dump(self.m, f, indent=2, sort_keys=True)


def main():
  parser = argparse.ArgumentParser(
    usage="""%(prog)s --input-image <dir> --input-seg <name> --output-sr <name>\n\n

Warning: This is a \"pyradiomics labs\" script, which means it is an experimental feature in development!
The intent of this helper script is to enable pyradiomics feature extraction directly from/to DICOM data.
The segmentation defining the region of interest must be defined as a DICOM Segmentation image.
Support for DICOM Radiotherapy Structure Sets for defining region of interest may be added in the future.""")
  parser.add_argument(
    '--input-image-dir',
    dest="inputDICOMImageDir",
    metavar="<folder>",
    help="Path to the directory with the input DICOM series."
         " It is expected that a single series is corresponding to a single scalar volume.",
    required=True)
  parser.add_argument(
    '--input-seg-file',
    dest="inputSEG",
    metavar="<file>",
    help="Path to the input segmentation defined as a DICOM Segmentation object.",
    required=True)
  parser.add_argument(
    '--output-dir',
    dest="outputDir",
    metavar="<folder>",
    help="Path to the directory for saving the resulting DICOM file.",
    required=True)
  parser.add_argument(
    '--parameters',
    dest="parameters",
    metavar="<parameters>",
    help="Pyradiomics feature extractor positional arguments")
  parser.add_argument(
    '--temp-dir',
    dest="tempDir",
    metavar="<folder>",
    help="Path to the directory to store intermediate results")
  parser.add_argument(
    '--features-dict',
    dest="featuresDict",
    metavar="<file>",
    help="Path to the dictionary mapping pyradiomics feature names to the IBSI defined features.")
  parser.add_argument(
    '--volume-reconstructor',
    dest="volumeReconstructor",
    metavar="<plastimatch or dcm2niix>",
    help="Choose the tool to be used for reconstructing image volume from the DICOM image series."
         " Allowed options are plastimatch or dcm2niix (should be installed on the system). plastimatch"
         " will be used by default.",
    choices=['plastimatch', 'dcm2niix'],
    default="plastimatch")
  parser.add_argument(
    '--geometry-tolerance',
    dest="geometryTolerance",
    metavar="<number>",
    help="Decimal number setting geometry tolerance for the extractor. Defaults to 1e-6.",
    default=1e-6)
  parser.add_argument(
    '--correct-mask',
    dest="correctMask",
    help="Boolean flag argument. If present, PyRadiomics will attempt to resample the mask to the image"
         " geometry if the mask check fails.",
    action='store_true',
    default=False)

  args = parser.parse_args()

  # with tempfile.mkdtemp() as tempDir:
  tempDir = args.tempDir
  if not tempDir:
    tempDir = tempfile.mkdtemp()

  scriptlogger.info("Temporary directory: " + tempDir)

  # convert input DICOM series into a scalar volume
  # plastimatch fails for prostate DWI Data! Need to report
  # Selection of the optimal volume reconstructor may depend
  # on the specific dataset!
  if args.volumeReconstructor == "plastimatch":
    scriptlogger.info(
      "Using Plastimatch for DICOM image volume reconstruction.")
    inputImage = dcmImageToNRRD(args.inputDICOMImageDir, tempDir)
  else:
    scriptlogger.info(
      "Using dcm2niix for DICOM image volume reconstruction.")
    inputImage = dcmImageToNIfTI(args.inputDICOMImageDir, tempDir)

  # convert segmentation into segments
  inputSegments = dcmSEGToNRRDs(args.inputSEG, tempDir)
  if len(inputSegments) == 0:
    scriptlogger.error("No segments found. Cannot compute features.")
    return -1

  featuresDir = os.path.join(tempDir, 'Features')
  if not os.path.isdir(featuresDir):
    os.mkdir(featuresDir)

  # initialize Metadata for the individual features
  # TODO: to be replaced with direct mapping in the pyradiomics feature functions
  # see https://github.com/Radiomics/pyradiomics/issues/435
  if args.featuresDict is not None:
    featuresDictPath = args.featuresDict
  else:
    featuresDictPath = "featuresDict.tsv"

  if not os.path.exists(featuresDictPath):
    scriptlogger.error(
      "Features dictionary file %s is not found!",
      featuresDictPath)
    return -1

  m = TID1500Metadata(featuresDictPath)

  # find a valid DICOM file in the input image DICOM directory
  dicomImage = None
  for f in os.listdir(args.inputDICOMImageDir):
    try:
      pydicom.read_file(os.path.join(args.inputDICOMImageDir, f))
      dicomImage = os.path.join(args.inputDICOMImageDir, f)
      break
    except BaseException:
      continue

  if dicomImage is None:
    scriptlogger.error(
      "Input DICOM image directory does not seem to contain any valid DICOM files!")
    return -1

  imageMetadataAccessor = DICOMMetadataAccessor(
    os.path.join(args.inputDICOMImageDir, f))
  segmentationMetadataAccessor = SEGMetadataAccessor(args.inputSEG)

  pyradiomicsVersion = None

  for inputSegment in inputSegments:
    scriptlogger.debug("Processing segmentation file %s", inputSegment)
    segmentNumber = os.path.split(inputSegment)[-1].split('.')[0]

    try:
      scriptlogger.debug("Initializing extractor")
      extractionSettings = {
        "geometryTolerance": float(args.geometryTolerance),
        "correctMask": True if args.correctMask else False
      }
      params = []
      if args.parameters is not None:
        params = [args.parameters]
      extractor = featureextractor.RadiomicsFeatureExtractor(*params, **extractionSettings)

    except Exception:
      scriptlogger.error(
        'Initialization of the pyradimics feature extraction failed.', exc_info=True)
      return -1

    featureVector = extractor.execute(
      inputImage, inputSegment, int(segmentNumber))

    if len(featureVector) == 0:
      scriptlogger.error("No features extracted!")
      return -1

    featuresFileName = os.path.join(featuresDir, segmentNumber + '.csv')
    scriptlogger.debug("Will save features as %s", featuresFileName)
    writer = csv.writer(open(featuresFileName, 'w'), lineterminator='\n')
    headers = list(featureVector.keys())
    writer.writerow(headers)

    row = []
    for h in headers:
      row.append(featureVector.get(h, ""))
    writer.writerow(row)

    scriptlogger.debug("Initializing TID 1500 Measurement groups.")
    m.addMeasurementGroup()

    includedFeatureVectorItems = 0
    for featureName in featureVector.keys():
      if featureName == 'diagnostics_Versions_PyRadiomics':
        pyradiomicsVersion = featureVector[featureName]
        continue

      featureValue = featureVector[featureName]
      featureNameSplit = featureName.split('_')
      if len(featureNameSplit) < 3:
        scriptlogger.warning(
          "Skipping unrecognized feature %s",
          featureName)
        continue
      includedFeatureVectorItems += 1
      m.addMeasurement(featureValue, featureName)
    scriptlogger.debug(
      "%d of %d total features included in the TID 1500 Measurement group.",
      len(featureVector), includedFeatureVectorItems)

    # initialize metadata common to all measurements
    scriptlogger.debug("Populating common metadata")
    m.m["Measurements"][-1]["SourceSeriesForImageSegmentation"] = imageMetadataAccessor.getSeriesInstanceUID()
    m.m["Measurements"][-1]["segmentationSOPInstanceUID"] = segmentationMetadataAccessor.getInstanceUID()

    # TODO: populate those from SEG SegmentationType / AnatomicLocation
    segmentationType = segmentationMetadataAccessor.getSegmentSegmentationTypeCode(
      int(segmentNumber) - 1)
    if segmentationType:
      m.m["Measurements"][-1]["Finding"] = CodedValue(segmentationType.CodeValue,
                                                      segmentationType.CodingSchemeDesignator,
                                                      segmentationType.CodeMeaning).getDict()

    segTrackingIdentifier = segmentationMetadataAccessor.getTrackingIdentifier(int(segmentNumber)-1)
    segTrackingUniqueIdentifier = segmentationMetadataAccessor.getTrackingUniqueIdentifier(int(segmentNumber)-1)

    if segTrackingIdentifier:
      m.m["Measurements"][-1]["TrackingIdentifier"] = segTrackingIdentifier
    else:
      m.m["Measurements"][-1]["TrackingIdentifier"] = segmentationType.CodeMeaning
      segmentDescription = segmentationMetadataAccessor.getSegmentDescription(int(segmentNumber)-1)
      # SegmentDescription is Type 3, and can be missing
      if segmentDescription is not None:
        m.m["Measurements"][-1]["TrackingIdentifier"] = segmentationType.CodeMeaning+" - "+segmentDescription

    if segTrackingUniqueIdentifier:
      m.m["Measurements"][-1]["TrackingUniqueIdentifier"] = segTrackingUniqueIdentifier

    segmentationLocation = segmentationMetadataAccessor.getSegmentAnatomicLocationCode(
      int(segmentNumber) - 1)
    if segmentationLocation:
      m.m["Measurements"][-1]["FindingSite"] = CodedValue(segmentationLocation.CodeValue,
                                                          segmentationLocation.CodingSchemeDesignator,
                                                          segmentationLocation.CodeMeaning).getDict()

    # AlgorithmIdentification
    m.m["Measurements"][-1]["measurementAlgorithmIdentification"] = {}
    m.m["Measurements"][-1]["measurementAlgorithmIdentification"]["AlgorithmName"] = "https://github.com/Radiomics/pyradiomics"
    m.m["Measurements"][-1]["measurementAlgorithmIdentification"]["AlgorithmVersion"] = pyradiomicsVersion
    m.m["Measurements"][-1]["measurementAlgorithmIdentification"]["AlgorithmParameters"] = [json.dumps(extractor.settings)]

  m.m["observerContext"] = {}
  m.m["observerContext"]["ObserverType"] = "DEVICE"
  m.m["observerContext"]["DeviceObserverName"] = "pyradiomics"
  m.m["observerContext"]["DeviceObserverModelName"] = pyradiomicsVersion

  m.m["compositeContext"] = [os.path.split(args.inputSEG)[-1]]
  m.m["imageLibrary"] = [os.path.split(f)[-1]
                         for f in os.listdir(args.inputDICOMImageDir)]

  m.m["SeriesDescription"] = segmentationMetadataAccessor.getSeriesDescription() + ' - pyradiomics features'

  scriptlogger.debug("Saving temporary files for DICOM SR writer.")
  dcmqiMetadataFile = os.path.join(featuresDir, "dcmqi_sr.json")
  outputSRTempFile = os.path.join(featuresDir, "sr.dcm")
  m.saveJSONToFile(dcmqiMetadataFile)

  scriptlogger.debug("Generating DICOM SR.")
  writeSR(
    args.inputSEG,
    dcmqiMetadataFile,
    args.inputDICOMImageDir,
    outputSRTempFile)

  # copy to the dest directory under UID as a name
  try:
    dcm = pydicom.read_file(outputSRTempFile)
    shutil.move(
      outputSRTempFile,
      os.path.join(args.outputDir, dcm.SOPInstanceUID + ".dcm"))
  except BaseException:
    scriptlogger.error("Failed to move output SR!")


if __name__ == "__main__":
  exeFound = {}
  for exe in ['tid1500writer', 'dcm2niix', 'plastimatch', 'segimage2itkimage']:
    if distutils.spawn.find_executable(exe) is None:
      exeFound[exe] = False
    else:
      exeFound[exe] = True
  if not (exeFound['tid1500writer'] and exeFound['segimage2itkimage']) or not (
    exeFound['plastimatch'] or exeFound['dcm2niix']):
    scriptlogger.error(
      "Dependency converter(s) not found in the path.")
    scriptlogger.error(
      "dcmqi (https://github.com/qiicr/dcmqi), and dcm2niix (https://github.com/rordenlab/dcm2niix/releases)")
    scriptlogger.error("or Plastimatch (http://plastimatch.org/)")
    scriptlogger.error(
      "need to be installed and available in the PATH for using this converter script.")
    sys.exit()
  main()
