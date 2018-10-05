import argparse, tempfile, os, sys, glob, json, csv, collections, pandas, pydicom, numpy, shutil
import distutils.spawn
from decimal import *
from subprocess import call
from pathlib import Path

def dcmImageToNRRD(inputDICOMImageDir, tempDir):
  scanNRRDFile = os.path.join(tempDir,"image.nrrd")
  if not os.path.isfile(scanNRRDFile):
    call(['plastimatch', 'convert','--input',inputDICOMImageDir,'--output-img',scanNRRDFile])
  return scanNRRDFile

def dcmImageToNIfTI(inputDICOMImageDir, tempDir):
  destScanNIfTIFile = os.path.join(tempDir,"volume.nii")
  scanNIfTIFile = os.path.join(inputDICOMImageDir,"volume.nii")
  scanJSONFile = os.path.join(inputDICOMImageDir,"volume.json")
  # will save to volume.nii
  if not os.path.isfile(destScanNIfTIFile):
    cmd = ['dcm2niix', "-m", "y", inputDICOMImageDir]
    call(cmd)
    shutil.move(scanNIfTIFile,destScanNIfTIFile)
    os.remove(scanJSONFile)
  return destScanNIfTIFile

def dcmSEGToNRRDs(inputSEG, tempDir):
  segmentsDir = os.path.join(tempDir,'Segments')
  if not Path(segmentsDir).is_dir():
    os.mkdir(segmentsDir)
  call(['segimage2itkimage', '--inputDICOM', inputSEG, '--outputDirectory', segmentsDir])
  return glob.glob(os.path.join(segmentsDir,"*nrrd"))

def writeSR(inputSEG,inputJSON,inputDICOMImageDir,outputSR):
  cmd = ['tid1500writer', '--inputImageLibraryDirectory', inputDICOMImageDir, '--inputCompositeContextDirectory', os.path.split(inputSEG)[0],'--inputMetadata',inputJSON,'--outputDICOM',outputSR]
  print("Writing SR with: "+str(cmd))
  call(cmd)

def getCTSeriesUID(imageDICOMDir):
  ctFile = os.listdir(imageDICOMDir)[0]
  dcm = pydicom.read_file(os.path.join(imageDICOMDir,ctFile))
  return dcm.SeriesInstanceUID

class DICOMMetadataAccessor:
  def __init__(self, dcmFileName):
    print("Calling parent constructor with "+dcmFileName)
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
      raise ValueError("SEGMetadataAccessor: DICOM object is not Segmentation!")

  def getSegmentSegmentationTypeCode(self,segmentNumber):
    try:
      return self.dcm.SegmentSequence[segmentNumber].SegmentedPropertyTypeCodeSequence[0]
    except:
      return None

  def getSegmentAnatomicLocationCode(self, segmentNumber):
    try:
      return self.dcm.SegmentSequence[segmentNumber].AnatomicRegionSequence[0]
    except:
      return None

class CodedValue:
  def __init__(self,codeValue):
    self.codeValue = codeValue
    self.codingSchemeDesignator = "99PYRADIOMICS"
    self.codeMeaning = codeValue

  def __init__(self,value,scheme="99PYRADIOMICS",meaning=None):
    self.codeValue = value
    self.codingSchemeDesignator = scheme
    if meaning is None:
      self.codeMeaning = value
    else:
      self.codeMeaning = meaning

  def getDict(self):
    return {"CodeValue": self.codeValue, "CodeMeaning": self.codeMeaning, "CodingSchemeDesignator": self.codingSchemeDesignator}

class Metadata:

  m = {}

  def __init__(self, featuresDictFile, seriesDescription="Radiomics features"):

    self.featuresDict = None
    self.readDictionary(featuresDictFile)

    self.m["@schema"] = "https://raw.githubusercontent.com/qiicr/dcmqi/master/doc/schemas/sr-tid1500-schema.json#"
    self.m["SeriesDescription"] = seriesDescription

    self.m["Measurements"] = []
    self.measurementGroupCount = 0

    imageTransformationMeaning = "Image transformation"
    LoGMeaning = "Laplacian of Gaussian"
    logarithmMeaning = "Logarithm transformation"
    waveletMeaning = "Wavelet transformation"
    sqrMeaning = "Square transformation"
    sqrtMeaning = "Square root transformation"
    kernelSizeMeaning = "Kernel size"
    waveletSubband = "Wavelet subband"
    gradientMeaning = "Gradient magnitude transformation"
    expMeaning = "Exponent transformation"

    self.LoGTransformation = { "modifier": CodedValue(self.makeHash(imageTransformationMeaning), "99PYRADIOMICS", imageTransformationMeaning).getDict(), \
                              "modifierValue": CodedValue(self.makeHash(LoGMeaning), "99PYRADIOMICS", LoGMeaning).getDict()}
    self.waveletTransformation = { "modifier": CodedValue(self.makeHash(imageTransformationMeaning), "99PYRADIOMICS", imageTransformationMeaning).getDict(), \
                              "modifierValue": CodedValue(self.makeHash(waveletMeaning), "99PYRADIOMICS", waveletMeaning).getDict()}
    self.gradient = { "modifier": CodedValue(self.makeHash(imageTransformationMeaning), "99PYRADIOMICS", imageTransformationMeaning).getDict(), \
                              "modifierValue": CodedValue(self.makeHash(gradientMeaning), "99PYRADIOMICS", gradientMeaning).getDict()}
    self.sqr = { "modifier": CodedValue(self.makeHash(imageTransformationMeaning), "99PYRADIOMICS", imageTransformationMeaning).getDict(), \
          "modifierValue": CodedValue(self.makeHash(sqrMeaning), "99PYRADIOMICS", sqrMeaning).getDict()}
    self.sqrt = { "modifier": CodedValue(self.makeHash(imageTransformationMeaning), "99PYRADIOMICS", imageTransformationMeaning).getDict(), \
                              "modifierValue": CodedValue(self.makeHash(sqrtMeaning), "99PYRADIOMICS", sqrtMeaning).getDict()}
    # note that some of those operations, such as log, include additional processing and are not pure application of mathematical
    # operations. For example, pyradiomics logarithm operation includes normalization and calculation of absolute values,
    # prior to calculating natural logarithm
    self.logarithm = { "modifier": CodedValue(self.makeHash(imageTransformationMeaning), "99PYRADIOMICS", imageTransformationMeaning).getDict(), \
          "modifierValue": CodedValue(self.makeHash(logarithmMeaning), "99PYRADIOMICS", logarithmMeaning).getDict()}
    self.exp = { "modifier": CodedValue(self.makeHash(imageTransformationMeaning), "99PYRADIOMICS", imageTransformationMeaning).getDict(), \
          "modifierValue": CodedValue(self.makeHash(expMeaning), "99PYRADIOMICS", expMeaning).getDict()}

    self.LoGParameter = { "derivationParameter": CodedValue(self.makeHash(kernelSizeMeaning), "99PYRADIOMICS", kernelSizeMeaning).getDict()}
    self.waveletSubband = { "modifier": CodedValue(self.makeHash(waveletSubband), "99PYRADIOMICS", waveletSubband).getDict(), \
                      "modifierValue": CodedValue("", "99PYRADIOMICS", "").getDict()}

  def addMeasurementGroup(self):
    self.measurementGroupCount = self.measurementGroupCount+1
    measurementsGroup = {}
    measurementsGroup["measurementItems"] = []
    measurementsGroup["ReferencedSegment"] = self.measurementGroupCount
    self.m["Measurements"].append(measurementsGroup)

  def readDictionary(self, featuresDictFile):
    # read features dictionary
    self.featuresDict = pandas.read_csv(featuresDictFile, sep='\t', low_memory=False)

  def makeHash(self,text, length=6):
    import base64, hashlib
    return base64.b64encode(hashlib.sha1(str.encode(text)).digest()).decode('ascii')[:length]

  # adds a single measurement to the last measurement group
  def addMeasurement(self,value,quantityCode,unitsCode=CodedValue("1","UCUM","no units")):
    if self.measurementGroupCount<1:
      print("ERROR: Cannot add measurement - no measurement groups initialized!")
      return

    (preprocessing,featureClass,featureName) = quantityCode.split('_')
    measurement = {}

    classSubset = self.featuresDict[self.featuresDict['pyradiomics_feature_class']==featureClass]
    featureTuple = classSubset[classSubset['pyradiomics_feature_name']==featureName]

    if featureTuple.empty:
        codeMeaning = featureClass+"_"+featureName
        import base64, hashlib
        # use first 16 characters of the hash as the code

        code = self.makeHash(codeMeaning)
        measurement["quantity"] = {"CodeValue": code, "CodeMeaning":   codeMeaning, "CodingSchemeDesignator": "99PYRADIOMICS"}
        if len(code)>16:
          print("Sorry, the code value is too long!")
          sys.exit()
    else:
      measurement["quantity"] = {"CodeValue": featureTuple["IBSI_code"].values[0], "CodeMeaning":   featureTuple["IBSI_meaning"].values[0], "CodingSchemeDesignator": "IBSI"}

    if numpy.isnan(value):
      return

    #print(str(featureTuple))

    measurement["value"] = '%E' % Decimal(value)

    measurement["units"] = unitsCode.getDict()

    #print(str(measurement))

    #print(measurement)
    measurementItems = self.m["Measurements"][-1]["measurementItems"]
    self.m["Measurements"][-1]["measurementItems"].append(measurement)

    # parse preprocessing parameters
    # Ignore for now!


    # encode preprocessing parameters as well, where applicable
    if preprocessing.startswith("log-"): # not log, since it could be logarithm
      # example: log-sigma-1-0-mm-3D
      try:
        LoGparameters = preprocessing.split("-")
        sigma = str(float(LoGparameters[2])+float("."+LoGparameters[3]))
      except IndexError:
        print("ERROR: Failed to parse LoG kernel from: %s!" % (preprocessing))
        sys.exit()
      except TypeError:
        print("ERROR: Failed to convert LoG kernel from: %s!" % (preprocessing))
        sys.exit()

      measurement["measurementModifiers"] = [self.LoGTransformation]
      measurement["measurementDerivationParameters"] = [self.LoGParameter]
      measurement["measurementDerivationParameters"][0]["derivationParameterValue"] = sigma
      measurement["measurementDerivationParameters"][0]["derivationParameterUnits"] = CodedValue("mm","UCUM","millimeters").getDict()

    elif preprocessing.startswith("wavelet"):
      try:
        subband = preprocessing.split("-")[1]
      except IndexError:
        print("ERROR: Failed to parse wavelet subband from %s!" % (preprocessing))
        sys.exit()

      measurement["measurementModifiers"] = [self.waveletTransformation, self.waveletSubband]
      measurement["measurementModifiers"][1]["modifierValue"]["CodeValue"] = self.makeHash(subband)
      measurement["measurementModifiers"][1]["modifierValue"]["CodeMeaning"] = subband

    elif preprocessing.startswith("logarithm"):
      measurement["measurementModifiers"] = [self.logarithm]

    elif preprocessing == "square":
      measurement["measurementModifiers"] = [self.sqr]

    elif preprocessing == "squareroot":
      measurement["measurementModifiers"] = [self.sqrt]

    elif preprocessing == "gradient":
      measurement["measurementModifiers"] = [self.gradient]

    elif preprocessing == "exponential":
      measurement["measurementModifiers"] = [self.exp]

    elif preprocessing == "original":
      # nothing to do here
      pass

    else:
      print('ERROR: Unrecognized pre-processing type: %s!' % (preprocessing))
      sys.exit()

    return

  def saveJSONToFile(self,fileName):
    import json
    with open(fileName,'w') as f:
      json.dump(self.m,f,indent=2, sort_keys=True)

def main():
  parser = argparse.ArgumentParser(usage="%(prog)s --input-image <dir> --input-seg <name> --output-sr <name>")
  parser.add_argument('--input-image', dest="inputDICOMImageDir", metavar='Input DICOM image directory',  help="Directory with the input DICOM series."
                      " It is expected that a single series is corresponding to a single scalar volume.", required=True)
  parser.add_argument('--input-seg', dest="inputSEG", metavar='Input DICOM SEG file', help="Input segmentation defined as a"
                      "DICOM Segmentation object.", required=True)
  parser.add_argument('--output-dir', dest="outputDir", metavar='Directory to store the output file', help="Directory for saving the resulting DICOM file.", required=True)
  parser.add_argument('--parameters', dest="parameters", metavar="pyradiomics extraction parameters")
  parser.add_argument('--temp-dir', dest="tempDir", metavar="Temporary directory")
  parser.add_argument('--features-dict', dest="featuresDict", metavar="Dictionary mapping pyradiomics feature names to the IBSI defined features.")
  parser.add_argument('--volume-reconstructor', dest="volumeReconstructor", metavar="Choose the tool to be used for reconstructing image volume from the DICOM image series. Allowed options are plastimatch or dcm2niix (should be installed on the system). plastimatch will be used by default.", default="plastimatch")

  args = parser.parse_args()

  #with tempfile.mkdtemp() as tempDir:
  tempDir = args.tempDir
  if not tempDir:
    tempDir = tempfile.mkdtemp()

  print("Temporary directory: "+tempDir)
  # convert input DICOM series into a scalar volume
  # plastimatch fails for prostate DWI Data! Need to report
  ##inputImage = dcmImageToNRRD(args.inputDICOMImageDir, tempDir)
  if args.volumeReconstructor == "plastimatch":
    inputImage = dcmImageToNRRD(args.inputDICOMImageDir, tempDir)
  elif args.volumeReconstructor == "dcm2niix":
    inputImage = dcmImageToNIfTI(args.inputDICOMImageDir, tempDir)
  else:
    print("Unknown option for volume reconstruction: %s. Allowed values are plastimatch and dcm2niix.")
    return -1

  # convert segmentation into segments
  inputSegments = dcmSEGToNRRDs(args.inputSEG, tempDir)

  # Initialize extractor
  import radiomics
  from radiomics import featureextractor

  featuresDir = os.path.join(tempDir,'Features')
  if not Path(featuresDir).is_dir():
    os.mkdir(featuresDir)

  # initialize Metadata for the individual features
  if args.featuresDict is not None:
    m = Metadata(args.featuresDict)
  else:
    m = Metadata("featuresDict.tsv")

  # find a valid DICOM file in the input image DICOM directory
  dicomImage = None
  for f in os.listdir(args.inputDICOMImageDir):
    try:
      pydicom.read_file(os.path.join(args.inputDICOMImageDir,f))
      dicomImage = os.path.join(args.inputDICOMImageDir,f)
      break
    except:
      continue

  imageMetadataAccessor = DICOMMetadataAccessor(os.path.join(args.inputDICOMImageDir,f))
  segmentationMetadataAccessor = SEGMetadataAccessor(args.inputSEG)

  for inputSegment in inputSegments:
    print("Processing "+inputSegment)
    segmentNumber = os.path.split(inputSegment)[-1].split('.')[0]

    try:
      print("Initializing extractor")
      parameters = collections.OrderedDict()
      #parameters["setting"] = {"geometryTolerance": 0.001, "correctMask": True}

      print("Will init extractor")
      extractor = featureextractor.RadiomicsFeaturesExtractor(parameters)
      if args.parameters is not None:
        print("Will init params")
        extractor.loadParams(args.parameters)

    except Exception:
      print('EXTRACTOR INITIALIZATION FAILED')
      exit(-1)
    featureVector = extractor.execute(inputImage, inputSegment, int(segmentNumber))

    if len(featureVector) == 0:
      print("No features extracted!")
      return

    #json.dump(featureVector, open(os.path.join(featuresDir,os.path.split(inputSegment+".json")[-1]),'w'))
    # Note - this one should start from 1!
    writer = csv.writer(open(os.path.join(featuresDir,segmentNumber+'.csv'),'w'), lineterminator='\n')
    headers = list(featureVector.keys())
    writer.writerow(headers)

    row = []
    for h in headers:
      row.append(featureVector.get(h, ""))
    writer.writerow(row)

    m.addMeasurementGroup()

    for featureName in featureVector.keys():
      featureValue = featureVector[featureName]
      featureNameSplit = featureName.split('_')
      if len(featureNameSplit)<3 or featureNameSplit[0] == "general":
        continue
      m.addMeasurement(featureValue, featureName)

    # initialize metadata common to all measurements
    m.m["Measurements"][-1]["SourceSeriesForImageSegmentation"] = imageMetadataAccessor.getSeriesInstanceUID()
    m.m["Measurements"][-1]["segmentationSOPInstanceUID"] = segmentationMetadataAccessor.getInstanceUID()
    m.m["Measurements"][-1]["TrackingIdentifier"] = "Region"

    # TODO: populate those from SEG SegmentationType / AnatomicLocation
    segmentationType = segmentationMetadataAccessor.getSegmentSegmentationTypeCode(int(segmentNumber)-1)
    if segmentationType:
      m.m["Measurements"][-1]["Finding"] = {  "CodeValue": segmentationType.CodeValue,
                                              "CodingSchemeDesignator": segmentationType.CodingSchemeDesignator,
                                              "CodeMeaning": segmentationType.CodeMeaning}
      m.m["Measurements"][-1]["TrackingIdentifier"] = segmentationType.CodeMeaning

    segmentationLocation = segmentationMetadataAccessor.getSegmentAnatomicLocationCode(int(segmentNumber)-1)
    if segmentationLocation:
      m.m["Measurements"][-1]["FindingSite"] = { "CodeValue": segmentationLocation.CodeValue,
                                            "CodingSchemeDesignator": segmentationLocation.CodingSchemeDesignator,
                                            "CodeMeaning": segmentationLocation.CodeMeaning}

  m.m["observerContext"] = {}
  m.m["observerContext"]["ObserverType"] = "PERSON"
  m.m["observerContext"]["PersonObserverName"] = "Reader1"
  m.m["compositeContext"] = [os.path.split(args.inputSEG)[-1]]
  m.m["imageLibrary"] = [os.path.split(f)[-1] for f in os.listdir(args.inputDICOMImageDir)]



  m.m["SeriesDescription"] = segmentationMetadataAccessor.getSeriesDescription()+' - pyradiomics features'

  dcmqiMetadataFile = os.path.join(featuresDir,"dcmqi_sr.json")
  outputSRTempFile = os.path.join(featuresDir,"sr.dcm")
  m.saveJSONToFile(dcmqiMetadataFile)

  writeSR(args.inputSEG,dcmqiMetadataFile,args.inputDICOMImageDir,outputSRTempFile)

  # copy to the dest directory under UID as a name
  try:
    dcm = pydicom.read_file(outputSRTempFile)
    shutil.move(outputSRTempFile,os.path.join(args.outputDir,dcm.SOPInstanceUID+".dcm"))
  except:
    print("Failed to read output SR")


if __name__ == "__main__":
  for exe in ['tid1500writer','dcm2niix','segimage2itkimage']:
    if distutils.spawn.find_executable(exe) is None:
      print("ERROR: Dependency converter not found in the path: " +exe)
      print("dcmqi (https://github.com/qiicr/dcmqi) and dcm2niix (https://github.com/rordenlab/dcm2niix/releases)")
      print("need to be installed and available in the PATH for using this converter script.")
      sys.exit()
  main()
