import collections
import logging
import sys

import numpy
import pywt
import SimpleITK as sitk

import radiomics


class GeneralInfo:
  def __init__(self):
    self.logger = logging.getLogger(self.__module__)

    self.generalInfo_prefix = 'general_info_'

    self.generalInfo = collections.OrderedDict()
    self.addStaticElements()

  def getGeneralInfo(self):
    """
    Return a dictionary containing all general info items. Format is <info_item>:<value>, where the type
    of the value is preserved. For CSV format, this will result in conversion to string and quotes where necessary, for
    JSON, the values will be interpreted and stored as JSON strings.
    """
    return self.generalInfo

  def addStaticElements(self):
    """
    Adds the following elements to the general info:

    - Version: current version of PyRadiomics
    - NumpyVersion: version of numpy used
    - SimpleITKVersion: version SimpleITK used
    - PyWaveletVersion: version of PyWavelet used
    - PythonVersion: version of the python interpreter running PyRadiomics
    """

    self.generalInfo[self.generalInfo_prefix + 'Version'] = radiomics.__version__
    self.generalInfo[self.generalInfo_prefix + 'NumpyVersion'] = numpy.__version__
    self.generalInfo[self.generalInfo_prefix + 'SimpleITKVersion'] = sitk.Version().VersionString()
    self.generalInfo[self.generalInfo_prefix + 'PyWaveletVersion'] = pywt.__version__
    self.generalInfo[self.generalInfo_prefix + 'PythonVersion'] = '%i.%i.%i' % sys.version_info[:3]

  def addImageElements(self, image, prefix='original'):
    """
    Calculates provenance info for the image

    Adds the following:
    - ImageHash: sha1 hash of the mask, which can be used to check if the same mask was used during reproducibility
      tests. (Only added when prefix is "original")
    - Spacing: Pixel spacing (x, y, z) in mm.
    - Size: Dimensions (x, y, z) of the image in number of voxels.
    - Mean: Mean intensity value over all voxels in the image.
    - Minimum: Minimum intensity value among all voxels in the image.
    - Maximum: Maximum intensity value among all voxels in the image.

    A prefix is added to indicate what type of image is described:

    - original: Image as loaded, without pre-processing.
    - interpolated: Image after it has been resampled to a new spacing (includes cropping).
    """
    if prefix == 'original':
      self.generalInfo[self.generalInfo_prefix + 'ImageHash'] = sitk.Hash(image)

    self.generalInfo[self.generalInfo_prefix + prefix + 'Spacing'] = image.GetSpacing()
    self.generalInfo[self.generalInfo_prefix + prefix + 'Size'] = image.GetSize()
    im_arr = sitk.GetArrayFromImage(image)
    self.generalInfo[self.generalInfo_prefix + prefix + 'Mean'] = numpy.mean(im_arr)
    self.generalInfo[self.generalInfo_prefix + prefix + 'Minimum'] = numpy.min(im_arr)
    self.generalInfo[self.generalInfo_prefix + prefix + 'Maximum'] = numpy.max(im_arr)

  def addMaskElements(self, image, mask, label, prefix='original'):
    """
    Calculates provenance info for the mask

    Adds the following:

    - MaskHash: sha1 hash of the mask, which can be used to check if the same mask was used during reproducibility
      tests. (Only added when prefix is "original")
    - BoundingBox: bounding box of the ROI defined by the specified label:
      Elements 0, 1 and 2 are the x, y and z coordinates of the lower bound, respectively.
      Elements 3, 4 and 5 are the size of the bounding box in x, y and z direction, respectively.
    - VoxelNum: Number of voxels included in the ROI defined by the specified label.
    - VolumeNum: Number of fully connected (26-connectivity) volumes in the ROI defined by the specified label.
    - CenterOfMassIndex: x, y and z coordinates of the center of mass of the ROI in terms of the image coordinate space (continuous index).
    - CenterOfMass: the real-world x, y and z coordinates of the center of mass of the ROI
    - ROIMean: Mean intensity value over all voxels in the ROI defined by the specified label.
    - ROIMinimum: Minimum intensity value among all voxels in the ROI defined by the specified label.
    - ROIMaximum: Maximum intensity value among all voxels in the ROI defined by the specified label.

    A prefix is added to indicate what type of mask is described:

    - original: Mask as loaded, without pre-processing.
    - corrected: Mask after it has been corrected by :py:func:`imageoperations.checkMask`.
    - interpolated: Mask after it has been resampled to a new spacing (includes cropping).
    - resegmented: Mask after resegmentation has been applied.
    """
    if mask is None:
      return

    if prefix == 'original':
      self.generalInfo[self.generalInfo_prefix + 'MaskHash'] = sitk.Hash(mask)

    self.generalInfo[self.generalInfo_prefix + prefix + 'ROISpacing'] = mask.GetSpacing()
    self.generalInfo[self.generalInfo_prefix + prefix + 'ROISize'] = mask.GetSize()

    lssif = sitk.LabelShapeStatisticsImageFilter()
    lssif.Execute(mask)

    self.generalInfo[self.generalInfo_prefix + prefix + 'BoundingBox'] = lssif.GetBoundingBox(label)
    self.generalInfo[self.generalInfo_prefix + prefix + 'VoxelNum'] = lssif.GetNumberOfPixels(label)

    labelMap = (mask == label)
    ccif = sitk.ConnectedComponentImageFilter()
    ccif.FullyConnectedOn()
    ccif.Execute(labelMap)
    self.generalInfo[self.generalInfo_prefix + prefix + 'VolumeNum'] = ccif.GetObjectCount()

    ma_arr = sitk.GetArrayFromImage(labelMap) == 1
    maskCoordinates = numpy.array(numpy.where(ma_arr))
    center_index = tuple(numpy.mean(maskCoordinates, axis=1)[::-1])  # also convert z, y, x to x, y, z order

    self.generalInfo[self.generalInfo_prefix + prefix + 'CenterOfMassIndex'] = center_index

    self.generalInfo[self.generalInfo_prefix + prefix + 'CenterOfMass'] = mask.TransformContinuousIndexToPhysicalPoint(center_index)

    if image is None:
      return

    im_arr = sitk.GetArrayFromImage(image)
    targetvoxels = im_arr[ma_arr]
    self.generalInfo[self.generalInfo_prefix + prefix + 'ROIMean'] = numpy.mean(targetvoxels)
    self.generalInfo[self.generalInfo_prefix + prefix + 'ROIMinimum'] = numpy.min(targetvoxels)
    self.generalInfo[self.generalInfo_prefix + prefix + 'ROIMaximum'] = numpy.max(targetvoxels)

  def addGeneralSettings(self, settings):
    """
    Add a string representation of the general settings.
    Format is {<settings_name>:<value>, ...}.
    """
    self.generalInfo[self.generalInfo_prefix + 'GeneralSettings'] = settings

  def addEnabledImageTypes(self, enabledImageTypes):
    """
    Add a string representation of the enabled image types and any custom settings for each image type.
    Format is {<imageType_name>:{<setting_name>:<value>, ...}, ...}.
    """
    self.generalInfo[self.generalInfo_prefix + 'EnabledImageTypes'] = enabledImageTypes
