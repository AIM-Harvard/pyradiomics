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

    self.generalInfo_prefix = 'diagnostics_'

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

    self.generalInfo[self.generalInfo_prefix + 'Versions_PyRadiomics'] = radiomics.__version__
    self.generalInfo[self.generalInfo_prefix + 'Versions_Numpy'] = numpy.__version__
    self.generalInfo[self.generalInfo_prefix + 'Versions_SimpleITK'] = sitk.Version().VersionString()
    self.generalInfo[self.generalInfo_prefix + 'Versions_PyWavelet'] = pywt.__version__
    self.generalInfo[self.generalInfo_prefix + 'Versions_Python'] = '%i.%i.%i' % sys.version_info[:3]

  def addImageElements(self, image, prefix='original'):
    """
    Calculates provenance info for the image

    Adds the following:

    - Hash: sha1 hash of the mask, which can be used to check if the same mask was used during reproducibility
      tests. (Only added when prefix is "original")
    - Dimensionality: Number of dimensions (e.g. 2D, 3D) in the image. (Only added when prefix is "original")
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
      self.generalInfo[self.generalInfo_prefix + 'Image-original_Hash'] = sitk.Hash(image)
      self.generalInfo[self.generalInfo_prefix + 'Image-original_Dimensionality'] = '%iD' % image.GetDimension()

    self.generalInfo[self.generalInfo_prefix + 'Image-' + prefix + '_Spacing'] = image.GetSpacing()
    self.generalInfo[self.generalInfo_prefix + 'Image-' + prefix + '_Size'] = image.GetSize()
    im_arr = sitk.GetArrayFromImage(image).astype('float')
    self.generalInfo[self.generalInfo_prefix + 'Image-' + prefix + '_Mean'] = numpy.mean(im_arr)
    self.generalInfo[self.generalInfo_prefix + 'Image-' + prefix + '_Minimum'] = numpy.min(im_arr)
    self.generalInfo[self.generalInfo_prefix + 'Image-' + prefix + '_Maximum'] = numpy.max(im_arr)

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
    - CenterOfMassIndex: x, y and z coordinates of the center of mass of the ROI in terms of the image coordinate space
      (continuous index).
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
      self.generalInfo[self.generalInfo_prefix + 'Mask-original_Hash'] = sitk.Hash(mask)

    self.generalInfo[self.generalInfo_prefix + 'Mask-' + prefix + '_Spacing'] = mask.GetSpacing()
    self.generalInfo[self.generalInfo_prefix + 'Mask-' + prefix + '_Size'] = mask.GetSize()

    lssif = sitk.LabelShapeStatisticsImageFilter()
    lssif.Execute(mask)

    self.generalInfo[self.generalInfo_prefix + 'Mask-' + prefix + '_BoundingBox'] = lssif.GetBoundingBox(label)
    self.generalInfo[self.generalInfo_prefix + 'Mask-' + prefix + '_VoxelNum'] = lssif.GetNumberOfPixels(label)

    labelMap = (mask == label)
    ccif = sitk.ConnectedComponentImageFilter()
    ccif.FullyConnectedOn()
    ccif.Execute(labelMap)
    self.generalInfo[self.generalInfo_prefix + 'Mask-' + prefix + '_VolumeNum'] = ccif.GetObjectCount()

    ma_arr = sitk.GetArrayFromImage(labelMap) == 1
    maskCoordinates = numpy.array(numpy.where(ma_arr))
    center_index = tuple(numpy.mean(maskCoordinates, axis=1)[::-1])  # also convert z, y, x to x, y, z order

    self.generalInfo[self.generalInfo_prefix + 'Mask-' + prefix + '_CenterOfMassIndex'] = center_index

    self.generalInfo[self.generalInfo_prefix + 'Mask-' + prefix + '_CenterOfMass'] = mask.TransformContinuousIndexToPhysicalPoint(center_index)

    if image is None:
      return

    im_arr = sitk.GetArrayFromImage(image)
    targetvoxels = im_arr[ma_arr].astype('float')
    self.generalInfo[self.generalInfo_prefix + 'Mask-' + prefix + '_Mean'] = numpy.mean(targetvoxels)
    self.generalInfo[self.generalInfo_prefix + 'Mask-' + prefix + '_Minimum'] = numpy.min(targetvoxels)
    self.generalInfo[self.generalInfo_prefix + 'Mask-' + prefix + '_Maximum'] = numpy.max(targetvoxels)

  def addGeneralSettings(self, settings):
    """
    Add a string representation of the general settings.
    Format is {<settings_name>:<value>, ...}.
    """
    self.generalInfo[self.generalInfo_prefix + 'Configuration_Settings'] = settings

  def addEnabledImageTypes(self, enabledImageTypes):
    """
    Add a string representation of the enabled image types and any custom settings for each image type.
    Format is {<imageType_name>:{<setting_name>:<value>, ...}, ...}.
    """
    self.generalInfo[self.generalInfo_prefix + 'Configuration_EnabledImageTypes'] = enabledImageTypes
