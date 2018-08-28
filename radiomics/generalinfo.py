import collections
import logging

import numpy
import pywt
import SimpleITK as sitk
import six

import radiomics


class GeneralInfo():
  def __init__(self, imagePath, maskPath, resampledMask, settings, enabledImageTypes):
    self.logger = logging.getLogger(self.__module__)

    self.elements = self._getElementNames()

    if isinstance(imagePath, six.string_types):
      self.image = sitk.ReadImage(imagePath)
    elif isinstance(imagePath, sitk.Image):
      self.image = imagePath
    else:
      self.logger.warning('Error reading image Filepath or SimpleITK object')
      self.image = None

    if isinstance(maskPath, six.string_types):
      self.mask = sitk.ReadImage(maskPath)
    elif isinstance(maskPath, sitk.Image):
      self.mask = maskPath
    else:
      self.logger.warning('Error reading mask Filepath or SimpleITK object')
      self.mask = None

    self.resampledMask = resampledMask

    self._settings = settings
    self._enabledImageTypes = enabledImageTypes

    self.label = self._settings.get('label', 1)

    if resampledMask is not None:
      self.lssif = sitk.LabelShapeStatisticsImageFilter()
      self.lssif.Execute(resampledMask)
    else:
      self.lssif = None

  def _getElementNames(self):
    return [member[3: -5] for member in dir(self) if member.startswith('get') and member.endswith('Value')]

  def execute(self):
    """
    Return a dictionary containing all general info items. Format is <info_item>:<value>, where the type
    of the value is preserved. For CSV format, this will result in conversion to string and quotes where necessary, for
    JSON, the values will be interpreted and stored as JSON strings.
    """
    generalInfo = collections.OrderedDict()
    for el in self.elements:
      generalInfo[el] = getattr(self, 'get%sValue' % el)()
    return generalInfo

  def getBoundingBoxValue(self):
    """
    Calculate and return the boundingbox extracted using the specified label.
    Elements 0, 1 and 2 are the x, y and z coordinates of the lower bound, respectively.
    Elements 3, 4 and 5 are the size of the bounding box in x, y and z direction, respectively.

    Values are based on the resampledMask.
    """
    if self.lssif is not None:
      return self.lssif.GetBoundingBox(self.label)
    else:
      return None

  def getGeneralSettingsValue(self):
    """
    Return a string representation of the general settings.
    Format is {<settings_name>:<value>, ...}.
    """
    return self._settings

  def getImageHashValue(self):
    """
    Returns the sha1 hash of the image. This enables checking whether two images are the same,
    regardless of the file location.

    If the reading of the image fails, an empty string is returned.
    """
    if self.image is not None:
      return sitk.Hash(self.image)
    else:
      return None

  def getImageSpacingValue(self):
    """
    Returns the original spacing (before any resampling) of the image.

    If the reading of the image fails, an empty string is returned.
    """
    if self.image is not None:
      return self.image.GetSpacing()
    else:
      return None

  def getCenterOfMassIndexValue(self):
    """
    Returns z, y and x coordinates of the center of mass of the ROI in terms of the image coordinate space (continuous index).

    Calculation is based on the original (non-resampled) mask.

    .. note::
      Because this represents the continuous index, the order of x, y and z is reversed, i.e. the first element is the z index, the second
      the y index and the last element is the x index.
    """
    if self.mask is not None:
      maskArray = sitk.GetArrayFromImage(self.mask)
      maskCoordinates = numpy.array(numpy.where(maskArray == self.label))
      center_index = numpy.mean(maskCoordinates, axis=1)
      return tuple(center_index)
    else:
      return None

  def getCenterOfMassValue(self):
    """
    Returns the real-world x, y and z coordinates of the center of mass of the ROI. This is the real-world transformation of
    :py:func:`~radiomics.generalinfo.getCenterOfMassIndexValue()`, taking into account the spacing, direction and origin of the mask.

    Calculation is based on the original (non-resampled) mask.
    """
    if self.mask is not None:
      return self.mask.TransformContinuousIndexToPhysicalPoint(self.getCenterOfMassIndexValue())
    else:
      return None

  def getEnabledImageTypesValue(self):
    """
    Return a string representation of the enabled image types and any custom settings for each image type.
    Format is {<imageType_name>:{<setting_name>:<value>, ...}, ...}.
    """
    return self._enabledImageTypes

  def getMaskHashValue(self):
    """
    Returns the sha1 hash of the mask. This enables checking whether two masks are the same,
    regardless of the file location.

    If the reading of the mask fails, an empty string is returned. Uses the original mask, specified in maskPath.
    """
    if self.mask is not None:
      return sitk.Hash(self.mask)
    else:
      return None

  @classmethod
  def getVersionValue(self):
    """
    Return the current version of this package.
    """
    return radiomics.__version__

  @classmethod
  def getNumpyVersionValue(self):
    """
    Return the current version of the numpy package, used for feature calculation.
    """
    return numpy.__version__

  @classmethod
  def getSimpleITKVersionValue(self):
    """
    Return the current version of the SimpleITK package, used for image processing.
    """
    return sitk.Version().VersionString()

  @classmethod
  def getPyWaveletVersionValue(self):
    """
    Return the current version of the PyWavelet package, used to apply the wavelet filter.
    """
    return pywt.__version__

  def getVolumeNumValue(self):
    """
    Calculate and return the number of zones within the mask for the specified label.
    A zone is defined as a group of connected neighbours that are segmented with the specified label, and a voxel is
    considered a neighbour using 26-connectedness for 3D and 8-connectedness for 2D.

    Values are based on the resampledMask.
    """
    if self.resampledMask is not None:
      labelMap = (self.resampledMask == self.label)
      ccif = sitk.ConnectedComponentImageFilter()
      ccif.FullyConnectedOn()
      ccif.Execute(labelMap)
      return ccif.GetObjectCount()
    else:
      return None

  def getVoxelNumValue(self):
    """
    Calculate and return the number of voxels that have been segmented using the specified label.

    Values are based on the resampledMask.
    """
    if self.lssif is not None:
      return self.lssif.GetNumberOfPixels(self.label)
    else:
      return None
