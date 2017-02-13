import collections
import logging
import SimpleITK as sitk
import radiomics


class GeneralInfo():
  def __init__(self, imagePath, maskPath, resampledMask, kwargs, inputImages):
    self.logger = logging.getLogger(self.__module__)

    if isinstance(imagePath, basestring):
      self.image = sitk.ReadImage(imagePath)
    elif isinstance(imagePath, sitk.Image):
      self.image = imagePath
    else:
      self.logger.warning('Error reading image Filepath or SimpleITK object')
      self.image = None

    if isinstance(maskPath, basestring):
      self.mask = sitk.ReadImage(maskPath)
    elif isinstance(maskPath, sitk.Image):
      self.mask = maskPath
    else:
      self.logger.warning('Error reading mask Filepath or SimpleITK object')
      self.mask = None

    self.resampledMask = resampledMask

    self.kwargs = kwargs
    self.inputImages = inputImages

    self.label = self.kwargs.get('label', 1)

    if resampledMask is not None:
      self.lssif = sitk.LabelShapeStatisticsImageFilter()
      self.lssif.Execute(resampledMask)
    else:
      self.lssif = None

  def execute(self):
    """
    Calculate and return a dictionary containing all general info items. Format is <info_item>:<value>, where any ',' in
    <value> are replaced by ';' to prevent column alignment errors in csv formatted output.
    """
    generalInfo = collections.OrderedDict()
    for mem in dir(self):
      if mem.startswith('get') and mem.endswith('Value'):
        generalInfo[mem[3:-5]] = str(eval('self.%s()' % (mem))).replace(',', ';')
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
    Return a string representation of the settings contained in kwargs.
    Format is {<settings_name>:<value>, ...}.
    """
    return self.kwargs

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
    Returns the original spacing of the image.

    If the reading of the image fails, an empty string is returned.
    """
    if self.image is not None:
      return self.image.GetSpacing()
    else:
      return None

  def getInputImagesValue(self):
    """
    Return a string representation of the enabled filters and any custom settings for the filter.
    Format is {<filter_name>:{<setting_name>:<value>, ...}, ...}.
    """
    return self.inputImages

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

  def getVersionValue(self):
    """
    Return the current version of this package.
    """
    return radiomics.__version__

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
