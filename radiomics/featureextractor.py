# -*- coding: utf-8 -*-
from __future__ import print_function

import collections
from itertools import chain
import logging
import os

import pykwalify.core
import SimpleITK as sitk
import six

import radiomics
from radiomics import generalinfo, getFeatureClasses, getInputImageTypes, imageoperations


class RadiomicsFeaturesExtractor:
  r"""
  Wrapper class for calculation of a radiomics signature.
  At and after initialisation various settings can be used to customize the resultant signature.
  This includes which classes and features to use, as well as what should be done in terms of preprocessing the image
  and what images (original and/or filtered) should be used as input.

  Then a call to :py:func:`execute` generates the radiomics
  signature specified by these settings for the passed image and labelmap combination. This function can be called
  repeatedly in a batch process to calculate the radiomics signature for all image and labelmap combinations.

  It initialisation, a parameters file can be provided containing all necessary settings. This is done by passing the
  location of the file as the single argument in the initialization call, without specifying it as a keyword argument.
  If such a file location is provided, any additional kwargs are ignored.
  Alternatively, at initialisation, the following general settings can be specified in the parameter file or ``kwargs``,
  with default values in brackets:

  - enableCExtensions [True]: Boolean, set to False to force calculation to full-python mode. See also
    :py:func:`~radiomics.enableCExtensions()`.
  - additionalInfo [True]: boolean, set to False to disable inclusion of additional information on the extraction in the
    output. See also :py:func:`~addProvenance()`.
  - binWidth [25]: Float, size of the bins when making a histogram and for discretization of the image gray level.
  - normalize [False]: Boolean, set to True to enable normalizing of the image before any resampling. See also
    :py:func:`~imageoperations.normalizeImage`.
  - normalizeScale [1]: Float, determines the scale after normalizing the image. If normalizing is disabled, this has
    no effect.
  - removeOutliers [None]: Float, defines the outliers to remove from the image. An outlier is defined as values that
    differ more than :math:`n\sigma_x` from the mean, where :math:`n>0` and equal to the value of this setting. If this
    parameter is omitted (providing it without a value (i.e. None) in the parameter file will throw an error), no
    outliers are removed. If normalizing is disabled, this has no effect. See also
    :py:func:`~imageoperations.normalizeImage`.
  - resampledPixelSpacing [None]: List of 3 floats, sets the size of the voxel in (x, y, z) plane when resampling.
  - interpolator [sitkBSpline]: Simple ITK constant or string name thereof, sets interpolator to use for resampling.
    Enumerated value, possible values:

    - sitkNearestNeighbor (= 1)
    - sitkLinear (= 2)
    - sitkBSpline (= 3)
    - sitkGaussian (= 4)
    - sitkLabelGaussian (= 5)
    - sitkHammingWindowedSinc (= 6)
    - sitkCosineWindowedSinc (= 7)
    - sitkWelchWindowedSinc (= 8)
    - sitkLanczosWindowedSinc (= 9)
    - sitkBlackmanWindowedSinc (= 10)

  - padDistance [5]: Integer, set the number of voxels pad cropped tumor volume with during resampling. Padding occurs
    in new feature space and is done on all faces, i.e. size increases in x, y and z direction by 2*padDistance.
    Padding is needed for some filters (e.g. LoG). Value of padded voxels are set to original gray level intensity,
    padding does not exceed original image boundaries. **N.B. After application of filters image is cropped again
    without padding.**

  .. note::

    Resampling is disabled when either `resampledPixelSpacing` or `interpolator` is set to `None`

  In addition to these general settings, filter or feature class specific settings can be defined here also.
  For more information on possible settings, see the respective filters and feature classes.

  By default, all features in all feature classes are enabled.
  By default, only `Original` input image is enabled (No filter applied).
  """

  def __init__(self, *args, **kwargs):
    self.logger = logging.getLogger(__name__)

    self.featureClasses = getFeatureClasses()

    self.kwargs = {}
    self.inputImages = {}
    self.enabledFeatures = {}

    if len(args) == 1 and isinstance(args[0], six.string_types):
      self.logger.info("Loading parameter file")
      self.loadParams(args[0])
    else:
      # Set default settings and update with and changed settings contained in kwargs
      self.kwargs = self._getDefaultSettings()
      if len(kwargs) > 0:
        self.logger.info('Applying custom settings')
        self.kwargs.update(kwargs)
      else:
        self.logger.info('No customized settings, applying defaults')

      self.logger.debug("Settings: %s", self.kwargs)

      self.inputImages = {'Original': {}}

      self.enabledFeatures = {}
      for featureClassName in self.getFeatureClassNames():
        self.enabledFeatures[featureClassName] = []

  @classmethod
  def _getDefaultSettings(cls):
    """
    Returns a dictionary containg the default settings specified in this class. These settings cover global toolbox
    settings, such as ``enableCExtensions``, as well as the image pre-processing settings (e.g. resampling). Feature
    class specific are defined in the respective feature classes and and not included here. Similarly, filter specific
    settings are defined in ``imageoperations.py`` and also not included here.
    """
    return {'normalize': False,
            'normalizeScale': 1,
            'removeOutliers': None,
            'resampledPixelSpacing': None,  # No resampling by default
            'interpolator': 'sitkBSpline',  # Alternative: sitk.sitkBSpline,
            'padDistance': 5,
            'label': 1,
            'enableCExtensions': True,
            'additionalInfo': True}

  def addProvenance(self, provenance_on=True):
    """
    Enable or disable reporting of additional information on the extraction. This information includes toolbox version,
    enabled input images and applied settings. Furthermore, additional information on the image and region of interest
    (ROI) is also provided, including original image spacing, total number of voxels in the ROI and total number of
    fully connected volumes in the ROI.

    To disable this, call ``addProvenance(False)``.
    """
    self.kwargs['additionalInfo'] = provenance_on

  def loadParams(self, paramsFile):
    """
    Parse specified parameters file and use it to update settings in kwargs, enabled feature(Classes) and input
    images:

    - settings not specified in parameters are set to their default value.
    - enabledFeatures are replaced by those in parameters. If no featureClass parameters were specified, all
      featureClasses and features are enabled.
    - inputImages are replaced by those in parameters. If no inputImage parameters were specified, only original
      image is used for feature extraction, with no additional custom settings

    The paramsFile is written according to the YAML-convention (www.yaml.org) and is checked by the code for
    consistency. Only one yaml document per file is allowed. Settings must be grouped by setting type as mentioned
    above are reflected in the structure of the document as follows::

        <Setting Type>:
          <Setting Name>: <value>
          ...
        <Setting Type>:
          ...

    Blank lines may be inserted to increase readability, the are ignored by the parser. Additional comments are also
    possible, these are preceded by an '#' and can be inserted on a blank line, or on a line containing settings::

        # This is a line containing only comments
        setting: # This is a comment placed after the declaration of the 'setting' group.

    Any keyword, such as a setting type or setting name may only be mentioned once. Multiple instances do not raise
    an error, but only the last encountered one is used.

    The three setting types are named as follows:

    - setting: Setting to use for preprocessing and class specific settings (``kwargs`` arguments). if no <value>
      is specified, the value for this setting is set to None.
    - featureClass: Feature class to enable, <value> is list of strings representing enabled features. If no
      <value> is specified or <value> is an empty list ('[]'), all features for this class are enabled.
    - inputImage: input image to calculate features on. <value> is custom kwarg settings (dictionary). if <value>
      is an empty dictionary ('{}'), no custom settings are added for this input image.

    If supplied params file does not match the requirements, a pykwalify error is raised.
    """
    dataDir = os.path.abspath(os.path.join(radiomics.__path__[0], 'schemas'))
    schemaFile = os.path.join(dataDir, 'paramSchema.yaml')
    schemaFuncs = os.path.join(dataDir, 'schemaFuncs.py')
    c = pykwalify.core.Core(source_file=paramsFile, schema_files=[schemaFile], extensions=[schemaFuncs])
    params = c.validate()

    inputImages = params.get('inputImage', {})
    enabledFeatures = params.get('featureClass', {})
    kwargs = params.get('setting', {})

    self.logger.debug("Parameter file parsed. Applying settings")

    if len(inputImages) == 0:
      self.inputImages = {'Original': {}}
    else:
      self.inputImages = inputImages

    self.logger.debug("Enabled input images: %s", self.inputImages)

    if len(enabledFeatures) == 0:
      self.enabledFeatures = {}
      for featureClassName in self.getFeatureClassNames():
        self.enabledFeatures[featureClassName] = []
    else:
      self.enabledFeatures = enabledFeatures

    self.logger.debug("Enabled features: %s", enabledFeatures)

    # Set default settings and update with and changed settings contained in kwargs
    self.kwargs = self._getDefaultSettings()
    self.kwargs.update(kwargs)

    self.logger.debug("Settings: %s", kwargs)

  def enableAllInputImages(self):
    """
    Enable all possible input images without any custom settings.
    """
    self.logger.debug('Enabling all input image types')
    for imageType in getInputImageTypes():
      self.inputImages[imageType] = {}
    self.logger.debug('Enabled input images types: %s', self.inputImages)

  def disableAllInputImages(self):
    """
    Disable all input images.
    """
    self.logger.debug('Disabling all input image types')
    self.inputImages = {}

  def enableInputImageByName(self, inputImage, enabled=True, customArgs=None):
    r"""
    Enable or disable specified input image. If enabling input image, optional custom settings can be specified in
    customArgs.

    Current possible input images are:

    - Original: No filter applied
    - Wavelet: Wavelet filtering, yields 8 decompositions per level (all possible combinations of applying either
      a High or a Low pass filter in each of the three dimensions.
      See also :py:func:`~radiomics.imageoperations.getWaveletImage`
    - LoG: Laplacian of Gaussian filter, edge enhancement filter. Emphasizes areas of gray level change, where sigma
      defines how coarse the emphasised texture should be. A low sigma emphasis on fine textures (change over a
      short distance), where a high sigma value emphasises coarse textures (gray level change over a large distance).
      See also :py:func:`~radiomics.imageoperations.getLoGImage`
    - Square: Takes the square of the image intensities and linearly scales them back to the original range.
      Negative values in the original image will be made negative again after application of filter.
    - SquareRoot: Takes the square root of the absolute image intensities and scales them back to original range.
      Negative values in the original image will be made negative again after application of filter.
    - Logarithm: Takes the logarithm of the absolute intensity + 1. Values are scaled to original range and
      negative original values are made negative again after application of filter.
    - Exponential: Takes the the exponential, where filtered intensity is e^(absolute intensity). Values are
      scaled to original range and negative original values are made negative again after application of filter.

    For the mathmetical formulas of square, squareroot, logarithm and exponential, see their respective functions in
    :ref:`imageoperations<radiomics-imageoperations-label>`
    (:py:func:`~radiomics.imageoperations.getSquareImage`,
    :py:func:`~radiomics.imageoperations.getSquareRootImage`,
    :py:func:`~radiomics.imageoperations.getLogarithmImage` and
    :py:func:`~radiomics.imageoperations.getExponentialImage`,
    respectively).
    """
    if inputImage not in getInputImageTypes():
      self.logger.warning('Input image type %s is not recognized', inputImage)
      return

    if enabled:
      if customArgs is None:
        customArgs = {}
        self.logger.debug('Enabling input image type %s (no additional custom settings)', inputImage)
      else:
        self.logger.debug('Enabling input image type %s (additional custom settings: %s)', inputImage, customArgs)
      self.inputImages[inputImage] = customArgs
    elif inputImage in self.inputImages:
      self.logger.debug('Disabling input image type %s', inputImage)
      del self.inputImages[inputImage]
    self.logger.debug('Enabled input images types: %s', self.inputImages)

  def enableInputImages(self, **inputImages):
    """
    Enable input images, with optionally custom settings, which are applied to the respective input image.
    Settings specified here override those in kwargs.
    The following settings are not customizable:

    - interpolator
    - resampledPixelSpacing
    - padDistance

    Updates current settings: If necessary, enables input image. Always overrides custom settings specified
    for input images passed in inputImages.
    To disable input images, use :py:func:`enableInputImageByName` or :py:func:`disableAllInputImages`
    instead.

    :param inputImages: dictionary, key is imagetype (original, wavelet or log) and value is custom settings (dictionary)
    """
    self.logger.debug('Updating enabled input images types with %s', inputImages)
    self.inputImages.update(inputImages)
    self.logger.debug('Enabled input images types: %s', self.inputImages)

  def enableAllFeatures(self):
    """
    Enable all classes and all features.
    """
    self.logger.debug('Enabling all features in all feature classes')
    for featureClassName in self.getFeatureClassNames():
      self.enabledFeatures[featureClassName] = []
    self.logger.debug('Enabled features: %s', self.enabledFeatures)

  def disableAllFeatures(self):
    """
    Disable all classes.
    """
    self.logger.debug('Disabling all feature classes')
    self.enabledFeatures = {}

  def enableFeatureClassByName(self, featureClass, enabled=True):
    """
    Enable or disable all features in given class.
    """
    if featureClass not in self.getFeatureClassNames():
      self.logger.warning('Feature class %s is not recognized', featureClass)
      return

    if enabled:
      self.logger.debug('Enabling all features in class %s', featureClass)
      self.enabledFeatures[featureClass] = []
    elif featureClass in self.enabledFeatures:
      self.logger.debug('Disabling feature class %s', featureClass)
      del self.enabledFeatures[featureClass]
    self.logger.debug('Enabled features: %s', self.enabledFeatures)

  def enableFeaturesByName(self, **enabledFeatures):
    """
    Specify which features to enable. Key is feature class name, value is a list of enabled feature names.

    To enable all features for a class, provide the class name with an empty list or None as value.
    Settings for feature classes specified in enabledFeatures.keys are updated, settings for feature classes
    not yet present in enabledFeatures.keys are added.
    To disable the entire class, use :py:func:`disableAllFeatures` or :py:func:`enableFeatureClassByName` instead.
    """
    self.logger.debug('Updating enabled features with %s', enabledFeatures)
    self.enabledFeatures.update(enabledFeatures)
    self.logger.debug('Enabled features: %s', self.enabledFeatures)

  def execute(self, imageFilepath, maskFilepath, label=None):
    """
    Compute radiomics signature for provide image and mask combination.
    First, image and mask are loaded and normalized/resampled if necessary. Second, if enabled, provenance information
    is calculated and stored as part of the result. Next, shape features are calculated on a cropped (no padding)
    version of the original image. Then other featureclasses are calculated using all specified input images in
    ``inputImages``. Images are cropped to tumor mask (no padding) after application of any filter and before being
    passed to the feature class. Finally, the dictionary containing all calculated features is returned.

    :param imageFilepath: SimpleITK Image, or string pointing to image file location
    :param maskFilepath: SimpleITK Image, or string pointing to labelmap file location
    :param label: Integer, value of the label for which to extract features. If not specified, last specified label
        is used. Default label is 1.
    :returns: dictionary containing calculated signature ("<filter>_<featureClass>_<featureName>":value).
    """
    # Enable or disable C extensions for high performance matrix calculation. Only logs a message (INFO) when setting is
    # successfully changed. If an error occurs, full-python mode is forced and a warning is logged.
    radiomics.enableCExtensions(self.kwargs['enableCExtensions'])

    if label is not None:
      self.kwargs.update({'label': label})

    self.logger.info('Calculating features with label: %d', self.kwargs['label'])
    self.logger.debug('Enabled input images types: %s', self.inputImages)
    self.logger.debug('Enabled features: %s', self.enabledFeatures)
    self.logger.debug('Current settings: %s', self.kwargs)

    featureVector = collections.OrderedDict()
    image, mask = self.loadImage(imageFilepath, maskFilepath)

    if image is not None and mask is not None:
      if self.kwargs['additionalInfo']:
        featureVector.update(self.getProvenance(imageFilepath, maskFilepath, mask))

      # Bounding box only needs to be calculated once after resampling, store the value, so it doesn't get calculated
      # after every filter
      boundingBox = None

      # If shape should be calculation, handle it separately here
      if 'shape' in self.enabledFeatures.keys():
        croppedImage, croppedMask, boundingBox = \
          imageoperations.cropToTumorMask(image, mask, self.kwargs['label'], boundingBox)
        enabledFeatures = self.enabledFeatures['shape']

        self.logger.info("Computing shape")
        shapeClass = self.featureClasses['shape'](croppedImage, croppedMask, **self.kwargs)
        if enabledFeatures is None or len(enabledFeatures) == 0:
          shapeClass.enableAllFeatures()
        else:
          for feature in enabledFeatures:
            shapeClass.enableFeatureByName(feature)

        shapeClass.calculateFeatures()
        for (featureName, featureValue) in six.iteritems(shapeClass.featureValues):
          newFeatureName = "original_shape_%s" % (featureName)
          featureVector[newFeatureName] = featureValue

      # Make generators for all enabled input image types
      self.logger.debug('Creating input image type iterator')
      imageGenerators = []
      for imageType, customKwargs in six.iteritems(self.inputImages):
        args = self.kwargs.copy()
        args.update(customKwargs)
        self.logger.info("Applying filter: '%s' with settings: %s" % (imageType, str(args)))
        imageGenerators = chain(imageGenerators, getattr(imageoperations, 'get%sImage' % imageType)(image, **args))

      self.logger.debug('Extracting features')
      # Calculate features for all (filtered) images in the generator
      for inputImage, inputImageName, inputKwargs in imageGenerators:
        self.logger.info('Calculating features for %s image, with settings: %s', inputImageName, str(inputKwargs))
        inputImage, inputMask, boundingBox = imageoperations.cropToTumorMask(inputImage,
                                                                             mask,
                                                                             self.kwargs['label'],
                                                                             boundingBox)
        featureVector.update(self.computeFeatures(inputImage, inputMask, inputImageName, **inputKwargs))

      self.logger.debug('Features extracted')

    return featureVector

  def loadImage(self, ImageFilePath, MaskFilePath):
    """
    Preprocess the image and labelmap.
    If ImageFilePath is a string, it is loaded as SimpleITK Image and assigned to image,
    if it already is a SimpleITK Image, it is just assigned to image.
    All other cases are ignored (nothing calculated).
    Equal approach is used for assignment of mask using MaskFilePath.

    If normalizing is enabled image is first normalized before any resampling is applied.

    If resampling is enabled, both image and mask are resampled and cropped to the tumor mask (with additional
    padding as specified in padDistance) after assignment of image and mask.
    """
    self.logger.info('Loading image and mask')
    if isinstance(ImageFilePath, six.string_types) and os.path.exists(ImageFilePath):
      image = sitk.ReadImage(ImageFilePath)
    elif isinstance(ImageFilePath, sitk.SimpleITK.Image):
      image = ImageFilePath
    else:
      self.logger.warning('Error reading image Filepath or SimpleITK object')
      return None, None  # this function is expected to always return a tuple of 2 elements

    if isinstance(MaskFilePath, six.string_types) and os.path.exists(MaskFilePath):
      mask = sitk.ReadImage(MaskFilePath)
    elif isinstance(MaskFilePath, sitk.SimpleITK.Image):
      mask = MaskFilePath
    else:
      self.logger.warning('Error reading mask Filepath or SimpleITK object')
      return None, None  # this function is expected to always return a tuple of 2 elements

    # This point is only reached if image and mask loaded correctly
    if self.kwargs['normalize']:
      image = imageoperations.normalizeImage(image, self.kwargs['normalizeScale'], self.kwargs['removeOutliers'])

    if self.kwargs['interpolator'] is not None and self.kwargs['resampledPixelSpacing'] is not None:
      image, mask = imageoperations.resampleImage(image, mask,
                                                  self.kwargs['resampledPixelSpacing'],
                                                  self.kwargs['interpolator'],
                                                  self.kwargs['label'],
                                                  self.kwargs['padDistance'])

    return image, mask

  def getProvenance(self, imageFilepath, maskFilepath, mask):
    """
    Generates provenance information for reproducibility. Takes the original image & mask filepath, as well as the
    resampled mask which is passed to the feature classes. Returns a dictionary with keynames coded as
    "general_info_<item>". For more information on generated items, see :ref:`generalinfo<radiomics-generalinfo-label>`
    """
    self.logger.info('Adding additional extraction information')

    provenanceVector = collections.OrderedDict()
    generalinfoClass = generalinfo.GeneralInfo(imageFilepath, maskFilepath, mask, self.kwargs, self.inputImages)
    for k, v in six.iteritems(generalinfoClass.execute()):
      provenanceVector['general_info_%s' % (k)] = v
    return provenanceVector

  def computeFeatures(self, image, mask, inputImageName, **kwargs):
    """
    Compute signature using image, mask, \*\*kwargs settings.

    This function computes the signature for just the passed image (original or derived), it does not preprocess or
    apply a filter to the passed image. Features / Classes to use for calculation of signature are defined in
    self.enabledFeatures. See also :py:func:`enableFeaturesByName`.

    .. note::

      shape descriptors are independent of gray level and therefore calculated separately (handeled in `execute`). In
      this function, no shape functions are calculated.
    """
    featureVector = collections.OrderedDict()

    # Calculate feature classes
    for featureClassName, enabledFeatures in six.iteritems(self.enabledFeatures):
      # Handle calculation of shape features separately
      if featureClassName == 'shape':
        continue

      if featureClassName in self.getFeatureClassNames():
        self.logger.info('Computing %s', featureClassName)

        featureClass = self.featureClasses[featureClassName](image, mask, **kwargs)

        if enabledFeatures is None or len(enabledFeatures) == 0:
          featureClass.enableAllFeatures()
        else:
          for feature in enabledFeatures:
            featureClass.enableFeatureByName(feature)

        featureClass.calculateFeatures()
        for (featureName, featureValue) in six.iteritems(featureClass.featureValues):
          newFeatureName = "%s_%s_%s" % (inputImageName, featureClassName, featureName)
          featureVector[newFeatureName] = featureValue

    return featureVector

  def getFeatureClassNames(self):
    """
    Returns a list of all possible feature classes.
    """
    return self.featureClasses.keys()

  def getFeatureNames(self, featureClassName):
    """
    Returns a list of all possible features in provided featureClass
    """
    return self.featureClasses[featureClassName].getFeatureNames()
