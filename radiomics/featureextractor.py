# -*- coding: utf-8 -*-
from __future__ import print_function

import collections
from itertools import chain
import json
import logging
import os
import pathlib

import pykwalify.core
import SimpleITK as sitk
import six

from radiomics import generalinfo, getFeatureClasses, getImageTypes, getParameterValidationFiles, imageoperations


logger = logging.getLogger(__name__)
geometryTolerance = None


class RadiomicsFeatureExtractor:
  r"""
  Wrapper class for calculation of a radiomics signature.
  At and after initialisation various settings can be used to customize the resultant signature.
  This includes which classes and features to use, as well as what should be done in terms of preprocessing the image
  and what images (original and/or filtered) should be used as input.

  Then a call to :py:func:`execute` generates the radiomics
  signature specified by these settings for the passed image and labelmap combination. This function can be called
  repeatedly in a batch process to calculate the radiomics signature for all image and labelmap combinations.

  At initialization, a parameters file (string pointing to yaml or json structured file) or dictionary can be provided
  containing all necessary settings (top level containing keys "setting", "imageType" and/or "featureClass). This is
  done by passing it as the first positional argument. If no positional argument is supplied, or the argument is not
  either a dictionary or a string pointing to a valid file, defaults will be applied.
  Moreover, at initialisation, custom settings (*NOT enabled image types and/or feature classes*) can be provided
  as keyword arguments, with the setting name as key and its value as the argument value (e.g. ``binWidth=25``).
  Settings specified here will override those in the parameter file/dict/default settings.
  For more information on possible settings and customization, see
  :ref:`Customizing the Extraction <radiomics-customization-label>`.

  By default, all features in all feature classes are enabled.
  By default, only `Original` input image is enabled (No filter applied).
  """

  def __init__(self, *args, **kwargs):
    global logger

    self.settings = {}
    self.enabledImagetypes = {}
    self.enabledFeatures = {}

    self.featureClassNames = list(getFeatureClasses().keys())

    if len(args) == 1 and isinstance(args[0], dict):
      logger.info("Loading parameter dictionary")
      self._applyParams(paramsDict=args[0])
    elif len(args) == 1 and (isinstance(args[0], six.string_types) or isinstance(args[0], pathlib.PurePath)):
      if not os.path.isfile(args[0]):
        raise IOError("Parameter file %s does not exist." % args[0])
      logger.info("Loading parameter file %s", str(args[0]))
      self._applyParams(paramsFile=args[0])
    else:
      # Set default settings and update with and changed settings contained in kwargs
      self.settings = self._getDefaultSettings()
      logger.info('No valid config parameter, using defaults: %s', self.settings)

      self.enabledImagetypes = {'Original': {}}
      logger.info('Enabled image types: %s', self.enabledImagetypes)

      for featureClassName in self.featureClassNames:
        if featureClassName == 'shape2D':  # Do not enable shape2D by default
          continue
        self.enabledFeatures[featureClassName] = []
      logger.info('Enabled features: %s', self.enabledFeatures)

    if len(kwargs) > 0:
      logger.info('Applying custom setting overrides: %s', kwargs)
      self.settings.update(kwargs)
      logger.debug("Settings: %s", self.settings)

    if self.settings.get('binCount', None) is not None:
      logger.warning('Fixed bin Count enabled! However, we recommend using a fixed bin Width. See '
                     'http://pyradiomics.readthedocs.io/en/latest/faq.html#radiomics-fixed-bin-width for more '
                     'details')

    self._setTolerance()

  def _setTolerance(self):
    global geometryTolerance, logger
    geometryTolerance = self.settings.get('geometryTolerance')
    if geometryTolerance is not None:
      logger.debug('Setting SimpleITK tolerance to %s', geometryTolerance)
      sitk.ProcessObject.SetGlobalDefaultCoordinateTolerance(geometryTolerance)
      sitk.ProcessObject.SetGlobalDefaultDirectionTolerance(geometryTolerance)

  def addProvenance(self, provenance_on=True):
    """
    Enable or disable reporting of additional information on the extraction. This information includes toolbox version,
    enabled input images and applied settings. Furthermore, additional information on the image and region of interest
    (ROI) is also provided, including original image spacing, total number of voxels in the ROI and total number of
    fully connected volumes in the ROI.

    To disable this, call ``addProvenance(False)``.
    """
    self.settings['additionalInfo'] = provenance_on

  @staticmethod
  def _getDefaultSettings():
    """
    Returns a dictionary containg the default settings specified in this class. These settings cover global settings,
    such as ``additionalInfo``, as well as the image pre-processing settings (e.g. resampling). Feature class specific
    are defined in the respective feature classes and and not included here. Similarly, filter specific settings are
    defined in ``imageoperations.py`` and also not included here.
    """
    return {'minimumROIDimensions': 2,
            'minimumROISize': None,  # Skip testing the ROI size by default
            'normalize': False,
            'normalizeScale': 1,
            'removeOutliers': None,
            'resampledPixelSpacing': None,  # No resampling by default
            'interpolator': 'sitkBSpline',  # Alternative: sitk.sitkBSpline
            'preCrop': False,
            'padDistance': 5,
            'distances': [1],
            'force2D': False,
            'force2Ddimension': 0,
            'resegmentRange': None,  # No resegmentation by default
            'label': 1,
            'additionalInfo': True}

  def loadParams(self, paramsFile):
    """
    Parse specified parameters file and use it to update settings, enabled feature(Classes) and image types. For more
    information on the structure of the parameter file, see
    :ref:`Customizing the extraction <radiomics-customization-label>`.

    If supplied file does not match the requirements (i.e. unrecognized names or invalid values for a setting), a
    pykwalify error is raised.
    """
    self._applyParams(paramsFile=paramsFile)

  def loadJSONParams(self, JSON_configuration):
    """
    Pars JSON structured configuration string and use it to update settings, enabled feature(Classes) and image types.
    For more information on the structure of the parameter file, see
    :ref:`Customizing the extraction <radiomics-customization-label>`.

    If supplied string does not match the requirements (i.e. unrecognized names or invalid values for a setting), a
    pykwalify error is raised.
    """
    parameter_data = json.loads(JSON_configuration)
    self._applyParams(paramsDict=parameter_data)

  def _applyParams(self, paramsFile=None, paramsDict=None):
    """
    Validates and applies a parameter dictionary. See :py:func:`loadParams` and :py:func:`loadJSONParams` for more info.
    """
    global logger

    # Ensure pykwalify.core has a log handler (needed when parameter validation fails)
    if len(pykwalify.core.log.handlers) == 0 and len(logging.getLogger().handlers) == 0:
      # No handler available for either pykwalify or root logger, provide first radiomics handler (outputs to stderr)
      pykwalify.core.log.addHandler(logging.getLogger('radiomics').handlers[0])

    schemaFile, schemaFuncs = getParameterValidationFiles()
    c = pykwalify.core.Core(source_file=paramsFile, source_data=paramsDict,
                            schema_files=[schemaFile], extensions=[schemaFuncs])
    params = c.validate()
    logger.debug('Parameters parsed, input is valid.')

    enabledImageTypes = params.get('imageType', {})
    enabledFeatures = params.get('featureClass', {})
    settings = params.get('setting', {})
    voxelSettings = params.get('voxelSetting', {})

    logger.debug("Applying settings")

    if len(enabledImageTypes) == 0:
      self.enabledImagetypes = {'Original': {}}
    else:
      self.enabledImagetypes = enabledImageTypes

    logger.debug("Enabled image types: %s", self.enabledImagetypes)

    if len(enabledFeatures) == 0:
      self.enabledFeatures = {}
      for featureClassName in self.featureClassNames:
        self.enabledFeatures[featureClassName] = []
    else:
      self.enabledFeatures = enabledFeatures

    logger.debug("Enabled features: %s", self.enabledFeatures)

    # Set default settings and update with and changed settings contained in kwargs
    self.settings = self._getDefaultSettings()
    self.settings.update(settings)
    self.settings.update(voxelSettings)

    logger.debug("Settings: %s", settings)

  def execute(self, imageFilepath, maskFilepath, label=None, label_channel=None, voxelBased=False):
    """
    Compute radiomics signature for provide image and mask combination. It comprises of the following steps:

    1. Image and mask are loaded and normalized/resampled if necessary.
    2. Validity of ROI is checked using :py:func:`~imageoperations.checkMask`, which also computes and returns the
       bounding box.
    3. If enabled, provenance information is calculated and stored as part of the result. (Not available in voxel-based
       extraction)
    4. Shape features are calculated on a cropped (no padding) version of the original image. (Not available in
       voxel-based extraction)
    5. If enabled, resegment the mask based upon the range specified in ``resegmentRange`` (default None: resegmentation
       disabled).
    6. Other enabled feature classes are calculated using all specified image types in ``_enabledImageTypes``. Images
       are cropped to tumor mask (no padding) after application of any filter and before being passed to the feature
       class.
    7. The calculated features is returned as ``collections.OrderedDict``.

    :param imageFilepath: SimpleITK Image, or string pointing to image file location
    :param maskFilepath: SimpleITK Image, or string pointing to labelmap file location
    :param label: Integer, value of the label for which to extract features. If not specified, last specified label
        is used. Default label is 1.
    :param label_channel: Integer, index of the channel to use when maskFilepath yields a SimpleITK.Image with a vector
        pixel type. Default index is 0.
    :param voxelBased: Boolean, default False. If set to true, a voxel-based extraction is performed, segment-based
        otherwise.
    :returns: dictionary containing calculated signature ("<imageType>_<featureClass>_<featureName>":value).
        In case of segment-based extraction, value type for features is float, if voxel-based, type is SimpleITK.Image.
        Type of diagnostic features differs, but can always be represented as a string.
    """
    global geometryTolerance, logger
    _settings = self.settings.copy()

    tolerance = _settings.get('geometryTolerance')
    additionalInfo = _settings.get('additionalInfo', False)
    resegmentShape = _settings.get('resegmentShape', False)

    if label is not None:
      _settings['label'] = label
    else:
      label = _settings.get('label', 1)

    if label_channel is not None:
      _settings['label_channel'] = label_channel

    if geometryTolerance != tolerance:
      self._setTolerance()

    if additionalInfo:
      generalInfo = generalinfo.GeneralInfo()
      generalInfo.addGeneralSettings(_settings)
      generalInfo.addEnabledImageTypes(self.enabledImagetypes)
    else:
      generalInfo = None

    if voxelBased:
      _settings['voxelBased'] = True
      kernelRadius = _settings.get('kernelRadius', 1)
      logger.info('Starting voxel based extraction')
    else:
      kernelRadius = 0

    logger.info('Calculating features with label: %d', label)
    logger.debug('Enabled images types: %s', self.enabledImagetypes)
    logger.debug('Enabled features: %s', self.enabledFeatures)
    logger.debug('Current settings: %s', _settings)

    # 1. Load the image and mask
    featureVector = collections.OrderedDict()
    image, mask = self.loadImage(imageFilepath, maskFilepath, generalInfo, **_settings)

    # 2. Check whether loaded mask contains a valid ROI for feature extraction and get bounding box
    # Raises a ValueError if the ROI is invalid
    boundingBox, correctedMask = imageoperations.checkMask(image, mask, **_settings)

    # Update the mask if it had to be resampled
    if correctedMask is not None:
      if generalInfo is not None:
        generalInfo.addMaskElements(image, correctedMask, label, 'corrected')
      mask = correctedMask

    logger.debug('Image and Mask loaded and valid, starting extraction')

    # 5. Resegment the mask if enabled (parameter regsegmentMask is not None)
    resegmentedMask = None
    if _settings.get('resegmentRange', None) is not None:
      resegmentedMask = imageoperations.resegmentMask(image, mask, **_settings)

      # Recheck to see if the mask is still valid, raises a ValueError if not
      boundingBox, correctedMask = imageoperations.checkMask(image, resegmentedMask, **_settings)

      if generalInfo is not None:
        generalInfo.addMaskElements(image, resegmentedMask, label, 'resegmented')

    # 3. Add the additional information if enabled
    if generalInfo is not None:
      featureVector.update(generalInfo.getGeneralInfo())

    # if resegmentShape is True and resegmentation has been enabled, update the mask here to also use the
    # resegmented mask for shape calculation (e.g. PET resegmentation)
    if resegmentShape and resegmentedMask is not None:
      mask = resegmentedMask

    if not voxelBased:
      # 4. If shape descriptors should be calculated, handle it separately here
      featureVector.update(self.computeShape(image, mask, boundingBox, **_settings))

    # (Default) Only use resegemented mask for feature classes other than shape
    # can be overridden by specifying `resegmentShape` = True
    if not resegmentShape and resegmentedMask is not None:
      mask = resegmentedMask

    # 6. Calculate other enabled feature classes using enabled image types
    # Make generators for all enabled image types
    logger.debug('Creating image type iterator')
    imageGenerators = []
    for imageType, customKwargs in six.iteritems(self.enabledImagetypes):
      args = _settings.copy()
      args.update(customKwargs)
      logger.info('Adding image type "%s" with custom settings: %s' % (imageType, str(customKwargs)))
      imageGenerators = chain(imageGenerators, getattr(imageoperations, 'get%sImage' % imageType)(image, mask, **args))

    logger.debug('Extracting features')
    # Calculate features for all (filtered) images in the generator
    for inputImage, imageTypeName, inputKwargs in imageGenerators:
      logger.info('Calculating features for %s image', imageTypeName)
      inputImage, inputMask = imageoperations.cropToTumorMask(inputImage, mask, boundingBox, padDistance=kernelRadius)
      featureVector.update(self.computeFeatures(inputImage, inputMask, imageTypeName, **inputKwargs))

    logger.debug('Features extracted')

    return featureVector

  @staticmethod
  def loadImage(ImageFilePath, MaskFilePath, generalInfo=None, **kwargs):
    """
    Load and pre-process the image and labelmap.
    If ImageFilePath is a string, it is loaded as SimpleITK Image and assigned to ``image``,
    if it already is a SimpleITK Image, it is just assigned to ``image``.
    All other cases are ignored (nothing calculated).
    Equal approach is used for assignment of ``mask`` using MaskFilePath. If necessary, a segmentation object (i.e. mask
    volume with vector-image type) is then converted to a labelmap (=scalar image type). Data type is forced to UInt32.
    See also :py:func:`~imageoperations.getMask()`.

    If normalizing is enabled image is first normalized before any resampling is applied.

    If resampling is enabled, both image and mask are resampled and cropped to the tumor mask (with additional
    padding as specified in padDistance) after assignment of image and mask.

    :param ImageFilePath: SimpleITK.Image object or string pointing to SimpleITK readable file representing the image
                          to use.
    :param MaskFilePath: SimpleITK.Image object or string pointing to SimpleITK readable file representing the mask
                         to use.
    :param generalInfo: GeneralInfo Object. If provided, it is used to store diagnostic information of the
                        pre-processing.
    :param kwargs: Dictionary containing the settings to use for this particular image type.
    :return: 2 SimpleITK.Image objects representing the loaded image and mask, respectively.
    """
    global logger

    normalize = kwargs.get('normalize', False)
    interpolator = kwargs.get('interpolator')
    resampledPixelSpacing = kwargs.get('resampledPixelSpacing')
    preCrop = kwargs.get('preCrop', False)
    label = kwargs.get('label', 1)

    logger.info('Loading image and mask')
    if isinstance(ImageFilePath, six.string_types) and os.path.isfile(ImageFilePath):
      image = sitk.ReadImage(ImageFilePath)
    elif isinstance(ImageFilePath, sitk.SimpleITK.Image):
      image = ImageFilePath
    else:
      raise ValueError('Error reading image Filepath or SimpleITK object')

    if isinstance(MaskFilePath, six.string_types) and os.path.isfile(MaskFilePath):
      mask = sitk.ReadImage(MaskFilePath)
    elif isinstance(MaskFilePath, sitk.SimpleITK.Image):
      mask = MaskFilePath
    else:
      raise ValueError('Error reading mask Filepath or SimpleITK object')

    # process the mask
    mask = imageoperations.getMask(mask, **kwargs)

    if generalInfo is not None:
      generalInfo.addImageElements(image)
      # Do not include the image here, as the overlap between image and mask have not been checked
      # It is therefore possible that image and mask do not align, or even have different sizes.
      generalInfo.addMaskElements(None, mask, label)

    # This point is only reached if image and mask loaded correctly
    if normalize:
      image = imageoperations.normalizeImage(image, **kwargs)

    if interpolator is not None and resampledPixelSpacing is not None:
      image, mask = imageoperations.resampleImage(image, mask, **kwargs)
      if generalInfo is not None:
        generalInfo.addImageElements(image, 'interpolated')
        generalInfo.addMaskElements(image, mask, label, 'interpolated')

    elif preCrop:
      bb, correctedMask = imageoperations.checkMask(image, mask, **kwargs)
      if correctedMask is not None:
        # Update the mask if it had to be resampled
        mask = correctedMask
      if bb is None:
        # Mask checks failed
        raise ValueError('Mask checks failed during pre-crop')

      image, mask = imageoperations.cropToTumorMask(image, mask, bb, **kwargs)

    return image, mask

  def computeShape(self, image, mask, boundingBox, **kwargs):
    """
    Calculate the shape (2D and/or 3D) features for the passed image and mask.

    :param image: SimpleITK.Image object representing the image used
    :param mask: SimpleITK.Image object representing the mask used
    :param boundingBox: The boundingBox calculated by :py:func:`~imageoperations.checkMask()`, i.e. a tuple with lower
      (even indices) and upper (odd indices) bound of the bounding box for each dimension.
    :param kwargs: Dictionary containing the settings to use.
    :return: collections.OrderedDict containing the calculated shape features. If no features are calculated, an empty
      OrderedDict will be returned.
    """
    global logger
    featureVector = collections.OrderedDict()

    enabledFeatures = self.enabledFeatures

    croppedImage, croppedMask = imageoperations.cropToTumorMask(image, mask, boundingBox)

    # Define temporary function to compute shape features
    def compute(shape_type):
      logger.info('Computing %s', shape_type)
      featureNames = enabledFeatures[shape_type]
      shapeClass = getFeatureClasses()[shape_type](croppedImage, croppedMask, **kwargs)

      if featureNames is not None:
        for feature in featureNames:
          shapeClass.enableFeatureByName(feature)

      for (featureName, featureValue) in six.iteritems(shapeClass.execute()):
        newFeatureName = 'original_%s_%s' % (shape_type, featureName)
        featureVector[newFeatureName] = featureValue

    Nd = mask.GetDimension()
    if 'shape' in enabledFeatures.keys():
      if Nd == 3:
        compute('shape')
      else:
        logger.warning('Shape features are only available 3D input (for 2D input, use shape2D). Found %iD input',
                       Nd)

    if 'shape2D' in enabledFeatures.keys():
      if Nd == 3:
        force2D = kwargs.get('force2D', False)
        force2Ddimension = kwargs.get('force2Ddimension', 0)
        if not force2D:
          logger.warning('parameter force2D must be set to True to enable shape2D extraction')
        elif not (boundingBox[1::2] - boundingBox[0::2] + 1)[force2Ddimension] > 1:
          logger.warning('Size in specified 2D dimension (%i) is greater than 1, cannot calculate 2D shape',
                         force2Ddimension)
        else:
          compute('shape2D')
      elif Nd == 2:
        compute('shape2D')
      else:
        logger.warning('Shape2D features are only available for 2D and 3D (with force2D=True) input. '
                       'Found %iD input', Nd)

    return featureVector

  def computeFeatures(self, image, mask, imageTypeName, **kwargs):
    r"""
    Compute signature using image, mask and \*\*kwargs settings.

    This function computes the signature for just the passed image (original or derived), it does not pre-process or
    apply a filter to the passed image. Features / Classes to use for calculation of signature are defined in
    ``self.enabledFeatures``. See also :py:func:`enableFeaturesByName`.

    :param image: The cropped (and optionally filtered) SimpleITK.Image object representing the image used
    :param mask: The cropped SimpleITK.Image object representing the mask used
    :param imageTypeName: String specifying the filter applied to the image, or "original" if no filter was applied.
    :param kwargs: Dictionary containing the settings to use for this particular image type.
    :return: collections.OrderedDict containing the calculated features for all enabled classes.
      If no features are calculated, an empty OrderedDict will be returned.

    .. note::

      shape descriptors are independent of gray level and therefore calculated separately (handled in `execute`). In
      this function, no shape features are calculated.
    """
    global logger
    featureVector = collections.OrderedDict()
    featureClasses = getFeatureClasses()

    enabledFeatures = self.enabledFeatures

    # Calculate feature classes
    for featureClassName, featureNames in six.iteritems(enabledFeatures):
      # Handle calculation of shape features separately
      if featureClassName.startswith('shape'):
        continue

      if featureClassName in featureClasses:
        logger.info('Computing %s', featureClassName)

        featureClass = featureClasses[featureClassName](image, mask, **kwargs)

        if featureNames is not None:
          for feature in featureNames:
            featureClass.enableFeatureByName(feature)

        for (featureName, featureValue) in six.iteritems(featureClass.execute()):
          newFeatureName = '%s_%s_%s' % (imageTypeName, featureClassName, featureName)
          featureVector[newFeatureName] = featureValue

    return featureVector

  def enableAllImageTypes(self):
    """
    Enable all possible image types without any custom settings.
    """
    global logger

    logger.debug('Enabling all image types')
    for imageType in getImageTypes():
      self.enabledImagetypes[imageType] = {}
    logger.debug('Enabled images types: %s', self.enabledImagetypes)

  def disableAllImageTypes(self):
    """
    Disable all image types.
    """
    global logger

    logger.debug('Disabling all image types')
    self.enabledImagetypes = {}

  def enableImageTypeByName(self, imageType, enabled=True, customArgs=None):
    r"""
    Enable or disable specified image type. If enabling image type, optional custom settings can be specified in
    customArgs.

    Current possible image types are:

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
    - Gradient: Returns the gradient magnitude.
    - LBP2D: Calculates and returns a local binary pattern applied in 2D.
    - LBP3D: Calculates and returns local binary pattern maps applied in 3D using spherical harmonics. Last returned
      image is the corresponding kurtosis map.

    For the mathmetical formulas of square, squareroot, logarithm and exponential, see their respective functions in
    :ref:`imageoperations<radiomics-imageoperations-label>`
    (:py:func:`~radiomics.imageoperations.getSquareImage`,
    :py:func:`~radiomics.imageoperations.getSquareRootImage`,
    :py:func:`~radiomics.imageoperations.getLogarithmImage`,
    :py:func:`~radiomics.imageoperations.getExponentialImage`,
    :py:func:`~radiomics.imageoperations.getGradientImage`,
    :py:func:`~radiomics.imageoperations.getLBP2DImage` and
    :py:func:`~radiomics.imageoperations.getLBP3DImage`,
    respectively).
    """
    global logger

    if imageType not in getImageTypes():
      logger.warning('Image type %s is not recognized', imageType)
      return

    if enabled:
      if customArgs is None:
        customArgs = {}
        logger.debug('Enabling image type %s (no additional custom settings)', imageType)
      else:
        logger.debug('Enabling image type %s (additional custom settings: %s)', imageType, customArgs)
      self.enabledImagetypes[imageType] = customArgs
    elif imageType in self.enabledImagetypes:
      logger.debug('Disabling image type %s', imageType)
      del self.enabledImagetypes[imageType]
    logger.debug('Enabled images types: %s', self.enabledImagetypes)

  def enableImageTypes(self, **enabledImagetypes):
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

    :param enabledImagetypes: dictionary, key is imagetype (original, wavelet or log) and value is custom settings
      (dictionary)
    """
    global logger

    logger.debug('Updating enabled images types with %s', enabledImagetypes)
    self.enabledImagetypes.update(enabledImagetypes)
    logger.debug('Enabled images types: %s', self.enabledImagetypes)

  def enableAllFeatures(self):
    """
    Enable all classes and all features.

    .. note::
      Individual features that have been marked "deprecated" are not enabled by this function. They can still be enabled
      manually by a call to :py:func:`~radiomics.base.RadiomicsBase.enableFeatureByName()`,
      :py:func:`~radiomics.featureextractor.RadiomicsFeaturesExtractor.enableFeaturesByName()`
      or in the parameter file (by specifying the feature by name, not when enabling all features).
      However, in most cases this will still result only in a deprecation warning.
    """
    global logger

    logger.debug('Enabling all features in all feature classes')
    for featureClassName in self.featureClassNames:
      self.enabledFeatures[featureClassName] = []
    logger.debug('Enabled features: %s', self.enabledFeatures)

  def disableAllFeatures(self):
    """
    Disable all classes.
    """
    global logger

    logger.debug('Disabling all feature classes')
    self.enabledFeatures = {}

  def enableFeatureClassByName(self, featureClass, enabled=True):
    """
    Enable or disable all features in given class.

    .. note::
      Individual features that have been marked "deprecated" are not enabled by this function. They can still be enabled
      manually by a call to :py:func:`~radiomics.base.RadiomicsBase.enableFeatureByName()`,
      :py:func:`~radiomics.featureextractor.RadiomicsFeaturesExtractor.enableFeaturesByName()`
      or in the parameter file (by specifying the feature by name, not when enabling all features).
      However, in most cases this will still result only in a deprecation warning.
    """
    global logger

    if featureClass not in self.featureClassNames:
      logger.warning('Feature class %s is not recognized', featureClass)
      return

    if enabled:
      logger.debug('Enabling all features in class %s', featureClass)
      self.enabledFeatures[featureClass] = []
    elif featureClass in self.enabledFeatures:
      logger.debug('Disabling feature class %s', featureClass)
      del self.enabledFeatures[featureClass]
    logger.debug('Enabled features: %s', self.enabledFeatures)

  def enableFeaturesByName(self, **enabledFeatures):
    """
    Specify which features to enable. Key is feature class name, value is a list of enabled feature names.

    To enable all features for a class, provide the class name with an empty list or None as value.
    Settings for feature classes specified in enabledFeatures.keys are updated, settings for feature classes
    not yet present in enabledFeatures.keys are added.
    To disable the entire class, use :py:func:`disableAllFeatures` or :py:func:`enableFeatureClassByName` instead.
    """
    global logger

    logger.debug('Updating enabled features with %s', enabledFeatures)
    self.enabledFeatures.update(enabledFeatures)
    logger.debug('Enabled features: %s', self.enabledFeatures)
