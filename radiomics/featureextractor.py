# -*- coding: utf-8 -*-
from __future__ import print_function

import collections
from itertools import chain
import json
import logging
import os

import pykwalify.core
import SimpleITK as sitk
import six

from radiomics import generalinfo, getFeatureClasses, getImageTypes, getParameterValidationFiles, \
  imageoperations


class RadiomicsFeaturesExtractor:
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
    self.logger = logging.getLogger(__name__)

    self.featureClasses = getFeatureClasses()

    self.generalInfo = None

    self.settings = {}
    self._enabledImagetypes = {}
    self._enabledFeatures = {}

    if len(args) == 1 and isinstance(args[0], six.string_types) and os.path.isfile(args[0]):
      self.logger.info("Loading parameter file")
      self._applyParams(paramsFile=args[0])
    elif len(args) == 1 and isinstance(args[0], dict):
        self.logger.info("Loading parameter dictionary")
        self._applyParams(paramsDict=args[0])
    else:
      # Set default settings and update with and changed settings contained in kwargs
      self.settings = self._getDefaultSettings()
      self.logger.info('No valid config parameter, applying defaults: %s', self.settings)

      self._enabledImagetypes = {'Original': {}}
      self.logger.info('Enabled image types: %s', self._enabledImagetypes)
      self._enabledFeatures = {}

      for featureClassName in self.getFeatureClassNames():
        self._enabledFeatures[featureClassName] = []
      self.logger.info('Enabled features: %s', self._enabledFeatures)

    if len(kwargs) > 0:
      self.logger.info('Applying custom setting overrides')
      self.settings.update(kwargs)
      self.logger.debug("Settings: %s", self.settings)

    self._setTolerance()

  @classmethod
  def _getDefaultSettings(cls):
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

  def _setTolerance(self):
    self.geometryTolerance = self.settings.get('geometryTolerance')
    if self.geometryTolerance is not None:
      sitk.ProcessObject.SetGlobalDefaultCoordinateTolerance(self.geometryTolerance)
      sitk.ProcessObject.SetGlobalDefaultDirectionTolerance(self.geometryTolerance)

  def addProvenance(self, provenance_on=True):
    """
    Enable or disable reporting of additional information on the extraction. This information includes toolbox version,
    enabled input images and applied settings. Furthermore, additional information on the image and region of interest
    (ROI) is also provided, including original image spacing, total number of voxels in the ROI and total number of
    fully connected volumes in the ROI.

    To disable this, call ``addProvenance(False)``.
    """
    self.settings['additionalInfo'] = provenance_on

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
    # Ensure pykwalify.core has a log handler (needed when parameter validation fails)
    if len(pykwalify.core.log.handlers) == 0 and len(logging.getLogger().handlers) == 0:
      # No handler available for either pykwalify or root logger, provide first radiomics handler (outputs to stderr)
      pykwalify.core.log.addHandler(logging.getLogger('radiomics').handlers[0])

    schemaFile, schemaFuncs = getParameterValidationFiles()
    c = pykwalify.core.Core(source_file=paramsFile, source_data=paramsDict,
                            schema_files=[schemaFile], extensions=[schemaFuncs])
    params = c.validate()
    self.logger.debug('Parameters parsed, input is valid.')

    enabledImageTypes = params.get('imageType', {})
    enabledFeatures = params.get('featureClass', {})
    settings = params.get('setting', {})
    voxelSettings = params.get('voxelSetting', {})

    self.logger.debug("Applying settings")

    if len(enabledImageTypes) == 0:
      self._enabledImagetypes = {'Original': {}}
    else:
      self._enabledImagetypes = enabledImageTypes

    self.logger.debug("Enabled image types: %s", self._enabledImagetypes)

    if len(enabledFeatures) == 0:
      self._enabledFeatures = {}
      for featureClassName in self.getFeatureClassNames():
        self._enabledFeatures[featureClassName] = []
    else:
      self._enabledFeatures = enabledFeatures

    self.logger.debug("Enabled features: %s", self._enabledFeatures)

    # Set default settings and update with and changed settings contained in kwargs
    self.settings = self._getDefaultSettings()
    self.settings.update(settings)
    self.settings.update(voxelSettings)

    self.logger.debug("Settings: %s", settings)

  def enableAllImageTypes(self):
    """
    Enable all possible image types without any custom settings.
    """
    self.logger.debug('Enabling all image types')
    for imageType in getImageTypes():
      self._enabledImagetypes[imageType] = {}
    self.logger.debug('Enabled images types: %s', self._enabledImagetypes)

  def disableAllImageTypes(self):
    """
    Disable all image types.
    """
    self.logger.debug('Disabling all image types')
    self._enabledImagetypes = {}

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

    For the mathmetical formulas of square, squareroot, logarithm and exponential, see their respective functions in
    :ref:`imageoperations<radiomics-imageoperations-label>`
    (:py:func:`~radiomics.imageoperations.getSquareImage`,
    :py:func:`~radiomics.imageoperations.getSquareRootImage`,
    :py:func:`~radiomics.imageoperations.getLogarithmImage` and
    :py:func:`~radiomics.imageoperations.getExponentialImage`,
    respectively).
    """
    if imageType not in getImageTypes():
      self.logger.warning('Image type %s is not recognized', imageType)
      return

    if enabled:
      if customArgs is None:
        customArgs = {}
        self.logger.debug('Enabling image type %s (no additional custom settings)', imageType)
      else:
        self.logger.debug('Enabling image type %s (additional custom settings: %s)', imageType, customArgs)
      self._enabledImagetypes[imageType] = customArgs
    elif imageType in self._enabledImagetypes:
      self.logger.debug('Disabling image type %s', imageType)
      del self._enabledImagetypes[imageType]
    self.logger.debug('Enabled images types: %s', self._enabledImagetypes)

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
    self.logger.debug('Updating enabled images types with %s', enabledImagetypes)
    self._enabledImagetypes.update(enabledImagetypes)
    self.logger.debug('Enabled images types: %s', self._enabledImagetypes)

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
    self.logger.debug('Enabling all features in all feature classes')
    for featureClassName in self.getFeatureClassNames():
      self._enabledFeatures[featureClassName] = []
    self.logger.debug('Enabled features: %s', self._enabledFeatures)

  def disableAllFeatures(self):
    """
    Disable all classes.
    """
    self.logger.debug('Disabling all feature classes')
    self._enabledFeatures = {}

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
    if featureClass not in self.getFeatureClassNames():
      self.logger.warning('Feature class %s is not recognized', featureClass)
      return

    if enabled:
      self.logger.debug('Enabling all features in class %s', featureClass)
      self._enabledFeatures[featureClass] = []
    elif featureClass in self._enabledFeatures:
      self.logger.debug('Disabling feature class %s', featureClass)
      del self._enabledFeatures[featureClass]
    self.logger.debug('Enabled features: %s', self._enabledFeatures)

  def enableFeaturesByName(self, **enabledFeatures):
    """
    Specify which features to enable. Key is feature class name, value is a list of enabled feature names.

    To enable all features for a class, provide the class name with an empty list or None as value.
    Settings for feature classes specified in enabledFeatures.keys are updated, settings for feature classes
    not yet present in enabledFeatures.keys are added.
    To disable the entire class, use :py:func:`disableAllFeatures` or :py:func:`enableFeatureClassByName` instead.
    """
    self.logger.debug('Updating enabled features with %s', enabledFeatures)
    self._enabledFeatures.update(enabledFeatures)
    self.logger.debug('Enabled features: %s', self._enabledFeatures)

  def execute(self, imageFilepath, maskFilepath, label=None, voxelBased=False):
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
    :returns: dictionary containing calculated signature ("<imageType>_<featureClass>_<featureName>":value).
    """
    geometryTolerance = self.settings.get('geometryTolerance')
    additionalInfo = self.settings.get('additionalInfo', False)
    resegmentShape = self.settings.get('resegmentShape', False)

    if label is not None:
      self.settings['label'] = label
    else:
      label = self.settings.get('label', 1)

    if self.geometryTolerance != geometryTolerance:
      self._setTolerance()

    if additionalInfo:
      self.generalInfo = generalinfo.GeneralInfo()
      self.generalInfo.addGeneralSettings(self.settings)
      self.generalInfo.addEnabledImageTypes(self._enabledImagetypes)
    else:
      self.generalInfo = None

    if voxelBased:
      self.settings['voxelBased'] = True
      self.logger.info('Starting voxel based extraction')

    self.logger.info('Calculating features with label: %d', label)
    self.logger.debug('Enabled images types: %s', self._enabledImagetypes)
    self.logger.debug('Enabled features: %s', self._enabledFeatures)
    self.logger.debug('Current settings: %s', self.settings)

    if self.settings.get('binCount', None) is not None:
      self.logger.warning('Fixed bin Count enabled! However, we recommend using a fixed bin Width. See '
                          'http://pyradiomics.readthedocs.io/en/latest/faq.html#radiomics-fixed-bin-width for more '
                          'details')

    # 1. Load the image and mask
    featureVector = collections.OrderedDict()
    image, mask = self.loadImage(imageFilepath, maskFilepath)

    # 2. Check whether loaded mask contains a valid ROI for feature extraction and get bounding box
    # Raises a ValueError if the ROI is invalid
    boundingBox, correctedMask = imageoperations.checkMask(image, mask, **self.settings)

    # Update the mask if it had to be resampled
    if correctedMask is not None:
      if self.generalInfo is not None:
        self.generalInfo.addMaskElements(image, correctedMask, label, 'corrected')
      mask = correctedMask

    self.logger.debug('Image and Mask loaded and valid, starting extraction')

    # 5. Resegment the mask if enabled (parameter regsegmentMask is not None)
    resegmentedMask = None
    if self.settings.get('resegmentRange', None) is not None:
      resegmentedMask = imageoperations.resegmentMask(image, mask, **self.settings)

      # Recheck to see if the mask is still valid, raises a ValueError if not
      boundingBox, correctedMask = imageoperations.checkMask(image, resegmentedMask, **self.settings)

      if self.generalInfo is not None:
        self.generalInfo.addMaskElements(image, resegmentedMask, label, 'resegmented')

    # if resegmentShape is True and resegmentation has been enabled, update the mask here to also use the
    # resegmented mask for shape calculation (e.g. PET resegmentation)
    if resegmentShape and resegmentedMask is not None:
      mask = resegmentedMask

    if not voxelBased:
      # 3. Add the additional information if enabled
      if self.generalInfo is not None:
        featureVector.update(self.generalInfo.getGeneralInfo())

      # 4. If shape descriptors should be calculated, handle it separately here
      if 'shape' in self._enabledFeatures.keys():
        croppedImage, croppedMask = imageoperations.cropToTumorMask(image, mask, boundingBox)
        enabledFeatures = self._enabledFeatures['shape']

        self.logger.info('Computing shape')
        shapeClass = self.featureClasses['shape'](croppedImage, croppedMask, **self.settings)
        if enabledFeatures is None or len(enabledFeatures) == 0:
          shapeClass.enableAllFeatures()
        else:
          for feature in enabledFeatures:
            shapeClass.enableFeatureByName(feature)

        for (featureName, featureValue) in six.iteritems(shapeClass.execute()):
          newFeatureName = 'original_shape_%s' % featureName
          featureVector[newFeatureName] = featureValue

    # (Default) Only use resegemented mask for feature classes other than shape
    # can be overridden by specifying `resegmentShape` = True
    if not resegmentShape and resegmentedMask is not None:
      mask = resegmentedMask

    # 6. Calculate other enabled feature classes using enabled image types
    # Make generators for all enabled image types
    self.logger.debug('Creating image type iterator')
    imageGenerators = []
    for imageType, customKwargs in six.iteritems(self._enabledImagetypes):
      args = self.settings.copy()
      args.update(customKwargs)
      self.logger.info('Adding image type "%s" with custom settings: %s' % (imageType, str(customKwargs)))
      imageGenerators = chain(imageGenerators, getattr(imageoperations, 'get%sImage' % imageType)(image, mask, **args))

    self.logger.debug('Extracting features')
    # Calculate features for all (filtered) images in the generator
    for inputImage, imageTypeName, inputKwargs in imageGenerators:
      self.logger.info('Calculating features for %s image', imageTypeName)
      inputImage, inputMask = imageoperations.cropToTumorMask(inputImage, mask, boundingBox)
      featureVector.update(self.computeFeatures(inputImage, inputMask, imageTypeName, **inputKwargs))

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
    normalize = self.settings.get('normalize', False)
    interpolator = self.settings.get('interpolator')
    resampledPixelSpacing = self.settings.get('resampledPixelSpacing')
    preCrop = self.settings.get('preCrop', False)
    label = self.settings.get('label', 1)

    self.logger.info('Loading image and mask')
    if isinstance(ImageFilePath, six.string_types) and os.path.isfile(ImageFilePath):
      image = sitk.ReadImage(ImageFilePath)
    elif isinstance(ImageFilePath, sitk.SimpleITK.Image):
      image = ImageFilePath
    else:
      raise ValueError('Error reading image Filepath or SimpleITK object')

    if isinstance(MaskFilePath, six.string_types) and os.path.isfile(MaskFilePath):
      mask = sitk.ReadImage(MaskFilePath, sitk.sitkUInt32)
    elif isinstance(MaskFilePath, sitk.SimpleITK.Image):
      mask = sitk.Cast(MaskFilePath, sitk.sitkUInt32)
    else:
      raise ValueError('Error reading mask Filepath or SimpleITK object')

    if self.generalInfo is not None:
      self.generalInfo.addImageElements(image)
      # Do not include the image here, as the overlap between image and mask have not been checked
      # It is therefore possible that image and mask do not align, or even have different sizes.
      self.generalInfo.addMaskElements(None, mask, label)

    # This point is only reached if image and mask loaded correctly
    if normalize:
      image = imageoperations.normalizeImage(image, **self.settings)

    if interpolator is not None and resampledPixelSpacing is not None:
      image, mask = imageoperations.resampleImage(image, mask, **self.settings)
      if self.generalInfo is not None:
        self.generalInfo.addImageElements(image, 'interpolated')
        self.generalInfo.addMaskElements(image, mask, self.settings.get('label', 1), 'interpolated')

    elif preCrop:
      bb, correctedMask = imageoperations.checkMask(image, mask, **self.settings)
      if correctedMask is not None:
        # Update the mask if it had to be resampled
        mask = correctedMask
      if bb is None:
        # Mask checks failed
        raise ValueError('Mask checks failed during pre-crop')

      image, mask = imageoperations.cropToTumorMask(image, mask, bb, **self.settings)

    return image, mask

  def computeFeatures(self, image, mask, imageTypeName, **kwargs):
    r"""
    Compute signature using image, mask, \*\*kwargs settings.

    This function computes the signature for just the passed image (original or derived), it does not preprocess or
    apply a filter to the passed image. Features / Classes to use for calculation of signature are defined in
    ``self._enabledFeatures``. See also :py:func:`enableFeaturesByName`.

    .. note::

      shape descriptors are independent of gray level and therefore calculated separately (handled in `execute`). In
      this function, no shape functions are calculated.
    """
    featureVector = collections.OrderedDict()

    # Calculate feature classes
    for featureClassName, enabledFeatures in six.iteritems(self._enabledFeatures):
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

        for (featureName, featureValue) in six.iteritems(featureClass.execute()):
          newFeatureName = '%s_%s_%s' % (imageTypeName, featureClassName, featureName)
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
