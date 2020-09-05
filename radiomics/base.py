import inspect
import logging
import traceback

import numpy
import SimpleITK as sitk
import six

from radiomics import getProgressReporter, imageoperations


class RadiomicsFeaturesBase(object):
  """
  This is the abstract class, which defines the common interface for the feature classes. All feature classes inherit
  (directly of indirectly) from this class.

  At initialization, image and labelmap are passed as SimpleITK image objects (``inputImage`` and ``inputMask``,
  respectively.) The motivation for using SimpleITK images as input is to keep the possibility of reusing the
  optimized feature calculators implemented in SimpleITK in the future. If either the image or the mask is None,
  initialization fails and a warning is logged (does not raise an error).

  Logging is set up using a child logger from the parent 'radiomics' logger. This retains the toolbox structure in
  the generated log. The child logger is named after the module containing the feature class (e.g. 'radiomics.glcm').

  Any pre calculations needed before the feature functions are called can be added by overriding the
  ``_initSegmentBasedCalculation`` function, which prepares the input for feature extraction. If image discretization is
  needed, this can be implemented by adding a call to ``_applyBinning`` to this initialization function, which also
  instantiates coefficients holding the maximum ('Ng') and unique ('GrayLevels') that can be found inside the ROI after
  binning. This function also instantiates the `matrix` variable, which holds the discretized image (the `imageArray`
  variable will hold only original gray levels).

  The following variables are instantiated at initialization:

  - kwargs: dictionary holding all customized settings passed to this feature class.
  - label: label value of Region of Interest (ROI) in labelmap. If key is not present, a default value of 1 is used.
  - featureNames: list containing the names of features defined in the feature class. See :py:func:`getFeatureNames`
  - inputImage: SimpleITK image object of the input image (dimensions x, y, z)

  The following variables are instantiated by the ``_initSegmentBasedCalculation`` function:

  - inputMask: SimpleITK image object of the input labelmap (dimensions x, y, z)
  - imageArray: numpy array of the gray values in the input image (dimensions z, y, x)
  - maskArray: numpy boolean array with elements set to ``True`` where labelmap = label, ``False`` otherwise,
    (dimensions z, y, x).
  - labelledVoxelCoordinates: tuple of 3 numpy arrays containing the z, x and y coordinates of the voxels included in
    the ROI, respectively. Length of each array is equal to total number of voxels inside ROI.
  - matrix: copy of the imageArray variable, with gray values inside ROI discretized using the specified binWidth.
    This variable is only instantiated if a call to ``_applyBinning`` is added to an override of
    ``_initSegmentBasedCalculation`` in the feature class.

  .. note::
    Although some variables listed here have similar names to customization settings, they do *not* represent all the
    possible settings on the feature class level. These variables are listed here to help developers develop new feature
    classes, which make use of these variables. For more information on customization, see
    :ref:`radiomics-customization-label`, which includes a comprehensive list of all possible settings, including
    default values and explanation of usage.
  """

  def __init__(self, inputImage, inputMask, **kwargs):
    self.logger = logging.getLogger(self.__module__)
    self.logger.debug('Initializing feature class')

    if inputImage is None or inputMask is None:
      raise ValueError('Missing input image or mask')

    self.progressReporter = getProgressReporter

    self.settings = kwargs

    self.label = kwargs.get('label', 1)
    self.voxelBased = kwargs.get('voxelBased', False)

    self.coefficients = {}

    # all features are disabled by default
    self.enabledFeatures = {}
    self.featureValues = {}

    self.featureNames = self.getFeatureNames()

    self.inputImage = inputImage
    self.inputMask = inputMask

    self.imageArray = sitk.GetArrayFromImage(self.inputImage)

    if self.voxelBased:
      self._initVoxelBasedCalculation()
    else:
      self._initSegmentBasedCalculation()

  def _initSegmentBasedCalculation(self):
    self.maskArray = (sitk.GetArrayFromImage(self.inputMask) == self.label)  # boolean array

  def _initVoxelBasedCalculation(self):
    self.masked = self.settings.get('maskedKernel', True)

    maskArray = sitk.GetArrayFromImage(self.inputMask) == self.label  # boolean array
    self.labelledVoxelCoordinates = numpy.array(numpy.where(maskArray))

    # Set up the mask array for the gray value discretization
    if self.masked:
      self.maskArray = maskArray
    else:
      # This will cause the discretization to use the entire image
      self.maskArray = numpy.ones(self.imageArray.shape, dtype='bool')

  def _initCalculation(self, voxelCoordinates=None):
    """
    Last steps to prepare the class for extraction. This function calculates the texture matrices and coefficients in
    the respective feature classes
    """
    pass

  def _applyBinning(self, matrix):
    matrix, _ = imageoperations.binImage(matrix, self.maskArray, **self.settings)
    self.coefficients['grayLevels'] = numpy.unique(matrix[self.maskArray])
    self.coefficients['Ng'] = int(numpy.max(self.coefficients['grayLevels']))  # max gray level in the ROI
    return matrix

  def enableFeatureByName(self, featureName, enable=True):
    """
    Enables or disables feature specified by ``featureName``. If feature is not present in this class, a lookup error is
    raised. ``enable`` specifies whether to enable or disable the feature.
    """
    if featureName not in self.featureNames:
      raise LookupError('Feature not found: ' + featureName)
    if self.featureNames[featureName]:
      self.logger.warning('Feature %s is deprecated, use with caution!', featureName)
    self.enabledFeatures[featureName] = enable

  def enableAllFeatures(self):
    """
    Enables all features found in this class for calculation.

    .. note::
      Features that have been marked "deprecated" are not enabled by this function. They can still be enabled manually by
      a call to :py:func:`~radiomics.base.RadiomicsBase.enableFeatureByName()`,
      :py:func:`~radiomics.featureextractor.RadiomicsFeaturesExtractor.enableFeaturesByName()`
      or in the parameter file (by specifying the feature by name, not when enabling all features).
      However, in most cases this will still result only in a deprecation warning.
    """
    for featureName, is_deprecated in six.iteritems(self.featureNames):
      # only enable non-deprecated features here
      if not is_deprecated:
        self.enableFeatureByName(featureName, True)

  def disableAllFeatures(self):
    """
    Disables all features. Additionally resets any calculated features.
    """
    self.enabledFeatures = {}
    self.featureValues = {}

  @classmethod
  def getFeatureNames(cls):
    """
    Dynamically enumerates features defined in the feature class. Features are identified by the
    ``get<Feature>FeatureValue`` signature, where <Feature> is the name of the feature (unique on the class level).

    Found features are returned as a dictionary of the feature names, where the value ``True`` if the
    feature is deprecated, ``False`` otherwise (``{<Feature1>:<deprecated>, <Feature2>:<deprecated>, ...}``).

    This function is called at initialization, found features are stored in the ``featureNames`` variable.
    """
    attributes = inspect.getmembers(cls)
    features = {a[0][3:-12]: getattr(a[1], '_is_deprecated', False) for a in attributes
                if a[0].startswith('get') and a[0].endswith('FeatureValue')}
    return features

  def execute(self):
    """
    Calculates all features enabled in  ``enabledFeatures``. A feature is enabled if it's key is present in this
    dictionary and it's value is True.

    Calculated values are stored in the ``featureValues`` dictionary, with feature name as key and the calculated
    feature value as value. If an exception is thrown during calculation, the error is logged, and the value is set to
    NaN.
    """
    if len(self.enabledFeatures) == 0:
      self.enableAllFeatures()

    if self.voxelBased:
      self._calculateVoxels()
    else:
      self._calculateSegment()

    return self.featureValues

  def _calculateVoxels(self):
    initValue = self.settings.get('initValue', 0)
    voxelBatch = self.settings.get('voxelBatch', -1)

    # Initialize the output with empty numpy arrays
    for feature, enabled in six.iteritems(self.enabledFeatures):
      if enabled:
        self.featureValues[feature] = numpy.full(list(self.inputImage.GetSize())[::-1], initValue, dtype='float')

    # Calculate the feature values for all enabled features
    voxel_count = self.labelledVoxelCoordinates.shape[1]
    voxel_batch_idx = 0
    if voxelBatch < 0:
      voxelBatch = voxel_count
    n_batches = numpy.ceil(float(voxel_count) / voxelBatch)
    with self.progressReporter(total=n_batches, desc='batch') as pbar:
      while voxel_batch_idx < voxel_count:
        self.logger.debug('Calculating voxel batch no. %i/%i', int(voxel_batch_idx / voxelBatch) + 1, n_batches)
        voxelCoords = self.labelledVoxelCoordinates[:, voxel_batch_idx:voxel_batch_idx + voxelBatch]
        # Calculate the feature values for the current kernel
        for success, featureName, featureValue in self._calculateFeatures(voxelCoords):
          if success:
            self.featureValues[featureName][tuple(voxelCoords)] = featureValue

        voxel_batch_idx += voxelBatch
        pbar.update(1)  # Update progress bar

    # Convert the output to simple ITK image objects
    for feature, enabled in six.iteritems(self.enabledFeatures):
      if enabled:
        self.featureValues[feature] = sitk.GetImageFromArray(self.featureValues[feature])
        self.featureValues[feature].CopyInformation(self.inputImage)

  def _calculateSegment(self):
    # Get the feature values using the current segment.
    for success, featureName, featureValue in self._calculateFeatures():
      # Always store the result. In case of an error, featureValue will be NaN
      self.featureValues[featureName] = numpy.squeeze(featureValue)

  def _calculateFeatures(self, voxelCoordinates=None):
    # Initialize the calculation
    # This function serves to calculate the texture matrices where applicable
    self._initCalculation(voxelCoordinates)

    self.logger.debug('Calculating features')
    for feature, enabled in six.iteritems(self.enabledFeatures):
      if enabled:
        try:
          # Use getattr to get the feature calculation methods, then use '()' to evaluate those methods
          yield True, feature, getattr(self, 'get%sFeatureValue' % feature)()
        except DeprecationWarning as deprecatedFeature:
          # Add a debug log message, as a warning is usually shown and would entail a too verbose output
          self.logger.debug('Feature %s is deprecated: %s', feature, deprecatedFeature.args[0])
        except Exception:
          self.logger.error('FAILED: %s', traceback.format_exc())
          yield False, feature, numpy.nan
