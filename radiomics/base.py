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
  - binWidth: bin width, as specified in ``**kwargs``. If key is not present, a default value of 25 is used.
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
  - boundingBoxSize: tuple of 3 integers containing the z, x and y sizes of the ROI bounding box, respectively.
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
    self.progressReporter = getProgressReporter

    self.kwargs = kwargs
    self.binWidth = kwargs.get('binWidth', 25)
    self.label = kwargs.get('label', 1)

    self.coefficients = {}

    # all features are disabled by default
    self.disableAllFeatures()

    self.featureNames = self.getFeatureNames()

    self.inputImage = inputImage
    self.inputMask = inputMask

    if inputImage is None or inputMask is None:
      self.logger.warning('Missing input image or mask')
      return

  def _initSegmentBasedCalculation(self):
    self.imageArray = sitk.GetArrayFromImage(self.inputImage)
    self.maskArray = (sitk.GetArrayFromImage(self.inputMask) == self.label)  # boolean array

    self.labelledVoxelCoordinates = numpy.where(self.maskArray != 0)
    self.boundingBoxSize = numpy.max(self.labelledVoxelCoordinates, 1) - numpy.min(self.labelledVoxelCoordinates, 1) + 1

  def _applyBinning(self):
    self.matrix, self.binEdges = imageoperations.binImage(self.binWidth,
                                                          self.imageArray,
                                                          self.maskArray)
    self.coefficients['grayLevels'] = numpy.unique(self.matrix[self.maskArray])
    self.coefficients['Ng'] = int(numpy.max(self.coefficients['grayLevels']))  # max gray level in the ROI

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
    for featureName, deprecated in six.iteritems(self.featureNames):
      # only enable non-deprecated features here
      if not deprecated:
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

  def calculateFeatures(self):
    """
    Calculates all features enabled in  ``enabledFeatures``. A feature is enabled if it's key is present in this
    dictionary and it's value is True.

    Calculated values are stored in the ``featureValues`` dictionary, with feature name as key and the calculated
    feature value as value. If an exception is thrown during calculation, the error is logged, and the value is set to
    NaN.
    """
    self.logger.debug('Calculating features')
    for feature, enabled in six.iteritems(self.enabledFeatures):
      if enabled:
        try:
          # Use getattr to get the feature calculation methods, then use '()' to evaluate those methods
          self.featureValues[feature] = getattr(self, 'get%sFeatureValue' % feature)()
        except DeprecationWarning as deprecatedFeature:
          self.logger.warning('Feature %s is deprecated: %s', feature, deprecatedFeature.message)
        except Exception:
          self.featureValues[feature] = numpy.nan
          self.logger.error('FAILED: %s', traceback.format_exc())
