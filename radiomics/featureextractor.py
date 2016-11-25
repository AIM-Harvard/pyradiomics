# -*- coding: utf-8 -*-
import os
import logging
import collections
from itertools import chain
import numpy
import SimpleITK as sitk
import pkgutil
import inspect
import pykwalify.core
import radiomics
from radiomics import base, imageoperations, generalinfo


class RadiomicsFeaturesExtractor:
  """
  Wrapper class for calculation of a radiomics signature.
  At and after initialisation various settings can be used to customize the resultant signature.
  This includes which classes and features to use, as well as what should be done in terms of preprocessing the image
  and what images (original and/or filtered) should be used as input.

  Then a call to computeSignature generates the radiomics signature specified by these settings for the passed
  image and labelmap combination. This function can be called repeatedly in a batch process to calculate the radiomics
  signature for all image and labelmap combinations.

  It initialisation, a parameters file can be provided containing all necessary settings. This is done by passing the
  location of the file as the single argument in the initialization call, without specifying it as a keyword argument.
  If such a file location is provided, any additional kwargs are ignored.
  Alternatively, at initialisation, the following general settings can be specified in kwargs:

  - verbose [True]: Boolean, set to False to disable status update printing.
  - binWidth [25]: Float, size of the bins when making a histogram and for discretization of the image gray level.
  - resampledPixelSpacing [None]: List of 3 floats, sets the size of the voxel in (x, y, z) plane when resampling.
  - interpolator [sitkBSpline]: Simple ITK constant or string name thereof, sets interpolator to use for resampling.
    Enumerated value, possible values:

    - sitkNearestNeighbor
    - sitkLinear
    - sitkBSpline
    - sitkGaussian
    - sitkLabelGaussian
    - sitkHammingWindowedSinc
    - sitkCosineWindowedSinc
    - sitkWelchWindowedSinc
    - sitkLanczosWindowedSinc
    - sitkBlackmanWindowedSinc

  - padDistance [5]: Integer, set the number of voxels pad cropped tumor volume with during resampling. Padding occurs
    in new feature space and is done on all faces, i.e. size increases in x, y and z direction by 2*padDistance.
    Padding is needed for some filters (e.g. LoG). After application of filters image is cropped again without
    padding. Value of padded voxels are set to original gray level intensity, padding does not exceed original image
    boundaries.

  N.B. Resampling is disabled when either `resampledPixelSpacing` or `interpolator` is set to `None`

  In addition to these general settings, filter or featureclass specific settings can be defined here also.
  For more information on possible settings, see the respective filters and feature classes.

  By default, all features in all feature classes are enabled.
  By default, only original input image is enabled
  N.B. for log, the sigma is set to range 0.5-5.0, step size 0.5
  """

  def __init__(self, *args, **kwargs):
    self.logger = logging.getLogger(__name__)

    self.featureClasses = self.getFeatureClasses()

    self.kwargs = {}
    self.provenance_on = True
    self.inputImages = {}
    self.enabledFeatures = {}

    if len(args) == 1 and isinstance(args[0], basestring):
      self.loadParams(args[0])
    else:
      # Set default settings and update with and changed settings contained in kwargs
      self.kwargs = {'resampledPixelSpacing': None,  # No resampling by default
                     'interpolator': sitk.sitkBSpline,
                     'padDistance': 5,
                     'label': 1,
                     'verbose': False}
      self.kwargs.update(kwargs)

      self.inputImages = {'original': {}}

      self.enabledFeatures = {}
      for featureClassName in self.getFeatureClassNames():
        self.enabledFeatures[featureClassName] = []

  def addProvenance(self, provenance_on=True):
    """
    Enable or disable reporting of settings used for calculated filters.
    By default, settings used are added to the dictionary of calculated features as {"settings_<filter>":<settings>}
    To disable this, call ``addProvenance(False)``
    """
    self.provenance_on = provenance_on

  def loadParams(self, paramsFile):
    """
    Parse specified parameters file and use it to update settings in kwargs, enabled feature(Classes) and input
    images:
    - kwarg settings not specified in parameters are set to their default value.
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
    dataDir = os.path.abspath(os.path.join(radiomics.__path__[0], '..', 'data'))
    schemaFile = os.path.join(dataDir, 'paramSchema.yaml')
    schemaFuncs = os.path.join(dataDir, 'schemaFuncs.py')
    c = pykwalify.core.Core(source_file=paramsFile, schema_files=[schemaFile], extensions=[schemaFuncs])
    params = c.validate()

    inputImages = params.get('inputImage', {})
    enabledFeatures = params.get('featureClass', {})
    kwargs = params.get('setting', {})

    if len(inputImages) == 0:
      self.inputImages = {'original': {}}
    else:
      self.inputImages = inputImages

    if len(enabledFeatures) == 0:
      self.enabledFeatures = {}
      for featureClassName in self.getFeatureClassNames():
        self.enabledFeatures[featureClassName] = []
    else:
      self.enabledFeatures = enabledFeatures

    # Set default settings and update with and changed settings contained in kwargs
    self.kwargs = {'resampledPixelSpacing': None,  # No resampling by default
                   'interpolator': sitk.sitkBSpline,
                   'padDistance': 5,
                   'label': 1,
                   'verbose': False}
    self.kwargs.update(kwargs)

  def enableAllInputImages(self):
    """
    Enable all possible input images without any custom settings.
    """
    for imageType in self.getInputImageTypes():
      self.inputImages[imageType] = {}

  def disableAllInputImages(self):
    """
    Disable all input images.
    """
    self.inputImages = {}

  def enableInputImageByName(self, inputImage, enabled=True, customArgs=None):
    r"""
    Enable or disable specified input image. If enabling input image, optional custom settings can be specified in
    customArgs.

    Current possible filters are:

    - original: No filter applied
    - wavelet: Wavelet filtering, yields 8 decompositions per level (all possible combinations of applying either
      a High or a Low pass filter in each of the three dimensions.
    - log: Laplacian of Gaussian filter, edge enhancement filter. Emphasizes areas of gray level change, where sigma
      defines how coarse the emphasised texture should be. A low sigma emphasis on fine textures (change over a
      short distance), where a high sigma value emphasises coarse textures (gray level change over a large distance)
    - square: Takes the square of the image intensities and linearly scales them back to the original range.
      Negative values in the original image will be made negative again after application of filter.
    - squareroot: Takes the square root of the absolute image intensities and scales them back to original range.
      Negative values in the original image will be made negative again after application of filter.
    - logarithm: Takes the logarithm of the absolute intensity + 1. Values are scaled to original range and
      negative original values are made negative again after application of filter.
    - exponential: Takes the the exponential, where filtered intensity is e^(absolute intensity). Values are
      scaled to original range and negative original values are made negative again after application of filter.

    For the mathmetical formulas of square, squareroot, logarithm and exponential, see their respective functions in
    :ref:`imageoperations<radiomics-imageoperations-label>`
    (:py:func:`~radiomics.imageoperations.applySquare`, :py:func:`~radiomics.imageoperations.applySquareRoot`,
    :py:func:`~radiomics.imageoperations.applyLogarithm` and :py:func:`~radiomics.imageoperations.applyExponential`,
    respectively).
    """
    if enabled:
      if customArgs is None:
        customArgs = {}
      self.inputImages[inputImage] = customArgs
    elif inputImage in self.inputImages:
      del self.inputImages[inputImage]

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
    To disable input images, use enableInputImageByName or disableAllInputImages instead.

    :param inputImages: dictionary, key is imagetype (original, wavelet or log) and value is custom settings (dictionary)
    """
    self.inputImages.update(inputImages)

  def enableAllFeatures(self):
    """
    Enable all classes and all features.
    """
    for featureClassName in self.getFeatureClassNames():
      self.enabledFeatures[featureClassName] = []

  def disableAllFeatures(self):
    """
    Disable all classes.
    """
    self.enabledFeatures = {}

  def enableFeatureClassByName(self, featureClass, enabled=True):
    """
    Enable or disable all features in given class.
    """
    if enabled:
      self.enabledFeatures[featureClass] = []
    elif featureClass in self.enabledFeatures:
      del self.enabledFeatures[featureClass]

  def enableFeaturesByName(self, **enabledFeatures):
    """
    Specify which features to enable. Key is feature class name, value is a list of enabled feature names.

    To enable all features for a class, provide the class name with an empty list or None as value.
    Settings for feature classes specified in enabledFeatures.keys are updated, settings for feature classes
    not yet present in enabledFeatures.keys are added.
    To disable the entire class, use disableAllFeatures or enableFeatureClassByName instead.
    """
    self.enabledFeatures.update(enabledFeatures)

  def execute(self, imageFilepath, maskFilepath, label=None):
    """
    Compute radiomics signature for provide image and mask combination.
    First, image and mask are loaded and resampled if necessary. Next shape features are calculated on a
    cropped (no padding) version of the original image. Then other featureclasses are calculated on using all
    specified filters in inputImages. Images are cropped to tumor mask (no padding) after application of filter and
    before being passed to the feature class.
    Finally, a dictionary containing all calculated features is returned.

    :param imageFilepath: SimpleITK Image, or string pointing to image file location
    :param maskFilepath: SimpleITK Image, or string pointing to labelmap file location
    :param label: Integer, value of the label for which to extract features. If not specified, last specified label
        is used. Default label is 1.
    :returns: dictionary containing calculated signature ("<filter>_<featureClass>_<featureName>":value).
    """
    if label is not None:
      self.kwargs.update({'label': label})

    self.logger.info('Calculating features with label: %d', self.kwargs['label'])

    featureVector = collections.OrderedDict()
    image, mask = self.loadImage(imageFilepath, maskFilepath)

    if self.provenance_on:
      # Add Provenance information
      generalinfoClass = generalinfo.GeneralInfo(imageFilepath, maskFilepath, mask, self.kwargs, self.inputImages)
      for k, v in generalinfoClass.execute().iteritems():
        featureVector['general_info_%s' % (k)] = v

    # If shape should be calculation, handle it separately here
    if 'shape' in self.enabledFeatures.keys():
      croppedImage, croppedMask = imageoperations.cropToTumorMask(image, mask, self.kwargs['label'])
      enabledFeatures = self.enabledFeatures['shape']
      shapeClass = self.featureClasses['shape'](croppedImage, croppedMask, **self.kwargs)
      if enabledFeatures is None or len(enabledFeatures) == 0:
        shapeClass.enableAllFeatures()
      else:
        for feature in enabledFeatures:
          shapeClass.enableFeatureByName(feature)

      if self.kwargs['verbose']: print "\t\tComputing shape"
      shapeClass.calculateFeatures()
      for (featureName, featureValue) in shapeClass.featureValues.iteritems():
        newFeatureName = "original_shape_%s" % (featureName)
        featureVector[newFeatureName] = featureValue

    # Make generators for all enabled input image types
    imageGenerators = []
    for imageType, customKwargs in self.inputImages.iteritems():
      args = self.kwargs.copy()
      args.update(customKwargs)
      self.logger.info("Applying filter: '%s' with settings: %s" % (imageType, str(args)))
      imageGenerators = chain(imageGenerators, eval('self.generate_%s(image, mask, **args)' % (imageType)))

    # Calculate features for all (filtered) images in the generator
    for inputImage, inputMask, inputImageName, inputKwargs in imageGenerators:
      featureVector.update(self.computeFeatures(inputImage, inputMask, inputImageName, **inputKwargs))

    return featureVector

  def loadImage(self, ImageFilePath, MaskFilePath):
    """
    Preprocess the image and labelmap.
    If ImageFilePath is a string, it is loaded as SimpleITK Image and assigned to image,
    if it already is a SimpleITK Image, it is just assigned to image.
    All other cases are ignored (nothing calculated).
    Equal approach is used for assignment of mask using MaskFilePath.

    If resampling is enabled, both image and mask are resampled and cropped to the tumormask (with additional
    padding as specified in padDistance) after assignment of image and mask.
    """
    if isinstance(ImageFilePath, basestring) and os.path.exists(ImageFilePath):
      image = sitk.ReadImage(ImageFilePath)
    elif isinstance(ImageFilePath, sitk.SimpleITK.Image):
      image = ImageFilePath
    else:
      self.logger.warning('Error reading image Filepath or SimpleITK object')
      if self.kwargs['verbose']: print "Error reading image Filepath or SimpleITK object"
      image = None

    if isinstance(MaskFilePath, basestring) and os.path.exists(MaskFilePath):
      mask = sitk.ReadImage(MaskFilePath)
    elif isinstance(ImageFilePath, sitk.SimpleITK.Image):
      mask = MaskFilePath
    else:
      self.logger.warning('Error reading mask Filepath or SimpleITK object')
      if self.kwargs['verbose']: print "Error reading mask Filepath or SimpleITK object"
      mask = None

    if self.kwargs['interpolator'] is not None and self.kwargs['resampledPixelSpacing'] is not None:
      image, mask = imageoperations.resampleImage(image, mask,
                                                  self.kwargs['resampledPixelSpacing'],
                                                  self.kwargs['interpolator'],
                                                  self.kwargs['label'],
                                                  self.kwargs['padDistance'])

    return image, mask

  def computeFeatures(self, image, mask, inputImageName, **kwargs):
    """
    Compute signature using image, mask, \*\*kwargs settings
    This function computes the signature for just the passed image (original or filtered),
    does not preprocess or apply a filter to the passed image.
    Features / Classes to use for calculation of signature are defined in self.enabledFeatures.
    see also enableFeaturesByName.
    """
    image, mask = imageoperations.cropToTumorMask(image, mask, self.kwargs['label'])
    featureVector = collections.OrderedDict()

    # Calculate feature classes
    for featureClassName, enabledFeatures in self.enabledFeatures.iteritems():
      # Handle calculation of shape features separately
      if featureClassName == 'shape':
        continue

      if featureClassName in self.getFeatureClassNames():
        featureClass = self.featureClasses[featureClassName](image, mask, **kwargs)

        if enabledFeatures is None or len(enabledFeatures) == 0:
          featureClass.enableAllFeatures()
        else:
          for feature in enabledFeatures:
            featureClass.enableFeatureByName(feature)

        if self.kwargs['verbose']: print "\t\tComputing %s" % (featureClassName)
        featureClass.calculateFeatures()
        for (featureName, featureValue) in featureClass.featureValues.iteritems():
          newFeatureName = "%s_%s_%s" % (inputImageName, featureClassName, featureName)
          featureVector[newFeatureName] = featureValue

    return featureVector

  def generate_original(self, image, mask, **kwargs):
    """
    No filter is applied.

    :return: Yields original image, mask, 'original' and ``kwargs``
    """
    yield image, mask, 'original', kwargs

  def generate_log(self, image, mask, **kwargs):
    """
    Apply Laplacian of Gaussian filter to input image and compute signature for each filtered image.

    Following settings are possible:

    - sigma: List of floats or integers, must be greater than 0. Sigma values to
      use for the filter (determines coarseness).

    N.B. Setting for sigma must be provided. If omitted, no LoG image features are calculated and the function
    will return an empty dictionary.

    Returned filter name reflects LoG settings:
    log-sigma-<sigmaValue>-3D-<featureName>.

    :return: Yields log filtered image, mask, filter name and ``kwargs``
    """
    # Check if size of image is > 4 in all 3D directions (otherwise, LoG filter will fail)
    size = numpy.array(image.GetSize())
    if numpy.min(size) < 4:
      self.logger.warning('Image too small to apply LoG filter, size: %s', size)
      if self.kwargs['verbose']: print 'Image too small to apply LoG filter'
      return

    sigmaValues = kwargs.get('sigma', numpy.arange(5., 0., -.5))

    for sigma in sigmaValues:
      self.logger.debug('Computing LoG with sigma %g', sigma)
      if self.kwargs['verbose']: print "\tComputing LoG with sigma %g" % (sigma)
      logImage = imageoperations.applyLoG(image, sigmaValue=sigma)
      if logImage is not None:
        inputImageName = "log-sigma-%s-mm-3D" % (str(sigma).replace('.', '-'))
        yield logImage, mask, inputImageName, kwargs
      else:
        # No log record needed here, this is handled by logging in imageoperations
        if self.kwargs['verbose']: print 'Application of LoG filter failed (sigma: %g)' % (sigma)

  def generate_wavelet(self, image, mask, **kwargs):
    """
    Apply wavelet filter to image and compute signature for each filtered image.

    Following settings are possible:

    - start_level [0]: integer, 0 based level of wavelet which should be used as first set of decompositions
      from which a signature is calculated
    - level [1]: integer, number of levels of wavelet decompositions from which a signature is calculated.
    - wavelet ["coif1"]: string, type of wavelet decomposition. Enumerated value, validated against possible values
      present in the ``pyWavelet.wavelist()``. Current possible values (pywavelet version 0.4.0) (where an
      aditional number is needed, range of values is indicated in []):

      - haar
      - dmey
      - sym[2-20]
      - db[1-20]
      - coif[1-5]
      - bior[1.1, 1.3, 1.5, 2.2, 2.4, 2.6, 2.8, 3.1, 3.3, 3.5, 3.7, 3.9, 4.4, 5.5, 6.8]
      - rbio[1.1, 1.3, 1.5, 2.2, 2.4, 2.6, 2.8, 3.1, 3.3, 3.5, 3.7, 3.9, 4.4, 5.5, 6.8]

    Returned filter name reflects wavelet type:
    wavelet[level]-<decompositionName>-<featureName>

    N.B. only levels greater than the first level are entered into the name.

    :return: Yields wavelet filtered image, mask, filter name and ``kwargs``
    """
    waveletArgs = {}
    waveletArgs['wavelet'] = kwargs.get('wavelet', 'coif1')
    waveletArgs['level'] = kwargs.get('level', 1)
    waveletArgs['start_level'] = kwargs.get('start_level', 0)

    approx, ret = imageoperations.swt3(image, **waveletArgs)

    for idx, wl in enumerate(ret, start=1):
      for decompositionName, decompositionImage in wl.items():
        self.logger.debug('Computing Wavelet %s', decompositionName)
        if self.kwargs['verbose']: print "\tComputing Wavelet %s" % (decompositionName)

        if idx == 1:
          inputImageName = 'wavelet-%s' % (decompositionName)
        else:
          inputImageName = 'wavelet%s-%s' % (idx, decompositionName)
        yield decompositionImage, mask, inputImageName, kwargs

    if len(ret) == 1:
      inputImageName = 'wavelet-LLL'
    else:
      inputImageName = 'wavelet%s-LLL' % (len(ret))
    yield approx, mask, inputImageName, kwargs

  def generate_square(self, image, mask, **kwargs):
    """
    Computes the square of the image intensities.

    Resulting values are rescaled on the range of the initial original image.

    :return: Yields square filtered image, mask, 'square' and ``kwargs``
    """
    squareImage = imageoperations.applySquare(image)
    yield squareImage, mask, 'square', kwargs

  def generate_squareroot(self, image, mask, **kwargs):
    """
    Computes the square root of the absolute value of image intensities.

    Resulting values are rescaled on the range of the initial original image and negative intensities are made
    negative in resultant filtered image.

    :return: Yields square root filtered image, mask, 'squareroot' and ``kwargs``
    """
    squarerootImage = imageoperations.applySquareRoot(image)
    yield squarerootImage, mask, 'squareroot', kwargs

  def generate_logarithm(self, image, mask, **kwargs):
    """
    Computes the logarithm of the absolute value of the original image + 1.

    Resulting values are rescaled on the range of the initial original image and negative intensities are made
    negative in resultant filtered image.

    :return: Yields logarithm filtered image, mask, 'logarithm' and ``kwargs``
    """
    logarithmImage = imageoperations.applyLogarithm(image)
    yield logarithmImage, mask, 'logarithm', kwargs

  def generate_exponential(self, image, mask, **kwargs):
    """
    Computes the exponential of the original image.

    Resulting values are rescaled on the range of the initial original image.

    :return: Yields exponential filtered image, mask, 'exponential' and ``kwargs``
    """
    exponentialImage = imageoperations.applyExponential(image)
    yield exponentialImage, mask, 'exponential', kwargs

  def getInputImageTypes(self):
    """
    Returns a list of possible input image types.
    """
    return [member[9:] for member in dir(self) if member.startswith('generate_')]

  def getFeatureClasses(self):
    """
    Iterates over all modules of the radiomics package using pkgutil and subsequently imports those modules.

    Return a dictionary of all modules containing featureClasses, with modulename as key, abstract
    class object of the featureClass as value. Assumes only one featureClass per module

    This is achieved by inspect.getmembers. Modules are added if it contains a memeber that is a class,
    with name starting with 'Radiomics' and is inherited from :py:class:`radiomics.base.RadiomicsFeaturesBase`.
    """
    featureClasses = {}
    for _, mod, _ in pkgutil.iter_modules(radiomics.__path__):
      __import__('radiomics.' + mod)
      attributes = inspect.getmembers(eval('radiomics.' + mod), inspect.isclass)
      for a in attributes:
        if a[0].startswith('Radiomics'):
          if radiomics.base.RadiomicsFeaturesBase in inspect.getmro(a[1])[1:]:
            featureClasses[mod] = a[1]

    return featureClasses

  def getFeatureClassNames(self):
    return self.featureClasses.keys()

  def getFeaturesNames(self, featureClassName):
    """
    Returns a list of all possible features in provided featureClass
    """
    return self.featureClasses[featureClassName].getFeatureNames()
