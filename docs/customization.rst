.. _radiomics-customization-label:

==========================
Customizing the Extraction
==========================

----------------------
Types of Customization
----------------------

There are 3 ways in which the feature extraction can be customized in PyRadiomics:

1. Specifying which image types (original/derived) to use to extract features from
2. Specifying which feature(class) to extract
3. Specifying settings, which control the pre processing and customize the behaviour of enabled filters and feature
   classes.

.. warning::
    At initialization of the feature extractor or an individual feature class, settings can be provided as keyword
    arguments in ``**kwargs``. These consist *only* of type 3 parameters (setting). Parameters of type 1 (image type)
    and 2 (feature class) can only provided at initialization when using the parameter file. When the parameter file is
    not used, or when these parameters have to be changed after initialization, use the respective function calls.

Image Types
###########

These are the image types (either the original image or derived images using filters) that can be used to extract
features from. The image types that are available are determined dynamically (all are functions in
``imageoperations.py`` that fit the :ref:`signature <radiomics-developers-filter>` of an image type).

The enabled types are stored in the ``_enabledImageTypes`` dictionary in the feature extractor class instance and can be
changed using the functions :py:func:`~radiomics.featureextractor.RadiomicsFeaturesExtractor.enableAllImageTypes`,
:py:func:`~radiomics.featureextractor.RadiomicsFeaturesExtractor.disableAllImageTypes`,
:py:func:`~radiomics.featureextractor.RadiomicsFeaturesExtractor.enableImageTypeByName` and
:py:func:`~radiomics.featureextractor.RadiomicsFeaturesExtractor.enableImageTypes`. Moreover, custom settings can be
provided for each enabled image type, which will then only be applied during feature extraction for that image type.
Please note that this will only work for settings that are applied at or after any filter is applied (i.e. not at the
feature extractor level).

By default, only the "Original" image type is enabled.

Currently available image types are:

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

Enabled Features
################

These are the features that are extracted from each (original and/or derived) image type. The available features are
determined dynamically, and are ordered in feature classes. For more information on the signature used to identify
features and feature classes, see the :ref:`radiomics-developers` section.

The enable features are stored in the ``_enabledFeatures`` dictionary in the feature extractor class instance and can be
changed using the functions :py:func:`~radiomics.featureextractor.RadiomicsFeaturesExtractor.enableAllFeatures`,
:py:func:`~radiomics.featureextractor.RadiomicsFeaturesExtractor.disableAllFeatures`,
:py:func:`~radiomics.featureextractor.RadiomicsFeaturesExtractor.enableFeatureClassByName` and
:py:func:`~radiomics.featureextractor.RadiomicsFeaturesExtractor.enableFeaturesByName`. Each key-value pair in the
dictionary represents one enabled feature class with the feature class name as the key and a list of enabled feature
names as value. If the value is ``None`` or an empty list, all features in that class are enabled. Otherwise only the
features specified.

By default, all feature classes and all features are enabled.

Currently available feature classes are:

- firstorder
- shape
- glcm
- glrlm
- glszm

An individual feature can be enabled by submitting the feature name as defined in the unique part of the function
signature (e.g. the First Order feature defined by ``get10PercentileFeatureValue()`` is enabled by specifying
``{firstorder: ['10Percentile']}``). Function signatures for all features are available in the
:ref:`radiomics-features-label` section.

.. _radiomics-settings-label:

Settings
########

Besides customizing what to extract (image types, features), PyRadiomics exposes various settings customizing how the
features are extracted. These settings operate at different levels. E.g. resampling is done just after the images are
loaded (in the feature extractor), so settings controlling the resampling operate only on the feature extractor level.
Settings are stored in the ``setttings`` dictionary in the feature extractor class instance, where the key is the case
sensitive setting name. Custom settings are provided as keyword arguments at initialization of the feature extractor
(with the setting name as keyword and value as the argument value, e.g. ``binWidth=25``), or by interacting directly
with the ``settings`` dictionary.

.. note::
    When using the feature classes directly, feature class level settings can be customized by providing them as keyword
    arguments at initialization of the feature class.

Below are the settings that control the behaviour of the extraction, ordered per level and category. Each setting is
listed as it's unique, case sensitive name, followed by it's default value in brackets. After the default value is the
documentation on the type of the value and what the setting controls.


Feature Extractor Level
+++++++++++++++++++++++

*Image Normalization*

- normalize [False]: Boolean, set to True to enable normalizing of the image before any resampling. See also
  :py:func:`~radiomics.imageoperations.normalizeImage`.
- normalizeScale [1]: Float, > 0, determines the scale after normalizing the image. If normalizing is disabled, this
  has no effect.
- removeOutliers [None]: Float, > 0, defines the outliers to remove from the image. An outlier is defined as values
  that differ more than :math:`n\sigma_x` from the mean, where :math:`n>0` and equal to the value of this setting. If
  this parameter is omitted (providing it without a value (i.e. None) in the parameter file will throw an error), no
  outliers are removed. If normalizing is disabled, this has no effect. See also
  :py:func:`~radiomics.imageoperations.normalizeImage`.

*Resampling the image*

- resampledPixelSpacing [None]: List of 3 floats (> 0), sets the size of the voxel in (x, y, z) plane when resampling.
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

- padDistance [5]: Integer, :math:`\geq 0`, set the number of voxels pad cropped tumor volume with during resampling.
  Padding occurs in new feature space and is done on all faces, i.e. size increases in x, y and z direction by
  2*padDistance. Padding is needed for some filters (e.g. LoG). Value of padded voxels are set to original gray level
  intensity, padding does not exceed original image boundaries. **N.B. After application of filters image is cropped
  again without padding.**

.. note::
    Resampling is disabled when either `resampledPixelSpacing` or `interpolator` is set to `None`

*Mask validation*

- minimumROIDimensions [1]: Integer, range 1-3, specifies the minimum dimensions (1D, 2D or 3D, respectively).
  Single-voxel segmentations are always excluded.
- minimumROISize [None]: Integer, > 0, specifies the minimum number of voxels required. Test is skipped
  if this parameter is omitted (specifying it as None in the parameter file will throw an error).
- geometryTolerance [None]: Float, determines the tolarance used by SimpleITK to compare origin, direction and spacing
  between image and mask. Affects the fist step in :py:func:`~radiomics.imageoperations.checkMask`. If set to ``None``,
  PyRadiomics will use SimpleITK default (1e-16).
- correctMask [False]: Boolean, if set to true, PyRadiomics will attempt to resample the mask to the image geometry when
  the first step in :py:func:`~radiomics.imageoperations.checkMask` fails. This uses a nearest neighbor interpolator.
  Mask check will still fail if the ROI defined in the mask includes areas outside of the image physical space.

*Miscellaneous*

- enableCExtensions [True]: Boolean, set to False to force calculation to full-python mode. See also
  :py:func:`~radiomics.enableCExtensions()`.
- additionalInfo [True]: boolean, set to False to disable inclusion of additional information on the extraction in the
  output. See also :py:func:`~radiomics.featureextractor.RadiomicsFeaturesExtractor.addProvenance()`.

Filter Level
++++++++++++

*Laplacian of Gaussian settings*

- sigma: List of floats or integers, must be greater than 0. Sigma values to use for the filter (determines coarseness).

.. warning::
    Setting for sigma must be provided if LoG filter is enabled. If omitted, no LoG image features are calculated and
    the function will return an empty dictionary.

*Wavelet settings*

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

Feature Class Level
+++++++++++++++++++

*Image discretization*

- binWidth [25]: Float, > 0, size of the bins when making a histogram and for discretization of the image gray level.

*Forced 2D extraction*

- force2D [False]: Boolean, set to true to force a by slice texture calculation. Dimension that identifies
  the 'slice' can be defined in ``force2Ddimension``. If input ROI is already a 2D ROI, features are automatically
  extracted in 2D. See also :py:func:`~radiomics.imageoperations.generateAngles`
- force2Ddimension [0]: int, range 0-2. Specifies the 'slice' dimension for a by-slice feature extraction. Value 0
  identifies the 'z' dimension (axial plane feature extraction), and features will be extracted from the xy plane.
  Similarly, 1 identifies the y dimension (coronal plane) and 2 the x dimension (saggital plane). if
  ``force2Dextraction`` is set to False, this parameter has no effect. See also
  :py:func:`~radiomics.imageoperations.generateAngles`

*Texture matrix weighting*

- weightingNorm [None]: string, indicates which norm should be used when applying distance weighting.
  Enumerated setting, possible values:

    - 'manhattan': first order norm
    - 'euclidean': second order norm
    - 'infinity': infinity norm.
    - 'no_weighting': GLCMs are weighted by factor 1 and summed
    - None: Applies no weighting, mean of values calculated on separate matrices is returned.

  In case of other values, an warning is logged and option 'no_weighting' is used.

.. note::
    This only affects the GLCM and GLRLM feature classes. Moreover, weighting is applied differently in those classes.
    For more information on how weighting is applied, see the documentation on :ref:`GLCM <radiomics-glcm-label>` and
    :ref:`GLRLM <radiomics-glszm-label>`.

Feature Class Specific Settings
+++++++++++++++++++++++++++++++

*First Order*

- voxelArrayShift [0]: Integer, This amount is added to the gray level intensity in features Energy, Total Energy and
  RMS, this is to prevent negative values. *If using CT data, or data normalized with mean 0, consider setting this
  parameter to a fixed value (e.g. 2000) that ensures non-negative numbers in the image. Bear in mind however, that
  the larger the value, the larger the volume confounding effect will be.*

*GLCM*

- distances [[1]]: List of integers. This specifies the distances between the center voxel and the neighbor, for which
  angles should be generated. See also :py:func:`~radiomics.imageoperations.generateAngles`

.. _radiomics-parameter-file-label:

--------------
Parameter File
--------------

All 3 categories of customization can be provided in a single yaml-structured text file, which can be provided in an
optional argument (``--param``) when running pyradiomics from the command line. In interactive mode, it can be provided
during initialization of the :ref:`feature extractor <radiomics-featureextractor-label>`, or using
:py:func:`~radiomics.featureextractor.RadiomicsFeaturesExtractor.loadParams` after initialization. This removes the need
to hard code a customized extraction in a python script through use of functions described above. Additionally, this
also makes it more easy to share settings for customized extractions. We encourage users to share their parameter files
in the PyRadiomics repository. See :ref:`radiomics-submit-parameter-file-label` for more information on how to submit
your parameter file.

.. note::
    Examples of the parameter file are provided in the ``pyradiomics/examples/exampleSettings`` folder.

The paramsFile is written according to the YAML-convention (www.yaml.org) and is checked by the code for
consistency. Only one yaml document per file is allowed. Parameters must be grouped by customization category as mentioned
above. This is reflected in the structure of the document as follows::

    <Customization Category>:
      <Setting Name>: <value>
      ...
    <Customization Categort>:
      ...

Blank lines may be inserted to increase readability, these are ignored by the parser. Additional comments are also
possible, these are preceded by an '#' and can be inserted on a blank line, or on a line containing parameters::

    # This is a line containing only comments
    setting: # This is a comment placed after the declaration of the 'setting' category.

Any keyword, such as a customization category or setting name may only be mentioned once. Multiple instances do not
raise an error, but only the last one encountered is used.

The three setting types are named as follows:

1. **imageType:** image type to calculate features on. <value> is custom kwarg settings (dictionary). if <value>
   is an empty dictionary ('{}'), no custom settings are added for this input image.
2. **featureClass:** Feature class to enable, <value> is list of strings representing enabled features. If no
   <value> is specified or <value> is an empty list ('[]'), all features for this class are enabled.
3. **setting:** Setting to use for pre processing and class specific settings. if no <value> is specified, the value for
   this setting is set to None.

Example::

    # This is a non-active comment on a separate line
    imageType:
        Original: {}
        LoG: {'sigma' : [1.0, 3.0]}  # This is a non active comment on a line with active code preceding it.
        Wavelet:
            binWidth: 10

    featureClass:
        glcm:
        glrlm: []
        firstorder: ['Mean',
                     'StandardDeviation']
        shape:
            - Volume
            - SurfaceArea

    setting:
        binWidth: 25
        resampledPixelSpacing:

In this example, 3 image types are enabled ("Original", "LoG" (Laplacian of Gaussian) and "Wavelet"), with custom
settings specified for "LoG" ("sigma") and "Wavelet" ("binWidth"). Note that the manner of specifying the custom
settings for "LoG" and "Wavelet" is equivalent.

Next, 4 feature classes are defined. "glcm" and "glrlm" are both enabled with all possible features in the respective
class, whereas only "Mean" and "StandardDeviation" are enabled for "firstorder", and only "Volume" and "SurfaceArea" for
shape. Note that the manner of specifying individual features for "firstorder" and "shape" is equivalent.

Finally, 2 settings are specified: "binWidth", whose value has been set to 25 (but will be set to 10 during extraction
of "Wavelet" derived features), and "resampledPixelSpacing", where no value is provided, which is equivalent to a
python "None" value.

.. note::
    - settings not specified in parameters are set to their default value.
    - enabledFeatures are replaced by those in parameters (i.e. only specified features/classes are enabled. If the
      'featureClass' customization type is omitted, all feature classes and features are enabled.
    - ImageTypes are replaced by those in parameters (i.e. only specified types are used to extract features from. If
      the 'inputImage' customization type is omitted, only "Original" image type is used for feature extraction, with no
      additional custom settings.
