.. _radiomics-customization-label:

==========================
Customizing the Extraction
==========================

----------------------
Types of Customization
----------------------

There are 4 ways in which the feature extraction can be customized in PyRadiomics:

1. Specifying which image types (original/derived) to use to extract features from
2. Specifying which feature(class) to extract
3. Specifying settings, which control the pre processing and customize the behaviour of enabled filters and feature
   classes.
4. Specifying the voxel-based specific settings, which are only needed when using PyRadiomics to generate feature maps

.. warning::
    At initialization of the feature extractor or an individual feature class, settings can be provided as keyword
    arguments in ``**kwargs``. These consist *only* of type 3 parameters (setting). Parameters of type 1 (image type)
    and 2 (feature class) can only provided at initialization when using the parameter file. When the parameter file is
    not used, or when these parameters have to be changed after initialization, use the respective function calls.

.. _radiomics-image-types-label:

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

.. note::
    This type of customization can be included in the :ref:`radiomics-parameter-file-label` using key ``imageType``.

.. note::
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
- Gradient: Returns the magnitude of the local gradient. See also :py:func:`~radiomics.imageoperations.getGradientImage`
- LocalBinaryPattern2D: Computes the Local Binary Pattern in a by-slice operation (2D).
  See also :py:func:`~radiomics.imageoperations.getLBP2DImage`
- LocalBinaryPattern3D: Computes the Local Binary Pattern in 3D using spherical harmonics.
  See also :py:func:`~radiomics.imageoperations.getLBP3DImage`


.. _radiomics-feature-classes-label:

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

.. note::
    This type of customization can be included in the :ref:`radiomics-parameter-file-label` using key ``featureClass``.

.. note::
    By default, all feature classes and all features are enabled.

Currently available feature classes are:

- firstorder
- shape
- glcm
- glrlm
- glszm
- gldm
- ngtdm

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
    This type of customization can be included in the :ref:`radiomics-parameter-file-label` using key ``setting``.

.. note::
    When using the feature classes directly, feature class level settings can be customized by providing them as keyword
    arguments at initialization of the feature class.

Below are the settings that control the behaviour of the extraction, ordered per level and category. Each setting is
listed as it's unique, case sensitive name, followed by it's default value in brackets. After the default value is the
documentation on the type of the value and what the setting controls.


Feature Extractor Level
+++++++++++++++++++++++

*Image Normalization*

- ``normalize`` [False]: Boolean, set to True to enable normalizing of the image before any resampling. See also
  :py:func:`~radiomics.imageoperations.normalizeImage`.
- ``normalizeScale`` [1]: Float, > 0, determines the scale after normalizing the image. If normalizing is disabled, this
  has no effect.
- ``removeOutliers`` [None]: Float, > 0, defines the outliers to remove from the image. An outlier is defined as values
  that differ more than :math:`n\sigma_x` from the mean, where :math:`n>0` and equal to the value of this setting. If
  this parameter is omitted (providing it without a value (i.e. None) in the parameter file will throw an error), no
  outliers are removed. If normalizing is disabled, this has no effect. See also
  :py:func:`~radiomics.imageoperations.normalizeImage`.

*Resampling the image/mask*

- ``resampledPixelSpacing`` [None]: List of 3 floats (>= 0), sets the size of the voxel in (x, y, z) plane when resampling.
  A value of 0 is replaced with the spacing for that dimension as it is in the original (non-resampled) image or mask. For example, to perform only in-plane resampling, the x and y values alone should be edited (e.g.: [2,2,0]). In-plane resolution is always relative to image acquisition plane (i.e. axial, coronal or sagittal).
- ``interpolator`` [sitkBSpline]: SimpleITK constant or string name thereof, sets interpolator to use for resampling. The choice of ``interpolator`` is only applied to resampling images, while ``sitkNearestNeighbor`` is always used for resampling masks in order to preserve label values.
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

- ``padDistance`` [5]: Integer, :math:`\geq 0`, set the number of voxels pad cropped tumor volume with during resampling.
  Padding occurs in new feature space and is done on all faces, i.e. size increases in x, y and z direction by
  2*padDistance. Padding is needed for some filters (e.g. LoG). Value of padded voxels are set to original gray level
  intensity, padding does not exceed original image boundaries. **N.B. After application of filters image is cropped
  again without padding.**

.. note::
    Resampling is disabled when either `resampledPixelSpacing` or `interpolator` is set to `None`

*Pre-Cropping*

- ``preCrop`` [False]: Boolean, if true and resampling is disabled, crops the image onto the bounding box with additional
  padding as specified in ``padDistance``. Similar to padding after resampling, padding does not exceed original image
  bounds after pre-cropping. Setting ``preCrop`` to true speeds up extraction and makes it less memory intensive,
  especially in the case of large images with only small ROIs.

.. note::
  Because image and mask are also cropped onto the bounding box before they are passed to the feature classes,
  pre-crop is only beneficial when filters are enabled.

*Resegmentation*

- ``resegmentRange`` [None]: List of 1 or 2 floats, specifies the lower and and optionally upper threshold,
  respectively. Segmented voxels outside this range are removed from the mask prior to feature calculation. When the
  value is None (default), no resegmentation is performed. Resegemented size is checked (using parameter
  ``minimumROISize``, default 1) and upon fail, an error is logged and extraction is skipped for this case.
- ``resegmentMode`` ['absolute']: string, specifying the method to use for defining the resegmentation thresholds:

  - 'absolute': The resegmentRange values are treated as absolute values, i.e. used directly to perform resegmentation.
  - 'relative': The resegmentRange values are treated as relative to the maximum in the ROI, i.e. the actual threshold
    used is defined as :math:`\text{threshold} = \text{value} * X_{max}`.
  - 'sigma': The resegmentRange values indicate a distance from the mean of the ROI in standard deviations. E.g. to
    exclude outliers farther from the mean than 3 sigma, specify mode 'sigma' and range [-3, 3]. Threshold is defined as
    :math:`\text{threshold} = \mu + \text{value} * \sigma`.

- ``resegmentShape`` [False]: Boolean, if set to True, the resegmented mask is also used for shape calculation. If set
  to False (default), only first order and texture classes are calculated using the resegmented mask (known in IBSI as
  the intensity mask). Shape is then calculated using the mask after any optional resampling and corrections (known in
  IBSI as the morphologic mask).

*Mask validation*

- ``minimumROIDimensions`` [2]: Integer, range 1-3, specifies the minimum dimensions (1D, 2D or 3D, respectively).
  Single-voxel segmentations are always excluded.
- ``minimumROISize`` [None]: Integer, > 0, specifies the minimum number of voxels required. Test is skipped
  if this parameter is omitted (specifying it as None in the parameter file will throw an error).
- ``geometryTolerance`` [None]: Float, determines the tolarance used by SimpleITK to compare origin, direction and spacing
  between image and mask. Affects the fist step in :py:func:`~radiomics.imageoperations.checkMask`. If set to ``None``,
  PyRadiomics will use SimpleITK default (1e-16).
- ``correctMask`` [False]: Boolean, if set to true, PyRadiomics will attempt to resample the mask to the image geometry when
  the first step in :py:func:`~radiomics.imageoperations.checkMask` fails. This uses a nearest neighbor interpolator.
  Mask check will still fail if the ROI defined in the mask includes areas outside of the image physical space.

*Miscellaneous*

- ``additionalInfo`` [True]: boolean, set to False to disable inclusion of additional information on the extraction in the
  output. See also :py:func:`~radiomics.featureextractor.RadiomicsFeaturesExtractor.addProvenance()`.

Filter Level
++++++++++++

*Laplacian of Gaussian settings*

- ``sigma``: List of floats or integers, must be greater than 0. Sigma values to use for the filter (determines coarseness).

.. warning::
    Setting for sigma must be provided if LoG filter is enabled. If omitted, no LoG image features are calculated and
    the function will return an empty dictionary.

*Wavelet settings*

- ``start_level`` [0]: integer, 0 based level of wavelet which should be used as first set of decompositions
  from which a signature is calculated
- ``level`` [1]: integer, number of levels of wavelet decompositions from which a signature is calculated.
- ``wavelet`` ["coif1"]: string, type of wavelet decomposition. Enumerated value, validated against possible values
  present in the ``pyWavelet.wavelist()``. Current possible values (pywavelet version 0.4.0) (where an
  aditional number is needed, range of values is indicated in []):

    - haar
    - dmey
    - sym[2-20]
    - db[1-20]
    - coif[1-5]
    - bior[1.1, 1.3, 1.5, 2.2, 2.4, 2.6, 2.8, 3.1, 3.3, 3.5, 3.7, 3.9, 4.4, 5.5, 6.8]
    - rbio[1.1, 1.3, 1.5, 2.2, 2.4, 2.6, 2.8, 3.1, 3.3, 3.5, 3.7, 3.9, 4.4, 5.5, 6.8]

*Gradient settings*

- ``gradientUseSpacing`` [True]: Boolean, if true, image spacing is taken into account when computing the
  gradient magnitude in the image.

*Local Binary Pattern 2D*

- ``lbp2DRadius`` [1]: Float, > 0, specifies the radius in which the neighbours should be sampled
- ``lbp2DSamples`` [9]: Integer, :math:`\geq 1`, specifies the number of samples to use
- ``lbp2DMethod`` ['uniform']: String, specifies the method for computing the LBP to use.

.. warning::
  Requires package ``skimage`` to function.

*Local Binary Pattern 3D*

- ``lbp3DLevels`` [2]: integer, :math:`\geq 1`, specifies the the number of levels in spherical harmonics to use.
- ``lbp3DIcosphereRadius`` [1]: Float, > 0, specifies the radius in which the neighbours should be sampled
- ``lbp3DIcosphereSubdivision`` [1]: Integer, :math:`\geq 0`, specifies the number of subdivisions to apply in the
  icosphere

.. warning::
  Requires package ``scipy`` and ``trimesh`` to function.

Feature Class Level
+++++++++++++++++++

- ``Label`` [1]: Integer, label value of Region of Interest (ROI) in labelmap.

*Image discretization*

- ``binWidth`` [25]: Float, > 0, size of the bins when making a histogram and for discretization of the image gray level.
- ``binCount`` [None]: integer, > 0, specifies the number of bins to create. The width of the bin is
  then determined by the range in the ROI. No definitive evidence is available on which method of discretization is
  superior, we advise a fixed bin width. See more :ref:`here <radiomics_fixed_bin_width>`.

*Forced 2D extraction*

- ``force2D`` [False]: Boolean, set to true to force a by slice texture calculation. Dimension that identifies
  the 'slice' can be defined in ``force2Ddimension``. If input ROI is already a 2D ROI, features are automatically
  extracted in 2D.
- ``force2Ddimension`` [0]: int, range 0-2. Specifies the 'slice' dimension for a by-slice feature extraction. A value of 0 represents the native acquisition plane for the images (usually axial for CT and axial, coronal or sagittal for MRI).
  Similarly, 1 identifies the out-of plane y dimension (e.g. coronal plane for an axial image) and 2 the out-of-plane x dimension (e.g. sagittal plane for an acial image). if
  ``force2D`` is set to False, this parameter has no effect.

*Texture matrix weighting*

- ``weightingNorm`` [None]: string, indicates which norm should be used when applying distance weighting.
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

*Distance to neighbour*

- ``distances`` [[1]]: List of integers. This specifies the distances between the center voxel and the neighbor, for which
  angles should be generated.

.. note::

    This only affects the GLCM and NGTDM feature classes. The GLSZM and GLRLM feature classes use a fixed distance of 1
    (infinity norm) to define neighbours.

Feature Class Specific Settings
+++++++++++++++++++++++++++++++

*First Order*

- ``voxelArrayShift`` [0]: Integer, This amount is added to the gray level intensity in features Energy, Total Energy and
  RMS, this is to prevent negative values. *If using CT data, or data normalized with mean 0, consider setting this
  parameter to a fixed value (e.g. 2000) that ensures non-negative numbers in the image. Bear in mind however, that
  the larger the value, the larger the volume confounding effect will be.*

*GLCM*

- ``symmetricalGLCM`` [True]: boolean, indicates whether co-occurrences should be assessed in two directions per angle,
  which results in a symmetrical matrix, with equal distributions for :math:`i` and :math:`j`. A symmetrical matrix
  corresponds to the GLCM as defined by Haralick et al.

*GLDM*

- ``gldm_a`` [0]: float, :math:`\alpha` cutoff value for dependence. A neighbouring voxel with gray level :math:`j` is
  considered dependent on center voxel with gray level :math:`i` if :math:`|i-j|\le\alpha`

.. _radiomics-voxel-settings-label:

Voxel-based specific settings
#############################

When using PyRadiomics to generate feature maps, additional customization options exist. These control the neighborhood
around each voxel that is used for calculation (kernel) and what the background value should be, i.e. the value of
voxels for which there is no calculated value.

- ``kernelRadius`` [1]: integer, specifies the size of the kernel to use as the radius from the center voxel. Therefore
  the actual size is ``2 * kernelRadius + 1``. E.g. a value of 1 yields a 3x3x3 kernel, a value of 2 5x5x5, etc. In case
  of 2D extraction, the generated kernel will also be a 2D shape (square instead of cube).

- ``maskedKernel`` [True]: boolean, specifies whether to mask the kernel with the overall mask. If True, only voxels in
  the kernel that are also segmented in the mask are used for calculation. Otherwise, all voxels inside the kernel are
  used. Moreover, gray value discretization is performed over the ROI if the setting is set to True, and over the entire
  image if False.

- ``initValue`` [0]: float, value to use for voxels outside the ROI, or voxels where calculation failed. If set to
  ``nan``, 3D slicer will treat them as transparent voxels

- ``voxelBatch`` [-1]: integer > 0, this value controls the maximum number of voxels that are calculated in one batch.
  Larger batches mean less loops in Python and therefore a quicker extraction, but do require more memory. This setting
  allows the user to compromise between extraction speed and memory usage.
  When providing this setting, the value is constrained to be > 0, only by not providing it is the default value of -1
  used (which means: all voxels in 1 batch).

.. _radiomics-parameter-file-label:

--------------
Parameter File
--------------

All 4 categories of customization can be provided in a single yaml or JSON structured text file, which can be provided
in an optional argument (``--param``) when running pyradiomics from the command line. In interactive mode, it can be
provided during initialization of the :ref:`feature extractor <radiomics-featureextractor-label>`, or using
:py:func:`~radiomics.featureextractor.RadiomicsFeaturesExtractor.loadParams` after initialization. This removes the need
to hard code a customized extraction in a python script through use of functions described above. Additionally, this
also makes it more easy to share settings for customized extractions. We encourage users to share their parameter files
in the PyRadiomics repository. See :ref:`radiomics-submit-parameter-file-label` for more information on how to submit
your parameter file.

.. note::
    For an extensive list of possible settings, see :ref:`Image Types<radiomics-image-types-label>`,
    :ref:`Feature Classes<radiomics-feature-classes-label>` and :ref:`Settings<radiomics-settings-label>`,
    which can be provided in the parameter file using key ``imageType``, ``featureClass`` and ``setting``, respectively.

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
4. **voxelSetting:** Settings used to control the voxel-based specific settings. E.g. the size of the kernel used and
   the background value in the parameter maps.

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
