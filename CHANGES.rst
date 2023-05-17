=============
Release Notes
=============

------------
Next Release
------------

-----------------
PyRadiomics 3.1.0
-----------------

Bug Fixes
#########

- Fix bug to allow a full-mask. 
  (`#660 <https://github.com/AIM-Harvard/pyradiomics/pull/660>`_)
- Fix label assignment in batch-processing example.
  (`#756 <https://github.com/AIM-Harvard/pyradiomics/pull/756>`_)
- Fix bug in label_channel selection for voxel-based extraction.
  (`54a37822 <https://github.com/AIM-Harvard/pyradiomics/commit/54a37822>`_)
- Force label and label_channel datatype to int.
  (`efb9756e <https://github.com/AIM-Harvard/pyradiomics/commit/efb9756e>`_)

Testing / Continuous Integration
################################

- Remove Python < 3.7 support, update support to Python versions 3.7, 3.8 and 3.9.
  (`#825 <https://github.com/AIM-Harvard/pyradiomics/pull/825>`_)
- Switch to pytest instead of nose and nose parameterized.
  (`#825 <https://github.com/AIM-Harvard/pyradiomics/pull/825>`_)
- Move MacOS CI to CircleCI, remove TravisCI.
  (`#825 <https://github.com/AIM-Harvard/pyradiomics/pull/825>`_)


Documentation
#############

- Update community link.
  (`038f8a8b <https://github.com/AIM-Harvard/pyradiomics/commit/038f8a8b>`_)
- Update usage.rst to reflect new notebooks folder location.
  (`#734 <https://github.com/AIM-Harvard/pyradiomics/pull/734>`_)
- Update ReadMe to reflect updated CI.
  (`0c53d1dd <https://github.com/AIM-Harvard/pyradiomics/commit/0c53d1dd>`_)
- Refactor project description to rely mainly on pyproject.toml and setup.cfg.
  (`c341b31f <https://github.com/AIM-Harvard/pyradiomics/commit/c341b31f>`_)
- Fix typos

Labs
####

- Expose mask correction and geometry tolerance settings in pyradiomics-dcm CLI
  (`#724 <https://github.com/AIM-Harvard/pyradiomics/pull/724>`_)

-----------------
PyRadiomics 3.0.1
-----------------

Bug Fixes
#########

- Fix bug causing ``IndexError`` when no gray levels are 'empty'.
  (`#592 <https://github.com/AIM-Harvard/pyradiomics/pull/592>`_)
- Fail initialization of feature extractor when the passed parameter file path
  does not point to existing file. (`#587 <https://github.com/AIM-Harvard/pyradiomics/pull/587>`_)
- Fix out-of-range check in GLSZM C calculation.
  (`#635 <https://github.com/AIM-Harvard/pyradiomics/pull/635>`_)
- Fix bug in Travis CI testing (MacOS platform).
  (`#643 <https://github.com/AIM-Harvard/pyradiomics/pull/643>`_,
  `#646 <https://github.com/AIM-Harvard/pyradiomics/pull/646>`_)
- Fix cmake URL and remove python2 support from DockerFiles.
  (`#645 <https://github.com/AIM-Harvard/pyradiomics/pull/645>`_)

Examples
########

- Add example settings for forced-2D extraction in MR.
  (`#613 <https://github.com/AIM-Harvard/pyradiomics/pull/613>`_,
  `#644 <https://github.com/AIM-Harvard/pyradiomics/pull/644>`_)

Documentation
#############

- Fix typos in documentation.  (`9d26a6b8 <https://github.com/AIM-Harvard/pyradiomics/commit/9d26a6b8>`_,
  `896682d7 <https://github.com/AIM-Harvard/pyradiomics/commit/896682d7>`_,
  `e100f1d0 <https://github.com/AIM-Harvard/pyradiomics/commit/e100f1d0>`_,
  `#639 <https://github.com/AIM-Harvard/pyradiomics/pull/639>`_)
- Further clarify resampling. (`#599 <https://github.com/AIM-Harvard/pyradiomics/pull/599>`_)

Internal API
############

- Fail gracefully when grayvalues < 1 are encountered in the discretized image.
  (`#602 <https://github.com/AIM-Harvard/pyradiomics/pull/602>`_)
- Add optional progress reporting for voxel-based extraction.
  (`#636 <https://github.com/AIM-Harvard/pyradiomics/pull/636>`_)


---------------
PyRadiomics 3.0
---------------

.. warning::
  As of this release, Python 2.7 testing is removed. Compatibility code such as it is will be left in place, but
  future changes will not be checked for backwards compatibility. Moreover, no pre-built binaries for python 2.7
  will be distributed on PyPi or Conda.
  Finally, some deprecated code is removed (``commandlinebatch.py`` and ``calculateFeatures()``).

Bug Fixes
#########

- Fix broken Conda deployment (`51c5849 <https://github.com/AIM-Harvard/pyradiomics/commit/51c5849>`_)
- Fix error in IBSI mapping (labs/pyradiomics-dcm) (`54d6689 <https://github.com/AIM-Harvard/pyradiomics/commit/54d6689>`_)
- Fix resampling error when spacing is correct, but sizes are different (`ac7458e <https://github.com/AIM-Harvard/pyradiomics/commit/ac7458e>`_)
- Fix label channel selection (`54a3782 <https://github.com/AIM-Harvard/pyradiomics/commit/54a3782>`_)
- Use local scope of settings, preventing race conditions in parallel extraction (`43578f7 <https://github.com/AIM-Harvard/pyradiomics/commit/43578f7>`_)
- Fix resampling for 2D input (`#545 <https://github.com/AIM-Harvard/pyradiomics/pull/545>`_)

Internal API
############

- Update C API to use large datatype for index pointers (`#500 <https://github.com/AIM-Harvard/pyradiomics/pull/500>`_,
  `#501 <https://github.com/AIM-Harvard/pyradiomics/pull/501>`_)
- Update docker CLI to use python 3.6.9 and fix bugs to allow integration with pyradiomics-dcm lab (`#527 <https://github.com/AIM-Harvard/pyradiomics/pull/527>`_)
- Add option to force path to UNIX style paths, even on windows (`3c0708a <https://github.com/AIM-Harvard/pyradiomics/commit/3c0708a>`_)
- Removed deprecated code (`fedaa5e <https://github.com/AIM-Harvard/pyradiomics/commit/fedaa5e>`_)

Testing
#######

- Remove testing and deployment for python 2.7 (`a5a7e61 <https://github.com/AIM-Harvard/pyradiomics/commit/a5a7e61>`_)

Documentation
#############

- Refactor documentation (`#536 <https://github.com/AIM-Harvard/pyradiomics/pull/536>`_)
- Fix various typos/wording
- Clarify use of force2D, and add example settings file (`#558 <https://github.com/AIM-Harvard/pyradiomics/pull/558>`_)

-----------------
PyRadiomics 2.2.0
-----------------

.. warning::
  In this release, the main interface class, :py:mod:`RadiomicsFeaturesExtractor <radiomics.featureextractor>`, was
  renamed to :py:mod:`RadiomicsFeatureExtractor <radiomics.featureextractor>`
  (no 's' between 'Feature' and 'Extractor'). This was done to avoid confusion between the module and class name.
  (`#481 <https://github.com/AIM-Harvard/pyradiomics/pull/481>`_)

New Features
############

- Add 2D shape features (`#442 <https://github.com/AIM-Harvard/pyradiomics/pull/442>`_)
- Expose voxel-based feature extraction on the PyRadiomics command line interface.
  (`#457 <https://github.com/AIM-Harvard/pyradiomics/pull/457>`_)

Labs
####

- Add notebook investigating reproducibility between PyRadiomics and USF tool (ITK-based;
  `#458 <https://github.com/AIM-Harvard/pyradiomics/pull/458>`_)

Bug Fixes
#########

- Flatten array when applying gray value discretization of the entire image (voxel-based, full kernel;
  `f87abcf <https://github.com/AIM-Harvard/pyradiomics/commit/f87abcf>`_)
- Fix incorrect removal of 'empty gray levels' in GLDM and GLRLM (voxel-based;
  `4b18ce2 <https://github.com/AIM-Harvard/pyradiomics/commit/4b18ce2>`_)
- Fix incorrect instantiation of firstorder voxel-based extraction.
  (`81e713a <https://github.com/AIM-Harvard/pyradiomics/commit/81e713a>`_)
- Force cast coefficients to float. Prevents overflow and type errors in feature calculation.
  (`e9d60c7 <https://github.com/AIM-Harvard/pyradiomics/commit/e9d60c7>`_)

Tests
#####

- Removed support and continuous integration for Python 3.4 (not maintained since March 2019). Added support and CI for
  Python 3.7. (`#486 <https://github.com/AIM-Harvard/pyradiomics/pull/486>`_)

Internal API
############

- Update C-extensions:

  - Rewrite C code to work with N-Dimensional input. (`#463 <https://github.com/AIM-Harvard/pyradiomics/pull/463>`_)
  - Add batch-calculation of kernels and vectorized feature calculation to improve voxel-based extraction duration.
    (`#466 <https://github.com/AIM-Harvard/pyradiomics/pull/466>`_)

- Add support for segmentation objects (multi-layer labelmaps;
  `#445 <https://github.com/AIM-Harvard/pyradiomics/pull/445>`_)

- Refactor the commandline interface (`#481 <https://github.com/AIM-Harvard/pyradiomics/pull/481>`_)

  - Extractor instantiated once (resulting in only 1 validation of the parameter file, outside of paralellization loop)
  - Simplify construction of the python generator of the cases that are to be extracted
  - Remove now unnecessary functions

Documentation
#############

- Update documentation (`#446 <https://github.com/AIM-Harvard/pyradiomics/pull/446>`_,
  `690891d <https://github.com/AIM-Harvard/pyradiomics/commit/690891d>`_)
- Fix some rendering errors (`723d868 <https://github.com/AIM-Harvard/pyradiomics/commit/723d868>`_,
  `e3eb427 <https://github.com/AIM-Harvard/pyradiomics/commit/e3eb427>`_)

-----------------
PyRadiomics 2.1.2
-----------------

Labs
####

- Include algorithm details in dcm output. (`f03145b <https://github.com/AIM-Harvard/pyradiomics/commit/f03145b>`_)

-----------------
PyRadiomics 2.1.1
-----------------

New Features
############

- Implement validation of commandline input. (`#433 <https://github.com/AIM-Harvard/pyradiomics/pull/433>`_)
- Implement thread-safe logging for python >= 3.2 (`#441 <https://github.com/AIM-Harvard/pyradiomics/pull/441>`_,
  `d8db675 <https://github.com/AIM-Harvard/pyradiomics/commit/d8db675>`_)

Labs
####

- Add script for using PyRadiomics with DICOM input and output.
  (`#437 <https://github.com/AIM-Harvard/pyradiomics/pull/437>`_)

Bug Fixes
#########

- Fix memory error in calculation of GLCM-MCC. (`167888b <https://github.com/AIM-Harvard/pyradiomics/commit/167888b>`_)
- Fix error in serialization for JSON output. (`9d992fe <https://github.com/AIM-Harvard/pyradiomics/commit/9d992fe>`_)

Tests
#####

- Expand testing to include more parts of PyRadiomics. (`#410 <https://github.com/AIM-Harvard/pyradiomics/pull/410>`_)

Internal API
############

- Force cast the mask to an integer datatype on load. (`#431 <https://github.com/AIM-Harvard/pyradiomics/pull/431>`_)

Dependencies
############

- Fix PyWavelets version to > 0.4.0, <= 1.0.0, due to compilation issue in SlicerRadiomics.
  (`c828b99 <https://github.com/AIM-Harvard/pyradiomics/commit/c828b99>`_,
  `SlicerRadiomics#50 <https://github.com/AIM-Harvard/SlicerRadiomics/issues/50>`_)

-----------------
PyRadiomics 2.1.0
-----------------

Feature Calculation Changes
###########################

- Switch Shape - Volume calculation to a mesh-based instead of a voxel-based one. This also affects all features derived
  from Volume. Original Volume calculation is retained as ``VoxelVolume``. Also switch calculation of maximum diameter
  to mesh based. Only PCA-derived are not affected. (`#427 <https://github.com/AIM-Harvard/pyradiomics/pull/427>`_)

New Features
############

- Add GLCM - Maximal Correlation Coefficient. (`#411 <https://github.com/AIM-Harvard/pyradiomics/pull/411>`_)

New Parameters
##############

- Update resegmentation function, add support for single (lower) threshold and new modes ``relative`` and ``sigma``,
  customizable in parameter ``resegmentMode``. (`#420 <https://github.com/AIM-Harvard/pyradiomics/pull/420>`_)
- Add ``resegmentShape``. Default ``False``, if set to ``True``, the resegmented mask (intensity mask) will also be used
  for shape calculation. Otherwise, the non-resegmented mask (morphological mask) is used for shape.
  (`#428 <https://github.com/AIM-Harvard/pyradiomics/pull/428>`_)

Bug fixes
#########

- Fix bug in dimension checking in ``checkMask``. (`623b836 <https://github.com/AIM-Harvard/pyradiomics/commit/623b836>`_)
- Fix some errors in the testUtils and baseline generation script.
  (`c285c15 <https://github.com/AIM-Harvard/pyradiomics/commit/c285c15>`_)
- Prevent division by 0 in NGTDM - Coarseness. Return 0 instead.
  (`a59861e <https://github.com/AIM-Harvard/pyradiomics/commit/a59861e>`_)
- Remove duplicate key in settings file example. (`828a7ac <https://github.com/AIM-Harvard/pyradiomics/commit/828a7ac>`_)
- Prevent duplicate log entries in parallel batch extraction.
  (`8cedd8f <https://github.com/AIM-Harvard/pyradiomics/commit/8cedd8f>`_)
- Build PyWavelets from source for AppVeyor (Windows) python 3.4 testing. Requires pre-installation of numpy and cython.
  (`6223d35 <https://github.com/AIM-Harvard/pyradiomics/commit/6223d35>`_)

Tests
#####

- Integrate automatic distribution to conda upon release. (`#422 <https://github.com/AIM-Harvard/pyradiomics/pull/422>`_)

Documentation
#############

- Update README and Setup.py with additional classifiers, urls. Update section in README on Docker usage.
  (`0fe737e <https://github.com/AIM-Harvard/pyradiomics/commit/0fe737e>`_)

Internal API
############

- Use ``ValueError`` exceptions when feature extraction pipeline fails (exceptions of individual features)
  (`#420 <https://github.com/AIM-Harvard/pyradiomics/pull/420>`_)
- Update generation and names of general info features (provenance information)
  (`#420 <https://github.com/AIM-Harvard/pyradiomics/pull/420>`_,
  `#426 <https://github.com/AIM-Harvard/pyradiomics/pull/426>`_)
- Rewrite signatures of pre-processing functions to accept all customization arguments in 1 ``**kwargs`` dict.
  Necessary parameters are obtained using ``kwargs.get`` inside the function. Full settings are passed to the function.
  (`#425 <https://github.com/AIM-Harvard/pyradiomics/pull/425>`_)

-----------------
PyRadiomics 2.0.1
-----------------

New Features
############

- Add Center of Mass to general info output. (`#416 <https://github.com/AIM-Harvard/pyradiomics/pull/416>`_)

Bug fixes
#########

- Fix invocation of numpy.histogram when using a fixed bin count.
  (`2a9fd79 <https://github.com/AIM-Harvard/pyradiomics/commit/2a9fd79>`_)
- Fix assignment of x and y pixelspacing in shape (no changes in results).
  (`#404 <https://github.com/AIM-Harvard/pyradiomics/pull/404>`_)
- Fix generation of approximation name (LLL or LL) in wavelet.
  (`#405 <https://github.com/AIM-Harvard/pyradiomics/pull/405>`_)
- Add missing requirements for new filters in Docker CLI file.
  (`#409 <https://github.com/AIM-Harvard/pyradiomics/pull/409>`_)
- Fix memory leak in C extensions. (`#419 <https://github.com/AIM-Harvard/pyradiomics/pull/419>`_)
- Fix Label column parsing in batch processing. (`217a840 <https://github.com/AIM-Harvard/pyradiomics/commit/217a840>`_)

Documentation
#############

- Fix math rendering in GLCM. (`c6a1f21 <https://github.com/AIM-Harvard/pyradiomics/commit/c6a1f21>`_)
- Add reference to GLDM feature class. (`9f9361a <https://github.com/AIM-Harvard/pyradiomics/commit/9f9361a>`_)
- Correct typo in IMC1 and 2 formulas. (`4ba909a <https://github.com/AIM-Harvard/pyradiomics/commit/4ba909a>`_)
- Update warning message in ROI check.  (`1f16b9e <https://github.com/AIM-Harvard/pyradiomics/commit/1f16b9e>`_)
- Update usage section in documentation on command line usage.
  (`fe0e2c3 <https://github.com/AIM-Harvard/pyradiomics/commit/fe0e2c3>`_)

Internal API
############

- Simplify calculation of various GLCM features (no changes in results).
  (`#407 <https://github.com/AIM-Harvard/pyradiomics/pull/407>`_)

-----------------
PyRadiomics 2.0.0
-----------------

Feature Calculation Changes
###########################

- Change calculation of filter coefficients to reflect absolute maximum (take into account negative values).
  (`#319 <https://github.com/AIM-Harvard/pyradiomics/pull/319>`_)
- Mark duplicate features as 'deprecated' and document mathematical proof of the equality.
  (`#321 <https://github.com/AIM-Harvard/pyradiomics/pull/321>`_)
- Fix error in calculation of NGTDM's Complexity and Contrast features
  (`#351 <https://github.com/AIM-Harvard/pyradiomics/pull/351>`_)

New Features
############

- Add ``preCrop``, which crops the image onto the bounding box with an additional padding specified in ``padDistance``.
  This is similar to cropping as performed during resampling and serves to decrease memory consumption and computation
  time. N.B. To ensure calculated values are not changed, a sufficient padding is required when using filters which
  include values outside of ROI (e.g. Wavelet, LoG). (`#317 <https://github.com/AIM-Harvard/pyradiomics/pull/317>`_)
- Add ``skip-nans`` as a commandline argument. If specified, features that compute NaN are removed from the output. In
  batch mode, NaN is replaced by an empty string. (`#318 <https://github.com/AIM-Harvard/pyradiomics/pull/318>`_)
- Add support to configure the feature extractor using a JSON structured string.
  (`#334 <https://github.com/AIM-Harvard/pyradiomics/pull/334>`_)
- Add Gradient Magnitude Filter. (`#356 <https://github.com/AIM-Harvard/pyradiomics/pull/356>`_)
- Add Local Binary Pattern Filter (2D/3D). (`#357 <https://github.com/AIM-Harvard/pyradiomics/pull/357>`_)
- Add support for Gray Value discretization using a fixed bin count.
  (`#386 <https://github.com/AIM-Harvard/pyradiomics/pull/386>`_)

Bug fixes
#########

- Ensure PyKwalify has a log handler, which is needed when parameter file validation fails.
  (`#309 <https://github.com/AIM-Harvard/pyradiomics/pull/309>`_)
- Fix bug in error handling in :py:func:`~radiomics.imageoperations.checkMask` (compatibility issue between python 2 and 3).
- Fix bug in GLCM (incorrect use of ``self.maskArray``) (`#322 <https://github.com/AIM-Harvard/pyradiomics/pull/322>`_)
- Fix bug in error handling during geometry checks of image and mask.
  (`0257217 <https://github.com/AIM-Harvard/pyradiomics/commit/0257217>`_)
- Fix broken continuous testing integration due to unavailability of pip script.
  (`#333 <https://github.com/AIM-Harvard/pyradiomics/pull/333>`_)
- Fix incorrect path separator in example scripts. (`c7c5d2e <https://github.com/AIM-Harvard/pyradiomics/commit/c7c5d2e>`_)
- Fix bug in the calculation of Wavelet. (`#346 <https://github.com/AIM-Harvard/pyradiomics/pull/346>`_)
- Fix machine-precision errors in Eigenvalue calculation (Shape)
  (`#355 <https://github.com/AIM-Harvard/pyradiomics/pull/355>`_)
- Update validation rule for image filters (remove hardcoded filters by package-detected filters).
  (`#364 <https://github.com/AIM-Harvard/pyradiomics/pull/364>`_)
- Add missing requirements for LBP filters in the dockerfile.
  (`#389 <https://github.com/AIM-Harvard/pyradiomics/pull/389>`_)
- Fix deprecation error in feature extractor. (`da1fc16 <https://github.com/AIM-Harvard/pyradiomics/commit/da1fc16>`_)
- Fix axis definition in wavelet. (`4027a52 <https://github.com/AIM-Harvard/pyradiomics/commit/4027a52>`_)
- Fix erroneous double return of wavelet approximation.
  (`c8ceee2 <https://github.com/AIM-Harvard/pyradiomics/commit/c8ceee2>`_)

Tests
#####

- Improve testing badge layout. (`#312 <https://github.com/AIM-Harvard/pyradiomics/pull/312>`_)
- Remove unused testing configuration files. (`#313 <https://github.com/AIM-Harvard/pyradiomics/pull/313>`_)
- Add testing for wavelet output. (`#387 <https://github.com/AIM-Harvard/pyradiomics/pull/387>`_)
- Integrate publication to PyPi into the Continuous Integration, revise the CI workflow to test
  python 2.7, 3.4, 3.5 and 3.6 for all 3 platforms (Windows, Mac and Linux).
  **N.B. This makes PyRadiomics installable via pip**
  (`#394 <https://github.com/AIM-Harvard/pyradiomics/pull/394>`_)

Documentation
#############

- Update documentation of ``base.py`` (`#306 <https://github.com/AIM-Harvard/pyradiomics/pull/306>`_)
- Update notebooks to reflect most recent version of PyRadiomics.
  (`ac66e6c <https://github.com/AIM-Harvard/pyradiomics/commit/ac66e6c>`_)
- Add documentation detailing rationale of enforcing a fixed bin width.
  (`#320 <https://github.com/AIM-Harvard/pyradiomics/pull/320>`_)
- Update reference to official publication. (`b395904 <https://github.com/AIM-Harvard/pyradiomics/commit/b395904>`_)
- Update installation instructions for docker. (`#329 <https://github.com/AIM-Harvard/pyradiomics/pull/329>`_)
- Add version of NumPy, SimpleITK and PyWavelet to the additional information in the output.
  (`#342 <https://github.com/AIM-Harvard/pyradiomics/pull/342>`_)
- Add documentation for the calculation of Laplacian of Gaussian.
  (`#345 <https://github.com/AIM-Harvard/pyradiomics/pull/345>`_)
- Add refrences for the newly implemented filters
  (`4464d1c <https://github.com/AIM-Harvard/pyradiomics/commit/4464d1c>`_)
- Fix an error in the firstorder-Uniformity documentation.
  (`da7321d <https://github.com/AIM-Harvard/pyradiomics/commit/da7321d>`_)

Examples
########

- Add example for batchprocessing using a multithreaded approach.
  (`#305 <https://github.com/AIM-Harvard/pyradiomics/pull/305>`_)

Internal API
############

- Update batch script for the commandline interface. Ensures all required input is available and relative filepaths are
  relative to the input file, not the current working directory.
  (`#307 <https://github.com/AIM-Harvard/pyradiomics/pull/307>`_)
- Remove support for 32-bits python, as memory errors can arise when extracting from many or large images in 32-bits
  python. (`#310 <https://github.com/AIM-Harvard/pyradiomics/pull/310>`_)
- Simplify Calculation of Wavelet Filter. Does not change output.
  (`#323 <https://github.com/AIM-Harvard/pyradiomics/pull/323>`_)
- Refactor commandline interface to work with only 1 entry point (``pyradiomics``). Also add parallel-processing option
  for batch-processing (argument ``-j``, which specifies number of CPU cores to use).
  (`#347 <https://github.com/AIM-Harvard/pyradiomics/pull/347>`_)
- Reconfigur testing to allow the removal of testcases from the repository itself (still available as binary data
  attached to release 1.0.0) and store the baseline in a different format (allowing for easier change-tracking)
  (`#353 <https://github.com/AIM-Harvard/pyradiomics/pull/353>`_)
- Add a check for number of bins generated (preventing construction of too large matrices in C)
  (`#391 <https://github.com/AIM-Harvard/pyradiomics/pull/391>`_,
  `#393 <https://github.com/AIM-Harvard/pyradiomics/pull/393>`_)

-----------------
PyRadiomics 1.3.0
-----------------

Feature Calculation Changes
###########################

- Remove feature *Sum Variance*, as this is mathematically equal to *Cluster Tendency*.
  (`#300 <https://github.com/AIM-Harvard/pyradiomics/pull/300>`_)
- Fix feature formula error in NGTDM (incorrect use of square in *Complexity* and *Contrast*).
  (`#351 <https://github.com/AIM-Harvard/pyradiomics/pull/351>`_)

New Features
############

- Add a row by row customization of the extraction label in the batch processing command line script, as well as both
  batchprocessing examples.
  (`#262 <https://github.com/AIM-Harvard/pyradiomics/pull/262>`_)
- Allow value 0 for a resampled pixel spacing (per dimension). Values of 0 are replaced by the spacing for that
  dimension as it is in the original (non-resampled) mask. This allows resampling over a subset of dimension (e.g. only
  in-plane resampling when out-of-plane spacing is set to 0).
  (`#299 <https://github.com/AIM-Harvard/pyradiomics/pull/299>`_)
- Add optional resegmentation of mask based on customizable threshold.
  (`#302 <https://github.com/AIM-Harvard/pyradiomics/pull/302>`_)
- Add Neighbouring Gray Tone Difference Matrix (NGTDM) (`#296 <https://github.com/AIM-Harvard/pyradiomics/pull/296>`_)
- Add Add Gray Level Dependence Matrix (GLDM) (`#295 <https://github.com/AIM-Harvard/pyradiomics/pull/295>`_)
- Add a docker file that exposes the PyRadiomics commandline tools.
  (`#297 <https://github.com/AIM-Harvard/pyradiomics/pull/297>`_,
  `#301 <https://github.com/AIM-Harvard/pyradiomics/pull/301>`_)
- Add voxel-based calculation, allowing for extraction of feature maps (values per voxel instead of per ROI).
  (`#337 <https://github.com/AIM-Harvard/pyradiomics/pull/337>`_)

Bug fixes
#########

- In GLCM, the matrix is made symmetrical by adding the transposed matrix. However, ``numpy.transpose`` returns a view
  and not a copy of the array, causing erroneous results when adding it to the original array. use
  ``numpy.ndarray.copy`` to prevent this bug. **N.B. This affects the feature values calculated by GLCM when symmetrical
  matrix is enabled (as is the default setting).**
  (`#261 <https://github.com/AIM-Harvard/pyradiomics/pull/261>`_)
- Use a python implementation to compute eigenvalues for ``shape.py`` instead of SimpleITK. The implementation in
  SimpleITK assumes segmented voxels to be consecutive on the x-axis lines. Furthermore, it also assumes that all voxels
  on a given line of x have the same values for y and z (which is not necessarily the case).
  (`#264 <https://github.com/AIM-Harvard/pyradiomics/pull/264>`_)
- Removal of outliers was not applied to returned object in ``normalizeImage``.
  (`#277 <https://github.com/AIM-Harvard/pyradiomics/pull/277>`_)
- Fix python 3 incompatibility when using ``urllib``
  (`#285 <https://github.com/AIM-Harvard/pyradiomics/pull/285>`_)
- Fix broken URL link in feature visualization notebooks.
- Update docker manually install python2 support (since recently not supported by default in
  jupyter/datascience-notebook).
  (`#287 <https://github.com/AIM-Harvard/pyradiomics/pull/287>`_)
- For GLRLM and GLSZM, force2D keyword is passed manually, but was incorrectly named and therefore ignored. Fix name to
  enable forced 2D extraction for GLRLM and GLSZM. (`26b9ef3 <https://github.com/AIM-Harvard/pyradiomics/commit/26b9ef3>`_)
- Fix bug in the calculation of eigen values due to machine precision errors.
  (`#355 <https://github.com/AIM-Harvard/pyradiomics/pull/355>`_)

Tests
#####

- Update the C Matrices test, so that the C and python calculated matrices will have the same dimensions when compared
  (In the previous implementation, the ``_calculateCoefficients`` function was applied to the C calculated matrix, but
  not in the python calculated matrix, for some texture matrices, this function can change the dimension of the matrix).
  This update ensures that ``_calculateCoefficients`` is applied to neither matrix.
  (`#265 <https://github.com/AIM-Harvard/pyradiomics/pull/265>`_)
- Add a test to check validity of parameter files included in ``examples/exampleSettings``.
  (`#294 <https://github.com/AIM-Harvard/pyradiomics/pull/294>`_)

Documentation
#############

`version 1.3.0 docs <http://pyradiomics.readthedocs.io/en/1.3.0>`_

- Update reference. (`#271 <https://github.com/AIM-Harvard/pyradiomics/pull/271>`_)
- Move section "Customizing the Extraction" to the top level, to make it more visible.
  (`#271 <https://github.com/AIM-Harvard/pyradiomics/pull/271>`_)
- Change License to 3-clause BSD (`#272 <https://github.com/AIM-Harvard/pyradiomics/pull/272>`_
- Document the extend of compliance between PyRadiomics and the IBSI feature definitions
  (`#289 <https://github.com/AIM-Harvard/pyradiomics/pull/289>`_)
- Fix typos in documentation.
- Expand documentation on customizing the extraction
  (`#291 <https://github.com/AIM-Harvard/pyradiomics/pull/291>`_)
- Include contributing guidelines in sphinx-generated documentation and add a section on sharing parameter files.
  (`#294 <https://github.com/AIM-Harvard/pyradiomics/pull/294>`_)
- Insert missing line to enable all features in documentation on using the feature classes directly.
  (`5ce9f48 <https://github.com/AIM-Harvard/pyradiomics/commit/5ce9f48>`_)
- Fix typo in NGTDM documentation. (`ea9a6ce <https://github.com/AIM-Harvard/pyradiomics/commit/ea9a6ce>`_)
- Fix some typos in documentation of firstorder - std and gldm - GLN
  (`#369 <https://github.com/AIM-Harvard/pyradiomics/pull/369>`_)
- Add additional comments to the code of the Wavelet filter (``_swt3``).
  (`#375 <https://github.com/AIM-Harvard/pyradiomics/pull/375>`_)
- Add references to the new filter functions. (`4464d1c <https://github.com/AIM-Harvard/pyradiomics/commit/4464d1c>`_)

Examples
########
- Add example settings for CT, MR (3 scenarios).
  (`#273 <https://github.com/AIM-Harvard/pyradiomics/pull/273>`_)

Internal API
############

- Remove unnecessary rows and columns from texture matrices prior to feature calculation. This does not affect the value
  of the calculated features, as the i and j vectors are updated accordingly, but it does reduce both computation time
  and memory requirements. This is especially the case when calculating GLSZM on large segmentations, where there may be
  many 'empty' zone sizes (i.e. no zones of that size are present in the ROI). This reduces the size of the matrix,
  which therefore reduces the memory needed and the number of calculations performed in the vectorized operations.
  (`#265 <https://github.com/AIM-Harvard/pyradiomics/pull/265>`_)
- Remove circular import statement in ``__init__.py`` (circular with ``radiomics.base``)
  (`#270 <https://github.com/AIM-Harvard/pyradiomics/pull/270>`_)
- Revise initialization of the feature class.
  (`#274 <https://github.com/AIM-Harvard/pyradiomics/pull/274>`_)
- Rename parts of the customization variables and functions to better reflect their definition
  (`#291 <https://github.com/AIM-Harvard/pyradiomics/pull/291>`_)
- Update C extensions: Make python wrapping more similar for different feature classes, simplify calculation of surface
  area, remove deprecated Numpy C-API references and implement angle-generation in C.
  (`#360 <https://github.com/AIM-Harvard/pyradiomics/pull/360>`_)
- Remove Python equivalents of C extensions: Some, but not all C extensions had python equivalents, which calculated
  equal values but, by using a python-only implementation, are much slower than the C extension. Only advantage is that
  it would also work when compiling the code fails. Also update the tests to check consistency of the calculated
  matrices against a baseline file (binary numpy array file) instead of python calculated matrices.
  (`#373 <https://github.com/AIM-Harvard/pyradiomics/pull/373>`_)

License
#######
- Switch to 3-clause BSD license.
  (`#272 <https://github.com/AIM-Harvard/pyradiomics/pull/272>`_)

-----------------
PyRadiomics 1.2.0
-----------------

Feature Calculation Changes
###########################

- Remove feature *SumVariance*, rename *SumVariance2*  to *SumVariance*. *SumVariance* reflected the formula as is
  defined in the paper by Haralick et al [1]_. However, the variance is calculated by subtracting the entropy as opposed to
  subtracting the average, most likely due to a typo('f8' instead of 'f6'). *SumVariance2* reflected the formula where
  the average is subtracted and is retained as the only *SumVariance*.
  (`#233 <https://github.com/AIM-Harvard/pyradiomics/pull/233>`_)
- Redefine features *Elongation* and *Flatness* as the inverse of the original definition. This prevents a returned
  value of NaN when the shape is completely flat. (`#234 <https://github.com/AIM-Harvard/pyradiomics/pull/234>`_)
- In certain edge cases, the calculated maximum diameters may be too small when calculating using the python
  implementation. This is corrected by the C extension and a warning is now logged when calculating these features in
  python. **N.B. As of this change, maximum diameter is not available for calculation in full-python mode**
  (`#257 <https://github.com/AIM-Harvard/pyradiomics/pull/257>`_)
- For certain formulas, a NaN value is returned in some edge cases. Catch this and return a predefined value instead.
  Document this behaviour in the docstrings of the features affected.
  (`#248 <https://github.com/AIM-Harvard/pyradiomics/pull/248>`_)

New Features
############

- Add Region of Interest checks. (`#223 <https://github.com/AIM-Harvard/pyradiomics/pull/223>`_,
  `#227 <https://github.com/AIM-Harvard/pyradiomics/pull/227>`_)
- Add variable column support for batch input file (`#228 <https://github.com/AIM-Harvard/pyradiomics/pull/228>`_)
- Add Docker support (`#236 <https://github.com/AIM-Harvard/pyradiomics/pull/236>`_)

Bug fixes
#########

- Instantiate output with input in ``commandlinebatch``
- Correct ``Np`` when weighting is applied in GLRLM (`#229 <https://github.com/AIM-Harvard/pyradiomics/pull/229>`_)
- Update CSV generators to reflect variable number of columns for input CSV in batch processing
  (`#246 <https://github.com/AIM-Harvard/pyradiomics/pull/246>`_)
- Return corrected mask when it had to be resampled due to geometry mismatch errors
  (`#260 <https://github.com/AIM-Harvard/pyradiomics/pull/260>`_)

Requirements
############

- Remove ``tqdm`` requirement (`#232 <https://github.com/AIM-Harvard/pyradiomics/pull/232>`_)
- Reorganize requirements, with requirements only needed during development moved to ``requirements-dev.txt``
  (`#231 <https://github.com/AIM-Harvard/pyradiomics/pull/231>`_)

Documentation
#############

`version 1.2.0 docs <http://pyradiomics.readthedocs.io/en/1.2.0>`_

- Update feature docstrings, making them more easily adaptable for article supplements
  (`#233 <https://github.com/AIM-Harvard/pyradiomics/pull/233>`_)
- Add FAQ concerning the cmatrices lib path (`#233 <https://github.com/AIM-Harvard/pyradiomics/pull/233>`_)
- Add developer install step to documentation (`#245 <https://github.com/AIM-Harvard/pyradiomics/pull/245>`_)
- Remove use of ``sudo`` (`#233 <https://github.com/AIM-Harvard/pyradiomics/pull/233>`_)
- Fix subclass name in feature class signature (section "Developers")
- Add subsection on customizing the extraction to the "Usage" section
  (`#252 <https://github.com/AIM-Harvard/pyradiomics/pull/252>`_)
- Remove SimpleITK installation workaround, this is no longer needed
  (`#249 <https://github.com/AIM-Harvard/pyradiomics/pull/249>`_)
- Add a changelog to keep track of changes and integrate this into the auto generated documentation
  (`#255 <https://github.com/AIM-Harvard/pyradiomics/pull/255>`_)

Examples
########

- Add ``pandas`` example, showing how to process PyRadiomics output/input using the ``pandas`` library
  (`#228 <https://github.com/AIM-Harvard/pyradiomics/pull/228>`_)

Internal API
############

- Add function to get or download test case (`#235 <https://github.com/AIM-Harvard/pyradiomics/pull/235>`_)
- Rewrite C Extension algorithm for GSLZM. Instead of searching over the image for the next voxel when
  growing a region, store all unprocessed voxels in a stack. This yields a significant increase in performance,
  especially in large ROIs. Requires slightly more memory (1 array, type integer, size equal to number of voxels in
  the ROI) (`#257 <https://github.com/AIM-Harvard/pyradiomics/pull/257>`_)
- Implement C extension for calculation of maximum diameters.
  (`#257 <https://github.com/AIM-Harvard/pyradiomics/pull/257>`_)

Cleanups
########

- Restructure repository (`#254 <https://github.com/AIM-Harvard/pyradiomics/pull/254>`_)

  - Move jupyter notebooks to separate root folder (``root/notebooks``)
  - Move example script to separate root folder (``root/examples``), with example settings in separate subfolder
    (``root/examples/exampleSettings``)
  - ``bin`` folder now only contains support scripts for the core code (i.e. generators for input files for batch
    processing and scripts to generate new baselines or to resample a mask to the image geometry)

-----------------
PyRadiomics 1.1.1
-----------------

Feature Calculation Changes
###########################

- Correct error in formula for *Compactness1*. **N.B. Baseline updated!**
  (`#218 <https://github.com/AIM-Harvard/pyradiomics/pull/218>`_)
- Remove feature *Roundness*, as this feature is identical to feature *Sphericity*, but uses different implementation
  for surface area calculation (all implemented in SimpleITK)
  (`#218 <https://github.com/AIM-Harvard/pyradiomics/pull/218>`_)
- Change handling of cases where ``max(X) mod binwidth = 0`` during image discretization. These used to be assigned to
  topmost bin, but this produces unexpected behaviour (i.e. in range 1, 2, 3, 4, 5 with binwidth 1, value 5 would be
  discretized to 4 in stead of 5). Value now assigned is topmost bin + 1 (in concordance with default behavior of
  ``numpy.digitize``) (`#219 <https://github.com/AIM-Harvard/pyradiomics/pull/219>`_)
- Change default value for ``voxelArrayShift`` (from 2000 to 0), this is to prevent unknowingly using a too large shift
  when not necessary. Document effect of this parameter in the first order formulas affected.
  (`#219 <https://github.com/AIM-Harvard/pyradiomics/pull/219>`_)

New features
############

- Add forced 2D extraction (as alternative to resampling for handling anisotropy in voxels spacing)
- Enable specification of distances between neighbors for GLCM matrix calculation

(`#215 <https://github.com/AIM-Harvard/pyradiomics/pull/215>`_)

New Parameters
##############

- ``force2D``, Boolean default ``False``. Set to ``True`` to force a by slice texture calculation. Dimension that
  identifies the 'slice' can be defined in ``force2Ddimension``. If input ROI is already a 2D ROI, features are
  automatically extracted in 2D.
- ``force2Ddimension``, int, range 0-2, default 0. Specifies the 'slice' dimension for a by-slice feature extraction.
  Value 0 identifies the 'z' dimension (axial plane feature extraction), and features will be extracted from the xy
  plane. Similarly, 1 identifies the y dimension (coronal plane) and 2 the x dimension (saggital plane).
- ``distances``, List of integers, default ``[1]``. This specifies the distances between the center voxel and the
  neighbor, for which angles should be generated.

(`#215 <https://github.com/AIM-Harvard/pyradiomics/pull/215>`_)

Bug fixes
#########

- Add some missing python 3 compatibility lines to the supporting script ``addClassToBaseline`` and command line script
  ``pyradiomicsbatch`` (`#210 <https://github.com/AIM-Harvard/pyradiomics/pull/210>`_,
  `#214 <https://github.com/AIM-Harvard/pyradiomics/pull/214>`_)
- Fix bug when loading image as file path and mask as SimpleITK object.
  (`#211 <https://github.com/AIM-Harvard/pyradiomics/pull/211>`_)
- Change location of parameter schema files. These files are otherwise not included in the wheel distribution.
  (`#221 <https://github.com/AIM-Harvard/pyradiomics/pull/221>`_)

Requirements
############

- Add sphinx_rtd_theme to requirements (needed to build documentation).
  (`#222 <https://github.com/AIM-Harvard/pyradiomics/pull/222>`_)

Documentation
#############

`version 1.1.1 docs <http://pyradiomics.readthedocs.io/en/1.1.1>`_

- Split package documentation into "Pipeline Modules" (all non-feature-class modules) and "Feature Definitions"
  (feature class modules)
- Add developers section with documentation on how to implement new filters, feature and feature classes.
- Add FAQ section with some trouble shooting tips
- Rename some GLSZM features, this is to make them more consistent with GLRLM features, which are similar, but
  calculated on a different matrix
- Add documentation for Elongation and Flatness
- Document mathematical correlation between various Shape features.

(`#216 <https://github.com/AIM-Harvard/pyradiomics/pull/216>`_)

Internal API
############

- Update logging with more extensive debug logging and more informative info log messages.
  (`#220 <https://github.com/AIM-Harvard/pyradiomics/pull/220>`_)
- Replace parameter verbose with output printing implemented in logging. Control verbosity level to output (stderr) by
  calling :py:func:`~radiomics.setVerbosity`, where level determines the verbosity level (as defined in python logging).
  This prints out the requested levels of the log messaging, where process reports with parameter verbose are now
  classified as INFO-level messages (i.e. specify INFO or DEBUG to enable these). **N.B. parameter verbose is not longer
  supported and will throw an error if passed in the parameter file**
  (`#220 <https://github.com/AIM-Harvard/pyradiomics/pull/220>`_)
- Add feature class and input image type checks in ``featureextractor`` when changing these settings.
  (`#213 <https://github.com/AIM-Harvard/pyradiomics/pull/213>`_)
- Remove usage of ``eval`` (replaced by implementations of ``getattr``), this is a more secure approach.
  (`#216 <https://github.com/AIM-Harvard/pyradiomics/pull/216>`_)
- Define default settings in featureextractor in a separate function. This is to ensure consistency in applied default
  settings, as well as make them easily available outside of featureextractor
  (`#216 <https://github.com/AIM-Harvard/pyradiomics/pull/216>`_)
- Update reference for citing PyRadiomics (`#224 <https://github.com/AIM-Harvard/pyradiomics/pull/224>`_)


Cleanups
########

- Remove unused variable (``self.provenance_on`` in ``featureextractor``, this value is now replaced by a customizable
  setting)

-----------------
PyRadiomics 1.1.0
-----------------

New features
############

- Image normalization. This feature enables the normalization of image intensity values prior to feeding them to the
  extraction pipeline (i.e. before any other preprocessing steps are performed). Normalization is based on the all gray
  values contained within the image, not just those defined by the ROI in the mask.
- C Extensions for texture matrix and surface area calculation. These extensions enhance performance of texture matrix
  calculation associated GLCM, GLRLM and GLSZM features and of surface area calculation. Below shows the decrease in
  computation time for the 5 test cases included in PyRadiomics.
  (`#158 <https://github.com/AIM-Harvard/pyradiomics/pull/158>`_,
  `#200 <https://github.com/AIM-Harvard/pyradiomics/pull/200>`_,
  `#202 <https://github.com/AIM-Harvard/pyradiomics/pull/202>`_)

  - GLCM 6913 ms -> 3 ms
  - GLRLM 1850 ms -> 10 ms
  - GLSZM 12064 ms -> 58 ms
  - Surface Area 3241 ms -> 1 ms

New Parameters
##############

- ``additionalInfo`` Boolean, default ``True``. Enables additional information in the output if set to ``True``.
  (`#190 <https://github.com/AIM-Harvard/pyradiomics/pull/190>`_)
- ``enableCExtensions`` Boolean, defailt ``True``. Enables enhanced performance for texture matrix calculation using C
  extensions if set to ``True``. (`#202 <https://github.com/AIM-Harvard/pyradiomics/pull/202>`_)
- ``normalize`` Boolean, default `` False``. If set to true, normalizes image before feeding it into the extraction
  pipeline. (`#209 <https://github.com/AIM-Harvard/pyradiomics/pull/209>`_)
- ``normalizeScale`` Float, > 0, default 1. Enables scaling of normalized intensities by specified value.
  (`#209 <https://github.com/AIM-Harvard/pyradiomics/pull/209>`_)
- ``removeOutliers`` Float, > 0, default ``None``. If set, outliers (defined by the value specified) are removed by
  setting them to the outlier value. Outlier value is defined on the non-scaled values.
  (`#209 <https://github.com/AIM-Harvard/pyradiomics/pull/209>`_)

Bug fixes
#########

- Unlink venv only when needed in Circle CI testing (`#199 <https://github.com/AIM-Harvard/pyradiomics/pull/199>`_)
- Fix datatype error when calling ``SimpleITK.ResampleImageFilter.SetSize()`` (only causes error in python 3,
  `#205 <https://github.com/AIM-Harvard/pyradiomics/pull/205>`_)

Requirements
############

- Add requirement for ``six>=1.10.0``, needed to make PyRadiomics compatible with both python 2 and 3.

Documentation
#############

`version 1.1.0 docs <http://pyradiomics.readthedocs.io/en/1.1.0>`_

- Documentation on installation and usage is upgraded, with the addition of an embedded instruction video (in section
  "Usage", cued at the section on usage examples). (`#187 <https://github.com/AIM-Harvard/pyradiomics/pull/187>`_)
- Updated contact information to point to the google groups.
- Updated the classifiers in the setup script to reflect the more advanced status of Pyradiomics.
  (`#193 <https://github.com/AIM-Harvard/pyradiomics/pull/193>`_)

Tests
#####

- Add support for multiple python versions and platforms, now including python 2.7, 3.4, 3.5 (32/64bits) for Linux,
  Windows and Mac. (`#183 <https://github.com/AIM-Harvard/pyradiomics/pull/183>`_,
  `#191 <https://github.com/AIM-Harvard/pyradiomics/pull/191>`_,
  `#199 <https://github.com/AIM-Harvard/pyradiomics/pull/199>`_)
- Testing output is upgraded to ensure unique feature names (`#195 <https://github.com/AIM-Harvard/pyradiomics/pull/195>`_,
  `#197 <https://github.com/AIM-Harvard/pyradiomics/pull/197>`_)
- Add ``test_cmatrices`` to assert conformity between output from Python and C based texture matrix calculation.

Internal API
############

- :py:func:`~radiomics.getFeatureClasses` and :py:func:`~radiomics.getInputImageTypes` are moved from
  `Feature Extractor <radiomics-featureextractor-label>` to the global radiomics namespace. This enumerates the possible
  feature classes and filters at initialization of the toolbox, and ensures feature classes are imported at
  initialization. (`#190 <https://github.com/AIM-Harvard/pyradiomics/pull/190>`_,
  `#198 <https://github.com/AIM-Harvard/pyradiomics/pull/198>`_)
- Python 3 Compatibility. Add support for compatibility with python 2.7 and python >= 3.4. This is achieved using
  package ``six``.
- Standardize function names for calculating matrices in python and with C extensions to ``_calculateMatrix`` and
  ``_calculateCMatrix``, respectively.
- Make C code consistent with C89 convention. All variables (pointers for python objects) are initialized at top of each
  block.
- Optimize GLSZM calculation (C extension)

  - Define temporary array for holding the calculated zones. During calculation, the matrix must be able to store all
    possible zones, ranging from zone size 1 to total number of voxels (Ns), for each gray level (Ng). In this case, the
    GLSZM would be initialized with size Ng * Ns, which is very memory intensive. Instead, use a temporary array of size
    (Ns * 2) + 1, which stores all calculated zones in pairs of 2 elements: the first element holds the gray level, the
    second the size of the calculated zone. The first element after the last zone is set to -1 to serve as a stop sign
    for the second function, which translates the temporary array into the final GLSZM, which can be directly
    initialized at optimum size.
  - Use ``calloc`` and ``free`` for the temporary array holding the calculated zones.
  - Use ``char`` datatype for mask. (signed char in GLSZM).
  - Uses ``while`` loops. This allows to reduce the memory usage. Additionally, we observed that with recursive
    functions it was 'unexpectedly' failing.
  - Optimized search that finds a new index to process in the region growing.

-----------------
PyRadiomics 1.0.1
-----------------

New features
############

- Added 2 commandline scripts ( pyradiomics and pyradiomicsbatch), which enable feature extraction directly from the
  commandline. For help on usage, run script with “-h” argument.
  (`#188 <https://github.com/AIM-Harvard/pyradiomics/pull/188>`_,
  `#194 <https://github.com/AIM-Harvard/pyradiomics/pull/194>`_,
  `#196 <https://github.com/AIM-Harvard/pyradiomics/pull/196>`_,
  `#205 <https://github.com/AIM-Harvard/pyradiomics/pull/205>`_)

Bug fixes
#########

- Fix hardcoded label in shape (`#175 <https://github.com/AIM-Harvard/pyradiomics/pull/175>`_)
- Fix incorrect axis when deleting empty angles in GLCM (`#176 <https://github.com/AIM-Harvard/pyradiomics/pull/176>`_)
- Numpy slicing error in application of wavelet filters. This error caused the derived image to be erroneously rotated
  and flipped, with misaligned mask as a result.(`#182 <https://github.com/AIM-Harvard/pyradiomics/pull/182>`_)

Requirements
############

- Revert numpy minimum requirement to ``1.9.2``. All operations in PyRadiomics are supported by this version, and it is
  the version used by Slicer. By reverting the minimum required version, installing PyRadiomics in the slicer extension
  does not cause an update of the numpy package distributed by slicer.
  (`#180 <https://github.com/AIM-Harvard/pyradiomics/pull/180>`_)

Documentation
#############

`version 1.0.1 docs <http://pyradiomics.readthedocs.io/en/v1.0.1>`_

- Update on the documentation, reflecting recent changes in the code.
- Add developers and affiliations to ReadMe and documentation
  (`#177 <https://github.com/AIM-Harvard/pyradiomics/pull/177>`_)
- Added additional references and updated installation and usage section.

Internal API
############

- Different implementation of the various filters. No changes to calculation, but has a changed signature.

  **N.B. This results in inputImages to be differently defined (different capitalization, e.g. "orginal" should now be
  "Original"). See documentation for definition of inputImages (featureextractor section).**

---------------
PyRadiomics 1.0
---------------

New features
############

- Initial Release of PyRadiomics

Work in progress
################

- Full python calculation (C matrices branch not stable and reserved for later release)

Documentation
#############

- Documentation published at `readthedocs <http://pyradiomics.readthedocs.io/en/v1.0>`_

.. [1] Haralick R, Shanmugan K, Dinstein I: Textural features for image classification. IEEE Trans Syst Man Cybern
       1973:610–621.
