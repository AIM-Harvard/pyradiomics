=============
Release Notes
=============

------------
Next Release
------------

New Features
############

- Add a row by row customization of the extraction label in the batch processing command line script, as well as both
  batchprocessing examples.
  (`#262 <https://github.com/Radiomics/pyradiomics/pull/262>`_)

Bug fixes
#########

- In GLCM, the matrix is made symmetrical by adding the transposed matrix. However, ``numpy.transpose`` returns a view
  and not a copy of the array, causing erroneous results when adding it to the original array. use
  ``numpy.ndarray.copy`` to prevent this bug. **N.B. This affects the feature values calculated by GLCM when symmetrical
  matrix is enabled (as is the default setting).**
  (`#261 <https://github.com/Radiomics/pyradiomics/pull/261>`_)
- Use a python implementation to compute eigenvalues for ``shape.py`` instead of SimpleITK. The implementation in
  SimpleITK assumes segmented voxels to be consecutive on the x-axis lines. Furthermore, it also assumes that all voxels
  on a given line of x have the same values for y and z (which is not necessarily the case).
  (`#264 <https://github.com/Radiomics/pyradiomics/pull/264>`_)
- Removal of outliers was not applied to returned object in ``normalizeImage``.
  (`#277 <https://github.com/Radiomics/pyradiomics/pull/277>`_)
- Fix python 3 incompatibility when using ``urllib``
  (`#285 <https://github.com/Radiomics/pyradiomics/pull/285>`_)
- Fix broken URL link in feature visualization notebooks.
- Update docker manually install python2 support (since recently not supported by default in
  jupyter/datascience-notebook).
  (`#287 <https://github.com/Radiomics/pyradiomics/pull/287>`_)
- For GLRLM and GLSZM, force2D keyword is passed manually, but was incorrectly named and therefore ignored. Fix name to
  enable forced 2D extraction for GLRLM and GLSZM.

Tests
#####

- Update the C Matrices test, so that the C and python calculated matrices will have the same dimensions when compared
  (In the previous implementation, the ``_calculateCoefficients`` function was applied to the C calculated matrix, but
  not in the python calculated matrix, for some texture matrices, this function can change the dimension of the matrix).
  This update ensures that ``_calculateCoefficients`` is applied to neither matrix.
  (`#265 <https://github.com/Radiomics/pyradiomics/pull/265>`_)
- Add a test to check validity of parameter files included in ``examples/exampleSettings``.
  (`#294 <https://github.com/Radiomics/pyradiomics/pull/294>`_)

Documentation
#############

- Update reference. (`#271 <https://github.com/Radiomics/pyradiomics/pull/271>`_)
- Move section "Customizing the Extraction" to the top level, to make it more visible.
  (`#271 <https://github.com/Radiomics/pyradiomics/pull/271>`_)
- Change License to 3-clause BSD (`#272 <https://github.com/Radiomics/pyradiomics/pull/272>`_
- Document the extend of compliance between PyRadiomics and the IBSI feature definitions
  (`#289 <https://github.com/Radiomics/pyradiomics/pull/289>`_)
- Fix typos in documentation.
- Expand documentation on customizing the extraction
  (`#291 <https://github.com/Radiomics/pyradiomics/pull/291>`_)
- Include contributing guidelines in sphinx-generated documentation and add a section on sharing parameter files.
  (`#294 <https://github.com/Radiomics/pyradiomics/pull/294>`_)
- Insert missing line to enable all features in documentation on using the feature classes directly.

Examples
########
- Add example settings for CT, MR (3 scenarios).
  (`#273 <https://github.com/Radiomics/pyradiomics/pull/273>`_)

Internal API
############

- Remove unnecessary rows and columns from texture matrices prior to feature calculation. This does not affect the value
  of the calculated features, as the i and j vectors are updated accordingly, but it does reduce both computation time
  and memory requirements. This is especially the case when calculating GLSZM on large segmentations, where there may be
  many 'empty' zone sizes (i.e. no zones of that size are present in the ROI). This reduces the size of the matrix,
  which therefore reduces the memory needed and the number of calculations performed in the vectorized operations.
  (`#265 <https://github.com/Radiomics/pyradiomics/pull/265>`_)
- Remove circular import statement in ``__init__.py` (circular with ``radiomics.base``)
  (`#270 <https://github.com/Radiomics/pyradiomics/pull/270>`_)
- Revise initialization of the feature class.
  (`#274 <https://github.com/Radiomics/pyradiomics/pull/274>`_)
- Rename parts of the customization variables and functions to better reflect their definition
  (`#291 <https://github.com/Radiomics/pyradiomics/pull/291>`_)

License
#######
- Switch to 3-clause BSD license.
  (`#272 <https://github.com/Radiomics/pyradiomics/pull/272>`_)

-----------------
PyRadiomics 1.2.0
-----------------

Feature Calculation Changes
###########################

- Remove feature *SumVariance*, rename *SumVariance2*  to *SumVariance*. *SumVariance* reflected the formula as is
  defined in the paper by Haralick et al [1]_. However, the variance is calculated by subtracting the entropy as opposed to
  subtracting the average, most likely due to a typo('f8' instead of 'f6'). *SumVariance2* reflected the formula where
  the average is subtracted and is retrained as the only *SumVariance*.
  (`#233 <https://github.com/Radiomics/pyradiomics/pull/233>`_)
- Redefine features *Elongation* and *Flatness* as the inverse of the original definition. This prevents a returned
  value of NaN when the shape is completely flat. (`#234 <https://github.com/Radiomics/pyradiomics/pull/234>`_)
- In certain edge cases, the calculated maximum diameters may be too small when calculating using the python
  implementation. This is corrected by the C extension and a warning is now logged when calculating these features in
  python. **N.B. As of this change, maximum diameter is not available for calculation in full-python mode**
  (`#257 <https://github.com/Radiomics/pyradiomics/pull/257>`_)
- For certain formulas, a NaN value is returned in some edge cases. Catch this and return a predefined value instead.
  Document this behaviour in the docstrings of the features affected.
  (`#248 <https://github.com/Radiomics/pyradiomics/pull/248>`_)

New Features
############

- Add Region of Interest checks. (`#223 <https://github.com/Radiomics/pyradiomics/pull/223>`_,
  `#227 <https://github.com/Radiomics/pyradiomics/pull/227>`_)
- Add variable column support for batch input file (`#228 <https://github.com/Radiomics/pyradiomics/pull/228>`_)
- Add Docker support (`#236 <https://github.com/Radiomics/pyradiomics/pull/236>`_)

Bug fixes
#########

- Instantiate output with input in ``commandlinebatch``
- Correct ``Np`` when weighting is applied in GLRLM (`#229 <https://github.com/Radiomics/pyradiomics/pull/229>`_)
- Update CSV generators to reflect variable number of columns for input CSV in batch processing
  (`#246 <https://github.com/Radiomics/pyradiomics/pull/246>`_)
- Return corrected mask when it had to be resampled due to geometry mismatch errors
  (`#260 <https://github.com/Radiomics/pyradiomics/pull/260>`_)

Requirements
############

- Remove ``tqdm`` requirement (`#232 <https://github.com/Radiomics/pyradiomics/pull/232>`_)
- Reorganize requirements, with requirements only needed during development moved to ``requirements-dev.txt``
  (`#231 <https://github.com/Radiomics/pyradiomics/pull/231>`_)

Documentation
#############

`version 1.2.0 docs <http://pyradiomics.readthedocs.io/en/1.2.0>`_

- Update feature docstrings, making them more easily adaptable for article supplements
  (`#233 <https://github.com/Radiomics/pyradiomics/pull/233>`_)
- Add FAQ concerning the cmatrices lib path (`#233 <https://github.com/Radiomics/pyradiomics/pull/233>`_)
- Add developer install step to documentation (`#245 <https://github.com/Radiomics/pyradiomics/pull/245>`_)
- Remove use of ``sudo`` (`#233 <https://github.com/Radiomics/pyradiomics/pull/233>`_)
- Fix subclass name in feature class signature (section "Developers")
- Add subsection on customizing the extraction to the "Usage" section
  (`#252 <https://github.com/Radiomics/pyradiomics/pull/252>`_)
- Remove SimpleITK installation workaround, this is no longer needed
  (`#249 <https://github.com/Radiomics/pyradiomics/pull/249>`_)
- Add a changelog to keep track of changes and integrate this into the auto generated documentation
  (`#255 <https://github.com/Radiomics/pyradiomics/pull/255>`_)

Examples
########

- Add ``pandas`` example, showing how to process PyRadiomics output/input using the ``pandas`` library
  (`#228 <https://github.com/Radiomics/pyradiomics/pull/228>`_)

Internal API
############

- Add function to get or download test case (`#235 <https://github.com/Radiomics/pyradiomics/pull/235>`_)
- Rewrite C Extension algorithm for GSLZM. Instead of searching over the image for the next voxel when
  growing a region, store all unprocessed voxels in a stack. This yields a significant increase in performance,
  especially in large ROIs. Requires slightly more memory (1 array, type integer, size equal to number of voxels in
  the ROI) (`#257 <https://github.com/Radiomics/pyradiomics/pull/257>`_)
- Implement C extension for calculation of maximum diameters.
  (`#257 <https://github.com/Radiomics/pyradiomics/pull/257>`_)

Cleanups
########

- Restructure repository (`#254 <https://github.com/Radiomics/pyradiomics/pull/254>`_)

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
  (`#218 <https://github.com/Radiomics/pyradiomics/pull/218>`_)
- Remove feature *Roundness*, as this feature is identical to feature *Sphericity*, but uses different implementation
  for surface area calculation (all implemented in SimpleITK)
  (`#218 <https://github.com/Radiomics/pyradiomics/pull/218>`_)
- Change handling of cases where ``max(X) mod binwidth = 0`` during image discretization. These used to be assigned to
  topmost bin, but this produces unexpected behaviour (i.e. in range 1, 2, 3, 4, 5 with binwidth 1, value 5 would be
  discretized to 4 in stead of 5). Value now assigned is topmost bin + 1 (in concordance with default behavior of
  ``numpy.digitize``) (`#219 <https://github.com/Radiomics/pyradiomics/pull/219>`_)
- Change default value for ``voxelArrayShift`` (from 2000 to 0), this is to prevent unknowingly using a too large shift
  when not necessary. Document effect of this parameter in the first order formulas affected.
  (`#219 <https://github.com/Radiomics/pyradiomics/pull/219>`_)

New features
############

- Add forced 2D extraction (as alternative to resampling for handling anisotropy in voxels spacing)
- Enable specification of distances between neighbors for GLCM matrix calculation

(`#215 <https://github.com/Radiomics/pyradiomics/pull/215>`_)

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

(`#215 <https://github.com/Radiomics/pyradiomics/pull/215>`_)

Bug fixes
#########

- Add some missing python 3 compatibility lines to the supporting script ``addClassToBaseline`` and command line script
  ``pyradiomicsbatch`` (`#210 <https://github.com/Radiomics/pyradiomics/pull/210>`_,
  `#214 <https://github.com/Radiomics/pyradiomics/pull/214>`_)
- Fix bug when loading image as file path and mask as SimpleITK object.
  (`#211 <https://github.com/Radiomics/pyradiomics/pull/211>`_)
- Change location of parameter schema files. These files are otherwise not included in the wheel distribution.
  (`#221 <https://github.com/Radiomics/pyradiomics/pull/221>`_)

Requirements
############

- Add sphinx_rtd_theme to requirements (needed to build documentation).
  (`#222 <https://github.com/Radiomics/pyradiomics/pull/222>`_)

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

(`#216 <https://github.com/Radiomics/pyradiomics/pull/216>`_)

Internal API
############

- Update logging with more extensive debug logging and more informative info log messages.
  (`#220 <https://github.com/Radiomics/pyradiomics/pull/220>`_)
- Replace parameter verbose with output printing implemented in logging. Control verbosity level to output (stderr) by
  calling :py:func:`~radiomics.setVerbosity`, where level determines the verbosity level (as defined in python logging).
  This prints out the requested levels of the log messaging, where process reports with parameter verbose are now
  classified as INFO-level messages (i.e. specify INFO or DEBUG to enable these). **N.B. parameter verbose is not longer
  supported and will throw an error if passed in the parameter file**
  (`#220 <https://github.com/Radiomics/pyradiomics/pull/220>`_)
- Add feature class and input image type checks in ``featureextractor`` when changing these settings.
  (`#213 <https://github.com/Radiomics/pyradiomics/pull/213>`_)
- Remove usage of ``eval`` (replaced by implementations of ``getattr``), this is a more secure approach.
  (`#216 <https://github.com/Radiomics/pyradiomics/pull/216>`_)
- Define default settings in featureextractor in a separate function. This is to ensure consistency in applied default
  settings, as well as make them easily available outside of featureextractor
  (`#216 <https://github.com/Radiomics/pyradiomics/pull/216>`_)
- Update reference for citing PyRadiomics (`#224 <https://github.com/Radiomics/pyradiomics/pull/224>`_)


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
  (`#158 <https://github.com/Radiomics/pyradiomics/pull/158>`_,
  `#200 <https://github.com/Radiomics/pyradiomics/pull/200>`_,
  `#202 <https://github.com/Radiomics/pyradiomics/pull/202>`_)

  - GLCM 6913 ms -> 3 ms
  - GLRLM 1850 ms -> 10 ms
  - GLSZM 12064 ms -> 58 ms
  - Surface Area 3241 ms -> 1 ms

New Parameters
##############

- ``additionalInfo`` Boolean, default ``True``. Enables additional information in the output if set to ``True``.
  (`#190 <https://github.com/Radiomics/pyradiomics/pull/190>`_)
- ``enableCExtensions`` Boolean, defailt ``True``. Enables enhanced performance for texture matrix calculation using C
  extensions if set to ``True``. (`#202 <https://github.com/Radiomics/pyradiomics/pull/202>`_)
- ``normalize`` Boolean, default `` False``. If set to true, normalizes image before feeding it into the extraction
  pipeline. (`#209 <https://github.com/Radiomics/pyradiomics/pull/209>`_)
- ``normalizeScale`` Float, > 0, default 1. Enables scaling of normalized intensities by specified value.
  (`#209 <https://github.com/Radiomics/pyradiomics/pull/209>`_)
- ``removeOutliers`` Float, > 0, default ``None``. If set, outliers (defined by the value specified) are removed by
  setting them to the outlier value. Outlier value is defined on the non-scaled values.
  (`#209 <https://github.com/Radiomics/pyradiomics/pull/209>`_)

Bug fixes
#########

- Unlink venv only when needed in Circle CI testing (`#199 <https://github.com/Radiomics/pyradiomics/pull/199>`_)
- Fix datatype error when calling ``SimpleITK.ResampleImageFilter.SetSize()`` (only causes error in python 3,
  `#205 <https://github.com/Radiomics/pyradiomics/pull/205>`_)

Requirements
############

- Add requirement for ``six>=1.10.0``, needed to make PyRadiomics compatible with both python 2 and 3.

Documentation
#############

`version 1.1.0 docs <http://pyradiomics.readthedocs.io/en/1.1.0>`_

- Documentation on installation and usage is upgraded, with the addition of an embedded instruction video (in section
  "Usage", cued at the section on usage examples). (`#187 <https://github.com/Radiomics/pyradiomics/pull/187>`_)
- Updated contact information to point to the google groups.
- Updated the classifiers in the setup script to reflect the more advanced status of Pyradiomics.
  (`#193 <https://github.com/Radiomics/pyradiomics/pull/193>`_)

Tests
#####

- Add support for multiple python versions and platforms, now including python 2.7, 3.4, 3.5 (32/64bits) for Linux,
  Windows and Mac. (`#183 <https://github.com/Radiomics/pyradiomics/pull/183>`_,
  `#191 <https://github.com/Radiomics/pyradiomics/pull/191>`_,
  `#199 <https://github.com/Radiomics/pyradiomics/pull/199>`_)
- Testing output is upgraded to ensure unique feature names (`#195 <https://github.com/Radiomics/pyradiomics/pull/195>`_,
  `#197 <https://github.com/Radiomics/pyradiomics/pull/197>`_)
- Add ``test_cmatrices`` to assert conformity between output from Python and C based texture matrix calculation.

Internal API
############

- :py:func:`~radiomics.getFeatureClasses` and :py:func:`~radiomics.getInputImageTypes` are moved from
  `Feature Extractor <radiomics-featureextractor-label>` to the global radiomics namespace. This enumerates the possible
  feature classes and filters at initialization of the toolbox, and ensures feature classes are imported at
  initialization. (`#190 <https://github.com/Radiomics/pyradiomics/pull/190>`_,
  `#198 <https://github.com/Radiomics/pyradiomics/pull/198>`_)
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
  (`#188 <https://github.com/Radiomics/pyradiomics/pull/188>`_,
  `#194 <https://github.com/Radiomics/pyradiomics/pull/194>`_,
  `#196 <https://github.com/Radiomics/pyradiomics/pull/196>`_,
  `#205 <https://github.com/Radiomics/pyradiomics/pull/205>`_)

Bug fixes
#########

- Fix hardcoded label in shape (`#175 <https://github.com/Radiomics/pyradiomics/pull/175>`_)
- Fix incorrect axis when deleting empty angles in GLCM (`#176 <https://github.com/Radiomics/pyradiomics/pull/176>`_)
- Numpy slicing error in application of wavelet filters. This error caused the derived image to be erroneously rotated
  and flipped, with misaligned mask as a result.(`#182 <https://github.com/Radiomics/pyradiomics/pull/182>`_)

Requirements
############

- Revert numpy minimum requirement to ``1.9.2``. All operations in PyRadiomics are supported by this version, and it is
  the version used by Slicer. By reverting the minimum required version, installing PyRadiomics in the slicer extension
  does not cause an update of the numpy package distributed by slicer.
  (`#180 <https://github.com/Radiomics/pyradiomics/pull/180>`_)

Documentation
#############

`version 1.0.1 docs <http://pyradiomics.readthedocs.io/en/v1.0.1>`_

- Update on the documentation, reflecting recent changes in the code.
- Add developers and affiliations to ReadMe and documentation
  (`#177 <https://github.com/Radiomics/pyradiomics/pull/177>`_)
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
