.. _developers:

==========
Developers
==========

This section contains information on how to add or customize the feature classes and filters available in PyRadiomics.
PyRadiomics enumerates the available feature classes and input image types at initialization of the toolbox. These are
available from the global ``radiomics`` namespace by use of the functions :py:func:`~radiomics.getFeatureClasses()` and
:py:func:`~radiomics.getInputImageTypes()`, respectively. Individual features in a feature class are enumerated at
initialization of the class. See also the
`contributing guidelines <https://github.com/Radiomics/pyradiomics/blob/master/CONTRIBUTING.md>`_

----------------------------
Signature of a feature class
----------------------------

Each feature class is defined in a separate module, the module name is used as the feature class name (e.g. if module
tex.py matches the feature class signature, it is available in the PyRadiomics toolbox as the 'tex' feature class). In
the module a class should be defined that fits the following signature::

    [required imports]
    from radiomics import base

    class Radiomics[Name](base.RadiomicsFeaturesBase):
        """
        Feature class docstring
        """

        def __init__(self, inputImage, inputMask, **kwargs):
            super(RadiomicsGLCM, self).__init__(inputImage, inputMask, **kwargs)

        def get[Feature]FeatureValue(self):
            """
            Feature docstring
            """

            return [value]

* At the top should be the import statements for packages required by the feature class. Unused import statements should
  be removed (flake8 will fail if unused import statements are encountered, or import statements are not structured as
  defined by appnexus).
* The class name should be 'Radiomics' followed by the name of the class (usually similar to the
  module name. However, this name is not used elsewhere and may be named arbitrarily).
* The class should inherit (directly or indirectly) from ``base.RadiomicsFeaturesBase``, which is an abstract class
  defining the common interface for the feature classes
* Additional initialization steps should be called in the ``__init__`` function. For default variables initialized, see
  :ref:`radiomics-base-label`.
* Documentation is required! Both at the class level (Feature class docstring) and at the level of the individual
  features (Feature docstring).
* If the feature class uses C extensions for matrix calculation enhancement, which should be tested using
  ``test_cmatrices``, matrix calculation should be implemented as follows:

  * The function calculating the matrix in python should be defined in a function called ``_calculateMatrix``.
  * The function calculating the matrix using the C extension should be defined in a function called ``_calculateCMatrix``.
  * The functions to calculate the matrix accept no additional input arguments other than the ``self`` argument, and
    return the fully processed matrix as a numpy array.
  * The fully processed matrix should be assigned to a variable in the feature class named ``P_[Name]``, where
    ``[Name]`` is identical to the feature class name (module name) (e.g. in feature class ``glcm``, matrix is stored in
    variable ``P_glcm``

--------------------------------
Signature of individual features
--------------------------------

Each individual feature is defined as a function in the feature class with the ``get[Name]FeatureValue(self)``
signature, where ``[Name]`` is the feature name (unique on the feature class level). It accepts no input arguments, and
should return a scalar value. The ``self`` argument represents the instantiated feature class that defines the function,
and identifies the feature function as non-static.

---------------------
Signature of a filter
---------------------

All filters are defined in the :ref:`imageoperations module<radiomics-imageoperations-label>`, and identified by the
signature ``get[Name]Image(inputImage, **kwargs)``. Here, ``[Name]`` represents the unique name for the image type,
which is also used to identify the filter during extraction. The input of a filter is fixed and consists of the
``inputImage``, a SimpleITK Image object of the original image and ``**kwargs``, which are the customized setting that
should be used for the extraction of features from the derived image.

One or more derived images are returned using the 'yield' statement: ``yield derivedImage, inputImageName, kwargs``.
Here, ``derivedImage`` is one SimpleITK image object representing the filtered imag
e, ``inputImageName`` is a unique
string identifying features calculated using this filter in the output and ``kwargs`` are the customized settings for
the extraction (``**kwargs`` passed as input, without the double asterisk). Multiple derived images can be
returned by multiple yield statements, or yield statements inside a loop. Please note that only one derived image should
be returned on each call to yield and that ``inputImageName`` is a unique name for each returned derived image. Derived
images must have the same dimensions and occupy the same physical space to ensure compatibility with the mask.

------------------------------
Addtional points for attention
------------------------------

Code style
**********

To keep the PyRadiomics code consistent and as readable as possible, some style rules are enforced. These are part of
the continuous testing and implemented using flake8. See also the ``.flake8`` configuration file in the root of the
repository. To aid in keeping a consistent code style, a ``.editorconfig`` file is provided in the root of the folder.

Module names should be lowercase, without underscores or spaces. Class names, function names and variables should be
declared using camelcase, with uppercase first letter for class names and lowercase first letter otherwise. Private
helper functions (which should not be included in the documentation) should be declared using a '_' prefix. This is
consistent with the python style for marking them as 'private', and will automatically exclude them from the generated
documentation.

Documentation
*************

The documentation of PyRadiomics is auto-generated from static files contained in the ``docs`` folder and the docstrings
of the Python code files. When a new feature class is added, this has to be added to the static file (``features.rst``)
describing the feature classes as well. If done so, sphinx will take care of the rest. A featureclass can be added as
follows::

    <Class Name> Features
    ---------------------

    .. automodule:: radiomics.<module name>
        :members:
        :undoc-members:
        :show-inheritance:
        :member-order: bysource

Documentation providing information of the feature class as a whole (e.g. how the feature matrix is calculated) should
be provided in the docstring of the class. Definition of individual features, including the mathematical formulas should
be provided in the docstrings of the feature functions. A docstring of the module is not required.

The presence of a docstring at the class level and at the level of each individual feature is required and checked
during testing. Missing docstrings will cause the test to fail.

Testing
*******

To ensure consistency in the extraction provided by PyRadiomics, continuous testing is used to test the PyRadiomics
source code after each commit. These tests are defined in the test folder and used to run tests for the following
environments:

    - Python 2.7 32 and 64 bits (Windows, Linux and Mac)
    - Python 3.4 32 and 64 bits (Windows and Linux)
    - Python 3.5 32 and 64 bits (Windows and Linux)

Python 3 testing for mac is currently disabled for Mac due to some issues with the SimpleITK package for python 3.

There are 3 testing scripts run for PyRadiomics. The first test is ``test_cmatrices``, which asserts if the matrices
calculated by the C extensions match those calculated by Python. A threshold of 1e-3 is used to allow for machine
precision errors. The second test is ``test_docstrings``, which asserts if there is missing documentation as described
above. The final and most important test is ``test_features``, which compares the features calculated by PyRadiomics
against a known baseline using 5 test cases. These test cases and the baseline are stored in the ``data`` folder of the
repository. This ensures that changes to the code do not silently change the calculated values of the features.

To add a new feature class to the baseline, run the ``addClassToBaseline.py`` script, contained in the ``bin`` folder.
This script detects if there are feature classes in PyRadiomics, for which there is no baseline available. If any are
found, a new baseline if calculated for these classes in the full-python mode and added to the baseline files. These new
baseline files then needed to be included in the repository and committed.
