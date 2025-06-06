.. _radiomics-developers:

==========
Developers
==========

This section contains information on how to add or customize the feature classes and filters available in PyRadiomics.
PyRadiomics enumerates the available feature classes and input image types at initialization of the toolbox. These are
available from the global ``radiomics`` namespace by use of the functions :py:func:`~radiomics.getFeatureClasses()` and
:py:func:`~radiomics.getImageTypes()`, respectively. Individual features in a feature class are enumerated at
initialization of the class. See also the :ref:`contributing guidelines <radiomics-contributing-label>`.

.. _radiomics-developers-featureclass:

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
            super(Radiomics[Name], self).__init__(inputImage, inputMask, **kwargs)
            # Feature class specific init

        def get[Feature]FeatureValue(self):
            """
            Feature docstring
            """
            # value = feature calculation using member variables of RadiomicsFeatureBase and this class.
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
  ``test_matrices``, matrix calculation should be implemented as follows:

  * The function calculating the matrix using the C extension should be defined in a function called ``_calculateMatrix``.
  * The functions to calculate the matrix accept no additional input arguments other than the ``self`` argument, and
    return the fully processed matrix as a numpy array.
  * The fully processed matrix should be assigned to a variable in the feature class named ``P_[Name]``, where
    ``[Name]`` is identical to the feature class name (module name) (e.g. in feature class ``glcm``, matrix is stored in
    variable ``P_glcm``

* A feature class specific logger is created by the base class, which will be a child logger (i.e. the 'radiomics.tex'
  logger in case of the feature class 'tex'). It is exposed in the feature class as ``self.logger``. Any log messages
  generated by the feature class should make use of this logger to ensure that the hierarchy of classes is correctly
  reflected in generated logs (i.e. ``self.logger.debug('message')`` to generate a debug log message).

-------------------
Adding the baseline
-------------------

During testing, calculated features are compared to a fixed baseline. If you implement a new class, it must be added to
the baseline, otherwise testing will fail. Fortunately, adding a new class to the baseline is fairly easy: just run the
``add_baseline.py`` script located in the ``tests`` folder. In case you change a feature and need to rebuild the
baseline for that class, run the script with the class name as an argument (e.g. ``python add_baseline.py glcm``).

.. _radiomics-developers-feature:

--------------------------------
Signature of individual features
--------------------------------

Each individual feature is defined as a function in the feature class with the ``get[Name]FeatureValue(self)``
signature, where ``[Name]`` is the feature name (unique on the feature class level). It accepts no input arguments, and
should return a scalar value. The ``self`` argument represents the instantiated feature class that defines the function,
and identifies the feature function as non-static.

.. _radiomics-developers-filter:

--------------------------
Signature of an image type
--------------------------

All image types are defined in the :ref:`imageoperations module<radiomics-imageoperations-label>`, and identified by the
signature ``get[Name]Image(inputImage, inputMask, **kwargs)``. Here, ``[Name]`` represents the unique name for the image
type, which is also used to identify the image type during extraction. The input of a image type function is fixed and
consists of the ``inputImage`` and ``inputMask``, SimpleITK Image objects of the original image and mask, respectively
and ``**kwargs``, which are the customized settings that should be used for the extraction of features from the derived
image.

One or more derived images are returned using the 'yield' statement: ``yield derivedImage, imageTypeName, kwargs``.
Here, ``derivedImage`` is one SimpleITK image object representing the filtered image, ``imageTypeName`` is a unique
string identifying features calculated using this filter in the output and ``kwargs`` are the customized settings for
the extraction (``**kwargs`` passed as input, without the double asterisk). Multiple derived images can be
returned by multiple yield statements, or yield statements inside a loop. Please note that only one derived image should
be returned on each call to yield and that ``imageTypeName`` is a unique name for *each* returned derived image. Derived
images must have the same dimensions and occupy the same physical space to ensure compatibility with the mask.

.. _radiomics-developers-progressreporting:

------------------
Progress Reporting
------------------

When operating in full-python mode, the calculation of the texture matrices can take some time. Therefore, PyRadiomics
provides the possibility to report the progress for calculation of GLCM and GLSZM.
This is only enabled in full-python mode when the verbosity (:py:func:`~radiomics.setVerbosity()`) is set to INFO or
DEBUG. By default, none is provided and no progress of matrix calculation will be reported.

To enable progress reporting, the ``radiomics.progressReporter`` variable should be set to a class object (NOT an
instance), which fits the following signature:

1. Accepts an iterable as the first positional argument and a keyword argument ('desc') specifying a label to display
2. Can be used in a 'with' statement (i.e. exposes a ``__enter__`` and ``__exit__`` function)
3. Is iterable (i.e. at least specifies an ``__iter__`` function, which iterates over the iterable passed at
   initialization)

It is also possible to create your own progress reporter. To achieve this, additionally specify a function ``__next__``,
and have the ``__iter__`` function return ``self``. The ``__next__`` function takes no arguments and returns a call to
the ``__next__`` function of the iterable (i.e. ``return self.iterable.__next__()``). Any prints/progress reporting
calls can then be inserted in this function prior to the return statement.

In ``radiomics\__init__.py`` a dummy progress reporter (``_DummyProgressReporter``) is defined, which is used when
calculating in full-python mode, but progress reporting is not enabled (verbosity > INFO) or the ``progressReporter``
variable is not set.

To design a custom progress reporter, the following code can be adapted and used as progressReporter::

    class MyProgressReporter(object):
        def __init__(self, iterable, desc=''):
            self.desc = desc  # A description is which describes the progress that is reported
            self.iterable = iterable  # Iterable is required

        # This function identifies the class as iterable and should return an object which exposes
        # the __next__ function that should be used to iterate over the object
        def __iter__(self):
            return self  # return self to 'intercept' the calls to __next__ to insert reporting code.

        def __next__(self):
            nextElement = self.iterable.__next__()
            # Insert custom progress reporting code here. This is called for every iteration in the loop
            # (once for each unique gray level in the ROI for GLCM and GLSZM)

            # By inserting after the call `self.iterable.__next__()` the function will exit before the
            # custom code is run when the stopIteration error is raised.
            return nextElement

        # This function is called when the 'with' statement is entered
        def __enter__(self):
            print (self.desc)  # Print out the description upon start of the loop
            return self  # The __enter__ function should return itself

        # This function is called when the 'with' statement is exited
        def __exit__(self, exc_type, exc_value, tb):
            pass  # If nothing needs to be closed or handled, so just specify 'pass'

------------------------------
Using feature classes directly
------------------------------

* This represents an example where feature classes are used directly, circumventing checks and preprocessing done by
  the radiomics feature extractor class, and is not intended as standard use.

* (LINUX) To run from source code, add pyradiomics to the environment variable PYTHONPATH (Not necessary when
  PyRadiomics is installed):

  *  ``setenv PYTHONPATH /path/to/pyradiomics/radiomics``

* Start the python interactive session:

  * ``python``

* Import the necessary classes::

     from radiomics import firstorder, glcm, imageoperations, shape, glrlm, glszm, getTestCase
     import SimpleITK as sitk
     import sys, os

* Set up a data directory variable::

    dataDir = '/path/to/pyradiomics/data'

* You will find sample data files brain1_image.nrrd and brain1_label.nrrd in that directory.

* Use SimpleITK to read a the brain image and mask::

     imageName, maskName = getTestCase('brain1', dataDir)
     image = sitk.ReadImage(imageName)
     mask = sitk.ReadImage(maskName)

* Calculate the first order features::

     firstOrderFeatures = firstorder.RadiomicsFirstOrder(image,mask)
     firstOrderFeatures.enableAllFeatures()  # On the feature class level, all features are disabled by default.
     firstOrderFeatures.calculateFeatures()
     for (key,val) in firstOrderFeatures.featureValues.items():
       print("\t%s: %s" % (key, val))

* See the :ref:`radiomics-features-label` section for more features that you can calculate.

--------------------------------
Additional points for attention
--------------------------------

Code style
***********

To keep the PyRadiomics code consistent and as readable as possible, some style rules are enforced. These are part of
the continuous testing and implemented using flake8. See also the ``.flake8`` configuration file in the root of the
repository. To aid in keeping a consistent code style, a ``.editorconfig`` file is provided in the root of the folder.

Module names should be lowercase, without underscores or spaces. Class names, function names and variables should be
declared using camelcase, with uppercase first letter for class names and lowercase first letter otherwise. Private
helper functions (which should not be included in the documentation) should be declared using a '_' prefix. This is
consistent with the python style for marking them as 'private', and will automatically exclude them from the generated
documentation.

Documentation
**************

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

If you make edits to the documentation, or if you want to test how documentation for the new classes you added is rendered, you can generate Sphinx documentation locally:

 1. ``pip install sphinx``
 2. ``pip install sphinx_rtd_theme``
 3. Run this command in the ``docs`` folder: ``make html``
 4. HTML version of the Sphinx documentation root will be in ``_build/html/index.html``

Testing
*******

To ensure consistency in the extraction provided by PyRadiomics, continuous testing is used to test the PyRadiomics
source code after each commit. These tests are defined in the test folder and used to run tests for the following
environments:

    - Python 2.7 64 bits (Windows, Linux and Mac)
    - Python 3.4 64 bits (Windows and Linux)
    - Python 3.5 64 bits (Windows and Linux)

.. note::

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
