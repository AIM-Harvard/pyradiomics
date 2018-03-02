.. _radiomics-usage-label:

=====
Usage
=====

-----------------
Instruction Video
-----------------

.. raw:: html

        <iframe width="560" height="315" src="https://www.youtube.com/embed/ZF1rTtRW3eQ" frameborder="0" allowfullscreen></iframe>

-------
Example
-------

* PyRadiomics example code and data is available in the `Github repository <https://github.com/Radiomics/pyradiomics>`_

* The sample sample data is provided in ``pyradiomics/data``

* Use `jupyter <http://jupyter.org/>`_ to run the helloRadiomics example, located in ``pyradiomics/examples/Notebooks``

* Jupyter can also be used to run the example notebook as shown in the instruction video

  * The example notebook can be found in ``pyradiomics/examples/Notebooks``

  * The parameter file used in the instruction video is available in ``pyradiomics/examples/exampleSettings``

* If jupyter is not installed, run the python script alternative (``pyradiomics/examples/helloRadiomics.py``):

  * ``python helloRadiomics.py``

----------------
Command Line Use
----------------

* PyRadiomics can be used directly from the commandline via the entry point ``pyradiomics``. Depending on the input
  provided, PyRadiomics is run in either single-extraction or batch-extraction mode.

* To extract features from a single image and segmentation run::

    pyradiomics <path/to/image> <path/to/segmentation>

* To extract features from a batch run::

    pyradiomics <path/to/input> <path/to/output>

* The input file for batch processing is a CSV file where the first row is contains headers and each subsequent row
  represents one combination of an image and a segmentation and contains at least 2 elements: 1) path/to/image,
  2) path/to/mask. The headers specify the column names and **must** be "Image" and "Mask" for image and mask location,
  respectively (capital sensitive). Additional columns may also be specified, all columns are copied to the output in
  the same order (with calculated features appended after last column). To specify custom label values for each
  combination, a column "Label" can optionally be added, which specifies the desired extraction label for each
  combination. Values specified in this column take precedence over label values specified in the parameter file or on
  the commandline. If a row contains no value, the default (or globally customized) value is used instead.

  .. note::

    All headers should be unique and different from headers provided by PyRadiomics (``<filter>_<class>_<feature>``).

* Extraction can be customized by specifying a `parameter file <radiomics-parameter-file-label>` in the ``--param``
  argument and/or by specifying override settings (only `type 3 customization <radiomics-settings-label>`) in the
  ``--setting`` argument. Multiple overrides can be used by specifying ``--setting`` multiple times.

* For more information on the possible command line arguments, run::

    pyradiomics -h


---------------
Interactive Use
---------------

* (LINUX) To run from source code, add pyradiomics to the environment variable PYTHONPATH (Not necessary when
  PyRadiomics is installed):

  *  ``setenv PYTHONPATH /path/to/pyradiomics/radiomics``

* Start the python interactive session:

  * ``python``

* Import the necessary classes::

     from radiomics import featureextractor, getTestCase
     import six
     import sys, os

* Set up a pyradiomics directory variable::

    dataDir = '/path/to/pyradiomics'

* You will find sample data files brain1_image.nrrd and brain1_label.nrrd in that directory.

* Store the path of your image and mask in two variables::

    imageName, maskName = getTestCase('brain1', dataDir)

* Also store the path to the file containing the extraction settings::

    params = os.path.join(dataDir, "examples", "exampleSettings", "Params.yaml")

* Instantiate the feature extractor class with the parameter file::

    extractor = featureextractor.RadiomicsFeaturesExtractor(params)

* Calculate the features::

    result = extractor.execute(imageName, maskName)
    for key, val in six.iteritems(result):
      print("\t%s: %s" %(key, val))

* See the :ref:`feature extractor class<radiomics-featureextractor-label>` for more information on using this core class.

------------------------
PyRadiomics in 3D Slicer
------------------------

A convenient front-end interface is provided as the 'Radiomics' extension for 3D Slicer. It is available
`here <https://github.com/Radiomics/SlicerRadiomics>`_.

------------------------------
Using feature classes directly
------------------------------

* This represents an example where feature classes are used directly, circumventing checks and preprocessing done by
  the radiomics feature extractor class, and is not intended as standard use example.

* (LINUX) To run from source code, add pyradiomics to the environment variable PYTHONPATH (Not necessary when
  PyRadiomics is installed):

  *  ``setenv PYTHONPATH /path/to/pyradiomics/radiomics``

* Start the python interactive session:

  * ``python``

* Import the necessary classes::

     from radiomics import firstorder, glcm, imageoperations, shape, glrlm, glszm, getTestCase
     import SimpleITK as sitk
     import six
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
     for (key,val) in six.iteritems(firstOrderFeatures.featureValues):
       print("\t%s: %s" % (key, val))

* See the :ref:`radiomics-features-label` section for more features that you can calculate.

.. _radiomics-logging-label:

------------------
Setting Up Logging
------------------

PyRadiomics features extensive logging to help track down any issues with the extraction of features.
By default PyRadiomics logging reports messages of level WARNING and up (reporting any warnings or errors that occur),
and prints this to the output (stderr). By default, PyRadiomics does not create a log file.

To change the amount of information that is printed to the output, use :py:func:`~radiomics.setVerbosity` in interactive
use and the optional ``--verbosity`` argument in commandline use.

When using PyRadiomics in interactive mode, enable storing the PyRadiomics logging in a file by adding an appropriate
handler to the pyradiomics logger::

    import radiomics

    log_file = 'path/to/log_file.txt'
    handler = logging.FileHandler(filename=log_file, mode='w')  # overwrites log_files from previous runs. Change mode to 'a' to append.
    formatter = logging.Formatter("%(levelname)s:%(name)s: %(message)s")  # format string for log messages
    handler.setFormatter(formatter)
    radiomics.logger.addHandler(handler)

    # Control the amount of logging stored by setting the level of the logger. N.B. if the level is higher than the
    # Verbositiy level, the logger level will also determine the amount of information printed to the output
    radiomics.logger.setLevel(logging.DEBUG)

To store a log file when running pyradiomics from the commandline, specify a file location in the optional
``--log-file`` argument. The amount of logging that is stored is controlled by the ``--log-level`` argument
(default level INFO and up).
