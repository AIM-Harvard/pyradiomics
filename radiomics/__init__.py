# For convenience, import collection and numpy
# into the "pyradiomics" namespace

import collections  # noqa: F401
import inspect
import logging
import os
import pkgutil
import sys
import tempfile

import numpy  # noqa: F401
from six.moves import urllib

from . import imageoperations


def deprecated(func):
  """
  Decorator function to mark functions as deprecated. This is used to ensure deprecated feature functions are not
  added to the enabled features list when enabling 'all' features.
  """
  func._is_deprecated = True
  return func


def setVerbosity(level):
  """
  Change the amount of information PyRadiomics should print out during extraction. The lower the level, the more
  information is printed to the output (stderr).

  Using the ``level`` (Python defined logging levels) argument, the following levels are possible:

  - 60: Quiet mode, no messages are printed to the stderr
  - 50: Only log messages of level "CRITICAL" are printed
  - 40: Log messages of level "ERROR" and up are printed
  - 30: Log messages of level "WARNING" and up are printed
  - 20: Log messages of level "INFO" and up are printed
  - 10: Log messages of level "DEBUG" and up are printed (i.e. all log messages)

  By default, the radiomics logger is set to level "INFO" and the stderr handler to level "WARNING". Therefore a log
  storing the extraction log messages from level "INFO" and up can be easily set up by adding an appropriate handler to
  the radiomics logger, while the output to stderr will still only contain warnings and errors.

  .. note::

    This function assumes the handler added to the radiomics logger at initialization of the toolbox is not removed from
    the logger handlers and therefore remains the first handler.

  .. note::

    This does not affect the level of the logger itself (e.g. if verbosity level = 3, log messages with DEBUG level can
    still be stored in a log file if an appropriate handler is added to the logger and the logging level of the logger
    has been set to the correct level. *Exception: In case the verbosity is set to DEBUG, the level of the logger is
    also lowered to DEBUG. If the verbosity level is then raised again, the logger level will remain DEBUG.*
  """
  global logger, handler
  if level < 10:  # Lowest level: DEBUG
    level = 10
  if level > 60:  # Highest level = 50 (CRITICAL), level 60 results in a 'quiet' mode
    level = 60

  handler.setLevel(level)
  if handler.level < logger.level:  # reduce level of logger if necessary
    logger.setLevel(level)


def enableCExtensions(enabled=True):
  """
  By default, calculation of GLCM, GLRLM and GLSZM is done in C, using extension ``_cmatrices.py``

  If an error occurs during loading of this extension, a warning is logged and the extension is disabled,
  matrices are then calculated in python.
  The C extension can be disabled by calling this function as ``enableCExtensions(False)``, which forces the calculation
  of the matrices to full-python mode.

  Re-enabling use of C implementation is also done by this function, but if the extension is not loaded correctly,
  a warning is logged and matrix calculation is forced to full-python mode.
  """
  global _cMatsState, logger
  if enabled:
    # If extensions are not yet enabled (_cMatsState == 2), check whether they are loaded (_cMatsState == 1) and if so,
    # enable them. Otherwise, log a warning.
    if _cMatsState == 1:  # Extension loaded but not enabled
      logger.info("Enabling C extensions")
      _cMatsState = 2  # Enables matrix calculation in C
    elif _cMatsState == 0:  # _Extension not loaded correctly, do not enable matrix calculation in C and log warning
      logger.warning("C Matrices not loaded correctly, cannot calculate matrices in C")
  elif _cMatsState == 2:  # enabled = False, _cMatsState = 2: extensions currently enabled, disable them
    logger.info("Disabling C extensions")
    _cMatsState = 1


def cMatsEnabled():
  """
  Returns a boolean indicating whether or not the C extensions are enabled. This function is called by the feature
  classes to switch between C-enhanced calculation and full python mode.
  """
  return _cMatsState == 2


def getFeatureClasses():
  """
  Iterates over all modules of the radiomics package using pkgutil and subsequently imports those modules.

  Return a dictionary of all modules containing featureClasses, with modulename as key, abstract
  class object of the featureClass as value. Assumes only one featureClass per module

  This is achieved by inspect.getmembers. Modules are added if it contains a member that is a class,
  with name starting with 'Radiomics' and is inherited from :py:class:`radiomics.base.RadiomicsFeaturesBase`.

  This iteration only runs once (at initialization of toolbox), subsequent calls return the dictionary created by the
  first call.
  """
  global _featureClasses
  if _featureClasses is None:  # On first call, enumerate possible feature classes and import PyRadiomics modules
    _featureClasses = {}
    for _, mod, _ in pkgutil.iter_modules([os.path.dirname(__file__)]):
      if str(mod).startswith('_'):  # Skip loading of 'private' classes, these don't contain feature classes
        continue
      __import__('radiomics.' + mod)
      module = sys.modules['radiomics.' + mod]
      attributes = inspect.getmembers(module, inspect.isclass)
      for a in attributes:
        if a[0].startswith('Radiomics'):
          for parentClass in inspect.getmro(a[1])[1:]:  # only include classes that inherit from RadiomicsFeaturesBase
            if parentClass.__name__ == 'RadiomicsFeaturesBase':
              _featureClasses[mod] = a[1]
              break

  return _featureClasses


def getImageTypes():
  """
  Returns a list of possible image types (i.e. the possible filters and the "Original", unfiltered image type). This
  function finds the image types dynamically by matching the signature ("get<imageType>Image") against functions defined
  in :ref:`imageoperations <radiomics-imageoperations-label>`. Returns a list containing available image type names
  (<imageType> part of the corresponding function name).

  This iteration only occurs once, at initialization of the toolbox. Found results are stored and returned on subsequent
  calls.
  """
  global _imageTypes
  if _imageTypes is None:  # On first cal, enumerate possible input image types (original and any filters)
    _imageTypes = [member[3:-5] for member in dir(imageoperations)
                   if member.startswith('get') and member.endswith("Image")]

  return _imageTypes


def getTestCase(testCase, repoDirectory=None):
  """
  This function provides an image and mask for testing PyRadiomics. One of five test cases can be selected:

   - brain1
   - brain2
   - breast1
   - lung1
   - lung2

  If the repository is available locally (including all five test cases, the path to the root folder of the repository
  can be specified in ``repoDirectory``, preventing unnecessary downloads. If the repository is not found, or the
  repository does not contain the requested test case, PyRadiomics checks if it is run in development mode (directly
  from the source code in the repository), and if so, if it can find the test case relative to it's own location.

  If the requested test case could not be found in the repository, PyRadiomics downloads the test case from the GitHub
  repository and stores it in temporary files. If the test case was already downloaded, this is returned instead.

  Returns a tuple of two strings: ``(path/to/image.nrrd, path/to/mask.nrrd)``
  """
  global logger
  if testCase not in ['brain1', 'brain2', 'breast1', 'lung1', 'lung2']:
    logger.error('Testcase "%s" not recognized!', testCase)
    return None, None

  logger.debug('Getting test case %s', testCase)

  # Use test cases included in the repository. If PyRadiomics is run from an installed version, the location of the
  # repository needs to be specified
  logger.debug('Looking for test case in repository')
  if repoDirectory is not None:
    dataDir = os.path.join(repoDirectory, 'data')
    imageFile = os.path.join(dataDir, '%s_image.nrrd' % testCase)
    maskFile = os.path.join(dataDir, '%s_label.nrrd' % testCase)
    if os.path.isfile(imageFile) and os.path.isfile(maskFile):
      logger.debug('Test case found in repository (repository specified)')
      return imageFile, maskFile

  # No repository directory specified, check if running in development mode (code run from repository)
  logger.debug('Repository not specified or test case not found, checking if running in development mode')
  # This folder exists if radiomics is run from the repository:
  dataDir = os.path.join(os.path.basename(__file__), '..', 'data')
  imageFile = os.path.join(dataDir, '%s_image.nrrd' % testCase)
  maskFile = os.path.join(dataDir, '%s_label.nrrd' % testCase)
  if os.path.isfile(imageFile) and os.path.isfile(maskFile):
    logger.debug('Test case found in repository (running development mode)')
    return imageFile, maskFile

  # Data folder not found (most likely running from installed version). Check if test case has been downloaded.
  logger.debug('Test case or repository not found, checking temporary data')
  dataDir = os.path.join(tempfile.gettempdir(), 'pyradiomics', 'data')
  imageFile = os.path.join(dataDir, '%s_image.nrrd' % testCase)
  maskFile = os.path.join(dataDir, '%s_label.nrrd' % testCase)
  if os.path.isfile(imageFile) and os.path.isfile(maskFile):
    logger.debug('Test case already downloaded')
    return imageFile, maskFile

  logger.info("Test case not available locally, downloading test case...")

  # Testcase not found in temporary files, download them. First check if the folder is available
  if not os.path.isdir(os.path.join(tempfile.gettempdir(), 'pyradiomics')):
    os.mkdir(os.path.join(tempfile.gettempdir(), 'pyradiomics'))
  if not os.path.isdir(dataDir):
    logger.debug('Creating temporary directory: %s', dataDir)
    os.mkdir(dataDir)

  # Download the test case files (image and label)
  url = r'https://github.com/Radiomics/pyradiomics/raw/master/data/%s_%s.nrrd'
  try:
    urllib.request.urlretrieve(url % (testCase, 'image'), imageFile)
    urllib.request.urlretrieve(url % (testCase, 'label'), maskFile)
  except Exception:
    logger.error('Download failed!', exc_info=True)
    return None, None

  logger.info('Test case %s downloaded', testCase)

  return imageFile, maskFile


def getParameterValidationFiles():
  """
  Returns file locations for the parameter schema and custom validation functions, which are needed when validating
  a parameter file using ``PyKwalify.core``.
  This functions returns a tuple with the file location of the schema as first and python script with custom validation
  functions as second element.
  """
  dataDir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'schemas'))
  schemaFile = os.path.join(dataDir, 'paramSchema.yaml')
  schemaFuncs = os.path.join(dataDir, 'schemaFuncs.py')
  return schemaFile, schemaFuncs


class _DummyProgressReporter(object):
  """
  This class represents the dummy Progress reporter and is used for where progress reporting is implemented, but not
  enabled (when the progressReporter is not set or verbosity level > INFO).

  PyRadiomics expects that the _getProgressReporter function returns an object that takes an iterable and 'desc' keyword
  argument at initialization. Furthermore, it should be iterable, where it iterates over the iterable passed at
  initialization and it should be used in a 'with' statement.

  In this class, the __iter__ function redirects to the __iter__ function of the iterable passed at initialization.
  The __enter__ and __exit__ functions enable usage in a 'with' statement
  """
  def __init__(self, iterable, desc=''):
    self.desc = desc  # A description is not required, but is provided by PyRadiomics
    self.iterable = iterable  # Iterable is required

  def __iter__(self):
    return self.iterable.__iter__()  # Just iterate over the iterable passed at initialization

  def __enter__(self):
    return self  # The __enter__ function should return itself

  def __exit__(self, exc_type, exc_value, tb):
    pass  # Nothing needs to be closed or handled, so just specify 'pass'


def getProgressReporter(*args, **kwargs):
  """
  This function returns an instance of the progressReporter, if it is set and the logging level is defined at level INFO
  or DEBUG. In all other cases a dummy progress reporter is returned.

  To enable progress reporting, the progressReporter variable should be set to a class object (NOT an instance), which
  fits the following signature:

  1. Accepts an iterable as the first positional argument and a keyword argument ('desc') specifying a label to display
  2. Can be used in a 'with' statement (i.e. exposes a __enter__ and __exit__ function)
  3. Is iterable (i.e. at least specifies an __iter__ function, which iterates over the iterable passed at
     initialization).

  It is also possible to create your own progress reporter. To achieve this, additionally specify a function `__next__`,
  and have the `__iter__` function return `self`. The `__next__` function takes no arguments and returns a call to the
  `__next__` function of the iterable (i.e. `return self.iterable.__next__()`). Any prints/progress reporting calls can
  then be inserted in this function prior to the return statement.
  """
  global handler, progressReporter
  if progressReporter is not None and logging.NOTSET < handler.level <= logging.INFO:
    return progressReporter(*args, **kwargs)
  else:
    return _DummyProgressReporter(*args, **kwargs)

progressReporter = None

# 1. Set up logging
debugging = True
logger = logging.getLogger(__name__)  # 'radiomics'
logger.setLevel(logging.INFO)  # Set default level of logger to INFO to reflect most common setting for a log file

# Set up a handler to print out to stderr (controlled by setVerbosity())
handler = logging.StreamHandler()
# formatter = logging.Formatter("%(asctime)s %(levelname)s: %(message)s", "%Y-%m-%d %H:%M")  # Alternative format
formatter = logging.Formatter("%(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)
# force level=WARNING for stderr handler, in case logging default is set differently (issue 102)
setVerbosity(logging.WARNING)

# 2. Attempt to load and enable the C extensions. If this fails, revert to full-python mode
_cMatsState = 0  # Indicates status of C extensions: 0 = not loaded, 1 = loaded but not enabled, 2 = enabled
try:
  from radiomics import _cmatrices as cMatrices
  from radiomics import _cshape as cShape
  _cMatsState = 1
  enableCExtensions()
except Exception:
  logger.warning("Error loading C extensions, switching to python calculation:", exc_info=True)
  cMatrices = None  # set cMatrices to None to prevent an import error in the feature classes.
  cShape = None

# 3. Enumerate implemented feature classes and input image types available in PyRadiomics
_featureClasses = None
_imageTypes = None
getFeatureClasses()
getImageTypes()

# 4. Set the version using the versioneer scripts
from ._version import get_versions  # noqa: I202

__version__ = get_versions()['version']
del get_versions
