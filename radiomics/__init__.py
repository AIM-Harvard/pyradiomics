# For convenience, import collection and numpy
# into the "pyradiomics" namespace

import collections  # noqa: F401
import inspect
import logging
import os
import pkgutil
import sys
import traceback

import numpy  # noqa: F401

from . import base, imageoperations

if sys.version_info < (2, 6, 0):
  raise ImportError("pyradiomics > 0.9.7 requires python 2.6 or later")


def debug(debug_on=True):
  """
  Control level of logger and stderr output of the toolbox. By default, this output reflects module hierarchy, as child
  loggers are created by module. This is achieved by the following line in base.py:
  ``self.logger = logging.getLogger(self.__module__)``. To use same instance in each module, set
  ``self.logger=logging.getLogger('radiomics')``.

  At command line, turn on debugging output to stderr for all pyradiomics functions with:

  ``import radiomics``\n
  ``radiomics.debug()``

  This set the level of both the logger and the handler for output to stderr to level = DEBUG. If level of logger is
  already at level "DEBUG" or "NOTSET", level of logger is not changed.

  Turn off debugging with (only changes the level of the handler for output to stderr):

  ``radiomics.debug(False)``

  By default, the radiomics logger is set to level "INFO" and the stderr handler to level "WARNING". Therefore a log
  storing the extraction log messages from level "INFO" and up can be easily set up by adding an appropriate handler to
  the radiomics logger.
  """
  global logger, debugging
  if debug_on:
    setVerbosity(logging.DEBUG)  # set the output to stderr to DEBUG
    debugging = True
  else:
    setVerbosity(logging.WARNING)  # set the output to stderr to WARNING
    debugging = False


def setVerbosity(level):
  """
  Assumes the handler added to the radiomics logger at initialization of the toolbox is not removed from the logger
  handlers.

  Using the ``level`` (Python defined logging levels) argument, determine how much PyRadiomics should print out to the
  stderr, the following levels are possible:

  - 60: Quiet mode, no messages are printed to the stderr
  - 50: Only log messages of level "CRITICAL" are printed
  - 40: Log messages of level "ERROR" and up are printed
  - 30: Log messages of level "WARNING" and up are printed
  - 20: Log messages of level "INFO" and up are printed
  - 10: Log messages of level "DEBUG" and up are printed (i.e. all log messages)

  **N.B. This does not affect the level of the logger itself (e.g. if verbosity level = 3, log messages with DEBUG level
  can still be stored in a log file if an appropriate handler is added to the logger and the logging level of the logger
  has been set to the correct level**
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
    for _, mod, _ in pkgutil.iter_modules([os.path.dirname(base.__file__)]):
      if str(mod).startswith('_'):  # Skip loading of 'private' classes, these don't contain feature classes
        continue
      __import__('radiomics.' + mod)
      module = sys.modules['radiomics.' + mod]
      attributes = inspect.getmembers(module, inspect.isclass)
      for a in attributes:
        if a[0].startswith('Radiomics'):
          if base.RadiomicsFeaturesBase in inspect.getmro(a[1])[1:]:
            _featureClasses[mod] = a[1]

  return _featureClasses


def getInputImageTypes():
  """
  Returns a list of possible input image types. This function finds the image types dynamically by matching the
  signature ("get<inputImage>Image") against functions defined in :ref:`imageoperations
  <radiomics-imageoperations-label>`. Returns a list containing available input image names (<inputImage> part of the
  corresponding function name).

  This iteration only occurs once, at initialization of the toolbox. Found results are stored and returned on subsequent
  calls.
  """
  global _inputImages
  if _inputImages is None:  # On first cal, enumerate possible input image types (original and any filters)
    _inputImages = [member[3:-5] for member in dir(imageoperations)
                    if member.startswith('get') and member.endswith("Image")]

  return _inputImages


def getProgressFunctions():
  """
  WIP
  """
  global initProgress, reportProgress, closeProgress
  return initProgress, reportProgress, closeProgress

# Instantiate 3 variables to hold the functions needed to report progress
initProgress = None
reportProgress = None
closeProgress = None

debugging = True
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # Set default level of logger to INFO to reflect most common setting for a log file

handler = logging.StreamHandler()
# formatter = logging.Formatter("%(asctime)s %(levelname)s: %(message)s", "%Y-%m-%d %H:%M")
formatter = logging.Formatter("%(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)
debug(False)  # force level=WARNING for stderr handler, in case logging default is set differently (issue 102)

_featureClasses = None
_inputImages = None

# Indicates status of C extensions: 0 = not loaded, 1 = loaded but not enabled, 2 = enabled
_cMatsState = 0

try:
  from radiomics import _cmatrices as cMatrices
  from radiomics import _cshape as cShape
  _cMatsState = 1
  enableCExtensions()
except Exception:
  logger.warning("Error loading C extensions, switching to python calculation:\n%s", traceback.format_exc())
  cMatrices = None  # set cMatrices to None to prevent an import error in the feature classes.
  cShape = None

from ._version import get_versions

__version__ = get_versions()['version']
del get_versions

getFeatureClasses()
getInputImageTypes()
