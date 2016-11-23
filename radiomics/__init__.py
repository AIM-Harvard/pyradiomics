import sys
import traceback

if sys.version_info < (2, 6, 0):
  raise ImportError("pyradiomics > 0.9.7 requires python 2.6 or later")
in_py3 = sys.version_info[0] > 2

import logging

def pythonMatrixCalculation(pyMatEnabled = False):
  """
  By default, calculation of matrices is done in C, using extension ``_cmatrices.py``

  If an error occurs during loading of this extension, a warning is logged and the extension is disabled,
  matrices are then calculated in python.
  Calculation in python can be forced by calling this function, specifying ``pyMatEnabled = True``

  Re-enabling use of C implementation is also done by this function, but if extension is not loaded correctly,
  a warning is logged and matrix calculation uses python.
  """
  global cMatsEnabled
  if pyMatEnabled:
    cMatsEnabled = False
  elif _cMatLoaded:
    cMatsEnabled = True
  else:
    logger.warning("C Matrices not loaded correctly, cannot calculate matrices in C")


def debug(debug_on=True):
  """
  Set up logging system for the whole package.
  By default, module hierarchy is reflected in log, as child loggers are created by module
  This is achieved by the following line in base.py: ``self.logger = logging.getLogger(self.__module__)``
  To use same instance in each module, set ``self.logger=logging.getLogger('radiomics')``.

  At command line, turn on debugging for all pyradiomics functions with:

  ``import radiomics``\n
  ``radiomics.debug()``

  Turn off debugging with:

  ``radiomics.debug(False)``
  """
  global logger, debugging
  if debug_on:
    logger.setLevel(logging.DEBUG)
    debugging = True
  else:
    logger.setLevel(logging.WARNING)
    debugging = False


debugging = True
logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
# formatter = logging.Formatter("%(asctime)s %(levelname)s: %(message)s", "%Y-%m-%d %H:%M")
formatter = logging.Formatter("%(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)
debug(False)  # force level=WARNING, in case logging default is set differently (issue 102)

try:
  import _cmatrices as cMatrices
  cMatsEnabled = True
  _cMatLoaded = True
except Exception as e:
  logger.warning("Error loading C Matrices, switching to python calculation\n%s", traceback.format_exc())
  cMatsEnabled = False
  _cMatLoaded = False

# For convenience, import the most used packages into the "pyradiomics" namespace
import collections, numpy

from ._version import get_versions

__version__ = get_versions()['version']
del get_versions
