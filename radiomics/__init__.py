from __future__ import print_function, unicode_literals, division, absolute_import
import sys

if sys.version_info < (2, 6, 0):
  raise ImportError("pyradiomics > 0.9.7 requires python 2.6 or later")
in_py3 = sys.version_info[0] > 2

if in_py3:
  c_str_type = str
  safe_xrange = lambda *x, **kwargs: iter(range(*x, **kwargs))
else:
  c_str_type = basestring
  safe_xrange = xrange


import pkgutil
import inspect
import os
import logging

from . import base, imageoperations


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


debugging = True
logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
# formatter = logging.Formatter("%(asctime)s %(levelname)s: %(message)s", "%Y-%m-%d %H:%M")
formatter = logging.Formatter("%(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)
debug(False)  # force level=WARNING, in case logging default is set differently (issue 102)

_featureClasses = None
_inputImages = None

# For convenience, import the most used packages into the "pyradiomics" namespace
import collections, numpy  # noqa: F401

from ._version import get_versions

__version__ = get_versions()['version']
del get_versions

getFeatureClasses()
getInputImageTypes()
