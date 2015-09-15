import sys
if sys.version_info < (2, 6, 0):
    raise ImportError("pydicom > 0.9.7 requires python 2.6 or later")
in_py3 = sys.version_info[0] > 2

# Set up logging system for the whole package.
# In each module, set logger=logging.getLogger('pyradiomics')  and the same instance
#     will be used by all
# At command line, turn on debugging for all pyradiomics functions with:
#        import radiomics
#        radiomics.debug()
#  Turn off debugging with
#       radiomics.debug(False)
import logging


def debug(debug_on=True):
    global logger, debugging
    if debug_on:
        logger.setLevel(logging.DEBUG)
        debugging = True
    else:
        logger.setLevel(logging.WARNING)
        debugging = False

logger = logging.getLogger('pyradiomics')
handler = logging.StreamHandler()
# formatter = logging.Formatter("%(asctime)s %(levelname)s: %(message)s", "%Y-%m-%d %H:%M") #'%(asctime)s %(levelname)s %(message)s'
formatter = logging.Formatter("%(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)
debug(False)  # force level=WARNING, in case logging default is set differently (issue 102)

# For convenience, import the most used packages into the "pyradiomics" namespace
import collections, numpy

__version__ = "0.0.1"
__version_info__ = (0, 0, 1)

print 'HELLO'
