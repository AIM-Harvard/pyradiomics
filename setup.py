#!/usr/bin/env python

from distutils import sysconfig
import platform

import numpy
from setuptools import Extension, setup
import versioneer

if platform.architecture()[0].startswith('32'):
  raise Exception('PyRadiomics requires 64 bits python')

commands = versioneer.get_cmdclass()
incDirs = [sysconfig.get_python_inc(), numpy.get_include()]

ext = [Extension("radiomics._cmatrices", ["radiomics/src/_cmatrices.c", "radiomics/src/cmatrices.c"],
                 include_dirs=incDirs),
       Extension("radiomics._cshape", ["radiomics/src/_cshape.c", "radiomics/src/cshape.c"],
                 include_dirs=incDirs)]

setup(
  name='pyradiomics',

  version=versioneer.get_version(),
  cmdclass=commands,

  packages=['radiomics', 'radiomics.scripts'],
  ext_modules=ext,
  zip_safe=False
)
