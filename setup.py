#!/usr/bin/env python

import numpy
import versioneer

from distutils import sysconfig
from setuptools import setup, Extension


incDirs = [sysconfig.get_python_inc(), numpy.get_include()]
setup(
    name='pyradiomics',

    url='http://github.com/Radiomics/pyradiomics#readme',

    author='pyradiomics community',
    author_email='pyradiomics@googlegroups.com',

    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),

    packages=['radiomics'],
    zip_safe=False,
    data_files=[('data', ['data/paramSchema.yaml', 'data/schemaFuncs.py', 'data/brain1_image.nrrd', 'data/brain1_label.nrrd'])],

    description='Radiomics features library for python',
    ext_modules=[Extension("_cmatrices",
                           ["radiomics/src/_cmatrices.c", "radiomics/src/cmatrices.c"],
                           include_dirs=incDirs),
                 Extension("_cshape", ["radiomics/src/_cshape.c", "radiomics/src/cshape.c"],
                           include_dirs=incDirs)],


    license='Slicer',

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'License :: Slicer',
        'Operating System :: OS Independent',
        'Programming Language :: C',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],

    keywords='radiomics cancerimaging medicalresearch',

    install_requires=requirements,
    setup_requires=setup_requires,
    tests_require=dev_requirements,
)
