#!/usr/bin/env python

import versioneer

from setuptools.command.test import test as TestCommand
from skbuild import setup

with open('requirements.txt', 'r') as fp:
    requirements = list(filter(bool, (line.strip() for line in fp)))

with open('requirements-dev.txt', 'r') as fp:
    dev_requirements = list(filter(bool, (line.strip() for line in fp)))

class NoseTestCommand(TestCommand):
    """Command to run unit tests using nose driver after in-place build"""

    user_options = TestCommand.user_options + [
        ("args=", None, "Arguments to pass to nose"),
    ]

    def initialize_options(self):
        self.args = []
        TestCommand.initialize_options(self)

    def finalize_options(self):
        TestCommand.finalize_options(self)
        if self.args:
            self.args = __import__('shlex').split(self.args)

    def run_tests(self):
        # Run nose ensuring that argv simulates running nosetests directly
        nose_args = ['nosetests']
        nose_args.extend(self.args)
        __import__('nose').run_exit(argv=nose_args)

commands = versioneer.get_cmdclass()
commands['test'] = NoseTestCommand

setup_requires = []

setup(
    name='pyradiomics',

    url='http://github.com/Radiomics/pyradiomics#readme',

    author='pyradiomics community',
    author_email='pyradiomics@googlegroups.com',

    version=versioneer.get_version(),
    cmdclass=commands,

    packages=['radiomics'],
    zip_safe=False,
    data_files=[('data', ['data/paramSchema.yaml', 'data/schemaFuncs.py', 'data/brain1_image.nrrd', 'data/brain1_label.nrrd'])],

    description='Radiomics features library for python',

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
    test_suite='nose.collector',
    tests_require=dev_requirements,
)
