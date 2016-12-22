import versioneer

from setuptools import setup

with open('requirements.txt', 'r') as fp:
  requirements = list(filter(bool, (line.strip() for line in fp)))

setup(name='pyradiomics',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      packages=['radiomics'],
      data_files=[
        ('data', ['data/paramSchema.yaml', 'data/schemaFuncs.py', 'data/brain1_image.nrrd', 'data/brain1_label.nrrd'])],
      description='Radiomics features library for python',
      url='http://github.com/Radiomics/pyradiomics',
      author='pyradiomics community',
      author_email='pyradiomics@googlegroups.com',
      license='Slicer',
      zip_safe=False,
      install_requires=requirements)
