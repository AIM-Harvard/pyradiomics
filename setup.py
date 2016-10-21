from setuptools import setup

setup(name='pyradiomics',
      version='0.0.1dev',
      packages=['radiomics'],
      install_requires=['numpy==1.11.0',
                        'SimpleITK==0.9.1',
                        'nose-parameterized==0.5.0',
                        'tqdm==4.7.1',
                        'PyWavelets==0.4.0'],
      description='Radiomics features library for python',
      url='http://github.com/Radiomics/pyradiomics',
      author='pyradiomics community',
      author_email='pyradiomics@googlegroups.com',
      license='Slicer',
      zip_safe=False)
