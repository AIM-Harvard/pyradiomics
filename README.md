
[![Appveyor](https://ci.appveyor.com/api/projects/status/tw69xbbeyluk7fl7/branch/master?svg=true)](https://ci.appveyor.com/project/Radiomics/pyradiomics/branch/master)

[![Circle CI](https://circleci.com/gh/Radiomics/pyradiomics.svg?style=svg&circle-token=a4748cf0de5fad2c12bc93a485282378551c3584)](https://circleci.com/gh/Radiomics/pyradiomics)

[![Travis CI](https://travis-ci.org/Radiomics/pyradiomics.svg?branch=master)](https://travis-ci.org/Radiomics/pyradiomics)

# pyradiomics v1.1.0

## Radiomics feature extraction in Python

This is an open-source python package for the extraction of Radiomics features from 2D and 3D images and 
binary masks.

Image loading and preprocessing (e.g. resampling and cropping) are first done using `SimpleITK`. 
Then, loaded data are converted into numpy arrays for further calculation using feature classes
outlined below.

### Feature Classes
Currently supports the following feature classes:

 - First Order Statistics
 - Shape-based
 - [Gray Level Cooccurence Matrix](https://en.wikipedia.org/wiki/Co-occurrence_matrix) (GLCM)
 - [Gray Level Run Length Matrix](http://www.insight-journal.org/browse/publication/231) (GLRLM)
 - [Gray Level Size Zone Matrix](https://en.wikipedia.org/wiki/Gray_level_size_zone_matrix) (GLSZM)

### Filter Classes
Aside from the feature classes, there are also some built-in optional filters:

- Laplacian of Gaussian (LoG, based on SimpleITK functionality)
- Wavelet (using the PyWavelets package)
- Square
- Square Root
- Logarithm
- Exponential

### Supporting reproducible extraction
Aside from calculating features, the pyradiomics package includes provenance information in the
output. This information contains information on used image and mask, as well as applied settings
and filters, thereby enabling fully reproducible feature extraction.

### Documentation

For more information, see the sphinx generated documentation available [here](http://pyradiomics.readthedocs.io/).

Alternatively, you can generate the documentation by checking out the master branch and running from the root directory:

    python setup.py build_sphinx

The documentation can then be viewed in a browser by opening `PACKAGE_ROOT\build\sphinx\html\index.html`. 

Furthermore, an instruction video is available [here](http://radiomics.io/pyradiomics.html).

### Installation

PyRadiomics is OS independent and compatible with both Python 2.7 and Python >=3.4.
To install this package on unix like systems run the following commands from the root directory:

    sudo python -m pip install -r requirements.txt
    sudo python setup.py install

Detailed installation instructions, as well as instructions for installing PyRadiomics on Windows are available in the 
[documentation](http://pyradiomics.readthedocs.io/en/latest/installation.html).

### Usage

PyRadiomics can be easily used in a Python script through the `featureextractor`
module. Furthermore, PyRadiomics provides two commandline scripts, `pyradiomics`
and `pyradiomicsbatch`, for single image extraction and batchprocessing, respectively.
Finally, a convenient front-end interface is provided as the 'Radiomics'
extension for 3D Slicer, available [here](https://github.com/Radiomics/SlicerRadiomics).

### Citation 
If you publish any work which uses this package, please cite the following publication:

Joost J.M. van Griethuysen et al, “Computational Radiomics System to Decode the Radiographic Phenotype”; Submitted

### 3rd-party packages used in pyradiomics:

 - SimpleITK
 - numpy
 - PyWavelets (Wavelet filter)
 - pykwalify (Enabling yaml parameters file checking)
 - tqdm (Progressbar)
 - sphinx (Generating documentation)
 - sphinx_rtd_theme (Template for documentation)
 - nose-parameterized (Testing)

See also the [requirements file](requirements.txt).

### WIP
 - Implementation of this package as an [extension](https://github.com/Radiomics/SlicerRadiomics) to [3D Slicer](slicer.org)

### License
This package is covered by the [3D Slicer License](LICENSE.txt).

**This work was supported in part by the US National Cancer Institute grant 
5U24CA194354, QUANTITATIVE RADIOMICS SYSTEM DECODING THE TUMOR PHENOTYPE.**

### Developers
 - [Joost van Griethuysen](https://github.com/JoostJM)<sup>1,3,4</sup>
 - [Andriy Fedorov](https://github.com/fedorov)<sup>2</sup>
 - [Nicole Aucoin](https://github.com/naucoin)<sup>2</sup>
 - [Jean-Christophe Fillion-Robin](https://github.com/jcfr)<sup>5</sup>
 - [Ahmed Hosny](https://github.com/ahmedhosny)<sup>1</sup>
 - [Steve Pieper](https://github.com/pieper)<sup>6</sup>
 - [Hugo Aerts (PI)](https://github.com/hugoaerts)<sup>1,2</sup>
 
<sup>1</sup>Department of Radiation Oncology, Dana-Farber Cancer Institute, Brigham and Women's Hospital, Harvard Medical School, Boston, MA,
<sup>2</sup>Department of Radiology, Brigham and Women's Hospital, Harvard Medical School, Boston, MA,
<sup>3</sup>Department of Radiology, Netherlands Cancer Institute, Amsterdam, The Netherlands, 
<sup>4</sup>GROW-School for Oncology and Developmental Biology, Maastricht University Medical Center, Maastricht, The Netherlands,
<sup>5</sup>Kitware,
<sup>6</sup>Isomics

### Contact

We are happy to help you with any questions. Please contact us on the [pyradiomics email list](https://groups.google.com/forum/#!forum/pyradiomics).

We welcome contributions to PyRadiomics. Please read the [contributing guidelines](CONTRIBUTING.md) on how to contribute
to PyRadiomics.

