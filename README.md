
[![Appveyor](https://ci.appveyor.com/api/projects/status/tw69xbbeyluk7fl7/branch/master?svg=true)](https://ci.appveyor.com/project/Radiomics/pyradiomics/branch/master)

[![Circle CI](https://circleci.com/gh/Radiomics/pyradiomics.svg?style=svg&circle-token=a4748cf0de5fad2c12bc93a485282378551c3584)](https://circleci.com/gh/Radiomics/pyradiomics)

[![Travis CI](https://travis-ci.org/Radiomics/pyradiomics.svg?branch=master)](https://travis-ci.org/Radiomics/pyradiomics)

# pyradiomics v1.2.0

## Radiomics feature extraction in Python

This is an open-source python package for the extraction of Radiomics features from medical imaging.

With this package we aim to establish a reference standard for Radiomic Analysis, and provide a tested and maintained
open-source platform for easy and reproducible Radiomic Feature extraction. By doing so, we hope to increase awareness
of radiomic capabilities and expand the community.

The platform supports both the feature extraction in 2D and 3D. **Not intended for clinical use.**

**If you publish any work which uses this package, please cite the following publication:**
*Joost JM van Griethuysen, Andriy Fedorov, Chintan Parmar, Ahmed Hosny, Nicole Aucoin, Vivek Narayan, Regina GH 
Beets-Tan, Jean-Christophe Fillion-Robin, Steve Pieper, Hugo JWL Aerts, “Computational Radiomics System to Decode the
Radiographic Phenotype”; Accepted Cancer Research 2017*

### Feature Classes
Currently supports the following feature classes:

 - First Order Statistics
 - Shape-based
 - Gray Level Cooccurence Matrix (GLCM)
 - Gray Level Run Length Matrix (GLRLM)
 - Gray Level Size Zone Matrix (GLSZM)

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

    python -m pip install -r requirements.txt
    python setup.py install

Detailed installation instructions, as well as instructions for installing PyRadiomics on Windows are available in the 
[documentation](http://pyradiomics.readthedocs.io/en/latest/installation.html).

### Docker

PyRadiomics also supports [Dockers](https://www.docker.com/).  Currently, the only docker available is a [Jupyter notebook](http://jupyter.org/) with PyRadiomics pre-installed with example Notebooks. To build the Docker:

    docker build -t radiomics/notebook .

The `radiomics/notebook` Docker has an exposed volume (`/data`) that can be mapped to the host system directory.  For example, to mount the current directory:

    docker run --rm -it --publish 8888:8888 -v `pwd`:/data radiomics/notebook

or for a less secure notebook, skip the randomly generated token

    docker run --rm -it --publish 8888:8888 -v `pwd`:/data radiomics/notebook start-notebook.sh --NotebookApp.token=''

and open the local webpage at http://localhost:8888/ with the current directory at http://localhost:8888/tree/data.

### Usage

PyRadiomics can be easily used in a Python script through the `featureextractor`
module. Furthermore, PyRadiomics provides two commandline scripts, `pyradiomics`
and `pyradiomicsbatch`, for single image extraction and batchprocessing, respectively.
Finally, a convenient front-end interface is provided as the 'Radiomics'
extension for 3D Slicer, available [here](https://github.com/Radiomics/SlicerRadiomics).

### 3rd-party packages used in pyradiomics:

 - SimpleITK (Image loading and preprocessing)
 - numpy (Feature calculation)
 - PyWavelets (Wavelet filter)
 - pykwalify (Enabling yaml parameters file checking)
 - six (Python 3 Compatibility)

See also the [requirements file](requirements.txt).

### 3D Slicer
PyRadiomics is also available as an [extension](https://github.com/Radiomics/SlicerRadiomics) to [3D Slicer](slicer.org). 
Download and install the 3D slicer [nightly build](http://download.slicer.org/), the extension is then available in the
extension manager under "SlicerRadiomics".

### License
This package is covered by the open source [3-clause BSD License](LICENSE.txt).

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

**This work was supported in part by the US National Cancer Institute grant 
5U24CA194354, QUANTITATIVE RADIOMICS SYSTEM DECODING THE TUMOR PHENOTYPE.**
