# pyradiomics v3.0.1

## Build Status

| Linux                          | macOS                         | Windows                       |
|--------------------------------|-------------------------------|-------------------------------|
| [![][circleci]][circleci-lnk]  | [![][travisci]][travisci-lnk] | [![][appveyor]][appveyor-lnk] |


[appveyor]: https://ci.appveyor.com/api/projects/status/tw69xbbeyluk7fl7/branch/master?svg=true
[appveyor-lnk]: https://ci.appveyor.com/project/Radiomics/pyradiomics/branch/master

[circleci]: https://circleci.com/gh/Radiomics/pyradiomics.svg?style=svg&circle-token=a4748cf0de5fad2c12bc93a485282378551c3584
[circleci-lnk]: https://circleci.com/gh/Radiomics/pyradiomics

[travisci]: https://travis-ci.org/Radiomics/pyradiomics.svg?branch=master
[travisci-lnk]: https://travis-ci.org/Radiomics/pyradiomics

## Radiomics feature extraction in Python
This is an open-source python package for the extraction of Radiomics features from medical imaging.

With this package we aim to establish a reference standard for Radiomic Analysis, and provide a tested and maintained
open-source platform for easy and reproducible Radiomic Feature extraction. By doing so, we hope to increase awareness
of radiomic capabilities and expand the community.

The platform supports both the feature extraction in 2D and 3D and can be used to calculate single values per feature
for a region of interest ("segment-based") or to generate feature maps ("voxel-based"). 

**Not intended for clinical use.**

**If you publish any work which uses this package, please cite the following publication:**
*van Griethuysen, J. J. M., Fedorov, A., Parmar, C., Hosny, A., Aucoin, N., Narayan, V., Beets-Tan, R. G. H.,
Fillion-Robin, J. C., Pieper, S.,  Aerts, H. J. W. L. (2017). Computational Radiomics System to Decode the Radiographic
Phenotype. Cancer Research, 77(21), e104–e107. https://doi.org/10.1158/0008-5472.CAN-17-0339*

### Join the Community!
Join the PyRadiomics community on google groups [here](https://groups.google.com/forum/#!forum/pyradiomics).

### Feature Classes
Currently supports the following feature classes:

 - First Order Statistics
 - Shape-based (2D and 3D)
 - Gray Level Cooccurence Matrix (GLCM)
 - Gray Level Run Length Matrix (GLRLM)
 - Gray Level Size Zone Matrix (GLSZM)
 - Gray Level Dependece Matrix (GLDM)
 - Neighboring Gray Tone Difference Matrix (NGTDM)

### Filter Classes
Aside from the feature classes, there are also some built-in optional filters:

- Laplacian of Gaussian (LoG, based on SimpleITK functionality)
- Wavelet (using the PyWavelets package)
- Square
- Square Root
- Logarithm
- Exponential
- Gradient (Magnitude)
- Local Binary Pattern (LBP) 2D / 3D

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
PyRadiomics is OS independent and compatible with Python >= 3.5. Pre-built binaries are available on
PyPi and Conda. To install PyRadiomics, ensure you have python
installed and run:

    `python -m pip install pyradiomics`

Detailed installation instructions, as well as instructions for building PyRadiomics from source, are available in the 
[documentation](http://pyradiomics.readthedocs.io/en/latest/installation.html).

### Docker
PyRadiomics also supports [Dockers](https://www.docker.com/).  Currently, 2 dockers are available:

The first one is a [Jupyter notebook](http://jupyter.org/) with PyRadiomics pre-installed with example Notebooks. 

To get the Docker:

    docker pull radiomics/pyradiomics:latest

The `radiomics/notebook` Docker has an exposed volume (`/data`) that can be mapped to the host system directory.  For example, to mount the current directory:

    docker run --rm -it --publish 8888:8888 -v `pwd`:/data radiomics/notebook

or for a less secure notebook, skip the randomly generated token

    docker run --rm -it --publish 8888:8888 -v `pwd`:/data radiomics/notebook start-notebook.sh --NotebookApp.token=''

and open the local webpage at http://localhost:8888/ with the current directory at http://localhost:8888/tree/data.

The second is a docker which exposes the PyRadiomics CLI interface. To get the CLI-Docker:

    docker pull radiomics/pyradiomics:CLI

You can then use the PyRadiomics CLI as follows:

    docker run radiomics/pyradiomics:CLI --help

For more information on using docker, see
[here](https://pyradiomics.readthedocs.io/en/latest/installation.html#use-pyradiomics-docker)

### Usage
PyRadiomics can be easily used in a Python script through the `featureextractor`
module. Furthermore, PyRadiomics provides a commandline script, `pyradiomics`, for both single image extraction and 
batchprocessing. Finally, a convenient front-end interface is provided as the 'Radiomics'
extension for 3D Slicer, available [here](https://github.com/Radiomics/SlicerRadiomics).

### 3rd-party packages used in pyradiomics:
 - SimpleITK (Image loading and preprocessing)
 - numpy (Feature calculation)
 - PyWavelets (Wavelet filter)
 - pykwalify (Enabling yaml parameters file checking)
 - six (Python 3 Compatibility)
 - scipy (Only for LBP filter, install separately to enable this filter)
 - scikit-image (Only for LBP filter, install separately to enable this filter)
 - trimesh (Only for LBP filter, install separately to enable this filter)

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
We are happy to help you with any questions. Please contact us on the [Radiomics community section of the 3D Slicer Discourse](https://discourse.slicer.org/c/community/radiomics/23).

We welcome contributions to PyRadiomics. Please read the [contributing guidelines](CONTRIBUTING.rst) on how to
contribute to PyRadiomics.

**This work was supported in part by the US National Cancer Institute grant 
5U24CA194354, QUANTITATIVE RADIOMICS SYSTEM DECODING THE TUMOR PHENOTYPE.**
