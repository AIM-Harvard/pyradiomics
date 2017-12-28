============
Installation
============

There are three ways you can use pyradiomics:
1. Install from source
2. Use 3D Slicer Radiomics extension
3. Use pyradiomics Docker

------------
1. Install from source
------------

* Ensure you have the version control system ``git`` installed on your machine.

* Ensure that you have ``python`` installed on your machine, at least version 2.7 or 3.4.

* Clone the repository:

  * ``git clone git://github.com/Radiomics/pyradiomics``

.. _installation-label:

* For unix like systems (MacOSX, linux):

  * ``cd pyradiomics``
  * ``python -m pip install -r requirements.txt``
  * ``python setup.py install``
  
  To use your build for interactive use and development:
  * ``python setup.py develop``

  * If you don't have sudo/admin rights on your machine, you need to locally install numpy, nose, tqdm, PyWavelets, SimpleITK (specified in requirements.txt).
    In a bash shell::

      pip install --user --upgrade pip
      export PATH=$HOME/.local/bin:$PATH
      pip install --user -r requirements.txt
      export PYTHONPATH=$HOME/.local/lib64/python2.7/site-packages

* For Windows:

  * ``cd pyradiomics``
  * ``python -m pip install -r requirements.txt``
  * ``python setup.py install``

---
2. Use 3D Slicer Radiomics extension
---

3D Slicer is a free open source research platform for medical image computing. Learn more and download 3D Slicer binary for your platform here: http://slicer.org.

Once installed, you can use 3D Slicer ExtensionManager to install Radiomics extension, which provides a graphical user interface to the pyradiomics library. The advantage of
using pyradiomics from 3D Slicer is that you can view images and segmentations, you can import existing segmentations and confirm their quality, or you can use the variety
of tools in 3D Slicer to automate your segmentation tasks.

More detailed instructions about installing 3D Slicer Radiomics extension are available here: https://github.com/Radiomics/SlicerRadiomics

---
3. Use pyradiomics Docker
---

This approach may be preferred if you are interested in using pyradiomics from the command line, but have difficulties installing the library on your system. 

First, you will need to install Docker on your system, if it is not installed already. You can follow the instructions on this page to do this.

Once Docker is installed, you can issue ``docker pull radiomics/pyradiomics:CLI`` command in the shell to download the pyradiomics Docker image.
After that you can invoke pyradiomics tool as follows:

  ``docker run radiomics/pyradiomics:CLI --help``

Docker containers cannot directly access the filesystem of the host. In order to pass files as arguments to pyradiomics and to access files that converters create, 
an extra step is required to specify which directories will be used for file exchange using the -v argument:

  ``-v <HOST_DIR>:<CONTAINER_DIR>``

The argument above will make the ``HOST_DIR`` path available within the container at ``CONTAINER_DIR`` location. The files that will be read or written by the 
converter run from the docker container should be referred to via the ``CONTAINER_DIR`` path.
