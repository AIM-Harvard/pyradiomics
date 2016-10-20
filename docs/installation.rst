============
Installation
============

------------
Get the code
------------

* Ensure you have the version control system ``git`` installed on your machine.

* Ensure that you have ``python`` installed on your machine, at least version 2.7.

* Clone the repository:

  * ``git clone git://github.com/Radiomics/pyradiomics``

.. _installation-label:

---------------------------
Installation on your system
---------------------------

* For unix like systems (MacOSX, linux):

  * ``cd pyradiomics``
  * ``sudo python setup.py install``

  * If you don't have sudo/admin rights on your machine, you need to locally install numpy, nose, tqdm, PyWavelets, SimpleITK. In a bash shell::

      pip install --user --upgrade pip
      export PATH=$HOME/.local/bin:$PATH
      pip install --user numpy
      pip install --user nose
      pip install --user nose-parameterized
      pip install --user tqdm
      export PYTHONPATH=$HOME/.local/lib64/python2.7/site-packages
      git clone https://github.com/nigma/pywt.git && cd pywt && python setup.py install --prefix=$HOME/.local && cd ..
      pip install --user SimpleITK

    * If the installation of SimpleITK fails (newer versions not available on the default servers, you can get it manually from `sourceforge <https://sourceforge.net/projects/simpleitk/files/SimpleITK/>`_

      * For linux::

          wget 'https://sourceforge.net/projects/simpleitk/files/SimpleITK/0.10.0/Python/SimpleITK-0.10.0-1-cp27-cp27m-manylinux1_x86_64.whl'
          pip install --user 'SimpleITK-0.10.0-1-cp27-cp27m-manylinux1_x86_64.whl'

      * For Mac::

          wget 'https://sourceforge.net/projects/simpleitk/files/SimpleITK/0.10.0/Python/SimpleITK-0.10.0-cp27-cp27m-macosx_10_6_intel.whl'
          pip install --user 'SimpleITK-0.10.0-cp27-cp27m-macosx_10_6_intel.whl'
