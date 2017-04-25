============
Installation
============

------------
Get the code
------------

* Ensure you have the version control system ``git`` installed on your machine.

* Ensure that you have ``python`` installed on your machine, at least version 2.7 or 3.4.

* Clone the repository:

  * ``git clone git://github.com/Radiomics/pyradiomics``

.. _installation-label:

---------------------------
Installation on your system
---------------------------

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
