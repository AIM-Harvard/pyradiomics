.. pyradiomics documentation master file, created by
   sphinx-quickstart on Tue Jun 28 13:27:28 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pyradiomics documentation!
=====================================

This is an open-source python package for the extraction of Radiomics features from 2D and 3D images and
segmentations.

Image loading and preprocessing (e.g. resampling and cropping) are first done using ``SimpleITK``.
Then, loaded data are converted into numpy arrays for further calculation using feature classes
outlined below.

With this package we aim to establish a reference standard for Radiomic Analysis, and provide a tested and mantained
open-source platform for easy and reproducible Radiomic Feature extraction.

By doing so, we hope to increase awareness of radiomic capabilities and expand the community.

**If you publish any work which uses this package, please cite the following publication:**
*van Griethuysen et al, “Computational Radiomics System to Decode the Radiographic Phenotype”; Submitted*

Table of Contents
-----------------

.. toctree::
   :hidden:

   Home <self>

.. toctree::
   :maxdepth: 2

   installation
   usage
   radiomics
   features
   developers
   FAQs <faq>

Feature Classes
---------------

Currently supports the following feature classes:

* :ref:`First Order Statistics <radiomics-firstorder-label>` (19 features)
* :ref:`Shape-based <radiomics-shape-label>` (13 features)
* :ref:`Gray Level Cooccurence Matrix <radiomics-glcm-label>` (GLCM) (28 features)
* :ref:`Gray Level Run Length Matrix <radiomics-glrlm-label>` (GLRLM) (16 features)
* :ref:`Gray Level Size Zone Matrix <radiomics-glszm-label>` (GLSZM) (16 features)

On average, Pyradiomics extracts 961 features per image, which consist of the 13 shape descriptors and features
extracted from original and derived images (LoG with 3 sigma levels, 1 level of Wavelet decomposistions yielding 8
derived images).

Detailed description on feature classes and individual features is provided in section :ref:`radiomics-features-label`.

Filter Classes
--------------

Aside from the feature classes, there are also some built-in optional filters:

* :py:func:`Laplacian of Gaussian <radiomics.imageoperations.getLoGImage>` (LoG, based on SimpleITK functionality)
* :py:func:`Wavelet <radiomics.imageoperations.getWaveletImage>` (using the PyWavelets package)
* :py:func:`Square <radiomics.imageoperations.getSquareImage>`
* :py:func:`Square Root <radiomics.imageoperations.getSquareRootImage>`
* :py:func:`Logarithm <radiomics.imageoperations.getLogarithmImage>`
* :py:func:`Exponential <radiomics.imageoperations.getExponentialImage>`

For more information, see also :ref:`radiomics-imageoperations-label`.

Supporting reproducible extraction
----------------------------------

Aside from calculating features, the pyradiomics package includes additional information in the
output. This information contains information on used image and mask, as well as applied settings
and filters, thereby enabling fully reproducible feature extraction. For more information, see
:ref:`radiomics-generalinfo-label`.

3rd-party packages used in pyradiomics
--------------------------------------

* SimpleITK
* numpy
* PyWavelets (Wavelet filter)
* pykwalify (Enabling yaml parameters file checking)
* tqdm (Progressbar)
* six (Python 3 Compatibility)
* sphinx (Generating documentation)
* sphinx_rtd_theme (Template for documentation)
* nose-parameterized (Testing)

See also the `requirements file <https://github.com/Radiomics/pyradiomics/blob/master/requirements.txt>`_.

Installation
------------

PyRadiomics is OS independent and compatible with both Python 2.7 and Python >=3.4.

* Clone the repository

  * ``git clone git://github.com/Radiomics/pyradiomics``

* Install on your system (Linux, Mac OSX), with prerequisites:

  * ``cd pyradiomics``
  * ``sudo python -m pip install -r requirements.txt``
  * ``sudo python setup.py install``

* For more detailed installation instructions and installation on Windows see
  :ref:`Installation Details<installation-label>`

Pyradiomics Indices and Tables
------------------------------

* :ref:`modindex`
* :ref:`genindex`
* :ref:`search`

License
-------

This package is covered by the open source `3D Slicer License <https://github.com/Radiomics/pyradiomics/blob/master/LICENSE.txt>`_.

Developers
----------

 - `Joost van Griethuysen <https://github.com/JoostJM>`_:sup:`1,3,4`
 - `Andriy Fedorov <https://github.com/fedorov>`_:sup:`2`
 - `Nicole Aucoin <https://github.com/naucoin>`_:sup:`2`
 - `Jean-Christophe Fillion-Robin <https://github.com/jcfr>`_:sup:`5`
 - `Ahmed Hosny <https://github.com/ahmedhosny>`_:sup:`1`
 - `Steve Pieper <https://github.com/pieper>`_:sup:`6`
 - `Hugo Aerts (PI) <https://github.com/hugoaerts>`_:sup:`1,2`
 
:sup:`1`\ Department of Radiation Oncology, Dana-Farber Cancer Institute, Brigham and Women's Hospital, Harvard Medical School, Boston, MA,
:sup:`2`\ Department of Radiology, Brigham and Women's Hospital, Harvard Medical School, Boston, MA
:sup:`3`\ Department of Radiology, Netherlands Cancer Institute, Amsterdam, The Netherlands,
:sup:`4`\ GROW-School for Oncology and Developmental Biology, Maastricht University Medical Center, Maastricht, The Netherlands,
:sup:`5`\ Kitware,
:sup:`6`\ Isomics

Contact
-------
We are happy to help you with any questions. Please contact us on the `pyradiomics email list <https://groups.google.com/forum/#!forum/pyradiomics>`_.

We'd welcome your contributions to PyRadiomics. Please read the
`contributing guidelines <https://github.com/Radiomics/pyradiomics/blob/master/CONTRIBUTING.md>`_ on how to contribute
to PyRadiomics. Information on adding / customizing feature classes and filters can be found in the :ref:`developers`
section.

.. note::

   This work was supported in part by the US National Cancer Institute grant
   5U24CA194354, QUANTITATIVE RADIOMICS SYSTEM DECODING THE TUMOR PHENOTYPE.
