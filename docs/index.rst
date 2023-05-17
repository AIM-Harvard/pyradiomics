.. pyradiomics documentation master file, created by
   sphinx-quickstart on Tue Jun 28 13:27:28 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pyradiomics documentation!
=====================================

This is an open-source python package for the extraction of Radiomics features from medical imaging.
With this package we aim to establish a reference standard for Radiomic Analysis, and provide a tested and maintained
open-source platform for easy and reproducible Radiomic Feature extraction. By doing so, we hope to increase awareness
of radiomic capabilities and expand the community. The platform supports both the feature extraction in 2D and 3D and
can be used to calculate single values per feature for a region of interest ("segment-based") or to generate feature
maps ("voxel-based").

**If you publish any work which uses this package, please cite the following publication:**
*van Griethuysen, J. J. M., Fedorov, A., Parmar, C., Hosny, A., Aucoin, N., Narayan, V., Beets-Tan, R. G. H.,
Fillon-Robin, J. C., Pieper, S.,  Aerts, H. J. W. L. (2017). Computational Radiomics System to Decode the Radiographic
Phenotype. Cancer Research, 77(21), e104–e107. `https://doi.org/10.1158/0008-5472.CAN-17-0339 <https://doi.org/10.1158/0008-5472.CAN-17-0339>`_*

.. note::

   This work was supported in part by the US National Cancer Institute grant
   5U24CA194354, QUANTITATIVE RADIOMICS SYSTEM DECODING THE TUMOR PHENOTYPE.

.. warning::

   Not intended for clinical use.

Join the Community!
-------------------
Join the PyRadiomics community on google groups `here <https://groups.google.com/forum/#!forum/pyradiomics>`_.


Table of Contents
-----------------

.. toctree::
   :hidden:

   Home <self>

.. toctree::
   :maxdepth: 2

   installation
   usage
   customization
   radiomics
   features
   removedfeatures
   contributing
   developers
   labs
   FAQs <faq>
   changes

Feature Classes
---------------

Currently supports the following feature classes:

* :py:class:`First Order Statistics <radiomics.firstorder.RadiomicsFirstOrder>`
* :py:class:`Shape-based (3D) <radiomics.shape.RadiomicsShape>`
* :py:class:`Shape-based (2D) <radiomics.shape2D.RadiomicsShape2D>`
* :py:class:`Gray Level Co-occurrence Matrix <radiomics.glcm.RadiomicsGLCM>` (GLCM)
* :py:class:`Gray Level Run Length Matrix <radiomics.glrlm.RadiomicsGLRLM>` (GLRLM)
* :py:class:`Gray Level Size Zone Matrix <radiomics.glszm.RadiomicsGLSZM>` (GLSZM)
* :py:class:`Neigbouring Gray Tone Difference Matrix <radiomics.ngtdm.RadiomicsNGTDM>` (NGTDM)
* :py:class:`Gray Level Dependence Matrix <radiomics.gldm.RadiomicsGLDM>` (GLDM)

On average, Pyradiomics extracts :math:`\approx 1500` features per image, which consist of the 16 shape descriptors and
features extracted from original and derived images (LoG with 5 sigma levels, 1 level of Wavelet decomposistions
yielding 8 derived images and images derived using Square, Square Root, Logarithm and Exponential filters).

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
* :py:func:`Gradient <radiomics.imageoperations.getGradientImage>`
* :py:func:`Local Binary Pattern (2D) <radiomics.imageoperations.getLBP2DImage>`
* :py:func:`Local Binary Pattern (3D) <radiomics.imageoperations.getLBP3DImage>`

For more information, see also :ref:`radiomics-imageoperations-label`.

Supporting reproducible extraction
----------------------------------

Aside from calculating features, the pyradiomics package includes additional information in the
output. This information contains information on used image and mask, as well as applied settings
and filters, thereby enabling fully reproducible feature extraction. For more information, see
:ref:`radiomics-generalinfo-label`.

3rd-party packages used in pyradiomics
--------------------------------------

* SimpleITK (Image loading and preprocessing)
* numpy (Feature calculation)
* PyWavelets (Wavelet filter)
* pykwalify (Enabling yaml parameters file checking)
* six (Python 3 Compatibility)

See also the `requirements file <https://github.com/Radiomics/pyradiomics/blob/master/requirements.txt>`_.

Installation
------------

PyRadiomics is OS independent and compatible with and Python >=3.5. Pre-built binaries are available on
PyPi and Conda. To install PyRadiomics, ensure you have python
installed and run:

*  ``python -m pip install pyradiomics``

For more detailed installation instructions and building from source,
see :ref:`radiomics-installation-label` section.

Pyradiomics Indices and Tables
------------------------------

* :ref:`modindex`
* :ref:`genindex`
* :ref:`search`

License
-------

This package is covered by the open source `3-clause BSD License <https://github.com/Radiomics/pyradiomics/blob/master/LICENSE.txt>`_.

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
We are happy to help you with any questions. Please contact us on the `Radiomics community section of the 3D Slicer Discourse <https://discourse.slicer.org/c/community/radiomics/23>`_.

We'd welcome your contributions to PyRadiomics. Please read the
:ref:`contributing guidelines <radiomics-contributing-label>` on how to contribute to PyRadiomics. Information on
adding / customizing feature classes and filters can be found in the :ref:`radiomics-developers` section.
