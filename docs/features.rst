.. _radiomics-features-label:

=================
Radiomic Features
=================

This section contains the definitions of the various features that can be extracted using PyRadiomics. They are
subdivided into the following classes:

* :py:class:`First Order Statistics <radiomics.firstorder.RadiomicsFirstOrder>` (19 features)
* :py:class:`Shape-based (3D) <radiomics.shape.RadiomicsShape>` (16 features)
* :py:class:`Shape-based (2D) <radiomics.shape2D.RadiomicsShape2D>` (10 features)
* :py:class:`Gray Level Co-occurrence Matrix <radiomics.glcm.RadiomicsGLCM>` (24 features)
* :py:class:`Gray Level Run Length Matrix <radiomics.glrlm.RadiomicsGLRLM>` (16 features)
* :py:class:`Gray Level Size Zone Matrix <radiomics.glszm.RadiomicsGLSZM>` (16 features)
* :py:class:`Neighbouring Gray Tone Difference Matrix <radiomics.ngtdm.RadiomicsNGTDM>` (5 features)
* :py:class:`Gray Level Dependence Matrix <radiomics.gldm.RadiomicsGLDM>` (14 features)

All feature classes, with the exception of shape can be calculated on either the original image and/or a derived image,
obtained by applying one of several filters. The shape descriptors are independent of gray value, and are extracted from
the label mask. If enabled, they are calculated separately of enabled input image types, and listed in the result as if
calculated on the original image.

Most features defined below are in compliance with feature definitions as described by the Imaging Biomarker
Standardization Initiative (IBSI), which are available in a separate document by Zwanenburg et al. (2016) [1]_.
Where features differ, a note has been added specifying the difference.

.. _radiomics-firstorder-label:

First Order Features
--------------------

.. automodule:: radiomics.firstorder
    :members:
    :undoc-members:
    :show-inheritance:
    :member-order: bysource

.. _radiomics-shape-label:

Shape Features (3D)
-------------------

.. automodule:: radiomics.shape
    :members:
    :undoc-members:
    :show-inheritance:
    :member-order: bysource

.. _radiomics-shape2D-label:

Shape Features (2D)
-------------------

.. automodule:: radiomics.shape2D
    :members:
    :undoc-members:
    :show-inheritance:
    :member-order: bysource

.. _radiomics-glcm-label:

Gray Level Co-occurrence Matrix (GLCM) Features
-----------------------------------------------

.. automodule:: radiomics.glcm
    :members:
    :undoc-members:
    :show-inheritance:
    :member-order: bysource

.. _radiomics-glszm-label:

Gray Level Size Zone Matrix (GLSZM) Features
--------------------------------------------

.. automodule:: radiomics.glszm
    :members:
    :undoc-members:
    :show-inheritance:
    :member-order: bysource

.. _radiomics-glrlm-label:

Gray Level Run Length Matrix (GLRLM) Features
---------------------------------------------

.. automodule:: radiomics.glrlm
    :members:
    :undoc-members:
    :show-inheritance:
    :member-order: bysource

.. _radiomics-ngtdm-label:

Neighbouring Gray Tone Difference Matrix (NGTDM) Features
---------------------------------------------------------

.. automodule:: radiomics.ngtdm
    :members:
    :undoc-members:
    :show-inheritance:
    :member-order: bysource

.. _radiomics-gldm-label:

Gray Level Dependence Matrix (GLDM) Features
--------------------------------------------

.. automodule:: radiomics.gldm
    :members:
    :undoc-members:
    :show-inheritance:
    :member-order: bysource

.. [1] Zwanenburg, A., Leger, S., Vallières, M., and Löck, S. (2016). Image biomarker
    standardisation initiative - feature definitions. In eprint arXiv:1612.07003 [cs.CV]
