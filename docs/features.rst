.. _radiomics-features-label:

=================
Radiomic Features
=================

This section contains the definitions of the various features that can be extracted using PyRadiomics. They are
subdivided into the following classes:

* :ref:`radiomics-firstorder-label` (19 features)
* :ref:`radiomics-shape-label` (13 features)
* :ref:`radiomics-glcm-label` (28 features)
* :ref:`radiomics-glszm-label` (16 features)
* :ref:`radiomics-glrlm-label` (16 features)

All feature classes, with the exception of shape can be calculated on either the original image and/or a derived image,
obtained by applying one of several filters. The shape descriptors are independent of gray value, and are extracted from
the label mask. If enabled, they are calculated separately of enabled input image types, and listed in the result as if
calculated on the original image.

.. _radiomics-firstorder-label:

First Order Features
--------------------

.. automodule:: radiomics.firstorder
    :members:
    :undoc-members:
    :show-inheritance:
    :member-order: bysource

.. _radiomics-shape-label:

Shape Features
--------------

.. automodule:: radiomics.shape
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

