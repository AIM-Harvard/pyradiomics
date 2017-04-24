==========================
Frequently Asked Questions
==========================

Installation
------------

**During setup, python is unable to compile the C extensions.**

This can occur when no compiler is available for python. If you're installing on Windows, you can find free compilers
for python `here <https://wiki.python.org/moin/WindowsCompilers>`_.

**I installed PyRadiomics successfully from the repository, but when I run the notebook, I get** ``Error loading C
extensions, switching to python calculation``

When PyRadiomics is installed, the C extensions are compiled and copied to the installation folder, by default the
``site-packages`` folder. However, when the notebook is run form the repository, it is possible that PyRadiomics uses
the source code directly (i.e. runs in development mode). You can check this by checking the ``radiomics.__path__``
field, which will be something like ``['radiomics']`` when it is running in development mode and
``['path/to/python/Lib/site-packages']`` when running from the installed folder. If running in development mode, the C
extensions are not available by default. To make them available in development mode, run ``python setup.py develop``
from the commandline, which is similar to the ``install`` command, but installs pyradiomics to the source folder
instead (i.e. does nothing to the python files, but compiles the C extensions and copies them to the source folder).

**Which python versions is PyRadiomics compatible with?**

PyRadiomics is compatible with both python 2 and python 3. The automated testing uses python versions 2.7, 3.4 and 3.5
(both 32 and 64 bits). Python < 2.6 is not supported. Other python versions may be compatible with PyRadiomics, but this
is not actively tested and therefore not guaranteed to work.

Input / Customization
---------------------

**I want to customize my extraction. How do I do that?**

See also the :ref:`radiomics-usage-label` section. PyRadiomics can be customized in various ways, but it's most easy to
do this by providing a parameter file. In this `yaml structured <http://yaml.org/>`_ text file you can define your
custom settings and which features and input image types to enable. Possible general settings, as well as possible image
types are documented :ref:`here <radiomics-featureextractor-label>`, with additional feature class specific settings
documented at the start of the :ref:`feature class <radiomics-features-label>`.

**When I try to load my own parameter file, I get error:"CoreError: Unable to load any data from source yaml file"**

This error is thrown by PyKwalify when it is unable to read the parameter file. The most common cause is when the file
is indented using tabs, which throws a "token ('\t') not recognized error". Instead, ensure that all indentations are
made using 2 or 4 spaces.

**What file types are supported by PyRadiomics for input image and mask?**

PyRadiomics uses SimpleITK for image loading and handling. Therefore, all image types supported by SimpleITK can be
used as input for PyRadiomics. Please note that only one file location can be provided for image/mask. If you want to
provide the image in DICOM format, load the DICOM images using SimpleITK functionality and pass the resultant image
object instead.

**My mask was generated using a another image than my input image, can I still extract features?**

For various reasons, both image and mask must have the same geometry (i.e. same spacing, size, direction and origin)
when passed the feature classes. To this end PyRadiomics includes checks in the pipeline to ensure this is the case.
If resampling is enabled, the mask is allowed to be of a different geometry, so long as the bounding box of the Region
of Interest (ROI) represents a physical space that is within the image. If resampling is disabled, mask and image are
required to be of the same geometry. For more information on resampling, see :py:func:`~imageoperations.resampleImage`.
For more information on the mask checks, see :py:func:`~imageoperations.checkMask`. Alternatively, you could also
resample the mask to image space. An example of this is provided in ``bin\resampleMask.py``.

**What modalities does PyRadiomics support?**

PyRadiomics is not developed for one specific modality. Multiple modalities can be processed by PyRadiomics, although
the optimal settings may differ between modalities. There are some constraints on the input however:

1. Gray scale volume: PyRadiomics currently does not provide extraction from color images or images with complex values
2. 3D or slice: Although PyRadiomics supports single slice (2D) feature extraction, the input is still required to have
   3 dimensions (where in case of 2D, a dimension may be of size 1).

**Can I use DICOM-RT struct for the input mask?**

PyRadiomics does not support DICOM-RT struct as input directly. We recommend to convert these using for example
`SlicerRT <http://slicerrt.github.io/>`_. We are working on providing support for DICOM-RT in the `Slicer extension
<https://github.com/Radiomics/SlicerRadiomics>`_, but this is not thoroughly tested yet.


Usage
-----

**How should the input file for** ``pyradiomicsbatch`` **be structured?**

Currently, the input file for ``pyradiomicsbatch`` is a csv file specifying the combinations of images and masks for
which to extract features. It must contain a header line, where at least header "Image" and "Mask" should be specified
(capital sensitive). These identify the columns that contain the file location of the image and the mask, respectively.
Each subsequent line represents one combination of an image and a mask. Additional columns are also allowed, these are
copied to the output in the same order as the input, with the additional columns of the calculated features appended
at the end. *N.B. All header names should be unique and not match any of the produced header names by pyradiomics.*

**I installed PyRadiomics, but when I run the jupyter notebook, I get** ``ImportError: No module named radiomics``

This can have two possible causes: 1) When installing PyRadiomics from the repository, your python path variable will be
updated to enable python to find the package. However, this value is only updated in commandline windows when they are
restarted. If your jupyter notebook was running during installation, you first need to restart it. 2) Multiple versions
of python can be installed on your machine simultaneously. Ensure PyRadiomics is installed on the same version you are
using in your Jupyter notebook.

**When I try to extract features, I get an error stating that image and mask do not occupy the same space.**

During extraction, SimpleITK checks whether the mask matches the image dimensions and physical space and raises an
error when this is not the case. This is to prevent attempts to extract features from an image using a mask that does
not match the image. If it is the correct mask, you can change the tolerance for differences in direction and origin
between the image and mask by setting ``SimpleITK.ProcessObject.SetGlobalDefaultDirectionTolerance`` and
``SimpleITK.ProcessObject.SetGlobalDefaultDirectionTolerance``, respectively. By default, these are set to 1e-6 mm.
Alternatively, you can resample your mask with SimpleITK.ResampleImageFilter(), using the image as a reference image and
Nearest Neighbour as interpolator. See the SimpleITK `documentation <https://itk.org/SimpleITKDoxygen/html/>`_ for more
information.

**I'm missing features from my output. How can I see what went wrong?**

If calculation of features or application of filters fails, a warning is logged. If you want to know exactly what
happens inside the toolbox, PyRadiomics provides extensive debug logging. You can enable this to be printed to the
out, or stored in a separate log file. The output is regulated by :py:func:`radiomics.setVerbosity` and the PyRadiomics
logger can be accessed via ``radiomics.logger``. See the examples included in the repository on how to set up logging.

**I'm able to extract features, but many are NaN, 0 or 1. What happend?**

It is possible that the segmentation was too small to extract a valid texture. Check the value of ``VoxelNum``, which is
part of the additional information in the output. This is the number of voxels in the ROI after pre processing and
therefore the number of voxels that are used for feature calculation.

Another problem can be that you have to many or too few gray values after discretization. You can check this by
comparing the range of gray values in the ROI (a First Order feature) with the value for your ``binWidth`` parameter.
More bins capture smaller differences in gray values, but too many bins (compared to number of voxels) will yield low
probabilities in the texture matrices, resulting in non-informative features. There is no definitive answer for the
ideal number of discretized gray values, and this may differ between modalities.
One study [1]_ assessed the number of bins in PET and found that in the range of 16 - 128 bins, texture features did not differ
significantly.

**Does PyRadiomics support voxel-wise feature extraction (for the generation of colormaps)?**

No, currently PyRadiomics only supports lesion-based feature extraction. However, voxel-based feature extraction may be
a good addition in the future. If you have thoughts or ideas on how to implement this, we'd welcome your input on the
`pyradiomics email list <https://groups.google.com/forum/#!forum/pyradiomics>`_.

Miscellaneous
-------------

**A new version of PyRadiomics is available! Where can I find out what changed?**

When a new version is released, a changelog is included in the
`release statement <https://github.com/Radiomics/pyradiomics/releases>`_. Between releases, changes are not explicitly
documented, but all significant changes are implemented using pull requests. Check the
`merged pull request <https://github.com/Radiomics/pyradiomics/pulls?utf8=%E2%9C%93&q=is%3Apr%20is%3Amerged>`_ for the
latest changes.

**I have some ideas for PyRadiomics. How can I contribute?**

We welcome suggestions and contributions to PyRadiomics. Check our
`guidelines <https://github.com/Radiomics/pyradiomics/blob/master/CONTRIBUTING.md>`_ to see how you can contribute to
PyRadiomics. Signatures and code styles used in PyRadiomics are documented in the :ref:`developers` section.

**I found a bug! Where do I report it?**

We strive to keep PyRadiomics as bug free as possible by thoroughly testing new additions before including them in the
stable version. However, nothing is perfect, and some bugs may therefore exist. Report yours by
`opening an issue <https://github.com/Radiomics/pyradiomics/issues>`_ on the GitHub or contact us at the
`pyradiomics email list <https://groups.google.com/forum/#!forum/pyradiomics>`_. If you want to help in fixing it, we'd
welcome you to open up a `pull request <https://github.com/Radiomics/pyradiomics/pulls>`_ with your suggested fix.

**My question is not listed here...**

If you have a question that is not listed here, check the
`pyradiomics email list <https://groups.google.com/forum/#!forum/pyradiomics>`_ or the
`issues on GitHub <https://github.com/Radiomics/pyradiomics/issues>`_. Feel free to post a new question or issue and
we'll try to get back to you ASAP.

.. [1] Tixier F, Cheze-Le Rest C, Hatt M, Albarghach NM, Pradier O, Metges J-P, et al. *Intratumor
        Heterogeneity Characterized by Textural Features on Baseline 18F-FDG PET Images Predicts Response to Concomitant
        Radiochemotherapy in Esophageal Cancer.* J Nucl Med. 2011;52:369–78.
