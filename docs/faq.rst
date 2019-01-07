.. _radiomics-faq-label:

==========================
Frequently Asked Questions
==========================

Installation
------------

During setup, python is unable to compile the C extensions.
###########################################################

This can occur when no compiler is available for python. If you're installing on Windows, you can find free compilers
for python `here <https://wiki.python.org/moin/WindowsCompilers>`_.

Error loading C extensions.
###########################

When I try to run PyRadiomics, I get ``Error loading C extensions``

When PyRadiomics is installed, the C extensions are compiled and copied to the installation folder, by default the
``site-packages`` folder. However, when the notebook is run form the repository, it is possible that PyRadiomics uses
the source code directly (i.e. runs in development mode). You can check this by checking the ``radiomics.__path__``
field, which will be something like ``['radiomics']`` when it is running in development mode and
``['path/to/python/Lib/site-packages']`` when running from the installed folder. If running in development mode, the C
extensions are not available by default. To make them available in development mode, run
``python setup.py build_ext --inplace`` from the commandline, which is similar to the ``install`` command, but just
compiles the C extensions end copies them to the source folder (making them available when running from the source tree).

Which python versions is PyRadiomics compatible with?
#####################################################

PyRadiomics is compatible with both python 2 and python 3. The automated testing uses python versions 2.7, 3.4, 3.5 and
3.6 (only 64 bits architecture). Python < 2.6 is not supported. Other python versions may be compatible with
PyRadiomics, but this is not actively tested and therefore not guaranteed to work.

Input / Customization
---------------------

I want to customize my extraction. How do I do that?
####################################################

See also :ref:`radiomics-customization-label`. PyRadiomics can be customized in various ways, but it's most easy to
do this by providing a :ref:`parameter file <radiomics-parameter-file-label>`. In this
`yaml structured <http://yaml.org/>`_ text file you can define your custom settings and which features and input image
types to enable.

.. _radiomics_fixed_bin_width:

What about gray value discretization? Fixed bin width? Fixed bin count?
#######################################################################

Currently, although many studies favour a fixed bin count over a fixed bin width, there is no hard evidence favouring
either a fixed bin width or a fixed bin count in all cases.
Therefore PyRadiomics implements both the option for setting a fixed bin count (``binCount``) and a fixed bin width
(``binWidth``, default).

The reason the a fixed bin width has been chosen as the default parameter is based in part on studies in PET that show
a better reproducibility of features when implementing a fixed bin width [1]_.
Furthermore, our reasoning is best illustrated by the following example:
Given an input with 2 images with 2 ROIs, with the range of gray values in the first being {0-100} and in the second
{0-10}. If you use a fixed bin count, the “meaning” of 1 (discretized) gray value difference is different (in the first
it means 10 gray values different, in the second just 1). This means you are looking at texture based on very different
contrasts.

This example does assume that the original gray values mean the same thing in both images, and in case of images with
definite/absolute gray values (e.g. HU in CT, SUV in PET imaging), this holds true. However, in case of
arbitrary/relative gray values (e.g. signal intensity in MR), this is not necessarily the case.
In this latter case, we still recommend a fixed bin width, but with additional pre-processing (e.g. normalization) to
ensure better comparability of gray values. Use of a fixed bin count would be possible here, but then the calculated
features may still be very influenced by the range of gray values seen in the image, as well as noise caused by the fact
that the original gray values are less comparable. Moreover, regardless of type of gray value discretization, steps must
be taken to ensure good comparability, as the first order features largely use the original gray values
(without discretization).

Finally, there is the issue of what value to use for the width of the bin. Again, there are currently no specific
guidelines from literature as to what constitutes an optimal bin width. We try to choose a bin width in such a way, that
the resulting amount of bins is somewhere between 30 and 130 bins, which shows good reproducibility and performance in
literature for a fixed bin count [2]_. This allows for differing ranges of intensity in
ROIs, while still keeping the texture features informative (and comparable inter lesion!).

Error loading parameter file
############################

When I try to load my own parameter file, I get error:"CoreError: Unable to load any data from source yaml file"

This error is thrown by PyKwalify when it is unable to read the parameter file. The most common cause is when the file
is indented using tabs, which throws a "token ('\t') not recognized error". Instead, ensure that all indentations are
made using 2 or 4 spaces.

What file types are supported by PyRadiomics for input image and mask?
######################################################################

PyRadiomics uses SimpleITK for image loading and handling. Therefore, all image types supported by SimpleITK can be
used as input for PyRadiomics. Please note that only one file location can be provided for image/mask. If you want to
provide the image in DICOM format, load the DICOM images using SimpleITK functionality and pass the resultant image
object instead.

If your input images are DICOM, you should first confirm the DICOM files you have correspond to a single image series. If you are not sure, you can sort the data such that you have a single directory per series using, for example, `dicomsort https://github.com/pieper/dicomsort`_. You can then convert the DICOM series into an ITK-readable volumetric format using `plastimatch convert http://plastimatch.org/plastimatch.html#plastimatch-convert`_ or `dcm2niix https://github.com/rordenlab/dcm2niix`_. We also provide a "labs" (experimental) script `pyradiomics-dcm https://github.com/Radiomics/pyradiomics/tree/master/labs/pyradiomics-dcm` that can do those conversions automatically and also saves the resulting features as DICOM SR.

.. _radiomics_geometry_mismatch:

Geometry mismatch between image and mask
########################################

My mask was generated using a another image than my input image, can I still extract features?

For various reasons, both image and mask must have the same geometry (i.e. same spacing, size, direction and origin)
when passed the feature classes. To this end PyRadiomics includes checks in the pipeline to ensure this is the case.
For more information on the mask checks, see :py:func:`~imageoperations.checkMask`. If the geometry error is small
difference in origin, spacing or direction, you can increase the tolerance by setting ``geometryTolerance``.
If the error is large, or the dimensions do not match, you could also resample the mask to image space. An example of
this is provided in ``bin\resampleMask.py`` and can be enabled in PyRadiomics by setting ``correctMask`` to ``True``,
which will only perform this correction in case of a geometery mismatch where the mask contains a valid ROI (i.e. the
mask contains the label value which does not include areas outside image physical bounds).

What modalities does PyRadiomics support?
#########################################

PyRadiomics is not developed for one specific modality. Multiple modalities can be processed by PyRadiomics, although
the optimal settings may differ between modalities. There are some constraints on the input however:

1. Gray scale volume: PyRadiomics currently does not provide extraction from color images or images with complex values
2. 3D or slice: Although PyRadiomics supports single slice (2D) feature extraction, the input is still required to have
   3 dimensions (where in case of 2D, a dimension may be of size 1).

If you want to use 2D, color and/or 4D volumes, additional preprocessing is required to convert the images.
See `this thread <https://groups.google.com/forum/#!topic/pyradiomics/QLdD_qEw3PY>`_ for some tips and tricks on how to achieve this.

Can I use DICOM-RT struct for the input mask?
#############################################

PyRadiomics does not support DICOM-RT struct as input directly. We recommend to convert these using for example `plastimatch convert http://plastimatch.org/plastimatch.html#plastimatch-convert`_. You can also
load DICOM RT in 3D Slicer after installing the
`SlicerRT <http://slicerrt.github.io/>`_ extension. DICOM RT loaded into 3D Slicer using SlicerRT extension
can then be passed as input to the `SlicerRadiomics extension
<https://github.com/Radiomics/SlicerRadiomics>`_.


Usage
-----

How should the input file for ``pyradiomics`` in batch-mode be structured?
##########################################################################

Currently, the batch input file for ``pyradiomics`` is a csv file specifying the combinations of images and masks for
which to extract features. It must contain a header line, where at least header "Image" and "Mask" should be specified
(capital sensitive). These identify the columns that contain the file location of the image and the mask, respectively.
Each subsequent line represents one combination of an image and a mask. Additional columns are also allowed, these are
copied to the output in the same order as the input, with the additional columns of the calculated features appended
at the end. *N.B. All header names should be unique and not match any of the produced header names by pyradiomics.*

Radiomics module not found in jupyter notebook
##############################################

I installed PyRadiomics, but when I run the jupyter notebook, I get ``ImportError: No module named radiomics``

This can have two possible causes:

1) When installing PyRadiomics from the repository, your python path variable will be updated to enable python to find
   the package. However, this value is only updated in commandline windows when they are restarted. If your jupyter
   notebook was running during installation, you first need to restart it.
2) Multiple versions of python can be installed on your machine simultaneously. Ensure PyRadiomics is installed on the
   same version you are using in your Jupyter notebook.

I'm missing features from my output. How can I see what went wrong?
###################################################################

If calculation of features or application of filters fails, a warning is logged. If you want to know exactly what
happens inside the toolbox, PyRadiomics provides extensive debug logging. You can enable this to be printed to the
out, or stored in a separate log file. The output is regulated by :py:func:`radiomics.setVerbosity` and the PyRadiomics
logger can be accessed via ``radiomics.logger``. See also :ref:`here <radiomics-logging-label>` and the examples
included in the repository on how to set up logging.

I'm unable to calculate texture matrices and getting a RunTimeError instead
###########################################################################

This error means that something went wrong during the calculation of the matrices in the C extensions.
There are several potential causes for this error:

- "Error parsing array arguments."

This error is thrown when either the Image or the Mask provided to the function could not be interpreted as a numpy array.

- "Expected a 3D array for image and mask."

Thrown when either the Image or Mask Provided did not have 3 dimensions (in case of a single slice calculation, the
input arrays should still have 3 dimensions, although one of them will then have a size of 1).

- "Dimensions of image and mask do not match."

This means that the size of the mask array does not match the size of the image array. Because numpy arrays do not
contain information on the transformation to the physical world, input arrays of differing sizes cannot be matched.
You can solve this error by resampling the SimplITK-Image object of the Mask to the geometry of the Image before
converting them to their respective numpy arrays for feature calculation. See also :ref:`radiomics_geometry_mismatch`.

- "Error parsing distances array."

This error is shown if the C extension was not able to interpret the distances argument that was provided. In the
settings, the ``distances`` parameter should be either a tuple or a list of values.

- "Expecting distances array to be 1-dimensional."

Again an error in the provided distances. The list provided should be 1 dimensional (i.e. no nested lists).

- "Error calculating angles."

This error means there was an issue in generating the angles based on the distances provided. Currently, this only
occurs when distances < 1 are provided.

- "Number of elements in <Matrix> would overflow index variable! (...)"

This error is shown when the size of the (flattened) output array would be larger than the maximum integer value
(~2 mln). This is generally caused by a too large number of bins after discretization, resulting in a too large range of
gray values in the discretized image used for texture calculation. We generally advise to chose a bin width so, that the
number of bins after discretization does not exceed 150-200. Running the code with DEBUG logging enabled shows the
number of bins that are generated and may help to give an indication as to how large your matrices are.

- "Failed to initialize output array for <Matrix>"

This means that the computer was unable to allocate memory for the output. This is most likely due to a too large output
size or too little free memory being available. Similar as above, run with DEBUG logging to see how many bins are
generated (giving an indication on how large the output matrices are).

- "Calculation of <Matrix> Failed."

This error means there was a problem in the calculation of the matrix itself. It is generally thrown if the code tries
to set an element in the output array that is out-of-range. This can happen if there are voxels inside the ROI that
have gray values that are larger than the ``Ng`` parameter that is provided when calling the C function from Python.

I'm able to extract features, but many are NaN, 0 or 1. What happened?
######################################################################

It is possible that the segmentation was too small to extract a valid texture. Check the value of ``VoxelNum``, which is
part of the additional information in the output. This is the number of voxels in the ROI after pre processing and
therefore the number of voxels that are used for feature calculation.

Another problem can be that you have to many or too few gray values after discretization. You can check this by
comparing the range of gray values in the ROI (a First Order feature) with the value for your ``binWidth`` parameter.
More bins capture smaller differences in gray values, but too many bins (compared to number of voxels) will yield low
probabilities in the texture matrices, resulting in non-informative features. There is no definitive answer for the
ideal number of discretized gray values, and this may differ between modalities.
One study [2]_ assessed the number of bins in PET and found that in the range of 16 - 128 bins, texture features did not
differ significantly.

Does PyRadiomics support voxel-wise feature extraction?
#######################################################

Yes, as of version 2.0, voxelwise calculation has been implemented. However, as this entails the calculations of
features for each voxel, performing a voxelwise extraction is much slower and as the output consists of a feature map
for each feature, output size is also much larger. See more on enabling a voxel-based extraction in the
:ref:`usage section<radiomics-usage-label>`.

Miscellaneous
-------------

A new version of PyRadiomics is available! Where can I find out what changed?
#############################################################################

When a new version is released, a changelog is included in the
`release statement <https://github.com/Radiomics/pyradiomics/releases>`_. Between releases, changes are not explicitly
documented, but all significant changes are implemented using pull requests. Check the
`merged pull request <https://github.com/Radiomics/pyradiomics/pulls?utf8=%E2%9C%93&q=is%3Apr%20is%3Amerged>`_ for the
latest changes.

I have some ideas for PyRadiomics. How can I contribute?
########################################################

We welcome suggestions and contributions to PyRadiomics. Check our
`guidelines <https://github.com/Radiomics/pyradiomics/blob/master/CONTRIBUTING.md>`_ to see how you can contribute to
PyRadiomics. Signatures and code styles used in PyRadiomics are documented in the :ref:`radiomics-developers` section.

I found a bug! Where do I report it?
####################################

We strive to keep PyRadiomics as bug free as possible by thoroughly testing new additions before including them in the
stable version. However, nothing is perfect, and some bugs may therefore exist. Report yours by
`opening an issue <https://github.com/Radiomics/pyradiomics/issues>`_ on the GitHub or contact us at the
`pyradiomics email list <https://groups.google.com/forum/#!forum/pyradiomics>`_. If you want to help in fixing it, we'd
welcome you to open up a `pull request <https://github.com/Radiomics/pyradiomics/pulls>`_ with your suggested fix.

My question is not listed here...
#################################

If you have a question that is not listed here, check the
`pyradiomics email list <https://groups.google.com/forum/#!forum/pyradiomics>`_ or the
`issues on GitHub <https://github.com/Radiomics/pyradiomics/issues>`_. Feel free to post a new question or issue and
we'll try to get back to you ASAP.

.. [1] Leijenaar RTH, Nalbantov G, Carvalho S, van Elmpt WJC, Troost EGC, Boellaard R, et al. ; The effect of SUV
        discretization in quantitative FDG-PET Radiomics: the need for standardized methodology in tumor texture
        analysis; Sci Rep. 2015;5(August):11075
.. [2] Tixier F, Cheze-Le Rest C, Hatt M, Albarghach NM, Pradier O, Metges J-P, et al. *Intratumor
        Heterogeneity Characterized by Textural Features on Baseline 18F-FDG PET Images Predicts Response to Concomitant
        Radiochemotherapy in Esophageal Cancer.* J Nucl Med. 2011;52:369–78.
