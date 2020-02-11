.. _radiomics-faq-label:

==========================
Frequently Asked Questions
==========================

Feature Extraction: Input, Customization and Reproducibility
------------------------------------------------------------

Does PyRadiomics adhere to IBSI definitions of features?
########################################################

For the most part, yes.

PyRadiomics development is also involved in the standardisation effort by the IBSI team.
Still, there are some differences between PyRadiomics and feature extraction as defined in the IBSI documents.
These arise in places where several equally valid alternatives exist. In some of these cases, PyRadiomics opted for the
one alternative, while the IBSI standard recommends another. For the sake of consistency in PyRadiomics development, we
have opted not to change the PyRadiomics implementation, but rather document the difference.

Most notably are the differences in gray value discretization (just for the fixed bin size type) and resampling. These
differences cannot be corrected by customization settings alone and require replacement by custom functions:

- **Binning:** When performing gray value discretization using a fixed bin width (AKA IBSI: FBS, Fixed Bins Size), and
  with Resegmentation set (most notable Case A and C), IBSI computes the bin edges as equally spaced from the minimum of
  the resegmentation range if absolute-resegmentation, and min(intensities) when sigma-resegmentation.

  In PyRadiomics, gray value discretization using a fixed bin width always utilizes bin edges which are equally spaced
  from 0 and where the lowest gray level is ensured to be discretized to the first bin. Regardless of any
  resegmentation, etc.
- **Resampling:**

  - Alignment of the grid: In IBSI, resampling grid is aligned by aligning the center of the image,
    whereas in Pyradiomics, we align the corner of the origin voxel. This can result in slightly different interpolated
    results, and even slightly different sizes of the resampled image and ROI, which in turn causes differences in
    extracted feature values.
  - **Gray value rounding:** In IBSI, they reason that if the original intensity values are from some lower precision
    datatype, resampled values (which are floating point numbers, usually 64-bit) should be resampled to a similar
    resolution. In the case of the IBSI phantoms, resampling to the nearest integer. PyRadiomics does not implement
    this, as differences are likely to be small and therefore serve more to add complexity than increase the meaning
    of the extracted values. Especially considering gray values are discretized anyway prior to calculation of most
    (all except firstorder) features. If some kind of normalization is performed, then meaning of the gray values
    changes as well. Differences arise here, because small rounding differences can cause a voxel to be assigned to a
    different bin, which can be a marked change in feature value results, especially in small ROIs.
  - Mask resampling: In IBSI, different interpolators can also be selected for resampling of the mask, with an
    additional thresholding the retrieve the binary mask. This only works if the mask is limited to zero and non-zero
    (i.e. 1) values. PyRadiomics also supports masks with different value labels, allowing extraction from different
    ROIs in the same mask file by indicating different label values. To prevent any incorrect re-assignment,
    PyRadiomics forces the mask resampling to be nearest neighbor.

Next, there are also some differences which can be solved by custom settings, in this case as only applies to
Configuration E, where both absolute AND sigma resegmentation are performed. In PyRadiomics, both types are implemented,
but only 1 can be selected at a time. To simulate applying both types, I calculated the absolute range after
resegmentation and used that as the absolute resegment range: [ -718, 400 ]

Finally, there is a difference between PyRadiomics and IBSI in the calculation of firstorder: Kurtosis. IBSI calculates
Excess Kurtosis, which is equal to Kurtosis - 3. PyRadiomics calculates Kurtosis, which is always +3 compared to IBSI.
The difference of 3 stems from the fact that a gaussian distribution has a kurtosis of 3.

So in summary, the cause of the difference between PyRadiomics results and IBSI benchmark, per case:

Configuration C: Due to differences in gray value discretization and resampling
Configuration D: Due to differences in resampling
Configuration E: Due to differences in resampling and resegmentation

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

What modalities does PyRadiomics support? 2D? 3D? Color?
########################################################

PyRadiomics is not developed for one specific modality. Multiple modalities can be processed by PyRadiomics, although
the optimal settings may differ between modalities. There are some constraints on the input however:

1. Gray scale volume: PyRadiomics currently does not provide extraction from color images or images with complex values.
   In those cases, each pixel/voxel has multiple values and PyRadiomics does not know how you'd want to combine those.
   It is possible to select a color channel and use that as input::

     import SimpleITK as sitk

     from radiomics import featureextractor

     # Instantiate extractor with parameter file
     extractor = featureextractor.RadiomicsFeatureExtractor(r'path/to/params.yml')

     # Set path to mask
     ma_path = 'path/to/maskfile.nrrd'
     label = 1  # Change this if the ROI in your mask is identified by a different value

     # Load the image and extract a color channel
     color_channel = 0
     im = sitk.ReadImage(r'path/to/image.jpg')
     selector = sitk.VectorIndexSelectionCastImageFilter()
     selector.SetIndex(color_channel)
     im = selector.Execute(im)

     # Run the extractor
     results = extractor.execute(im, ma_path, label=label)

2. File format: Currently, PyRadiomics requires the image and mask input to be either a string pointing to a single file
   containing the image/mask, or a SimpleITK.Image object (only possible in interactive mode). When e.g. using DICOMs,
   the separate files need to be first combined into 1 volume prior to extracting features by either converting to e.g.
   NRRD or NIFTII, or reading the DICOM in a python script and calling PyRadiomics from that script. See also
   :ref:`radiomics_input_formats`.
3. Dimensionality: PyRadiomics supports both 2D and 3D input images, but be aware that feature class ``shape`` only
   extracts 3D shape descriptors and ``shape2D`` only 2D shape descriptors. If you have a 3D volume, but a single-slice
   segmentation and want the results to include 2D shape descriptors, enable ``shape2D`` and set ``force2D=True``. This
   allows you to extract 2D shape features from a 3D volume with single-slice segmentation (but fails when segmentation
   represents a volume segmentation spanning multiple slices).

.. _radiomics_input_formats:

What file types are supported by PyRadiomics for input image and mask?
######################################################################

PyRadiomics uses SimpleITK for image loading and handling. Therefore,
`all image types supported by SimpleITK <https://simpleitk.readthedocs.io/en/master/Documentation/docs/source/IO.html>`_
can be used as input for PyRadiomics. Please note that only one file location can be provided for image/mask.

If your input images are DICOM, things become more complicated. If you want to process a single 2D image slice stored in
DICOM format, you can use it as any other format. However, if you are processing a volumetric dataset, you should first
confirm the DICOM files you have correspond to a single image series.
If you are not sure, you can sort the data such that you have a single directory per series using, for example,
`dicomsort <https://github.com/pieper/dicomsort>`_. You can then convert the DICOM series into an ITK-readable
volumetric format using `plastimatch convert <http://plastimatch.org/plastimatch.html#plastimatch-convert>`_ or
`dcm2niix <https://github.com/rordenlab/dcm2niix>`_.

If your label is defined in DICOM format, this can mean different things. First, check what is the modality of the
dataset with the label. You can check that by using `dcmdump <https://support.dcmtk.org/docs/dcmdump.html>`_, and then
checking the line that says "Modality". You can find the binary packages of this tool
`here <https://github.com/QIICR/atom-dicom-dump#install-dcmtk-andor-gdcm>`_ (you can also use
`dicom-dump package <https://github.com/QIICR/atom-dicom-dump>`_ if you want to look at DICOM files more conveniently
from the `Atom editor <https://atom.io>`_).

* If the modality is an image (e.g., CT or MR), use `plastimatch` or `dcm2niix` to convert the image into a 3D volume.
* If the modality is RT, use `plastimatch` to convert contours of the structure sets into 3D volumes.
* If the modality is SEG, use `dcmqi <https://github.com/QIICR/dcmqi>`_ to convert voxel segmentations into 3D volumes.

We also provide a "labs" (experimental) script
`pyradiomics-dcm <https://github.com/Radiomics/pyradiomics/tree/master/labs/pyradiomics-dcm>`_ that can do those
conversions automatically and also saves the resulting features as DICOM SR.

I want to customize my extraction. How do I do that?
####################################################

See also :ref:`radiomics-customization-label`. PyRadiomics can be customized in various ways, but it's most easy to
do this by providing a :ref:`parameter file <radiomics-parameter-file-label>`. In this
`yaml structured <http://yaml.org/>`_ text file you can define your custom settings and which features and input image
types to enable.
We strongly recommend to use this method for customization, as it contains all the customization in 1 human-readable
file, which can be shared for other users and, more importantly, is checked for validity prior to the extraction.

Does PyRadiomics support voxel-wise feature extraction?
#######################################################

Yes, as of version 2.0, voxelwise calculation has been implemented. However, as this entails the calculations of
features for each voxel, performing a voxelwise extraction is much slower and as the output consists of a feature map
for each feature, output size is also much larger. See more on enabling a voxel-based extraction in the
:ref:`usage section<radiomics-usage-label>`.

How should the input file for ``pyradiomics`` in batch-mode be structured?
##########################################################################

Currently, the batch input file for ``pyradiomics`` is a csv file specifying the combinations of images and masks for
which to extract features. It must contain a header line, where at least header "Image" and "Mask" should be specified
(capital sensitive). These identify the columns that contain the file location of the image and the mask, respectively.
Each subsequent line represents one combination of an image and a mask. Additional columns are also allowed, these are
copied to the output in the same order as the input, with the additional columns of the calculated features appended
at the end. *N.B. All header names should be unique and not match any of the produced header names by pyradiomics.*

Common Errors
-------------

Error loading parameter file
############################

When I try to load my own parameter file, I get error:"CoreError: Unable to load any data from source yaml file"

This error is thrown by PyKwalify when it is unable to read the parameter file. The most common cause is when the file
is indented using tabs, which throws a "token ('\t') not recognized error". Instead, ensure that all indentations are
made using 2 or 4 spaces.

.. _radiomics_geometry_mismatch:

Geometry mismatch between image and mask
########################################

My mask was generated using a another image than my input image, can I still extract features?

For various reasons, both image and mask must have the same geometry (i.e. same spacing, size, direction and origin)
when passed the feature classes. To this end PyRadiomics includes checks in the pipeline to ensure this is the case.
For more information on the mask checks, see :py:func:`~imageoperations.checkMask`. If the geometry error is due to a
small difference in origin, spacing or direction, you can increase the tolerance by setting ``geometryTolerance``.
If the error is large, or the dimensions do not match, you could also resample the mask to image reference space. An
example of this is provided in ``bin\resampleMask.py`` and can be enabled in PyRadiomics by setting ``correctMask`` to
``True``, which will only perform this correction in case of a geometery mismatch where the mask contains a valid ROI
(i.e. the mask contains the label value which does not include areas outside image physical bounds).

ValueError: ('Label (...) not present in mask. Choose from [...]'
#################################################################

This error indicates that the ROI identified by the value of ``label`` does not exist. I.e. there are no voxels in the
mask volume that have the value specified in ``label``. To help you fix the error, this error also lists the possible
alternative label values that have been found. If no values are listed at the end of this error, it means that there
are no segmented voxels in your mask.

.. note::

  It is possible that in the original mask there were voxels segmented, but were lost during resampling (when
  upsampling). In that case, the ROI is too small for the requested ``resampledPixelSpacing`` and should be treated as
  'nothing is segmented'.

I'm unable to calculate texture matrices and getting a RunTimeError instead
###########################################################################

This error means that something went wrong during the calculation of the matrices in the C extensions.
There are several potential causes for this error:

- "Error parsing array arguments."

This error is thrown when either the Image or the Mask provided to the function could not be interpreted as a numpy array.

- "Expected image and mask to have equal number of dimensions."

Thrown when the Image and Mask Provided did not have the same number of dimensions. N-Dimensional extraction is
supported, but the image and mask are still required to match in both the size and number of dimensions.

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

I'm missing features from my output. How can I see what went wrong?
###################################################################

If calculation of features or application of filters fails, a warning is logged. If you want to know exactly what
happens inside the toolbox, PyRadiomics provides extensive debug logging. You can enable this to be printed to the
out, or stored in a separate log file. The output is regulated by :py:func:`radiomics.setVerbosity` and the PyRadiomics
logger can be accessed via ``radiomics.logger``. See also :ref:`here <radiomics-logging-label>` and the examples
included in the repository on how to set up logging.

Radiomics module not found in jupyter notebook
##############################################

I installed PyRadiomics, but when I run the jupyter notebook, I get ``ImportError: No module named radiomics``

This can have two possible causes:

1) When installing PyRadiomics from the repository, your python path variable will be updated to enable python to find
   the package. However, this value is only updated in commandline windows when they are restarted. If your jupyter
   notebook was running during installation, you first need to restart it.
2) Multiple versions of python can be installed on your machine simultaneously. Ensure PyRadiomics is installed on the
   same version you are using in your Jupyter notebook.

Building PyRadiomics from source
--------------------------------

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

Miscellaneous
-------------

Which python versions is PyRadiomics compatible with?
#####################################################

PyRadiomics is compatible with python 3. Python 2 support was dropped in PyRadiomics version 3.0, though compatibility
code was retained. However, the automated testing only uses python versions 3.5, 3.6 and 3.7 (64 bits architecture).
Python < 2.6 is not supported. Other python versions may be compatible with PyRadiomics, but this is not actively tested
and therefore not guaranteed to work. Pre-built wheels are only available for the tested versions of python (3.5, 3.6
and 3.7)

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
`opening an issue <https://github.com/Radiomics/pyradiomics/issues>`_ on the GitHub. If you want to help in fixing it,
we'd welcome you to open up a `pull request <https://github.com/Radiomics/pyradiomics/pulls>`_ with your suggested fix.

My question is not listed here...
#################################

If you have a question that is not listed here, check the
`pyradiomics forum on 3D Slicer discourse <https://discourse.slicer.org/c/community/radiomics/23>`_ or the
`issues on GitHub <https://github.com/Radiomics/pyradiomics/issues>`_. Feel free to post a new question or issue and
we'll try to get back to you ASAP.

.. [1] Leijenaar RTH, Nalbantov G, Carvalho S, van Elmpt WJC, Troost EGC, Boellaard R, et al. ; The effect of SUV
        discretization in quantitative FDG-PET Radiomics: the need for standardized methodology in tumor texture
        analysis; Sci Rep. 2015;5(August):11075
.. [2] Tixier F, Cheze-Le Rest C, Hatt M, Albarghach NM, Pradier O, Metges J-P, et al. *Intratumor
        Heterogeneity Characterized by Textural Features on Baseline 18F-FDG PET Images Predicts Response to Concomitant
        Radiochemotherapy in Esophageal Cancer.* J Nucl Med. 2011;52:369–78.
