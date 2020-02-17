from __future__ import print_function

import logging

import numpy
import pywt
import SimpleITK as sitk
import six
from six.moves import range

logger = logging.getLogger(__name__)


def getMask(mask, **kwargs):
  """
  Function to get the correct mask. Includes enforcing a correct pixel data type (UInt32).

  Also supports extracting the mask for a segmentation (stored as SimpleITK Vector image) if necessary.
  In this case, the mask at index ``label_channel`` is extracted. The resulting 3D volume is then treated as it were a
  scalar input volume (i.e. with the region of interest defined by voxels with value matching ``label``).

  Finally, checks if the mask volume contains an ROI identified by ``label``. Raises a value error if the label is not
  present (including a list of valid labels found).

  :param mask: SimpleITK Image object representing the mask. Can be a vector image to allow for overlapping masks.
  :param kwargs: keyword arguments. If argument ``label_channel`` is present, this is used to select the channel.
    Otherwise label_channel ``0`` is assumed.
  :return: SimpleITK.Image with pixel type UInt32 representing the mask volume
  """
  global logger
  label = kwargs.get('label', 1)
  label_channel = kwargs.get('label_channel', 0)
  if 'vector' in mask.GetPixelIDTypeAsString().lower():
    logger.debug('Mask appears to be a segmentation object (=stored as vector image).')
    n_components = mask.GetNumberOfComponentsPerPixel()
    assert label_channel < n_components, \
        "Mask %i requested, but segmentation object only contains %i objects" % (label_channel, n_components)

    logger.info('Extracting mask at index %i', label_channel)
    selector = sitk.VectorIndexSelectionCastImageFilter()
    selector.SetIndex(label_channel)
    mask = selector.Execute(mask)

  logger.debug('Force casting mask to UInt32 to ensure correct datatype.')
  mask = sitk.Cast(mask, sitk.sitkUInt32)

  labels = numpy.unique(sitk.GetArrayFromImage(mask))
  if len(labels) == 1:
    raise ValueError('No labels found in this mask (i.e. nothing is segmented)!')
  if label not in labels:
    raise ValueError('Label (%g) not present in mask. Choose from %s' % (label, labels[labels != 0]))

  return mask


def getBinEdges(parameterValues, **kwargs):
  r"""
  Calculate and return the histogram using parameterValues (1D array of all segmented voxels in the image).

  **Fixed bin width:**

  Returns the bin edges, a list of the edges of the calculated bins, length is N(bins) + 1. Bins are defined such, that
  the bin edges are equally spaced from zero, and that the leftmost edge :math:`\leq \min(X_{gl})`. These bin edges
  represent the half-open ranges of each bin :math:`[\text{lower_edge}, \text{upper_edge})` and result in gray value
  discretization as follows:

  .. math::
    X_{b, i} = \lfloor \frac{X_{gl, i}}{W} \rfloor - \lfloor \frac {\min(X_{gl})}{W} \rfloor + 1

  Here, :math:`X_{gl, i}` and :math:`X_{b, i}` are gray level intensities before and after discretization, respectively.
  :math:`{W}` is the bin width value (specfied in ``binWidth`` parameter). The first part of the formula ensures that
  the bins are equally spaced from 0, whereas the second part ensures that the minimum gray level intensity inside the
  ROI after binning is always 1.

  In the case where the maximum gray level intensity is equally dividable by the binWidth, i.e.
  :math:`\max(X_{gl}) \mod W = 0`, this will result in that maximum gray level being assigned to bin
  :math:`[\max(X_{gl}), \max(X_{gl}) + W)`, which is consistent with numpy.digitize, but different from the behaviour
  of numpy.histogram, where the final bin has a closed range, including the maximum gray level, i.e.
  :math:`[\max(X_{gl}) - W, \max(X_{gl})]`.

  .. note::
    This method is slightly different from the fixed bin size discretization method described by IBSI. The two most
    notable differences are 1) that PyRadiomics uses a floor division (and adds 1), as opposed to a ceiling division and
    2) that in PyRadiomics, bins are always equally spaced from 0, as opposed to equally spaced from the minimum
    gray level intensity.

  *Example: for a ROI with values ranging from 54 to 166, and a bin width of 25, the bin edges will be [50, 75, 100,
  125, 150, 175].*

  This value can be directly passed to ``numpy.histogram`` to generate a histogram or ``numpy.digitize`` to discretize
  the ROI gray values. See also :py:func:`binImage()`.

  **Fixed bin Count:**

  .. math::
    X_{b, i} = \left\{ {\begin{array}{lcl}
    \lfloor N_b\frac{(X_{gl, i} - \min(X_{gl})}{\max(X_{gl}) - \min(X_{gl})} \rfloor + 1 &
    \mbox{for} & X_{gl, i} < \max(X_{gl}) \\
    N_b & \mbox{for} & X_{gl, i} = \max(X_{gl}) \end{array}} \right.

  Here, :math:`N_b` is the number of bins to use, as defined in ``binCount``.

  References

  - Leijenaar RTH, Nalbantov G, Carvalho S, et al. The effect of SUV discretization in quantitative FDG-PET Radiomics:
    the need for standardized methodology in tumor texture analysis. Sci Rep. 2015;5(August):11075.
  """
  global logger
  binWidth = kwargs.get('binWidth', 25)
  binCount = kwargs.get('binCount')

  if binCount is not None:
    binEdges = numpy.histogram(parameterValues, binCount)[1]
    binEdges[-1] += 1  # Ensures that the maximum value is included in the topmost bin when using numpy.digitize
  else:
    minimum = min(parameterValues)
    maximum = max(parameterValues)

    # Start binning form the first value lesser than or equal to the minimum value and evenly dividable by binwidth
    lowBound = minimum - (minimum % binWidth)
    # Add + 2* binwidth to ensure the maximum value is included in the range generated by numpy.arange, and that values
    # equal to highbound are binned into a separate bin by numpy.histogram (This ensures ALL bins are half open, as
    # numpy.histogram treats the last bin as a closed interval. Moreover, this ensures consistency with numpy.digitize,
    # which will assign len(bins) + 1 to values equal to rightmost bin edge, treating all bins as half-open)
    highBound = maximum + 2 * binWidth

    binEdges = numpy.arange(lowBound, highBound, binWidth)

    # if min(parameterValues) % binWidth = 0 and min(parameterValues) = max(parameterValues), binEdges will only contain
    # 1 value. If this is the case (flat region) ensure that numpy.histogram creates 1 bin (requires 2 edges). For
    # numpy.histogram, a binCount (1) would also suffice, however, this is not accepted by numpy.digitize, which also uses
    # binEdges calculated by this function.
    if len(binEdges) == 1:  # Flat region, ensure that there is 1 bin
      binEdges = [binEdges[0] - .5, binEdges[0] + .5]  # Simulates binEdges returned by numpy.histogram if bins = 1

    logger.debug('Calculated %d bins for bin width %g with edges: %s)', len(binEdges) - 1, binWidth, binEdges)

  return binEdges  # numpy.histogram(parameterValues, bins=binedges)


def binImage(parameterMatrix, parameterMatrixCoordinates=None, **kwargs):
  r"""
  Discretizes the parameterMatrix (matrix representation of the gray levels in the ROI) using the binEdges calculated
  using :py:func:`getBinEdges`. Only voxels defined by parameterMatrixCoordinates (defining the segmentation) are used
  for calculation of histogram and subsequently discretized. Voxels outside segmentation are left unchanged.
  """
  global logger
  logger.debug('Discretizing gray levels inside ROI')

  discretizedParameterMatrix = numpy.zeros(parameterMatrix.shape, dtype='int')
  if parameterMatrixCoordinates is None:
    binEdges = getBinEdges(parameterMatrix.flatten(), **kwargs)
    discretizedParameterMatrix = numpy.digitize(parameterMatrix, binEdges)
  else:
    binEdges = getBinEdges(parameterMatrix[parameterMatrixCoordinates], **kwargs)
    discretizedParameterMatrix[parameterMatrixCoordinates] = numpy.digitize(parameterMatrix[parameterMatrixCoordinates], binEdges)

  return discretizedParameterMatrix, binEdges


def checkMask(imageNode, maskNode, **kwargs):
  """
  Checks whether the Region of Interest (ROI) defined in the mask size and dimensions match constraints, specified in
  settings. The following checks are performed.

  1. Check whether the mask corresponds to the image (i.e. has a similar size, spacing, direction and origin). **N.B.
     This check is performed by SimpleITK, if it fails, an error is logged, with additional error information from
     SimpleITK logged with level DEBUG (i.e. logging-level has to be set to debug to store this information in the log
     file).** The tolerance can be increased using the ``geometryTolerance`` parameter. Alternatively, if the
     ``correctMask`` parameter is ``True``, PyRadiomics will check if the mask contains a valid ROI (inside image
     physical area) and if so, resample the mask to image geometry. See :ref:`radiomics-settings-label` for more info.

  2. Check if the label is present in the mask
  3. Count the number of dimensions in which the size of the ROI > 1 (i.e. does the ROI represent a single voxel (0), a
     line (1), a surface (2) or a volume (3)) and compare this to the minimum number of dimension required (specified in
     ``minimumROIDimensions``).
  4. Optional. Check if there are at least N voxels in the ROI. N is defined in ``minimumROISize``, this test is skipped
     if ``minimumROISize = None``.

  This function returns a tuple of two items. The first item is the bounding box of the mask. The second item is the
  mask that has been corrected by resampling to the input image geometry (if that resampling was successful).

  If a check fails, a ValueError is raised. No features will be extracted for this mask.
  If the mask passes all tests, this function returns the bounding box, which is used in the :py:func:`cropToTumorMask`
  function.

  The bounding box is calculated during (1.) and used for the subsequent checks. The bounding box is
  calculated by SimpleITK.LabelStatisticsImageFilter() and returned as a tuple of indices: (L_x, U_x, L_y, U_y, L_z,
  U_z), where 'L' and 'U' are lower and upper bound, respectively, and 'x', 'y' and 'z' the three image dimensions.

  By reusing the bounding box calculated here, calls to SimpleITK.LabelStatisticsImageFilter() are reduced, improving
  performance.

  Uses the following settings:

  - minimumROIDimensions [1]: Integer, range 1-3, specifies the minimum dimensions (1D, 2D or 3D, respectively).
    Single-voxel segmentations are always excluded.
  - minimumROISize [None]: Integer, > 0,  specifies the minimum number of voxels required. Test is skipped if
    this parameter is set to None.

  .. note::

    If the first check fails there are generally 2 possible causes:

     1. The image and mask are matched, but there is a slight difference in origin, direction or spacing. The exact
        cause, difference and used tolerance are stored with level DEBUG in a log (if enabled). For more information on
        setting up logging, see ":ref:`setting up logging <radiomics-logging-label>`" and the helloRadiomics examples
        (located in the ``pyradiomics/examples`` folder). This problem can be fixed by changing the global tolerance
        (``geometryTolerance`` parameter) or enabling mask correction (``correctMask`` parameter).
     2. The image and mask do not match, but the ROI contained within the mask does represent a physical volume
        contained within the image. If this is the case, resampling is needed to ensure matching geometry between image
        and mask before features can be extracted. This can be achieved by enabling mask correction using the
        ``correctMask`` parameter.
  """
  global logger

  correctedMask = None

  label = kwargs.get('label', 1)
  minDims = kwargs.get('minimumROIDimensions', 2)
  minSize = kwargs.get('minimumROISize', None)

  logger.debug('Checking mask with label %d', label)
  logger.debug('Calculating bounding box')
  # Determine bounds
  lsif = sitk.LabelStatisticsImageFilter()
  try:
    lsif.Execute(imageNode, maskNode)

    # If lsif fails, and mask is corrected, it includes a check whether the label is present. Therefore, perform
    # this test here only if lsif does not fail on the first attempt.
    if label not in lsif.GetLabels():
      raise ValueError('Label (%g) not present in mask' % label)
  except RuntimeError as e:
    # If correctMask = True, try to resample the mask to the image geometry, otherwise return None ("fail")
    if not kwargs.get('correctMask', False):
      if "Both images for LabelStatisticsImageFilter don't match type or dimension!" in e.args[0]:
        logger.debug('Additional information on error.', exc_info=True)
        raise ValueError('Image/Mask datatype or size mismatch. Potential fix: enable correctMask, see '
                         'Documentation:Usage:Customizing the Extraction:Settings:correctMask for more information')
      elif "Inputs do not occupy the same physical space!" in e.args[0]:
        logger.debug('Additional information on error.', exc_info=True)
        raise ValueError('Image/Mask geometry mismatch. Potential fix: increase tolerance using geometryTolerance, '
                         'see Documentation:Usage:Customizing the Extraction:Settings:geometryTolerance for more '
                         'information')
      else:
        raise e  # unhandled error

    logger.warning('Image/Mask geometry mismatch, attempting to correct Mask')

    correctedMask = _correctMask(imageNode, maskNode, **kwargs)  # Raises Value error if ROI outside image physical space

    # Resampling successful, try to calculate boundingbox
    try:
      lsif.Execute(imageNode, correctedMask)
    except RuntimeError:
      logger.debug('Bounding box calculation with resampled mask failed', exc_info=True)
      raise ValueError('Calculation of bounding box failed, for more information run with DEBUG logging and check log')

  # LBound and UBound of the bounding box, as (L_X, U_X, L_Y, U_Y, L_Z, U_Z)
  boundingBox = numpy.array(lsif.GetBoundingBox(label))

  logger.debug('Checking minimum number of dimensions requirements (%d)', minDims)
  ndims = numpy.sum((boundingBox[1::2] - boundingBox[0::2] + 1) > 1)  # UBound - LBound + 1 = Size
  if ndims == 0:
    raise ValueError('mask only contains 1 segmented voxel! Cannot extract features for a single voxel.')
  elif ndims < minDims:
    raise ValueError('mask has too few dimensions (number of dimensions %d, minimum required %d)' % (ndims, minDims))

  if minSize is not None:
    logger.debug('Checking minimum size requirements (minimum size: %d)', minSize)
    roiSize = lsif.GetCount(label)
    if roiSize <= minSize:
      raise ValueError('Size of the ROI is too small (minimum size: %g, ROI size: %g' % (minSize, roiSize))

  return boundingBox, correctedMask


def _correctMask(imageNode, maskNode, **kwargs):
  """
  If the mask geometry does not match the image geometry, this function can be used to resample the mask to the image
  physical space.

  First, the mask is checked for a valid ROI (i.e. maskNode contains an ROI with the given label value, which does not
  include areas outside of the physical image bounds).

  If the ROI is valid, the maskNode is resampled using the imageNode as a reference image and a nearest neighbor
  interpolation.

  If the ROI is valid, the resampled mask is returned, otherwise ``None`` is returned.
  """
  global logger
  logger.debug('Resampling mask to image geometry')

  _checkROI(imageNode, maskNode, **kwargs)  # Raises a value error if ROI is invalid

  rif = sitk.ResampleImageFilter()
  rif.SetReferenceImage(imageNode)
  rif.SetInterpolator(sitk.sitkNearestNeighbor)

  logger.debug('Resampling...')

  return rif.Execute(maskNode)


def _checkROI(imageNode, maskNode, **kwargs):
  """
  Check whether maskNode contains a valid ROI defined by label:

  1. Check whether the label value is present in the maskNode.
  2. Check whether the ROI defined by the label does not include an area outside the physical area of the image.

  For the second check, a tolerance of 1e-3 is allowed.

  If the ROI is valid, the bounding box (lower bounds, followd by size in all dimensions (X, Y, Z ordered)) is
  returned. Otherwise, a ValueError is raised.
  """
  global logger
  label = kwargs.get('label', 1)

  logger.debug('Checking ROI validity')

  # Determine bounds of cropped volume in terms of original Index coordinate space
  lssif = sitk.LabelShapeStatisticsImageFilter()
  lssif.Execute(maskNode)

  logger.debug('Checking if label %d is persent in the mask', label)
  if label not in lssif.GetLabels():
    raise ValueError('Label (%d) not present in mask', label)

  # LBound and size of the bounding box, as (L_X, L_Y, [L_Z], S_X, S_Y, [S_Z])
  bb = numpy.array(lssif.GetBoundingBox(label))
  Nd = maskNode.GetDimension()

  # Determine if the ROI is within the physical space of the image

  logger.debug('Comparing physical space of bounding box to physical space of image')
  # Step 1: Get the origin and UBound corners of the bounding box in physical space
  # The additional 0.5 represents the difference between the voxel center and the voxel corner
  # Upper bound index of ROI = bb[:Nd] + bb[Nd:] - 1 (LBound + Size - 1), .5 is added to get corner
  ROIBounds = (maskNode.TransformContinuousIndexToPhysicalPoint(bb[:Nd] - .5),  # Origin
               maskNode.TransformContinuousIndexToPhysicalPoint(bb[:Nd] + bb[Nd:] - 0.5))  # UBound
  # Step 2: Translate the ROI physical bounds to the image coordinate space
  ROIBounds = (imageNode.TransformPhysicalPointToContinuousIndex(ROIBounds[0]),  # Origin
               imageNode.TransformPhysicalPointToContinuousIndex(ROIBounds[1]))

  logger.debug('ROI bounds (image coordinate space): %s', ROIBounds)

  # Check if any of the ROI bounds are outside the image indices (i.e. -0.5 < ROI < Im.Size -0.5)
  # The additional 0.5 is to allow for different spacings (defines the edges, not the centers of the edge-voxels
  tolerance = 1e-3  # Define a tolerance to correct for machine precision errors
  if numpy.any(numpy.min(ROIBounds, axis=0) < (- .5 - tolerance)) or \
     numpy.any(numpy.max(ROIBounds, axis=0) > (numpy.array(imageNode.GetSize()) - .5 + tolerance)):
    raise ValueError('Bounding box of ROI is larger than image space:\n\t'
                     'ROI bounds (x, y, z image coordinate space) %s\n\tImage Size %s' %
                     (ROIBounds, imageNode.GetSize()))

  logger.debug('ROI valid, calculating resampling grid')

  return bb


def cropToTumorMask(imageNode, maskNode, boundingBox, **kwargs):
  """
  Create a sitkImage of the segmented region of the image based on the input label.

  Create a sitkImage of the labelled region of the image, cropped to have a
  cuboid shape equal to the ijk boundaries of the label.

  :param boundingBox: The bounding box used to crop the image. This is the bounding box as returned by
    :py:func:`checkMask`.
  :param label: [1], value of the label, onto which the image and mask must be cropped.
  :return: Cropped image and mask (SimpleITK image instances).

  """
  global logger
  padDistance = kwargs.get('padDistance', 0)

  size = numpy.array(maskNode.GetSize())

  ijkMinBounds = boundingBox[0::2] - padDistance
  ijkMaxBounds = size - boundingBox[1::2] - padDistance - 1

  # Ensure cropped area is not outside original image bounds
  ijkMinBounds = numpy.maximum(ijkMinBounds, 0)
  ijkMaxBounds = numpy.maximum(ijkMaxBounds, 0)

  # Crop Image
  logger.debug('Cropping to size %s', (boundingBox[1::2] - boundingBox[0::2]) + 1)
  cif = sitk.CropImageFilter()
  try:
    cif.SetLowerBoundaryCropSize(ijkMinBounds)
    cif.SetUpperBoundaryCropSize(ijkMaxBounds)
  except TypeError:
    # newer versions of SITK/python want a tuple or list
    cif.SetLowerBoundaryCropSize(ijkMinBounds.tolist())
    cif.SetUpperBoundaryCropSize(ijkMaxBounds.tolist())
  croppedImageNode = cif.Execute(imageNode)
  croppedMaskNode = cif.Execute(maskNode)

  return croppedImageNode, croppedMaskNode


def resampleImage(imageNode, maskNode, **kwargs):
  """
  Resamples image and mask to the specified pixel spacing (The default interpolator is Bspline).

  Resampling can be enabled using the settings 'interpolator' and 'resampledPixelSpacing' in the parameter file or as
  part of the settings passed to the feature extractor. See also
  :ref:`feature extractor <radiomics-featureextractor-label>`.

  'imageNode' and 'maskNode' are SimpleITK Objects, and 'resampledPixelSpacing' is the output pixel spacing (sequence of
  3 elements).

  If only in-plane resampling is required, set the output pixel spacing for the out-of-plane dimension (usually the last
  dimension) to 0. Spacings with a value of 0 are replaced by the spacing as it is in the original mask.

  Only part of the image and labelmap are resampled. The resampling grid is aligned to the input origin, but only voxels
  covering the area of the image ROI (defined by the bounding box) and the padDistance are resampled. This results in a
  resampled and partially cropped image and mask. Additional padding is required as some filters also sample voxels
  outside of segmentation boundaries. For feature calculation, image and mask are cropped to the bounding box without
  any additional padding, as the feature classes do not need the gray level values outside the segmentation.

  The resampling grid is calculated using only the input mask. Even when image and mask have different directions, both
  the cropped image and mask will have the same direction (equal to direction of the mask). Spacing and size are
  determined by settings and bounding box of the ROI.

  .. note::
    Before resampling the bounds of the non-padded ROI are compared to the bounds. If the ROI bounding box includes
    areas outside of the physical space of the image, an error is logged and (None, None) is returned. No features will
    be extracted. This enables the input image and mask to have different geometry, so long as the ROI defines an area
    within the image.

  .. note::
    The additional padding is adjusted, so that only the physical space within the mask is resampled. This is done to
    prevent resampling outside of the image. Please note that this assumes the image and mask to image the same physical
    space. If this is not the case, it is possible that voxels outside the image are included in the resampling grid,
    these will be assigned a value of 0. It is therefore recommended, but not enforced, to use an input mask which has
    the same or a smaller physical space than the image.
  """
  global logger
  resampledPixelSpacing = kwargs['resampledPixelSpacing']
  interpolator = kwargs.get('interpolator', sitk.sitkBSpline)
  padDistance = kwargs.get('padDistance', 5)
  label = kwargs.get('label', 1)

  logger.debug('Resampling image and mask')

  if imageNode is None or maskNode is None:
    raise ValueError('Requires both image and mask to resample')

  maskSpacing = numpy.array(maskNode.GetSpacing())
  imageSpacing = numpy.array(imageNode.GetSpacing())

  Nd_resampled = len(resampledPixelSpacing)
  Nd_mask = len(maskSpacing)
  assert Nd_resampled == Nd_mask, \
      'Wrong dimensionality (%i-D) of resampledPixelSpacing!, %i-D required' % (Nd_resampled, Nd_mask)

  # If spacing for a direction is set to 0, use the original spacing (enables "only in-slice" resampling)
  logger.debug('Where resampled spacing is set to 0, set it to the original spacing (mask)')
  resampledPixelSpacing = numpy.array(resampledPixelSpacing)
  resampledPixelSpacing = numpy.where(resampledPixelSpacing == 0, maskSpacing, resampledPixelSpacing)

  # Check if the maskNode contains a valid ROI. If ROI is valid, the bounding box needed to calculate the resampling
  # grid is returned.
  bb = _checkROI(imageNode, maskNode, **kwargs)

  # Do not resample in those directions where labelmap spans only one slice.
  maskSize = numpy.array(maskNode.GetSize())
  resampledPixelSpacing = numpy.where(bb[Nd_mask:] != 1, resampledPixelSpacing, maskSpacing)

  # If current spacing is equal to resampledPixelSpacing, no interpolation is needed
  # Tolerance = 1e-5 + 1e-8*abs(resampledSpacing)
  logger.debug('Comparing resampled spacing to original spacing (image')
  if numpy.allclose(imageSpacing, resampledPixelSpacing):
    logger.info('New spacing equal to original image spacing, just resampling the mask')

    # Ensure that image and mask geometry match
    rif = sitk.ResampleImageFilter()
    rif.SetReferenceImage(imageNode)
    rif.SetInterpolator(sitk.sitkNearestNeighbor)
    maskNode = rif.Execute(maskNode)

    # re-calculate the bounding box of the mask
    lssif = sitk.LabelShapeStatisticsImageFilter()
    lssif.Execute(maskNode)
    bb = numpy.array(lssif.GetBoundingBox(label))

    low_up_bb = numpy.empty(Nd_mask * 2, dtype=int)
    low_up_bb[::2] = bb[:Nd_mask]
    low_up_bb[1::2] = bb[:Nd_mask] + bb[Nd_mask:] - 1
    return cropToTumorMask(imageNode, maskNode, low_up_bb, **kwargs)

  spacingRatio = maskSpacing / resampledPixelSpacing

  # Determine bounds of cropped volume in terms of new Index coordinate space,
  # round down for lowerbound and up for upperbound to ensure entire segmentation is captured (prevent data loss)
  # Pad with an extra .5 to prevent data loss in case of upsampling. For Ubound this is (-1 + 0.5 = -0.5)
  bbNewLBound = numpy.floor((bb[:Nd_mask] - 0.5) * spacingRatio - padDistance)
  bbNewUBound = numpy.ceil((bb[:Nd_mask] + bb[Nd_mask:] - 0.5) * spacingRatio + padDistance)

  # Ensure resampling is not performed outside bounds of original image
  maxUbound = numpy.ceil(maskSize * spacingRatio) - 1
  bbNewLBound = numpy.where(bbNewLBound < 0, 0, bbNewLBound)
  bbNewUBound = numpy.where(bbNewUBound > maxUbound, maxUbound, bbNewUBound)

  # Calculate the new size. Cast to int to prevent error in sitk.
  newSize = numpy.array(bbNewUBound - bbNewLBound + 1, dtype='int').tolist()

  # Determine continuous index of bbNewLBound in terms of the original Index coordinate space
  bbOriginalLBound = bbNewLBound / spacingRatio

  # Origin is located in center of first voxel, e.g. 1/2 of the spacing
  # from Corner, which corresponds to 0 in the original Index coordinate space.
  # The new spacing will be in 0 the new Index coordinate space. Here we use continuous
  # index to calculate where the new 0 of the new Index coordinate space (of the original volume
  # in terms of the original spacing, and add the minimum bounds of the cropped area to
  # get the new Index coordinate space of the cropped volume in terms of the original Index coordinate space.
  # Then use the ITK functionality to bring the continuous index into the physical space (mm)
  newOriginIndex = numpy.array(.5 * (resampledPixelSpacing - maskSpacing) / maskSpacing)
  newCroppedOriginIndex = newOriginIndex + bbOriginalLBound
  newOrigin = maskNode.TransformContinuousIndexToPhysicalPoint(newCroppedOriginIndex)

  imagePixelType = imageNode.GetPixelID()
  maskPixelType = maskNode.GetPixelID()

  direction = numpy.array(maskNode.GetDirection())

  logger.info('Applying resampling from spacing %s and size %s to spacing %s and size %s',
              maskSpacing, maskSize, resampledPixelSpacing, newSize)

  try:
    if isinstance(interpolator, six.string_types):
      interpolator = getattr(sitk, interpolator)
  except Exception:
    logger.warning('interpolator "%s" not recognized, using sitkBSpline', interpolator)
    interpolator = sitk.sitkBSpline

  rif = sitk.ResampleImageFilter()

  rif.SetOutputSpacing(resampledPixelSpacing)
  rif.SetOutputDirection(direction)
  rif.SetSize(newSize)
  rif.SetOutputOrigin(newOrigin)

  logger.debug('Resampling image')
  rif.SetOutputPixelType(imagePixelType)
  rif.SetInterpolator(interpolator)
  resampledImageNode = rif.Execute(imageNode)

  logger.debug('Resampling mask')
  rif.SetOutputPixelType(maskPixelType)
  rif.SetInterpolator(sitk.sitkNearestNeighbor)
  resampledMaskNode = rif.Execute(maskNode)

  return resampledImageNode, resampledMaskNode


def normalizeImage(image, **kwargs):
  r"""
  Normalizes the image by centering it at the mean with standard deviation. Normalization is based on all gray values in
  the image, not just those inside the segmentation.

  :math:`f(x) = \frac{s(x - \mu_x)}{\sigma_x}`

  Where:

  - :math:`x` and :math:`f(x)` are the original and normalized intensity, respectively.
  - :math:`\mu_x` and :math:`\sigma_x` are the mean and standard deviation of the image instensity values.
  - :math:`s` is an optional scaling defined by ``scale``. By default, it is set to 1.

  Optionally, outliers can be removed, in which case values for which :math:`x > \mu_x + n\sigma_x` or
  :math:`x < \mu_x - n\sigma_x` are set to :math:`\mu_x + n\sigma_x` and :math:`\mu_x - n\sigma_x`, respectively.
  Here, :math:`n>0` and defined by ``outliers``. This, in turn, is controlled by the ``removeOutliers`` parameter.
  Removal of outliers is done after the values of the image are normalized, but before ``scale`` is applied.
  """
  global logger
  scale = kwargs.get('normalizeScale', 1)
  outliers = kwargs.get('removeOutliers')

  logger.debug('Normalizing image with scale %d', scale)
  image = sitk.Normalize(image)

  if outliers is not None:
    logger.debug('Removing outliers > %g standard deviations', outliers)
    imageArr = sitk.GetArrayFromImage(image)

    imageArr[imageArr > outliers] = outliers
    imageArr[imageArr < -outliers] = -outliers

    newImage = sitk.GetImageFromArray(imageArr)
    newImage.CopyInformation(image)
    image = newImage

  image *= scale

  return image


def resegmentMask(imageNode, maskNode, **kwargs):
  r"""
  Resegment the Mask based on the range specified by the threshold(s) in ``resegmentRange``. Either 1 or 2 thresholds
  can be defined. In case of 1 threshold, all values equal to or higher than that threshold are included. If there are
  2 thresholds, all voxels with a value inside the closed-range defined by these thresholds is included
  (i.e. a voxels is included if :math:`T_{lower} \leq X_gl \leq T_{upper}`).
  The resegmented mask is therefore always equal or smaller in size than the original mask.
  In the case where either resegmentRange or resegmentMode contains illigal values, a ValueError is raised.

  There are 3 modes for defining the threshold:

  1. absolute (default): The values in resegmentRange define  as absolute values (i.e. corresponding to the gray values
     in the image
  2. relative: The values in resegmentRange define the threshold as relative to the maximum value found in the ROI.
     (e.g. 0.5 indicates a threshold at 50% of maximum gray value)
  3. sigma: The threshold is defined as the number of sigma from the mean. (e.g. resegmentRange [-3, 3] will include
     all voxels that have a value that differs 3 or less standard deviations from the mean).

  """
  global logger
  resegmentRange = kwargs['resegmentRange']
  resegmentMode = kwargs.get('resegmentMode', 'absolute')
  label = kwargs.get('label', 1)

  if resegmentRange is None:
    raise ValueError('resegmentRange is None.')
  if len(resegmentRange) == 0 or len(resegmentRange) > 2:
    raise ValueError('Length %i is not allowed for resegmentRange' % len(resegmentRange))

  logger.debug('Resegmenting mask (range %s, mode %s)', resegmentRange, resegmentMode)

  im_arr = sitk.GetArrayFromImage(imageNode)
  ma_arr = (sitk.GetArrayFromImage(maskNode) == label)  # boolean array

  oldSize = numpy.sum(ma_arr)

  if resegmentMode == 'absolute':
    logger.debug('Resegmenting in absolute mode')
    thresholds = sorted(resegmentRange)
  elif resegmentMode == 'relative':
    max_gl = numpy.max(im_arr[ma_arr])
    logger.debug('Resegmenting in relative mode, max %g', max_gl)
    thresholds = [max_gl * th for th in sorted(resegmentRange)]
  elif resegmentMode == 'sigma':
    mean_gl = numpy.mean(im_arr[ma_arr])
    sd_gl = numpy.std(im_arr[ma_arr])
    logger.debug('Resegmenting in sigma mode, mean %g, std %g', mean_gl, sd_gl)
    thresholds = [mean_gl + sd_gl * th for th in sorted(resegmentRange)]
  else:
    raise ValueError('Resegment mode %s not recognized.' % resegmentMode)

  # Apply lower threshold
  logger.debug('Applying lower threshold (%g)', thresholds[0])
  ma_arr[ma_arr] = im_arr[ma_arr] >= thresholds[0]

  # If 2 thresholds are defined, also apply an upper threshold
  if len(thresholds) == 2:
    logger.debug('Applying upper threshold (%g)', thresholds[1])
    ma_arr[ma_arr] = im_arr[ma_arr] <= thresholds[1]

  roiSize = numpy.sum(ma_arr)

  if roiSize <= 1:
    raise ValueError("Resegmentation excluded too many voxels with label %i (retained %i voxel(s))! "
                     "Cannot extract features" % (label, roiSize))

  # Transform the boolean array back to an image with the correct voxels set to the label value
  newMask_arr = numpy.zeros(ma_arr.shape, dtype='int')
  newMask_arr[ma_arr] = label

  newMask = sitk.GetImageFromArray(newMask_arr)
  newMask.CopyInformation(maskNode)
  logger.debug('Resegmentation complete, new size: %d voxels (excluded %d voxels)', roiSize, oldSize - roiSize)

  return newMask


def getOriginalImage(inputImage, inputMask, **kwargs):
  """
  This function does not apply any filter, but returns the original image. This function is needed to
  dynamically expose the original image as a valid image type.

  :return: Yields original image, 'original' and ``kwargs``
  """
  global logger
  logger.debug('Yielding original image')
  yield inputImage, 'original', kwargs


def getLoGImage(inputImage, inputMask, **kwargs):
  r"""
  Applies a Laplacian of Gaussian filter to the input image and yields a derived image for each sigma value specified.

  A Laplacian of Gaussian image is obtained by convolving the image with the second derivative (Laplacian) of a Gaussian
  kernel.

  The Gaussian kernel is used to smooth the image and is defined as

  .. math::

    G(x, y, z, \sigma) = \frac{1}{(\sigma \sqrt{2 \pi})^3}e^{-\frac{x^2 + y^2 + z^2}{2\sigma^2}}

  The Gaussian kernel is convolved by the laplacian kernel :math:`\nabla^2G(x, y, z)`, which is sensitive to areas with
  rapidly changing intensities, enhancing edges. The width of the filter in the Gaussian kernel is determined by
  :math:`\sigma` and can be used to emphasize more fine (low :math:`\sigma` values) or coarse (high :math:`\sigma`
  values) textures.

  .. warning::

    The LoG filter implemented in PyRadiomics is a 3D LoG filter, and therefore requires 3D input. Features using a
    single slice (2D) segmentation can still be extracted, but the input image *must* be a 3D image, with a minimum size
    in all dimensions :math:`\geq \sigma`. If input image is too small, a warning is logged and :math:`\sigma` value is
    skipped. Moreover, the image size *must* be at least 4 voxels in each dimensions, if this constraint is not met, no
    LoG derived images can be generated.

  Following settings are possible:

  - sigma: List of floats or integers, must be greater than 0. Filter width (mm) to use for the Gaussian kernel
    (determines coarseness).

  .. warning::
    Setting for sigma must be provided. If omitted, no LoG image features are calculated and the function
    will return an empty dictionary.

  Returned filter name reflects LoG settings:
  log-sigma-<sigmaValue>-3D.

  References:

  - `SimpleITK Doxygen documentation
    <https://itk.org/SimpleITKDoxygen/html/classitk_1_1simple_1_1LaplacianRecursiveGaussianImageFilter.html>`_
  - `ITK Doxygen documentation <https://itk.org/Doxygen/html/classitk_1_1LaplacianRecursiveGaussianImageFilter.html>`_
  - `<https://en.wikipedia.org/wiki/Blob_detection#The_Laplacian_of_Gaussian>`_

  :return: Yields log filtered image for each specified sigma, corresponding image type name and ``kwargs`` (customized
    settings).
  """
  global logger

  logger.debug('Generating LoG images')

  # Check if size of image is > 4 in all 3D directions (otherwise, LoG filter will fail)
  size = numpy.array(inputImage.GetSize())
  spacing = numpy.array(inputImage.GetSpacing())

  if numpy.min(size) < 4:
    logger.warning('Image too small to apply LoG filter, size: %s', size)
    return

  sigmaValues = kwargs.get('sigma', [])

  for sigma in sigmaValues:
    logger.info('Computing LoG with sigma %g', sigma)

    if sigma > 0.0:
      if numpy.all(size >= numpy.ceil(sigma / spacing) + 1):
        lrgif = sitk.LaplacianRecursiveGaussianImageFilter()
        lrgif.SetNormalizeAcrossScale(True)
        lrgif.SetSigma(sigma)
        inputImageName = 'log-sigma-%s-mm-3D' % (str(sigma).replace('.', '-'))
        logger.debug('Yielding %s image', inputImageName)
        yield lrgif.Execute(inputImage), inputImageName, kwargs
      else:
        logger.warning('applyLoG: sigma(%g)/spacing(%s) + 1 must be greater than the size(%s) of the inputImage',
                       sigma,
                       spacing,
                       size)
    else:
      logger.warning('applyLoG: sigma must be greater than 0.0: %g', sigma)


def getWaveletImage(inputImage, inputMask, **kwargs):
  """
  Applies wavelet filter to the input image and yields the decompositions and the approximation.

  Following settings are possible:

  - start_level [0]: integer, 0 based level of wavelet which should be used as first set of decompositions
    from which a signature is calculated
  - level [1]: integer, number of levels of wavelet decompositions from which a signature is calculated.
  - wavelet ["coif1"]: string, type of wavelet decomposition. Enumerated value, validated against possible values
    present in the ``pyWavelet.wavelist()``. Current possible values (pywavelet version 0.4.0) (where an
    aditional number is needed, range of values is indicated in []):

    - haar
    - dmey
    - sym[2-20]
    - db[1-20]
    - coif[1-5]
    - bior[1.1, 1.3, 1.5, 2.2, 2.4, 2.6, 2.8, 3.1, 3.3, 3.5, 3.7, 3.9, 4.4, 5.5, 6.8]
    - rbio[1.1, 1.3, 1.5, 2.2, 2.4, 2.6, 2.8, 3.1, 3.3, 3.5, 3.7, 3.9, 4.4, 5.5, 6.8]

  Returned filter name reflects wavelet type:
  wavelet[level]-<decompositionName>

  N.B. only levels greater than the first level are entered into the name.

  :return: Yields each wavelet decomposition and final approximation, corresponding imaget type name and ``kwargs``
    (customized settings).
  """
  global logger

  logger.debug('Generating Wavelet images')

  Nd = inputImage.GetDimension()
  axes = list(range(Nd - 1, -1, -1))
  if kwargs.get('force2D', False):
    axes.remove(kwargs.get('force2Ddimension', 0))

  approx, ret = _swt3(inputImage, tuple(axes), **kwargs)

  for idx, wl in enumerate(ret, start=1):
    for decompositionName, decompositionImage in wl.items():
      logger.info('Computing Wavelet %s', decompositionName)

      if idx == 1:
        inputImageName = 'wavelet-%s' % (decompositionName)
      else:
        inputImageName = 'wavelet%s-%s' % (idx, decompositionName)
      logger.debug('Yielding %s image', inputImageName)
      yield decompositionImage, inputImageName, kwargs

  if len(ret) == 1:
    inputImageName = 'wavelet-%s' % ('L' * len(axes))
  else:
    inputImageName = 'wavelet%s-%s' % (len(ret), ('L' * len(axes)))
  logger.debug('Yielding approximation (%s) image', inputImageName)
  yield approx, inputImageName, kwargs


def _swt3(inputImage, axes, **kwargs):  # Stationary Wavelet Transform 3D
  wavelet = kwargs.get('wavelet', 'coif1')
  level = kwargs.get('level', 1)
  start_level = kwargs.get('start_level', 0)

  matrix = sitk.GetArrayFromImage(inputImage)  # This function gets a numpy array from the SimpleITK Image "inputImage"
  matrix = numpy.asarray(matrix) # The function np.asarray converts "matrix" (which could be also a tuple) into an array.

  original_shape = matrix.shape
  # original_shape becomes a tuple (?,?,?) containing the number of rows, columns, and slices of the image
  # this is of course dependent on the number of dimensions, but the same principle holds
  padding = tuple([(0, 1 if dim % 2 != 0 else 0) for dim in original_shape])
  # padding is necessary because of pywt.swtn (see function Notes)
  data = matrix.copy()  # creates a modifiable copy of "matrix" and we call it "data"
  data = numpy.pad(data, padding, 'wrap')  # padding the tuple "padding" previously computed

  if not isinstance(wavelet, pywt.Wavelet):
    wavelet = pywt.Wavelet(wavelet)

  for i in range(0, start_level):  # if start_level = 0 (default) this for loop never gets executed
    # compute all decompositions and saves them in "dec" dict
    dec = pywt.swtn(data, wavelet, level=1, start_level=0, axes=axes)[0]
    # copies in "data" just the "aaa" decomposition (i.e. approximation; No of consecutive 'a's = len(axes))
    data = dec['a' * len(axes)].copy()

  ret = []  # initialize empty list
  for i in range(start_level, start_level + level):
    # compute the n-dimensional stationary wavelet transform
    dec = pywt.swtn(data, wavelet, level=1, start_level=0, axes=axes)[0]
    # Copy the approximation into data (approximation in output / input for next levels)
    data = dec['a' * len(axes)].copy()

    dec_im = {}  # initialize empty dict
    for decName, decImage in six.iteritems(dec):
      # Returning the approximiation is done only for the last loop,
      # and is handled separately below (by building it from `data`)
      # There for, skip it here
      if decName == 'a' * len(axes):
        continue
      decTemp = decImage.copy()
      decTemp = decTemp[tuple(slice(None, -1 if dim % 2 != 0 else None) for dim in original_shape)]
      sitkImage = sitk.GetImageFromArray(decTemp)
      sitkImage.CopyInformation(inputImage)
      dec_im[str(decName).replace('a', 'L').replace('d', 'H')] = sitkImage
      # modifies 'a' with 'L' (Low-pass filter) and 'd' with 'H' (High-pass filter)

    ret.append(dec_im)  # appending all the filtered sitk images (stored in "dec_im") to the "ret" list

  data = data[tuple(slice(None, -1 if dim % 2 != 0 else None) for dim in original_shape)]
  approximation = sitk.GetImageFromArray(data)
  approximation.CopyInformation(inputImage)

  return approximation, ret  # returns the approximation and the detail (ret) coefficients of the stationary wavelet decomposition


def getSquareImage(inputImage, inputMask, **kwargs):
  r"""
  Computes the square of the image intensities.

  Resulting values are rescaled on the range of the initial original image and negative intensities are made
  negative in resultant filtered image.

  :math:`f(x) = (cx)^2,\text{ where } c=\displaystyle\frac{1}{\sqrt{\max(|x|)}}`

  Where :math:`x` and :math:`f(x)` are the original and filtered intensity, respectively.

  :return: Yields square filtered image, 'square' and ``kwargs`` (customized settings).
  """
  global logger

  im = sitk.GetArrayFromImage(inputImage)
  im = im.astype('float64')
  coeff = 1 / numpy.sqrt(numpy.max(numpy.abs(im)))
  im = (coeff * im) ** 2
  im = sitk.GetImageFromArray(im)
  im.CopyInformation(inputImage)

  logger.debug('Yielding square image')
  yield im, 'square', kwargs


def getSquareRootImage(inputImage, inputMask, **kwargs):
  r"""
  Computes the square root of the absolute value of image intensities.

  Resulting values are rescaled on the range of the initial original image and negative intensities are made
  negative in resultant filtered image.

  :math:`f(x) = \left\{ {\begin{array}{lcl}
  \sqrt{cx} & \mbox{for} & x \ge 0 \\
  -\sqrt{-cx} & \mbox{for} & x < 0\end{array}} \right.,\text{ where } c=\max(|x|)`

  Where :math:`x` and :math:`f(x)` are the original and filtered intensity, respectively.

  :return: Yields square root filtered image, 'squareroot' and ``kwargs`` (customized settings).
  """
  global logger

  im = sitk.GetArrayFromImage(inputImage)
  im = im.astype('float64')
  coeff = numpy.max(numpy.abs(im))
  im[im > 0] = numpy.sqrt(im[im > 0] * coeff)
  im[im < 0] = - numpy.sqrt(-im[im < 0] * coeff)
  im = sitk.GetImageFromArray(im)
  im.CopyInformation(inputImage)

  logger.debug('Yielding squareroot image')
  yield im, 'squareroot', kwargs


def getLogarithmImage(inputImage, inputMask, **kwargs):
  r"""
  Computes the logarithm of the absolute value of the original image + 1.

  Resulting values are rescaled on the range of the initial original image and negative intensities are made
  negative in resultant filtered image.

  :math:`f(x) = \left\{ {\begin{array}{lcl}
  c\log{(x + 1)} & \mbox{for} & x \ge 0 \\
  -c\log{(-x + 1)} & \mbox{for} & x < 0\end{array}} \right. \text{, where } c=\frac{\max(|x|)}{\log(\max(|x|) + 1)}`

  Where :math:`x` and :math:`f(x)` are the original and filtered intensity, respectively.

  :return: Yields logarithm filtered image, 'logarithm' and ``kwargs`` (customized settings)
  """
  global logger

  im = sitk.GetArrayFromImage(inputImage)
  im = im.astype('float64')
  im_max = numpy.max(numpy.abs(im))
  im[im > 0] = numpy.log(im[im > 0] + 1)
  im[im < 0] = - numpy.log(- (im[im < 0] - 1))
  im = im * (im_max / numpy.max(numpy.abs(im)))
  im = sitk.GetImageFromArray(im)
  im.CopyInformation(inputImage)

  logger.debug('Yielding logarithm image')
  yield im, 'logarithm', kwargs


def getExponentialImage(inputImage, inputMask, **kwargs):
  r"""
  Computes the exponential of the original image.

  Resulting values are rescaled on the range of the initial original image.

  :math:`f(x) = e^{cx},\text{ where } c=\displaystyle\frac{\log(\max(|x|))}{\max(|x|)}`

  Where :math:`x` and :math:`f(x)` are the original and filtered intensity, respectively.

  :return: Yields exponential filtered image, 'exponential' and ``kwargs`` (customized settings)
  """
  global logger

  im = sitk.GetArrayFromImage(inputImage)
  im = im.astype('float64')
  im_max = numpy.max(numpy.abs(im))
  coeff = numpy.log(im_max) / im_max
  im = numpy.exp(coeff * im)
  im = sitk.GetImageFromArray(im)
  im.CopyInformation(inputImage)

  logger.debug('Yielding exponential image')
  yield im, 'exponential', kwargs


def getGradientImage(inputImage, inputMask, **kwargs):
  r"""
  Compute and return the Gradient Magnitude in the image.
  By default, takes into account the image spacing, this can be switched off by specifying
  ``gradientUseSpacing = False``.

  References:

  - `SimpleITK documentation
    <https://itk.org/SimpleITKDoxygen/html/classitk_1_1simple_1_1GradientMagnitudeImageFilter.html>`_
  - `<https://en.wikipedia.org/wiki/Image_gradient>`_
  """
  gmif = sitk.GradientMagnitudeImageFilter()
  gmif.SetUseImageSpacing(kwargs.get('gradientUseSpacing', True))
  im = gmif.Execute(inputImage)
  yield im, 'gradient', kwargs


def getLBP2DImage(inputImage, inputMask, **kwargs):
  """
  Compute and return the Local Binary Pattern (LBP) in 2D. If ``force2D`` is set to false (= feature extraction in 3D) a
  warning is logged, as this filter processes the image in a by-slice operation. The plane in which the LBP is
  applied can be controlled by the ``force2Ddimension`` parameter (see also :py:func:`generateAngles`).

  Following settings are possible (in addition to ``force2Ddimension``):

    - ``lbp2DRadius`` [1]: Float, specifies the radius in which the neighbours should be sampled
    - ``lbp2DSamples`` [9]: Integer, specifies the number of samples to use
    - ``lbp2DMethod`` ['uniform']: String, specifies the method for computing the LBP to use.

  For more information see `scikit documentation
  <http://scikit-image.org/docs/dev/api/skimage.feature.html#skimage.feature.local_binary_pattern>`_

  :return: Yields LBP filtered image, 'lbp-2D' and ``kwargs`` (customized settings)

  .. note::
    LBP can often return only a very small number of different gray levels. A customized bin width is often needed.
  .. warning::
    Requires package ``scikit-image`` to function. If not available, this filter logs a warning and does not yield an image.

  References:

  - T. Ojala, M. Pietikainen, and D. Harwood (1994), "Performance evaluation of texture measures with classification
    based on Kullback discrimination of distributions", Proceedings of the 12th IAPR International Conference on Pattern
    Recognition (ICPR 1994), vol. 1, pp. 582 - 585.
  - T. Ojala, M. Pietikainen, and D. Harwood (1996), "A Comparative Study of Texture Measures with Classification Based
    on Feature Distributions", Pattern Recognition, vol. 29, pp. 51-59.
  """
  global logger
  try:
    from skimage.feature import local_binary_pattern
  except ImportError:
    logger.warning('Could not load required package "skimage", cannot implement filter LBP 2D')
    return

  lbp_radius = kwargs.get('lbp2DRadius', 1)
  lbp_samples = kwargs.get('lbp2DSamples', 8)
  lbp_method = kwargs.get('lbp2DMethod', 'uniform')

  im_arr = sitk.GetArrayFromImage(inputImage)

  Nd = inputImage.GetDimension()
  if Nd == 3:
    # Warn the user if features are extracted in 3D, as this function calculates LBP in 2D
    if not kwargs.get('force2D', False):
      logger.warning('Calculating Local Binary Pattern in 2D, but extracting features in 3D. Use with caution!')
    lbp_axis = kwargs.get('force2Ddimension', 0)

    im_arr = im_arr.swapaxes(0, lbp_axis)
    for idx in range(im_arr.shape[0]):
      im_arr[idx, ...] = local_binary_pattern(im_arr[idx, ...], P=lbp_samples, R=lbp_radius, method=lbp_method)
    im_arr = im_arr.swapaxes(0, lbp_axis)
  elif Nd == 2:
    im_arr = local_binary_pattern(im_arr, P=lbp_samples, R=lbp_radius, method=lbp_method)
  else:
    logger.warning('LBP 2D is only available for 2D or 3D with forced 2D extraction')
    return

  im = sitk.GetImageFromArray(im_arr)
  im.CopyInformation(inputImage)

  yield im, 'lbp-2D', kwargs


def getLBP3DImage(inputImage, inputMask, **kwargs):
  """
  Compute and return the Local Binary Pattern (LBP) in 3D using spherical harmonics.
  If ``force2D`` is set to true (= feature extraction in 2D) a warning is logged.

  LBP is only calculated for voxels segmented in the mask

  Following settings are possible:

    - ``lbp3DLevels`` [2]: integer, specifies the the number of levels in spherical harmonics to use.
    - ``lbp3DIcosphereRadius`` [1]: Float, specifies the radius in which the neighbours should be sampled
    - ``lbp3DIcosphereSubdivision`` [1]: Integer, specifies the number of subdivisions to apply in the icosphere

  :return: Yields LBP filtered image for each level, 'lbp-3D-m<level>' and ``kwargs`` (customized settings).
           Additionally yields the kurtosis image, 'lbp-3D-k' and ``kwargs``.

  .. note::
    LBP can often return only a very small number of different gray levels. A customized bin width is often needed.
  .. warning::
    Requires package ``scipy`` and ``trimesh`` to function. If not available, this filter logs a warning and does not
    yield an image.

  References:

  - Banerjee, J, Moelker, A, Niessen, W.J, & van Walsum, T.W. (2013), "3D LBP-based rotationally invariant region
    description." In: Park JI., Kim J. (eds) Computer Vision - ACCV 2012 Workshops. ACCV 2012. Lecture Notes in Computer
    Science, vol 7728. Springer, Berlin, Heidelberg. doi:10.1007/978-3-642-37410-4_3
  """
  global logger
  Nd = inputImage.GetDimension()
  if Nd != 3:
    logger.warning('LBP 3D only available for 3 dimensional images, found %i dimensions', Nd)
    return

  try:
    from scipy.stats import kurtosis
    from scipy.ndimage.interpolation import map_coordinates
    from scipy.special import sph_harm
    from trimesh.creation import icosphere
  except ImportError:
    logger.warning('Could not load required package "scipy" or "trimesh", cannot implement filter LBP 3D')
    return

  # Warn the user if features are extracted in 2D, as this function calculates LBP in 3D
  if kwargs.get('force2D', False):
    logger.warning('Calculating Local Binary Pattern in 3D, but extracting features in 2D. Use with caution!')

  label = kwargs.get('label', 1)

  lbp_levels = kwargs.get('lbp3DLevels', 2)
  lbp_icosphereRadius = kwargs.get('lbp3DIcosphereRadius', 1)
  lbp_icosphereSubdivision = kwargs.get('lbp3DIcosphereSubdivision', 1)

  im_arr = sitk.GetArrayFromImage(inputImage)
  ma_arr = sitk.GetArrayFromImage(inputMask)

  # Variables used in the shape comments:
  # Np Number of voxels
  # Nv Number of vertices

  # Vertices icosahedron for spherical sampling
  coords_icosahedron = numpy.array(icosphere(lbp_icosphereSubdivision, lbp_icosphereRadius).vertices)  # shape(Nv, 3)

  # Corresponding polar coordinates
  theta = numpy.arccos(numpy.true_divide(coords_icosahedron[:, 2], lbp_icosphereRadius))
  phi = numpy.arctan2(coords_icosahedron[:, 1], coords_icosahedron[:, 0])

  # Corresponding spherical harmonics coefficients Y_{m, n, theta, phi}
  Y = sph_harm(0, 0, theta, phi)  # shape(Nv,)
  n_ix = numpy.array(0)

  for n in range(1, lbp_levels):
    for m in range(-n, n + 1):
      n_ix = numpy.append(n_ix, n)
      Y = numpy.column_stack((Y, sph_harm(m, n, theta, phi)))
  # shape (Nv, x) where x is the number of iterations in the above loops + 1

  # Get labelled coordinates
  ROI_coords = numpy.where(ma_arr == label)  # shape(3, Np)

  # Interpolate f (samples on the spheres across the entire volume)
  coords = numpy.array(ROI_coords).T[None, :, :] + coords_icosahedron[:, None, :]  # shape(Nv, Np, 3)
  f = map_coordinates(im_arr, coords.T, order=3)  # Shape(Np, Nv)  Note that 'Np' and 'Nv' are swapped due to .T

  # Compute spherical Kurtosis
  k = kurtosis(f, axis=1)  # shape(Np,)

  # Apply sign function
  f_centroids = im_arr[ROI_coords]  # Shape(Np,)
  f = numpy.greater_equal(f, f_centroids[:, None]).astype(int)  # Shape(Np, Nv)

  # Compute c_{m,n} coefficients
  c = numpy.multiply(f[:, :, None], Y[None, :, :])  # Shape(Np, Nv, x)
  c = c.sum(axis=1)  # Shape(Np, x)

  # Integrate over m
  f = numpy.multiply(c[:, None, n_ix == 0], Y[None, :, n_ix == 0])  # Shape (Np, Nv, 1)
  for n in range(1, lbp_levels):
    f = numpy.concatenate((f,
                           numpy.sum(numpy.multiply(c[:, None, n_ix == n], Y[None, :, n_ix == n]),
                                     axis=2, keepdims=True)
                           ),
                          axis=2)
  # Shape f (Np, Nv, levels)

  # Compute L2-Norm
  f = numpy.sqrt(numpy.sum(f ** 2, axis=1))  # shape(Np, levels)

  # Keep only Real Part
  f = numpy.real(f)  # shape(Np, levels)
  k = numpy.real(k)  # shape(Np,)

  # Yield the derived images for each level
  result = numpy.ndarray(im_arr.shape)
  for l_idx in range(lbp_levels):
    result[ROI_coords] = f[:, l_idx]

    # Create a SimpleITK image
    im = sitk.GetImageFromArray(result)
    im.CopyInformation(inputImage)

    yield im, 'lbp-3D-m%d' % (l_idx + 1), kwargs

  # Yield Kurtosis
  result[ROI_coords] = k

  # Create a SimpleITK image
  im = sitk.GetImageFromArray(result)
  im.CopyInformation(inputImage)

  yield im, 'lbp-3D-k', kwargs
