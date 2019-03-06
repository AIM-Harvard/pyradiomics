import numpy
import SimpleITK as sitk

from radiomics import base, cShape, deprecated


class RadiomicsShape2D(base.RadiomicsFeaturesBase):
  r"""
  In this group of features we included descriptors of the two-dimensional size and shape of the ROI. These features
  are independent from the gray level intensity distribution in the ROI and are therefore only calculated on the
  non-derived image and mask.

  Unless otherwise specified, features are derived from the approximated shape defined by the circumference mesh. To
  build this mesh, vertices (points) are first defined as points halfway on an edge between a pixel included in the ROI
  and one outside the ROI. By connecting these vertices a mesh of connected lines is obtained, with each line
  defined by 2 adjacent vertices, which shares each a point with exactly one other line.

  This mesh is generated using an adapted version marching cubes algorithm. In this algorithm, a 2x2 square is moved
  through the mask space (2d). For each position, the corners of the square are then marked 'segmented' (1) or
  'not segmented' (0). Treating the corners as specific bits in a binary number, a unique square-index is obtained
  (0-15). This index is then used to determine which lines are present in the square, which are defined in a lookup
  table.

  These lines are defined in such a way, that the normal of the triangle defined by these points and the origin
  is always oriented in the a consistent direction. This results in signed values for the surface area of each triangle,
  so that when summed, the superfluous (postive) area included by triangles partly inside and outside the ROI is
  perfectly cancelled out by the (negative) area of triangles entirely outside the ROI.

  Let:

  - :math:`N_p` represent the number of pixels included in the ROI
  - :math:`N_f` represent the number of lines defining the circumference (perimeter) Mesh.
  - :math:`A` the surface area of the mesh in mm\ :sup:`2`, calculated by :py:func:`getMeshSurfaceFeatureValue`
  - :math:`P` the perimeter of the mesh in mm, calculated by :py:func:`getPerimeterFeatureValue`

  .. note::
    This class can **only** be calculated for truly 2D masks. To ensure correct processing, it is required that
    ``force2D`` is set to ``True`` and ``force2Ddimension`` to the dimension that is out-of plane (e.g. 0 (z-axis) for
    an axial slice). Furthermore, this dimension is required to have size 1. If not set correctly, a ValueError is
    raised.

  References:

  - Lorensen WE, Cline HE. Marching cubes: A high resolution 3D surface construction algorithm. ACM SIGGRAPH Comput
    Graph `Internet <http://portal.acm.org/citation.cfm?doid=37402.37422>`_. 1987;21:163-9.
  """

  def __init__(self, inputImage, inputMask, **kwargs):
    super(RadiomicsShape2D, self).__init__(inputImage, inputMask, **kwargs)

  def _initVoxelBasedCalculation(self):
    raise NotImplementedError('Shape features are not available in pixel-based mode')

  def _initSegmentBasedCalculation(self):
    self.maskArray = (sitk.GetArrayFromImage(self.inputMask) == self.label)  # boolean array

    Nd = self.inputMask.GetDimension()
    if Nd == 3:
      if not self.settings.get('force2D', False):
        raise ValueError('Shape2D is can only be calculated when input is 2D or 3D with `force2D=True`')

      force2DDimension = self.settings.get('force2Ddimension', 0)
      axes = [0, 1, 2]
      axes.remove(force2DDimension)

      self.pixelSpacing = numpy.array(self.inputImage.GetSpacing()[::-1])[(axes,)]

      if self.maskArray.shape[force2DDimension] > 1:
        raise ValueError('Size of the mask in dimension %i is more than 1, cannot compute 2D shape')

      # Drop the 2D axis, ensuring the input is truly 2D
      self.maskArray = numpy.squeeze(self.maskArray, axis=force2DDimension)
    elif Nd == 2:
      self.pixelSpacing = numpy.array(self.inputImage.GetSpacing()[::-1])
    else:
      raise ValueError('Shape2D is can only be calculated when input is 2D or 3D with `force2D=True`')

    # Pad maskArray to prevent index-out-of-range errors
    self.logger.debug('Padding the mask with 0s')
    self.maskArray = numpy.pad(self.maskArray, pad_width=1, mode='constant', constant_values=0)

    self.labelledPixelCoordinates = numpy.where(self.maskArray != 0)

    self.logger.debug('Pre-calculate surface, perimeter, diameter and eigenvalues')

    # Volume, Surface Area and eigenvalues are pre-calculated

    # Compute Surface Area and volume
    self.Perimeter, self.Surface, self.Diameter = cShape.calculate_coefficients2D(self.maskArray, self.pixelSpacing)

    # Compute eigenvalues and -vectors
    Np = len(self.labelledPixelCoordinates[0])
    coordinates = numpy.array(self.labelledPixelCoordinates, dtype='int').transpose((1, 0))  # Transpose equals zip(*a)
    physicalCoordinates = coordinates * self.pixelSpacing[None, :]
    physicalCoordinates -= numpy.mean(physicalCoordinates, axis=0)  # Centered at 0
    physicalCoordinates /= numpy.sqrt(Np)
    covariance = numpy.dot(physicalCoordinates.T.copy(), physicalCoordinates)
    self.eigenValues = numpy.linalg.eigvals(covariance)

    # Correct machine precision errors causing very small negative eigen values in case of some 2D segmentations
    machine_errors = numpy.bitwise_and(self.eigenValues < 0, self.eigenValues > -1e-10)
    if numpy.sum(machine_errors) > 0:
      self.logger.warning('Encountered %d eigenvalues < 0 and > -1e-10, rounding to 0', numpy.sum(machine_errors))
      self.eigenValues[machine_errors] = 0

    self.eigenValues.sort()  # Sort the eigenValues from small to large

    self.logger.debug('Shape feature class initialized')

  def getMeshSurfaceFeatureValue(self):
    r"""
    **1. Mesh Surface**

    .. math::
      A_i = \frac{1}{2}\text{Oa}_i \times \text{Ob}_i \text{ (1)}

      A = \displaystyle\sum^{N_f}_{i=1}{A_i} \text{ (2)}

    where:

    :math:`\text{O}_i\text{a}_i` and :math:`\text{O}_i\text{b}_i` are edges of the :math:`i^{\text{th}}` triangle in the
    mesh, formed by vertices :math:`\text{a}_i`, :math:`\text{b}_i` of the perimiter and the origin :math:`\text{O}`.

    To calculate the surface area, first the signed surface area :math:`A_i` of each triangle in the mesh is calculated
    (1). The total surface area is then obtained by taking the sum of all calculated sub-areas (2), where the sign will
    ensure correct surface area, as the negative area of triangles outside the ROI will cancel out the surplus area
    included by triangles partly inside and partly outside the ROI.
    """
    return self.Surface

  def getPixelSurfaceFeatureValue(self):
    r"""
    **2. Pixel Surface**

    .. math::
      A_{pixel} = \displaystyle\sum^{N_v}_{k=1}{A_k}

    The surface area of the ROI :math:`A_{pixel}` is approximated by multiplying the number of pixels in the ROI by the
    surface area of a single pixel :math:`A_k`. This is a less precise approximation of the surface area.
    This feature does not make use of the mesh and is not used in calculation of other 2D shape features.
    """
    y, x = self.pixelSpacing
    Np = len(self.labelledPixelCoordinates[0])
    return Np * (x * y)

  def getPerimeterFeatureValue(self):
    r"""
    **3. Perimeter**

    .. math::
      P_i = \sqrt{(\text{a}_i-\text{b}_i)^2} \text{ (1)}

      P = \displaystyle\sum^{N_f}_{i=1}{P_i} \text{ (2)}

    where:

    :math:`\text{a}_i` and :math:`\text{b}_i` are vertices of the :math:`i^{\text{th}}` line in the
    perimeter mesh.

    To calculate the perimeter, first the perimeter :math:`A_i` of each line in the mesh circumference is calculated
    (1). The total perimeter is then obtained by taking the sum of all calculated sub-areas (2).

    """
    return self.Perimeter

  def getPerimeterSurfaceRatioFeatureValue(self):
    r"""
    **4. Perimeter to Surface ratio**

    .. math::
      \textit{perimeter to surface ratio} = \frac{P}{A}

    Here, a lower value indicates a more compact (circle-like) shape. This feature is not dimensionless, and is
    therefore (partly) dependent on the surface area of the ROI.
    """
    return self.Perimeter / self.Surface

  def getSphericityFeatureValue(self):
    r"""
    **5. Sphericity**

    .. math::
      \textit{sphericity} = \frac{2\pi R}{P} = \frac{2\sqrt{\pi A}}{P}

    Where :math:`R` is the radius of a circle with the same surface as the ROI, and equal to
    :math:`\sqrt{\frac{A}{\pi}}`.

    Sphericity is the ratio of the perimeter of the tumor region to the perimeter of a circle with
    the same surface area as the tumor region and therefore a measure of the roundness of the shape of the tumor region
    relative to a circle. It is a dimensionless measure, independent of scale and orientation. The value range is
    :math:`0 < sphericity \leq 1`, where a value of 1 indicates a perfect circle (a circle has the smallest possible
    perimeter for a given surface area, compared to other shapes).

    .. note::
      This feature is correlated to Spherical Disproportion. Therefore, only this feature is enabled by default.
    """
    return (2 * numpy.sqrt(numpy.pi * self.Surface)) / self.Perimeter

  @deprecated
  def getSphericalDisproportionFeatureValue(self):
    r"""
    **6. Spherical Disproportion**

    .. math::
      \textit{spherical disproportion} = \frac{P}{2\sqrt{\pi A}}

    Spherical Disproportion is the ratio of the perimeter of the tumor region to the perimeter of a circle with
    the same surface area as the tumor region, and by definition, the inverse of Sphericity. Therefore, the value range
    is :math:`spherical\ disproportion \geq 1`, with a value of 1 indicating a perfect sphere.

    .. note::
      This feature is correlated to Sphericity.
      Therefore, this feature is marked, so it is not enabled by default (i.e. this feature will not be enabled if no
      individual features are specified (enabling 'all' features), but will be enabled when individual features are
      specified, including this feature). To include this feature in the extraction, specify it by name in the enabled
      features.
    """
    return 1.0 / self.getSphericityFeatureValue()

  def getMaximumDiameterFeatureValue(self):
    r"""
    **7. Maximum 2D diameter**

    Maximum diameter is defined as the largest pairwise Euclidean distance between tumor surface mesh
    vertices.
    """
    return self.Diameter

  def getMajorAxisLengthFeatureValue(self):
    r"""
    **8. Major Axis Length**

    .. math::
      \textit{major axis} = 4 \sqrt{\lambda_{major}}

    This feature yield the largest axis length of the ROI-enclosing ellipsoid and is calculated using the largest
    principal component :math:`\lambda_{major}`.

    The principal component analysis is performed using the physical coordinates of the pixel centers defining the ROI.
    It therefore takes spacing into account, but does not make use of the shape mesh.
    """
    if self.eigenValues[1] < 0:
      self.logger.warning('Major axis eigenvalue negative! (%g)', self.eigenValues[1])
      return numpy.nan
    return numpy.sqrt(self.eigenValues[1]) * 4

  def getMinorAxisLengthFeatureValue(self):
    r"""
    **9. Minor Axis Length**

    .. math::
      \textit{minor axis} = 4 \sqrt{\lambda_{minor}}

    This feature yield the second-largest axis length of the ROI-enclosing ellipsoid and is calculated using the largest
    principal component :math:`\lambda_{minor}`.

    The principal component analysis is performed using the physical coordinates of the pixel centers defining the ROI.
    It therefore takes spacing into account, but does not make use of the shape mesh.
    """
    if self.eigenValues[0] < 0:
      self.logger.warning('Minor axis eigenvalue negative! (%g)', self.eigenValues[0])
      return numpy.nan
    return numpy.sqrt(self.eigenValues[0]) * 4

  def getElongationFeatureValue(self):
    r"""
    **10. Elongation**

    Elongation shows the relationship between the two largest principal components in the ROI shape.
    For computational reasons, this feature is defined as the inverse of true elongation.

    .. math::
      \textit{elongation} = \sqrt{\frac{\lambda_{minor}}{\lambda_{major}}}

    Here, :math:`\lambda_{\text{major}}` and :math:`\lambda_{\text{minor}}` are the lengths of the largest and second
    largest principal component axes. The values range between 1 (where the cross section through the first and second
    largest principal moments is circle-like (non-elongated)) and 0 (where the object is a maximally elongated: i.e. a 1
    dimensional line).

    The principal component analysis is performed using the physical coordinates of the pixel centers defining the ROI.
    It therefore takes spacing into account, but does not make use of the shape mesh.
    """
    if self.eigenValues[0] < 0 or self.eigenValues[1] < 0:
      self.logger.warning('Elongation eigenvalue negative! (%g, %g)', self.eigenValues[0], self.eigenValues[1])
      return numpy.nan
    return numpy.sqrt(self.eigenValues[0] / self.eigenValues[1])
