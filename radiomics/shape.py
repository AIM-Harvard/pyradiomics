import numpy
import SimpleITK as sitk

from radiomics import base, cShape, deprecated


class RadiomicsShape(base.RadiomicsFeaturesBase):
  r"""
  In this group of features we included descriptors of the three-dimensional size and shape of the ROI. These features
  are independent from the gray level intensity distribution in the ROI and are therefore only calculated on the
  non-derived image and mask.

  Let:

  - :math:`V` the volume of the ROI in mm\ :sup:`3`
  - :math:`A` the surface area of the ROI in mm\ :sup:`2`
  """

  def __init__(self, inputImage, inputMask, **kwargs):
    super(RadiomicsShape, self).__init__(inputImage, inputMask, **kwargs)

  def _initVoxelBasedCalculation(self):
    raise NotImplementedError('Shape features are not available in voxel-based mode')

  def _initSegmentBasedCalculation(self):

    self.pixelSpacing = numpy.array(self.inputImage.GetSpacing()[::-1])

    # Pad inputMask to prevent index-out-of-range errors
    self.logger.debug('Padding the mask with 0s')

    cpif = sitk.ConstantPadImageFilter()

    padding = numpy.tile(1, 3)
    try:
      cpif.SetPadLowerBound(padding)
      cpif.SetPadUpperBound(padding)
    except TypeError:
      # newer versions of SITK/python want a tuple or list
      cpif.SetPadLowerBound(padding.tolist())
      cpif.SetPadUpperBound(padding.tolist())

    self.inputMask = cpif.Execute(self.inputMask)

    # Reassign self.maskArray using the now-padded self.inputMask
    self.maskArray = (sitk.GetArrayFromImage(self.inputMask) == self.label)
    self.labelledVoxelCoordinates = numpy.where(self.maskArray != 0)

    self.logger.debug('Pre-calculate Volume, Surface Area and Eigenvalues')

    # Volume, Surface Area and eigenvalues are pre-calculated

    # Compute Surface Area and volume
    self.SurfaceArea, self.Volume, self.diameters = cShape.calculate_coefficients(self.maskArray, self.pixelSpacing)

    # Compute eigenvalues and -vectors
    Np = len(self.labelledVoxelCoordinates[0])
    coordinates = numpy.array(self.labelledVoxelCoordinates, dtype='int').transpose((1, 0))  # Transpose equals zip(*a)
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

  def getVolumeFeatureValue(self):
    r"""
    **1. Volume**

    .. math::
      V_{f} = \displaystyle\frac{a \dot (b \times c)}{6} (1)

      V = \displaystyle\sum^{N_f}_{f=1}{V_f} (2)

    The volume of the ROI :math:`V` is calculated from the triangle mesh of the ROI. First the volume is calculated of
    the tetrahedron that is defined by face f and the origin of the image (1) for all :math:`N_f` faces that make up the
    triangle mesh. This mesh is also used for the calculation of surface area. By ensuring that the normals (obtained
    from the cross product) have a consistent orientation (outward or inward facing), the volumes obtained are signed.
    All sub-volumes are then summed (2), yielding the volume of the ROI.

    .. note::
      For more extensive documentation on how the volume is obtained using the surface mesh, see the IBSI document on
      the calculation of Volume.
    """
    return self.Volume

  def getApproximateVolumeFeatureValue(self):
    r"""
    **1. Volume**

    .. math::
      V_{approx} = \displaystyle\sum^{N}_{i=1}{V_i}

    The volume of the ROI :math:`V` is approximated by multiplying the number of voxels in the ROI by the volume of a
    single voxel :math:`V_i`. This is a less precise approximation of the volume and is not used in subsequent features.
    """
    z, y, x = self.pixelSpacing
    Np = len(self.labelledVoxelCoordinates[0])
    return Np * (z * x * y)

  def getSurfaceAreaFeatureValue(self):
    r"""
    **2. Surface Area**

    .. math::
      A = \displaystyle\sum^{N}_{i=1}{\frac{1}{2}|\text{a}_i\text{b}_i \times \text{a}_i\text{c}_i|}

    Where:

    :math:`N` is the number of triangles forming the surface mesh of the volume (ROI)

    :math:`\text{a}_i\text{b}_i` and :math:`\text{a}_i\text{c}_i` are the edges of the :math:`i^{\text{th}}` triangle
    formed by points :math:`\text{a}_i`, :math:`\text{b}_i` and :math:`\text{c}_i`

    Surface Area is an approximation of the surface of the ROI in mm2, calculated using a the triangle mesh generated
    using the marching cubes algorithm.

    References:

    - Lorensen WE, Cline HE. Marching cubes: A high resolution 3D surface construction algorithm. ACM SIGGRAPH Comput
      Graph `Internet <http://portal.acm.org/citation.cfm?doid=37402.37422>`_. 1987;21:163-9.
    """
    return self.SurfaceArea

  def getSurfaceVolumeRatioFeatureValue(self):
    r"""
    **3. Surface Area to Volume ratio**

    .. math::
      \textit{surface to volume ratio} = \frac{A}{V}

    Here, a lower value indicates a more compact (sphere-like) shape. This feature is not dimensionless, and is
    therefore (partly) dependent on the volume of the ROI.
    """
    return self.SurfaceArea / self.Volume

  def getSphericityFeatureValue(self):
    r"""
    **4. Sphericity**

    .. math::
      \textit{sphericity} = \frac{\sqrt[3]{36 \pi V^2}}{A}

    Sphericity is a measure of the roundness of the shape of the tumor region relative to a sphere. It is a
    dimensionless measure, independent of scale and orientation. The value range is :math:`0 < sphericity \leq 1`, where
    a value of 1 indicates a perfect sphere (a sphere has the smallest possible surface area for a given volume,
    compared to other solids).

    .. note::
      This feature is correlated to Compactness 1, Compactness 2 and Spherical Disproportion. In the default
      parameter file provided in the ``pyradiomics/examples/exampleSettings`` folder, Compactness 1 and Compactness 2
      are therefore disabled.
    """
    return (36 * numpy.pi * self.Volume ** 2) ** (1.0 / 3.0) / self.SurfaceArea

  @deprecated
  def getCompactness1FeatureValue(self):
    r"""
    **5. Compactness 1**

    .. math::
      \textit{compactness 1} = \frac{V}{\sqrt{\pi A^3}}

    Similar to Sphericity, Compactness 1 is a measure of how compact the shape of the tumor is relative to a sphere
    (most compact). It is therefore correlated to Sphericity and redundant. It is provided here for completeness.
    The value range is :math:`0 < compactness\ 1 \leq \frac{1}{6 \pi}`, where a value of :math:`\frac{1}{6 \pi}`
    indicates a perfect sphere.

    By definition, :math:`compactness\ 1 = \frac{1}{6 \pi}\sqrt{compactness\ 2} =
    \frac{1}{6 \pi}\sqrt{sphericity^3}`.

    .. note::
      This feature is correlated to Compactness 2, Sphericity and Spherical Disproportion.
      Therefore, this feature is marked, so it is not enabled by default (i.e. this feature will not be enabled if no
      individual features are specified (enabling 'all' features), but will be enabled when individual features are
      specified, including this feature). To include this feature in the extraction, specify it by name in the enabled features.
    """
    return self.Volume / (self.SurfaceArea ** (3.0 / 2.0) * numpy.sqrt(numpy.pi))

  @deprecated
  def getCompactness2FeatureValue(self):
    r"""
    **6. Compactness 2**

    .. math::
      \textit{compactness 2} = 36 \pi \frac{V^2}{A^3}

    Similar to Sphericity and Compactness 1, Compactness 2 is a measure of how compact the shape of the tumor is
    relative to a sphere (most compact). It is a dimensionless measure, independent of scale and orientation. The value
    range is :math:`0 < compactness\ 2 \leq 1`, where a value of 1 indicates a perfect sphere.

    By definition, :math:`compactness\ 2 = (sphericity)^3`

    .. note::
      This feature is correlated to Compactness 1, Sphericity and Spherical Disproportion.
      Therefore, this feature is marked, so it is not enabled by default (i.e. this feature will not be enabled if no
      individual features are specified (enabling 'all' features), but will be enabled when individual features are
      specified, including this feature). To include this feature in the extraction, specify it by name in the enabled features.
    """
    return (36.0 * numpy.pi) * (self.Volume ** 2.0) / (self.SurfaceArea ** 3.0)

  @deprecated
  def getSphericalDisproportionFeatureValue(self):
    r"""
    **7. Spherical Disproportion**

    .. math::
      \textit{spherical disproportion} = \frac{A}{4\pi R^2} = \frac{A}{\sqrt[3]{36 \pi V^2}}

    Where :math:`R` is the radius of a sphere with the same volume as the tumor, and equal to
    :math:`\sqrt[3]{\frac{3V}{4\pi}}`.

    Spherical Disproportion is the ratio of the surface area of the tumor region to the surface area of a sphere with
    the same volume as the tumor region, and by definition, the inverse of Sphericity. Therefore, the value range is
    :math:`spherical\ disproportion \geq 1`, with a value of 1 indicating a perfect sphere.

    .. note::
      This feature is correlated to Compactness 2, Compactness2 and Sphericity.
      Therefore, this feature is marked, so it is not enabled by default (i.e. this feature will not be enabled if no
      individual features are specified (enabling 'all' features), but will be enabled when individual features are
      specified, including this feature). To include this feature in the extraction, specify it by name in the enabled features.
    """
    return self.SurfaceArea / (36 * numpy.pi * self.Volume ** 2) ** (1.0 / 3.0)

  def getMaximum3DDiameterFeatureValue(self):
    r"""
    **8. Maximum 3D diameter**

    Maximum 3D diameter is defined as the largest pairwise Euclidean distance between tumor surface mesh
    vertices.

    Also known as Feret Diameter.
    """
    return self.diameters[3]

  def getMaximum2DDiameterSliceFeatureValue(self):
    r"""
    **9. Maximum 2D diameter (Slice)**

    Maximum 2D diameter (Slice) is defined as the largest pairwise Euclidean distance between tumor surface mesh
    vertices in the row-column (generally the axial) plane.
    """
    return self.diameters[0]

  def getMaximum2DDiameterColumnFeatureValue(self):
    r"""
    **10. Maximum 2D diameter (Column)**

    Maximum 2D diameter (Column) is defined as the largest pairwise Euclidean distance between tumor surface mesh
    vertices in the row-slice (usually the coronal) plane.
    """
    return self.diameters[1]

  def getMaximum2DDiameterRowFeatureValue(self):
    r"""
    **11. Maximum 2D diameter (Row)**

    Maximum 2D diameter (Row) is defined as the largest pairwise Euclidean distance between tumor surface mesh
    vertices in the column-slice (usually the sagittal) plane.
    """
    return self.diameters[2]

  def getMajorAxisFeatureValue(self):
    r"""
    **12. Major Axis**

    .. math::
      \textit{major axis} = 4 \sqrt{\lambda_{\text{major}}}

    """
    if self.eigenValues[2] < 0:
      self.logger.warning('Major axis eigenvalue negative! (%g)', self.eigenValues[2])
      return numpy.nan
    return numpy.sqrt(self.eigenValues[2]) * 4

  def getMinorAxisFeatureValue(self):
    r"""
    **13. Minor Axis**

    .. math::
      \textit{minor axis} = 4 \sqrt{\lambda_{\text{minor}}}

    """
    if self.eigenValues[1] < 0:
      self.logger.warning('Minor axis eigenvalue negative! (%g)', self.eigenValues[1])
      return numpy.nan
    return numpy.sqrt(self.eigenValues[1]) * 4

  def getLeastAxisFeatureValue(self):
    r"""
    **14. Least Axis**

    .. math::
      \textit{least axis} = 4 \sqrt{\lambda_{\text{least}}}

    """
    if self.eigenValues[0] < 0:
      self.logger.warning('Least axis eigenvalue negative! (%g)', self.eigenValues[0])
      return numpy.nan
    return numpy.sqrt(self.eigenValues[0]) * 4

  def getElongationFeatureValue(self):
    r"""
    **15. Elongation**

    Elongation is calculated using its implementation in SimpleITK, and is defined as:

    .. math::
      \textit{elongation} = \sqrt{\frac{\lambda_{\text{minor}}}{\lambda_{\text{major}}}}

    Here, :math:`\lambda_{\text{major}}` and :math:`\lambda_{\text{minor}}` are the lengths of the largest and second
    largest principal component axes. The values range between 1 (where the cross section through the first and second
    largest principal moments is circle-like (non-elongated)) and 0 (where the object is a single point or 1 dimensional
    line).
    """
    if self.eigenValues[1] < 0 or self.eigenValues[2] < 0:
      self.logger.warning('Elongation eigenvalue negative! (%g, %g)', self.eigenValues[1], self.eigenValues[2])
      return numpy.nan
    return numpy.sqrt(self.eigenValues[1] / self.eigenValues[2])

  def getFlatnessFeatureValue(self):
    r"""
    **16. Flatness**

    Flatness is calculated using its implementation in SimpleITK, and is defined as:

    .. math::
      \textit{flatness} = \sqrt{\frac{\lambda_{\text{least}}}{\lambda_{\text{major}}}}

    Here, :math:`\lambda_{\text{major}}` and :math:`\lambda_{\text{least}}` are the lengths of the largest and smallest
    principal component axes. The values range between 1 (non-flat, sphere-like) and 0 (a flat object).
    """
    if self.eigenValues[0] < 0 or self.eigenValues[2] < 0:
      self.logger.warning('Elongation eigenvalue negative! (%g, %g)', self.eigenValues[0], self.eigenValues[2])
      return numpy.nan
    return numpy.sqrt(self.eigenValues[0] / self.eigenValues[2])
