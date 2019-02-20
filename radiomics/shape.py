import numpy
import SimpleITK as sitk

from radiomics import base, cShape, deprecated


class RadiomicsShape(base.RadiomicsFeaturesBase):
  r"""
  In this group of features we included descriptors of the three-dimensional size and shape of the ROI. These features
  are independent from the gray level intensity distribution in the ROI and are therefore only calculated on the
  non-derived image and mask.

  Unless otherwise specified, features are derived from the approximated shape defined by the triangle mesh. To build
  this mesh, vertices (points) are first defined as points halfway on an edge between a voxel included in the ROI and
  one outside the ROI. By connecting these vertices a mesh of connected triangles is obtained, with each triangle
  defined by 3 adjacent vertices, which shares each side with exactly one other triangle.

  This mesh is generated using a marching cubes algorithm. In this algorithm, a 2x2 cube is moved through the mask
  space. For each position, the corners of the cube are then marked 'segmented' (1) or 'not segmented' (0). Treating the
  corners as specific bits in a binary number, a unique cube-index is obtained (0-255). This index is then used to
  determine which triangles are present in the cube, which are defined in a lookup table.

  These triangles are defined in such a way, that the normal (obtained from the cross product of vectors describing 2
  out of 3 edges) are always oriented in the same direction. For PyRadiomics, the calculated normals are always pointing
  outward. This is necessary to obtain the correct signed volume used in calculation of ``MeshVolume``.

  Let:

  - :math:`N_v` represent the number of voxels included in the ROI
  - :math:`N_f` represent the number of faces (triangles) defining the Mesh.
  - :math:`V` the volume of the mesh in mm\ :sup:`3`, calculated by :py:func:`getMeshVolumeFeatureValue`
  - :math:`A` the surface area of the mesh in mm\ :sup:`2`, calculated by :py:func:`getMeshSurfaceAreaFeatureValue`

  References:

  - Lorensen WE, Cline HE. Marching cubes: A high resolution 3D surface construction algorithm. ACM SIGGRAPH Comput
    Graph `Internet <http://portal.acm.org/citation.cfm?doid=37402.37422>`_. 1987;21:163-9.
  """

  def __init__(self, inputImage, inputMask, **kwargs):
    assert inputMask.GetDimension() == 3, 'Shape features are only available in 3D. If 2D, use shape2D instead'
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

  def getMeshVolumeFeatureValue(self):
    r"""
    **1. Mesh Volume**

    .. math::
      V_i = \displaystyle\frac{Oa_i \cdot (Ob_i \times Oc_i)}{6} \text{ (1)}

      V = \displaystyle\sum^{N_f}_{i=1}{V_i} \text{ (2)}

    The volume of the ROI :math:`V` is calculated from the triangle mesh of the ROI.
    For each face :math:`i` in the mesh, defined by points :math:`a_i, b_i` and :math:`c_i`, the (signed) volume
    :math:`V_f` of the tetrahedron defined by that face and the origin of the image (:math:`O`) is calculated. (1)
    The sign of the volume is determined by the sign of the normal, which must be consistently defined as either facing
    outward or inward of the ROI.

    Then taking the sum of all :math:`V_i`, the total volume of the ROI is obtained (2)

    .. note::
      For more extensive documentation on how the volume is obtained using the surface mesh, see the IBSI document,
      where this feature is defined as ``Volume``.
    """
    return self.Volume

  def getVoxelVolumeFeatureValue(self):
    r"""
    **2. Voxel Volume**

    .. math::
      V_{voxel} = \displaystyle\sum^{N_v}_{k=1}{V_k}

    The volume of the ROI :math:`V_{voxel}` is approximated by multiplying the number of voxels in the ROI by the volume
    of a single voxel :math:`V_k`. This is a less precise approximation of the volume and is not used in subsequent
    features. This feature does not make use of the mesh and is not used in calculation of other shape features.

    .. note::
      Defined in IBSI as ``Approximate Volume``.
    """
    z, y, x = self.pixelSpacing
    Np = len(self.labelledVoxelCoordinates[0])
    return Np * (z * x * y)

  def getSurfaceAreaFeatureValue(self):
    r"""
    **3. Surface Area**

    .. math::
      A_i = \frac{1}{2}|\text{a}_i\text{b}_i \times \text{a}_i\text{c}_i| \text{ (1)}

      A = \displaystyle\sum^{N_f}_{i=1}{A_i} \text{ (2)}

    where:

    :math:`\text{a}_i\text{b}_i` and :math:`\text{a}_i\text{c}_i` are edges of the :math:`i^{\text{th}}` triangle in the
    mesh, formed by vertices :math:`\text{a}_i`, :math:`\text{b}_i` and :math:`\text{c}_i`.

    To calculate the surface area, first the surface area :math:`A_i` of each triangle in the mesh is calculated (1).
    The total surface area is then obtained by taking the sum of all calculated sub-areas (2).

    .. note::
      Defined in IBSI as ``Surface Area``.
    """
    return self.SurfaceArea

  def getSurfaceVolumeRatioFeatureValue(self):
    r"""
    **4. Surface Area to Volume ratio**

    .. math::
      \textit{surface to volume ratio} = \frac{A}{V}

    Here, a lower value indicates a more compact (sphere-like) shape. This feature is not dimensionless, and is
    therefore (partly) dependent on the volume of the ROI.
    """
    return self.SurfaceArea / self.Volume

  def getSphericityFeatureValue(self):
    r"""
    **5. Sphericity**

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
    **6. Compactness 1**

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
      specified, including this feature). To include this feature in the extraction, specify it by name in the enabled
      features.
    """
    return self.Volume / (self.SurfaceArea ** (3.0 / 2.0) * numpy.sqrt(numpy.pi))

  @deprecated
  def getCompactness2FeatureValue(self):
    r"""
    **7. Compactness 2**

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
      specified, including this feature). To include this feature in the extraction, specify it by name in the enabled
      features.
    """
    return (36.0 * numpy.pi) * (self.Volume ** 2.0) / (self.SurfaceArea ** 3.0)

  @deprecated
  def getSphericalDisproportionFeatureValue(self):
    r"""
    **8. Spherical Disproportion**

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
      specified, including this feature). To include this feature in the extraction, specify it by name in the enabled
      features.
    """
    return self.SurfaceArea / (36 * numpy.pi * self.Volume ** 2) ** (1.0 / 3.0)

  def getMaximum3DDiameterFeatureValue(self):
    r"""
    **9. Maximum 3D diameter**

    Maximum 3D diameter is defined as the largest pairwise Euclidean distance between tumor surface mesh
    vertices.

    Also known as Feret Diameter.
    """
    return self.diameters[3]

  def getMaximum2DDiameterSliceFeatureValue(self):
    r"""
    **10. Maximum 2D diameter (Slice)**

    Maximum 2D diameter (Slice) is defined as the largest pairwise Euclidean distance between tumor surface mesh
    vertices in the row-column (generally the axial) plane.
    """
    return self.diameters[0]

  def getMaximum2DDiameterColumnFeatureValue(self):
    r"""
    **11. Maximum 2D diameter (Column)**

    Maximum 2D diameter (Column) is defined as the largest pairwise Euclidean distance between tumor surface mesh
    vertices in the row-slice (usually the coronal) plane.
    """
    return self.diameters[1]

  def getMaximum2DDiameterRowFeatureValue(self):
    r"""
    **12. Maximum 2D diameter (Row)**

    Maximum 2D diameter (Row) is defined as the largest pairwise Euclidean distance between tumor surface mesh
    vertices in the column-slice (usually the sagittal) plane.
    """
    return self.diameters[2]

  def getMajorAxisLengthFeatureValue(self):
    r"""
    **13. Major Axis Length**

    .. math::
      \textit{major axis} = 4 \sqrt{\lambda_{major}}

    This feature yield the largest axis length of the ROI-enclosing ellipsoid and is calculated using the largest
    principal component :math:`\lambda_{major}`.

    The principal component analysis is performed using the physical coordinates of the voxel centers defining the ROI.
    It therefore takes spacing into account, but does not make use of the shape mesh.
    """
    if self.eigenValues[2] < 0:
      self.logger.warning('Major axis eigenvalue negative! (%g)', self.eigenValues[2])
      return numpy.nan
    return numpy.sqrt(self.eigenValues[2]) * 4

  def getMinorAxisLengthFeatureValue(self):
    r"""
    **14. Minor Axis Length**

    .. math::
      \textit{minor axis} = 4 \sqrt{\lambda_{minor}}

    This feature yield the second-largest axis length of the ROI-enclosing ellipsoid and is calculated using the largest
    principal component :math:`\lambda_{minor}`.

    The principal component analysis is performed using the physical coordinates of the voxel centers defining the ROI.
    It therefore takes spacing into account, but does not make use of the shape mesh.
    """
    if self.eigenValues[1] < 0:
      self.logger.warning('Minor axis eigenvalue negative! (%g)', self.eigenValues[1])
      return numpy.nan
    return numpy.sqrt(self.eigenValues[1]) * 4

  def getLeastAxisLengthFeatureValue(self):
    r"""
    **15. Least Axis Length**

    .. math::
      \textit{least axis} = 4 \sqrt{\lambda_{least}}

    This feature yield the smallest axis length of the ROI-enclosing ellipsoid and is calculated using the largest
    principal component :math:`\lambda_{least}`. In case of a 2D segmentation, this value will be 0.

    The principal component analysis is performed using the physical coordinates of the voxel centers defining the ROI.
    It therefore takes spacing into account, but does not make use of the shape mesh.
    """
    if self.eigenValues[0] < 0:
      self.logger.warning('Least axis eigenvalue negative! (%g)', self.eigenValues[0])
      return numpy.nan
    return numpy.sqrt(self.eigenValues[0]) * 4

  def getElongationFeatureValue(self):
    r"""
    **16. Elongation**

    Elongation shows the relationship between the two largest principal components in the ROI shape.
    For computational reasons, this feature is defined as the inverse of true elongation.

    .. math::
      \textit{elongation} = \sqrt{\frac{\lambda_{minor}}{\lambda_{major}}}

    Here, :math:`\lambda_{\text{major}}` and :math:`\lambda_{\text{minor}}` are the lengths of the largest and second
    largest principal component axes. The values range between 1 (where the cross section through the first and second
    largest principal moments is circle-like (non-elongated)) and 0 (where the object is a maximally elongated: i.e. a 1
    dimensional line).

    The principal component analysis is performed using the physical coordinates of the voxel centers defining the ROI.
    It therefore takes spacing into account, but does not make use of the shape mesh.
    """
    if self.eigenValues[1] < 0 or self.eigenValues[2] < 0:
      self.logger.warning('Elongation eigenvalue negative! (%g, %g)', self.eigenValues[1], self.eigenValues[2])
      return numpy.nan
    return numpy.sqrt(self.eigenValues[1] / self.eigenValues[2])

  def getFlatnessFeatureValue(self):
    r"""
    **17. Flatness**

    Flatness shows the relationship between the largest and smallest principal components in the ROI shape.
    For computational reasons, this feature is defined as the inverse of true flatness.

    .. math::
      \textit{flatness} = \sqrt{\frac{\lambda_{least}}{\lambda_{major}}}

    Here, :math:`\lambda_{\text{major}}` and :math:`\lambda_{\text{least}}` are the lengths of the largest and smallest
    principal component axes. The values range between 1 (non-flat, sphere-like) and 0 (a flat object, or single-slice
    segmentation).

    The principal component analysis is performed using the physical coordinates of the voxel centers defining the ROI.
    It therefore takes spacing into account, but does not make use of the shape mesh.
    """
    if self.eigenValues[0] < 0 or self.eigenValues[2] < 0:
      self.logger.warning('Elongation eigenvalue negative! (%g, %g)', self.eigenValues[0], self.eigenValues[2])
      return numpy.nan
    return numpy.sqrt(self.eigenValues[0] / self.eigenValues[2])
