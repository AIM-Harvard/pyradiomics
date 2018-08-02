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

    # Reassign self.maskArray using the now-padded self.inputMask and force and interger datatype
    self.maskArray = (sitk.GetArrayFromImage(self.inputMask) == self.label).astype('int')
    self.labelledVoxelCoordinates = numpy.where(self.maskArray != 0)

    self.logger.debug('Pre-calculate Volume, Surface Area and Eigenvalues')

    # Volume, Surface Area and eigenvalues are pre-calculated
    # Compute volume
    z, y, x = self.pixelSpacing
    Np = len(self.labelledVoxelCoordinates[0])
    self.Volume = Np * (z * x * y)

    # Compute Surface Area
    self.SurfaceArea = self._calculateSurfaceArea()

    # Compute eigenvalues and -vectors
    coordinates = numpy.array(self.labelledVoxelCoordinates, dtype='int').transpose((1, 0))  # Transpose equals zip(*a)
    physicalCoordinates = coordinates * self.pixelSpacing[None, :]
    physicalCoordinates -= numpy.mean(physicalCoordinates, axis=0)  # Centered at 0
    physicalCoordinates /= numpy.sqrt(Np)
    covariance = numpy.dot(physicalCoordinates.T.copy(), physicalCoordinates)
    self.eigenValues, eigenVectors = numpy.linalg.eig(covariance)  # eigenVectors are not used

    # Correct machine precision errors causing very small negative eigen values in case of some 2D segmentations
    machine_errors = numpy.bitwise_and(self.eigenValues < 0, self.eigenValues > -1e-10)
    if numpy.sum(machine_errors) > 0:
      self.logger.warning('Encountered %d eigenvalues < 0 and > -1e-10, rounding to 0', numpy.sum(machine_errors))
      self.eigenValues[machine_errors] = 0

    self.eigenValues.sort()  # Sort the eigenValues from small to large

    self.diameters = None  # Do not precompute diameters, but instantiate the variable for lazy assignment

    self.logger.debug('Shape feature class initialized')

  def _calculateSurfaceArea(self):
    self.logger.debug('Calculating Surface Area in C')

    return cShape.calculate_surfacearea(self.maskArray, self.pixelSpacing)

  def calculateDiameters(self):
    """
    Calculate maximum diameters in 2D and 3D using C extension. Function returns a tuple with 4 elements:

    0. Maximum 2D diameter Slice (XY Plane, Axial)
    1. Maximum 2D diameter Column (ZX Plane, Coronal)
    2. Maximum 2D diameter Row (ZY Plane, Sagittal)
    3. Maximum 3D diameter
    """
    self.logger.debug('Calculating Maximum diameters in C')
    Ns = len(self.labelledVoxelCoordinates[0])
    return cShape.calculate_diameter(self.maskArray, self.pixelSpacing, Ns)

  def getVolumeFeatureValue(self):
    r"""
    **1. Volume**

    .. math::
      V = \displaystyle\sum^{N}_{i=1}{V_i}

    The volume of the ROI :math:`V` is approximated by multiplying the number of voxels in the ROI by the volume of a
    single voxel :math:`V_i`.

    .. note::
      In the IBSI feature definitions, a more precise approximation of the volume is used. That method uses tetrahedrons
      consisting of the origin and faces in the ROI. Although the method implemented here overestimates the volume,
      especially in small volumes, the difference will be negligible in large ROIs.
    """
    return self.Volume

  def getSurfaceAreaFeatureValue(self):
    r"""
    **2. Surface Area**

    .. math::
      A = \displaystyle\sum^{N}_{i=1}{\frac{1}{2}|\text{a}_i\text{b}_i \times \text{a}_i\text{c}_i|}

    Where:

    :math:`N` is the number of triangles forming the surface mesh of the volume (ROI)

    :math:`\text{a}_i\text{b}_i` and :math:`\text{a}_i\text{c}_i` are the edges of the :math:`i^{\text{th}}` triangle
    formed by points :math:`\text{a}_i`, :math:`\text{b}_i` and :math:`\text{c}_i`

    Surface Area is an approximation of the surface of the ROI in mm2, calculated using a marching cubes algorithm.

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

    Maximum 3D diameter is defined as the largest pairwise Euclidean distance between surface voxels in the ROI.

    Also known as Feret Diameter.

    .. warning::
      This feature is only available when C Extensions are enabled
    """
    if self.diameters is None:
      self.diameters = self.calculateDiameters()
    return self.diameters[3]

  def getMaximum2DDiameterSliceFeatureValue(self):
    r"""
    **9. Maximum 2D diameter (Slice)**

    Maximum 2D diameter (Slice) is defined as the largest pairwise Euclidean distance between tumor surface voxels in
    the row-column (generally the axial) plane.

    .. warning::
      This feature is only available when C Extensions are enabled
    """
    if self.diameters is None:
      self.diameters = self.calculateDiameters()
    return self.diameters[0]

  def getMaximum2DDiameterColumnFeatureValue(self):
    r"""
    **10. Maximum 2D diameter (Column)**

    Maximum 2D diameter (Column) is defined as the largest pairwise Euclidean distance between tumor surface voxels in
    the row-slice (usually the coronal) plane.

    .. warning::
      This feature is only available when C Extensions are enabled
    """
    if self.diameters is None:
      self.diameters = self.calculateDiameters()
    return self.diameters[1]

  def getMaximum2DDiameterRowFeatureValue(self):
    r"""
    **11. Maximum 2D diameter (Row)**

    Maximum 2D diameter (Row) is defined as the largest pairwise Euclidean distance between tumor surface voxels in the
    column-slice (usually the sagittal) plane.

    .. warning::
      This feature is only available when C Extensions are enabled
    """
    if self.diameters is None:
      self.diameters = self.calculateDiameters()
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

  def _interpolate(self, grid, p1, p2):
    diff = (.5 - self.maskArray[tuple(grid[p1])]) / (self.maskArray[tuple(grid[p2])] - self.maskArray[tuple(grid[p1])])
    return (grid[p1] + ((grid[p2] - grid[p1]) * diff)) * self.pixelSpacing

  def _getMarchingTables(self):
    vertList = numpy.array(((0, 0, 0.5), (0, 0.5, 1), (0, 1, 0.5), (0, 0.5, 0),
                            (1, 0, 0.5), (1, 0.5, 1), (1, 1, 0.5), (1, 0.5, 0),
                            (0.5, 0, 0), (0.5, 0, 1), (0.5, 1, 1), (0.5, 1, 0)), dtype='float64')
    vertList *= self.pixelSpacing[None, :]

    triTable = [[],
                [[0, 8, 3]],
                [[0, 1, 9]],
                [[1, 8, 3], [9, 8, 1]],
                [[1, 2, 10]],
                [[0, 8, 3], [1, 2, 10]],
                [[9, 2, 10], [0, 2, 9]],
                [[2, 8, 3], [2, 10, 8], [10, 9, 8]],
                [[3, 11, 2]],
                [[0, 11, 2], [8, 11, 0]],
                [[1, 9, 0], [2, 3, 11]],
                [[1, 11, 2], [1, 9, 11], [9, 8, 11]],
                [[3, 10, 1], [11, 10, 3]],
                [[0, 10, 1], [0, 8, 10], [8, 11, 10]],
                [[3, 9, 0], [3, 11, 9], [11, 10, 9]],
                [[9, 8, 10], [10, 8, 11]],
                [[4, 7, 8]],
                [[4, 3, 0], [7, 3, 4]],
                [[0, 1, 9], [8, 4, 7]],
                [[4, 1, 9], [4, 7, 1], [7, 3, 1]],
                [[1, 2, 10], [8, 4, 7]],
                [[3, 4, 7], [3, 0, 4], [1, 2, 10]],
                [[9, 2, 10], [9, 0, 2], [8, 4, 7]],
                [[2, 10, 9], [2, 9, 7], [2, 7, 3], [7, 9, 4]],
                [[8, 4, 7], [3, 11, 2]],
                [[11, 4, 7], [11, 2, 4], [2, 0, 4]],
                [[9, 0, 1], [8, 4, 7], [2, 3, 11]],
                [[4, 7, 11], [9, 4, 11], [9, 11, 2], [9, 2, 1]],
                [[3, 10, 1], [3, 11, 10], [7, 8, 4]],
                [[1, 11, 10], [1, 4, 11], [1, 0, 4], [7, 11, 4]],
                [[4, 7, 8], [9, 0, 11], [9, 11, 10], [11, 0, 3]],
                [[4, 7, 11], [4, 11, 9], [9, 11, 10]],
                [[9, 5, 4]],
                [[9, 5, 4], [0, 8, 3]],
                [[0, 5, 4], [1, 5, 0]],
                [[8, 5, 4], [8, 3, 5], [3, 1, 5]],
                [[1, 2, 10], [9, 5, 4]],
                [[3, 0, 8], [1, 2, 10], [4, 9, 5]],
                [[5, 2, 10], [5, 4, 2], [4, 0, 2]],
                [[2, 10, 5], [3, 2, 5], [3, 5, 4], [3, 4, 8]],
                [[9, 5, 4], [2, 3, 11]],
                [[0, 11, 2], [0, 8, 11], [4, 9, 5]],
                [[0, 5, 4], [0, 1, 5], [2, 3, 11]],
                [[2, 1, 5], [2, 5, 8], [2, 8, 11], [4, 8, 5]],
                [[10, 3, 11], [10, 1, 3], [9, 5, 4]],
                [[4, 9, 5], [0, 8, 1], [8, 10, 1], [8, 11, 10]],
                [[5, 4, 0], [5, 0, 11], [5, 11, 10], [11, 0, 3]],
                [[5, 4, 8], [5, 8, 10], [10, 8, 11]],
                [[9, 7, 8], [5, 7, 9]],
                [[9, 3, 0], [9, 5, 3], [5, 7, 3]],
                [[0, 7, 8], [0, 1, 7], [1, 5, 7]],
                [[1, 5, 3], [3, 5, 7]],
                [[9, 7, 8], [9, 5, 7], [10, 1, 2]],
                [[10, 1, 2], [9, 5, 0], [5, 3, 0], [5, 7, 3]],
                [[8, 0, 2], [8, 2, 5], [8, 5, 7], [10, 5, 2]],
                [[2, 10, 5], [2, 5, 3], [3, 5, 7]],
                [[7, 9, 5], [7, 8, 9], [3, 11, 2]],
                [[9, 5, 7], [9, 7, 2], [9, 2, 0], [2, 7, 11]],
                [[2, 3, 11], [0, 1, 8], [1, 7, 8], [1, 5, 7]],
                [[11, 2, 1], [11, 1, 7], [7, 1, 5]],
                [[9, 5, 8], [8, 5, 7], [10, 1, 3], [10, 3, 11]],
                [[5, 7, 0], [5, 0, 9], [7, 11, 0], [1, 0, 10], [11, 10, 0]],
                [[11, 10, 0], [11, 0, 3], [10, 5, 0], [8, 0, 7], [5, 7, 0]],
                [[11, 10, 5], [7, 11, 5]],
                [[10, 6, 5]],
                [[0, 8, 3], [5, 10, 6]],
                [[9, 0, 1], [5, 10, 6]],
                [[1, 8, 3], [1, 9, 8], [5, 10, 6]],
                [[1, 6, 5], [2, 6, 1]],
                [[1, 6, 5], [1, 2, 6], [3, 0, 8]],
                [[9, 6, 5], [9, 0, 6], [0, 2, 6]],
                [[5, 9, 8], [5, 8, 2], [5, 2, 6], [3, 2, 8]],
                [[2, 3, 11], [10, 6, 5]],
                [[11, 0, 8], [11, 2, 0], [10, 6, 5]],
                [[0, 1, 9], [2, 3, 11], [5, 10, 6]],
                [[5, 10, 6], [1, 9, 2], [9, 11, 2], [9, 8, 11]],
                [[6, 3, 11], [6, 5, 3], [5, 1, 3]],
                [[0, 8, 11], [0, 11, 5], [0, 5, 1], [5, 11, 6]],
                [[3, 11, 6], [0, 3, 6], [0, 6, 5], [0, 5, 9]],
                [[6, 5, 9], [6, 9, 11], [11, 9, 8]],
                [[5, 10, 6], [4, 7, 8]],
                [[4, 3, 0], [4, 7, 3], [6, 5, 10]],
                [[1, 9, 0], [5, 10, 6], [8, 4, 7]],
                [[10, 6, 5], [1, 9, 7], [1, 7, 3], [7, 9, 4]],
                [[6, 1, 2], [6, 5, 1], [4, 7, 8]],
                [[1, 2, 5], [5, 2, 6], [3, 0, 4], [3, 4, 7]],
                [[8, 4, 7], [9, 0, 5], [0, 6, 5], [0, 2, 6]],
                [[7, 3, 9], [7, 9, 4], [3, 2, 9], [5, 9, 6], [2, 6, 9]],
                [[3, 11, 2], [7, 8, 4], [10, 6, 5]],
                [[5, 10, 6], [4, 7, 2], [4, 2, 0], [2, 7, 11]],
                [[0, 1, 9], [4, 7, 8], [2, 3, 11], [5, 10, 6]],
                [[9, 2, 1], [9, 11, 2], [9, 4, 11], [7, 11, 4], [5, 10, 6]],
                [[8, 4, 7], [3, 11, 5], [3, 5, 1], [5, 11, 6]],
                [[5, 1, 11], [5, 11, 6], [1, 0, 11], [7, 11, 4], [0, 4, 11]],
                [[0, 5, 9], [0, 6, 5], [0, 3, 6], [11, 6, 3], [8, 4, 7]],
                [[6, 5, 9], [6, 9, 11], [4, 7, 9], [7, 11, 9]],
                [[10, 4, 9], [6, 4, 10]],
                [[4, 10, 6], [4, 9, 10], [0, 8, 3]],
                [[10, 0, 1], [10, 6, 0], [6, 4, 0]],
                [[8, 3, 1], [8, 1, 6], [8, 6, 4], [6, 1, 10]],
                [[1, 4, 9], [1, 2, 4], [2, 6, 4]],
                [[3, 0, 8], [1, 2, 9], [2, 4, 9], [2, 6, 4]],
                [[0, 2, 4], [4, 2, 6]],
                [[8, 3, 2], [8, 2, 4], [4, 2, 6]],
                [[10, 4, 9], [10, 6, 4], [11, 2, 3]],
                [[0, 8, 2], [2, 8, 11], [4, 9, 10], [4, 10, 6]],
                [[3, 11, 2], [0, 1, 6], [0, 6, 4], [6, 1, 10]],
                [[6, 4, 1], [6, 1, 10], [4, 8, 1], [2, 1, 11], [8, 11, 1]],
                [[9, 6, 4], [9, 3, 6], [9, 1, 3], [11, 6, 3]],
                [[8, 11, 1], [8, 1, 0], [11, 6, 1], [9, 1, 4], [6, 4, 1]],
                [[3, 11, 6], [3, 6, 0], [0, 6, 4]],
                [[6, 4, 8], [11, 6, 8]],
                [[7, 10, 6], [7, 8, 10], [8, 9, 10]],
                [[0, 7, 3], [0, 10, 7], [0, 9, 10], [6, 7, 10]],
                [[10, 6, 7], [1, 10, 7], [1, 7, 8], [1, 8, 0]],
                [[10, 6, 7], [10, 7, 1], [1, 7, 3]],
                [[1, 2, 6], [1, 6, 8], [1, 8, 9], [8, 6, 7]],
                [[2, 6, 9], [2, 9, 1], [6, 7, 9], [0, 9, 3], [7, 3, 9]],
                [[7, 8, 0], [7, 0, 6], [6, 0, 2]],
                [[7, 3, 2], [6, 7, 2]],
                [[2, 3, 11], [10, 6, 8], [10, 8, 9], [8, 6, 7]],
                [[2, 0, 7], [2, 7, 11], [0, 9, 7], [6, 7, 10], [9, 10, 7]],
                [[1, 8, 0], [1, 7, 8], [1, 10, 7], [6, 7, 10], [2, 3, 11]],
                [[11, 2, 1], [11, 1, 7], [10, 6, 1], [6, 7, 1]],
                [[8, 9, 6], [8, 6, 7], [9, 1, 6], [11, 6, 3], [1, 3, 6]],
                [[0, 9, 1], [11, 6, 7]],
                [[7, 8, 0], [7, 0, 6], [3, 11, 0], [11, 6, 0]],
                [[7, 11, 6]]]
    return vertList, triTable
