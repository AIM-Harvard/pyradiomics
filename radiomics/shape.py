import numpy
from radiomics import base
import SimpleITK as sitk


class RadiomicsShape(base.RadiomicsFeaturesBase):
  r"""
  In this group of features we included descriptors of the three-dimensional size and shape of the tumor region.
  Let in the following definitions denote :math:`V` the volume and :math:`A` the surface area of the volume of interest.
  """

  def __init__(self, inputImage, inputMask, **kwargs):
    super(RadiomicsShape, self).__init__(inputImage, inputMask, **kwargs)

    self.pixelSpacing = inputImage.GetSpacing()[::-1]
    z, x, y = self.pixelSpacing
    self.cubicMMPerVoxel = z * x * y

    # Use SimpleITK for some shape features
    self.lssif = sitk.LabelShapeStatisticsImageFilter()
    self.lssif.SetComputeFeretDiameter(True)
    self.lssif.Execute(inputMask)

    # Pad inputMask to prevent index-out-of-range errors
    cpif = sitk.ConstantPadImageFilter()

    padding = numpy.tile(1, 3)
    cpif.SetPadLowerBound(padding)
    cpif.SetPadUpperBound(padding)

    self.inputMask = cpif.Execute(self.inputMask)

    # Reassign self.maskArray using the now-padded self.inputMask and make it binary
    self.maskArray = (sitk.GetArrayFromImage(self.inputMask) == self.label).astype('int')
    self.matrixCoordinates = numpy.where(self.maskArray != 0)

    # Volume and Surface Area are pre-calculated
    self.Volume = self.lssif.GetPhysicalSize(self.label)
    self.SurfaceArea = self._calculateSurfaceArea()

  def _calculateSurfaceArea(self):
    # define relative locations of the 8 voxels of a sampling cube
    gridAngles = numpy.array([(0, 0, 0), (0, 0, 1), (0, 1, 1), (0, 1, 0),
                              (1, 0, 0), (1, 0, 1), (1, 1, 1), (1, 1, 0)])
    # instantiate lookup tables
    edgeTable, triTable = self._getMarchingTables()

    minBounds = numpy.array([numpy.min(self.matrixCoordinates[0]), numpy.min(self.matrixCoordinates[1]),
                             numpy.min(self.matrixCoordinates[2])])
    maxBounds = numpy.array([numpy.max(self.matrixCoordinates[0]), numpy.max(self.matrixCoordinates[1]),
                             numpy.max(self.matrixCoordinates[2])])
    minBounds = numpy.where(minBounds < 1, 1, minBounds)
    maxBounds = numpy.where(maxBounds > self.maskArray.shape, self.maskArray.shape, maxBounds)

    S_A = 0.0
    # iterate over all voxels which may border segmentation or are a part of it
    for v_z in xrange(minBounds[0] - 1, maxBounds[0] + 1):
      for v_y in xrange(minBounds[1] - 1, maxBounds[1] + 1):
        for v_x in xrange(minBounds[2] - 1, maxBounds[2] + 1):
          # indices to corners of current sampling cube
          gridCell = gridAngles + [v_z, v_y, v_x]

          # generate lookup index for current cube
          cube_idx = 0
          for p_idx, p in enumerate(gridCell):
            if self.maskArray[tuple(p)] == 0:
              cube_idx |= 1 << p_idx

          # full lookup tables are symmetrical, if cube_idx >= 128, take the XOR,
          # this allows for lookup tables to be half the size.
          if cube_idx & 128:
            cube_idx = cube_idx ^ 0xff

          # lookup which edges contain vertices and calculate vertices coordinates
          if edgeTable[cube_idx] == 0:
            continue
          vertList = numpy.zeros((12, 3), dtype='float64')
          if edgeTable[cube_idx] & 1:
            vertList[0] = self._interpolate(gridCell, 0, 1)
          if edgeTable[cube_idx] & 2:
            vertList[1] = self._interpolate(gridCell, 1, 2)
          if edgeTable[cube_idx] & 4:
            vertList[2] = self._interpolate(gridCell, 2, 3)
          if edgeTable[cube_idx] & 8:
            vertList[3] = self._interpolate(gridCell, 3, 0)
          if edgeTable[cube_idx] & 16:
            vertList[4] = self._interpolate(gridCell, 4, 5)
          if edgeTable[cube_idx] & 32:
            vertList[5] = self._interpolate(gridCell, 5, 6)
          if edgeTable[cube_idx] & 64:
            vertList[6] = self._interpolate(gridCell, 6, 7)
          if edgeTable[cube_idx] & 128:
            vertList[7] = self._interpolate(gridCell, 7, 4)
          if edgeTable[cube_idx] & 256:
            vertList[8] = self._interpolate(gridCell, 0, 4)
          if edgeTable[cube_idx] & 512:
            vertList[9] = self._interpolate(gridCell, 1, 5)
          if edgeTable[cube_idx] & 1024:
            vertList[10] = self._interpolate(gridCell, 2, 6)
          if edgeTable[cube_idx] & 2048:
            vertList[11] = self._interpolate(gridCell, 3, 7)

          # calculate triangles
          for triangle in triTable[cube_idx]:
            a = vertList[triangle[1]] - vertList[triangle[0]]
            b = vertList[triangle[2]] - vertList[triangle[0]]
            c = numpy.cross(a, b)
            S_A += .5 * numpy.sqrt(numpy.sum(c ** 2))

    return S_A

  def _getMaximum2Ddiameter(self, dim):
    otherDims = tuple(set([0, 1, 2]) - set([dim]))

    a = numpy.array(zip(*self.matrixCoordinates))

    maxDiameter = 0
    # Check maximum diameter in every slice, retain the overall maximum
    for i in numpy.unique(a[:, dim]):
      # Retrieve all indices of mask in current slice
      plane = a[numpy.where(a[:, dim] == i)]

      minBounds = numpy.min(plane, 0)
      maxBounds = numpy.max(plane, 0)

      # Generate 2 sets of indices: one set of indices in zSlice where at least the x or y component of the index is equal to the
      # minimum indices in the current slice, and one set of indices where at least one element it is equal to the maximum
      edgeVoxelsMinCoords = numpy.vstack(
        [plane[plane[:, otherDims[0]] == minBounds[otherDims[0]]],
         plane[plane[:, otherDims[1]] == minBounds[otherDims[1]]]]) * self.pixelSpacing
      edgeVoxelsMaxCoords = numpy.vstack(
        [plane[plane[:, otherDims[0]] == maxBounds[otherDims[0]]],
         plane[plane[:, otherDims[1]] == maxBounds[otherDims[1]]]]) * self.pixelSpacing

      # generate a matrix of distances for every combination of an index in edgeVoxelsMinCoords and edgeVoxelsMaxCoords
      # By subtraction the distance between the x, y and z components are obtained. The euclidean distance is then calculated:
      # Sum the squares of dx, dy and dz components and then take the square root; Sqrt( Sum( dx^2 + dy^2 + dz^2 ) )
      distances = numpy.sqrt(numpy.sum((edgeVoxelsMaxCoords[:, None] - edgeVoxelsMinCoords[None, :]) ** 2, 2))
      tempMax = numpy.max(distances)
      if tempMax > maxDiameter:
        maxDiameter = tempMax

    return maxDiameter

  def getVolumeFeatureValue(self):
    r"""
    Calculate the volume of the tumor region in cubic millimeters.
    """
    return self.Volume

  def getSurfaceAreaFeatureValue(self):
    r"""
    Calculate the surface area of the tumor region in square millimeters.

    :math:`A = \displaystyle\sum^{N}_{i=1}{\frac{1}{2}|\textbf{a}_i\textbf{b}_i \times \textbf{a}_i\textbf{c}_i|}`

    Where:

    :math:`N` is the number of triangles forming the surface of the volume

    :math:`a_ib_i` and :math:`a_ic_i` are the edges of the :math:`i`\ :sup:`th` triangle formed by points :math:`a_i`,
    :math:`b_i` and :math:`c_i`

    """
    return (self.SurfaceArea)

  def getSurfaceVolumeRatioFeatureValue(self):
    r"""
    Calculate the surface area to volume ratio of the tumor region

    :math:`surface\ to\ volume\ ratio = \frac{A}{V}`
    """
    return (self.SurfaceArea / self.Volume)

  def getCompactness1FeatureValue(self):
    r"""
    Calculate the compactness (1) of the tumor region.

    :math:`compactness\ 1 = \frac{V}{\sqrt{\pi}A^{\frac{2}{3}}}`

    Compactness 1 is a measure of how compact the shape of the tumor is
    relative to a sphere (most compact). It is a dimensionless measure,
    independent of scale and orientation. Compactness 1 is defined as the
    ratio of volume to the :math:`\sqrt{\text{surface area}^3}`. This is a measure of the
    compactness of the shape of the image ROI
    """
    return ((self.Volume) / ((self.SurfaceArea) ** (2.0 / 3.0) * numpy.sqrt(numpy.pi)))

  def getCompactness2FeatureValue(self):
    r"""
    Calculate the Compactness (2) of the tumor region.

    :math:`compactness\ 2 = 36\pi\frac{V^2}{A^3}`

    Compactness 2 is a measure of how compact the shape of the tumor is
    relative to a sphere (most compact). It is a dimensionless measure,
    independent of scale and orientation. This is a measure of the compactness
    of the shape of the image ROI.
    """
    return ((36.0 * numpy.pi) * ((self.Volume) ** 2.0) / ((self.SurfaceArea) ** 3.0))

  def getMaximum3DDiameterFeatureValue(self):
    r"""
    Calculate the largest pairwise euclidean distance between tumor surface voxels.
    Also known as Feret Diameter.
    """
    return self.lssif.GetFeretDiameter(self.label)

  def getMaximum2DDiameterSliceFeatureValue(self):
    r"""
    Calculate the largest pairwise euclidean distance between tumor surface voxels in the row-column plane.
    """

    return self._getMaximum2Ddiameter(0)

  def getMaximum2DDiameterColumnFeatureValue(self):
    r"""
    Calculate the largest pairwise euclidean distance between tumor surface voxels in the row-slice plane.
    """

    return self._getMaximum2Ddiameter(1)

  def getMaximum2DDiameterRowFeatureValue(self):
    r"""
    Calculate the largest pairwise euclidean distance between tumor surface voxels in the column-slice plane.
    """

    return self._getMaximum2Ddiameter(2)

  def getSphericalDisproportionFeatureValue(self):
    r"""
    Calculate the Spherical Disproportion of the tumor region.

    :math:`spherical\ disproportion = \frac{A}{4\pi R^2}`

    Where :math:`R` is the radius of a sphere with the same volume as the tumor.

    Spherical Disproportion is the ratio of the surface area of the
    tumor region to the surface area of a sphere with the same
    volume as the tumor region.
    """
    R = self.lssif.GetEquivalentSphericalRadius(self.label)
    return ((self.SurfaceArea) / (4.0 * numpy.pi * (R ** 2.0)))

  def getSphericityFeatureValue(self):
    r"""
    Calculate the Sphericity of the tumor region.

    :math:`sphericity = \frac{\pi^{\frac{1}{3}}(6V)^{\frac{2}{3}}}{A}`

    Sphericity is a measure of the roundness of the shape of the tumor region
    relative to a sphere. This is another measure of the compactness of a tumor.
    """
    return (((numpy.pi) ** (1.0 / 3.0) * (6.0 * self.Volume) ** (2.0 / 3.0)) / (self.SurfaceArea))

  def getElongationFeatureValue(self):
    """

    """
    return self.lssif.GetElongation(self.label)

  def getFlatnessFeatureValue(self):
    """

    """
    return self.lssif.GetFlatness(self.label)

  def getRoundnessFeatureValue(self):
    """

    """
    return self.lssif.GetRoundness(self.label)

  def _interpolate(self, grid, p1, p2):
    diff = (.5 - self.maskArray[tuple(grid[p1])]) / (self.maskArray[tuple(grid[p2])] - self.maskArray[tuple(grid[p1])])
    return (grid[p1] + ((grid[p2] - grid[p1]) * diff)) * self.pixelSpacing

  def _getMarchingTables(self):
    edgeTable = [0x0, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
                 0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
                 0x190, 0x99, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
                 0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
                 0x230, 0x339, 0x33, 0x13a, 0x636, 0x73f, 0x435, 0x53c,
                 0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
                 0x3a0, 0x2a9, 0x1a3, 0xaa, 0x7a6, 0x6af, 0x5a5, 0x4ac,
                 0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
                 0x460, 0x569, 0x663, 0x76a, 0x66, 0x16f, 0x265, 0x36c,
                 0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
                 0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff, 0x3f5, 0x2fc,
                 0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
                 0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55, 0x15c,
                 0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
                 0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc,
                 0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0]

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
    return edgeTable, triTable
