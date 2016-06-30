import numpy
import SimpleITK as sitk
from radiomics import base, imageoperations
import pdb
from tqdm import trange

class RadiomicsGLSZM(base.RadiomicsFeaturesBase):
  """GLSZM feature calculation."""
  def __init__(self, inputImage, inputMask, **kwargs):
    super(RadiomicsGLSZM,self).__init__(inputImage, inputMask, **kwargs)

    if inputImage == None or inputMask == None:
      if self.verbose: print('ERROR GLSZM: missing input image or mask')
      return

    self.imageArray = sitk.GetArrayFromImage(inputImage)
    self.maskArray = sitk.GetArrayFromImage(inputMask)

    (self.matrix, self.matrixCoordinates) = \
      imageoperations.padTumorMaskToCube(self.imageArray,self.maskArray)

    self.targetVoxelArray = self.matrix[self.matrixCoordinates]
    self.coefficients = {}
    self.P_glszm = {}

    # binning
    self.matrix, self.histogram = imageoperations.binImage(self.binWidth, self.targetVoxelArray, self.matrix, self.matrixCoordinates)
    self.coefficients['Ng'] = len(self.histogram[0])
    self.coefficients['grayLevels'] = numpy.linspace(1,self.coefficients['Ng'],num=self.coefficients['Ng'])
    self.coefficients['Np'] = self.targetVoxelArray.size

    self.calculateGLSZM()
    self.calculateCoefficients()


  def calculateGLSZM(self):
    """
    Number of times a 26-connected region with a
    gray level and voxel count occurs in an image. P_glszm[level, voxel_count] = # occurrences
    """
    angles = numpy.array([ (0, 0, 1),
                           (0, 1, 0),
                           (0, 1, 1),
                           (0, 1, -1),
                           (1, 0, 0),
                           (1, 0, 1),
                           (1, 0, -1),
                           (1, 1, 0),
                           (1, 1, 1),
                           (1, 1, -1),
                           (1, -1, 0),
                           (1, -1, 1),
                           (1, -1, -1) ])

    # Empty GLSZ matrix
    P_glszm = numpy.zeros((self.coefficients['grayLevels'].size, self.coefficients['Np']))

    # Iterate over all gray levels in the image
    numGrayLevels = self.coefficients['grayLevels'].size+1

    if self.verbose: bar = trange(numGrayLevels-1, desc= 'calculate GLSZM')

    for i in xrange(1, numGrayLevels):
      # give some progress
      if self.verbose: bar.update()

      dataTemp = numpy.where(self.matrix==i, 1, 0)
      ind = zip(*numpy.where(dataTemp==1))
      ind = list(set(ind).intersection(set(zip(*self.matrixCoordinates))))

      while ind: # check if ind is not empty
        # Current label number and first coordinate
        ind_node = ind[0]

        # get all coordinates in the 26-connected region
        region_full = [ind_node] + [tuple(sum(a) for a in zip(ind_node,angle_i)) for angle_i in angles]

        # get coordinates in 26-connected region with same grey level
        region_level = list(set(ind).intersection(set(region_full)))

        # Set already processed indices to zero
        for pos in region_level:
          dataTemp[pos] = 0

        # Size of the region (# of voxels in region)
        regionSize = len(region_level)

        # Update the gray level size zone matrix
        P_glszm[i-1,regionSize-1] += 1

        # Find unprocessed nonzero positions for current gray level
        ind = zip(*numpy.where(dataTemp==1))
        ind = list(set(ind).intersection(set(zip(*self.matrixCoordinates))))

    if self.verbose: bar.close()

    # Crop gray-level axis of GLSZM matrix to between minimum and maximum observed gray-levels
    # Crop size-zone area axis of GLSZM matrix up to maximum observed size-zone area
    P_glszm_bounds = numpy.argwhere(P_glszm)
    (xstart, ystart), (xstop, ystop) = P_glszm_bounds.min(0), P_glszm_bounds.max(0) + 1
    self.P_glszm = P_glszm[xstart:xstop,:ystop]

  def calculateCoefficients(self):
    sumP_glszm = numpy.sum(self.P_glszm, (0, 1) )

    # set sum to numpy.spacing(1) if sum is 0?
    if sumP_glszm == 0:
      sumP_glszm = 1

    pr = numpy.sum(self.P_glszm, 0)
    pg = numpy.sum(self.P_glszm, 1)

    ivector = numpy.arange(1, self.P_glszm.shape[0] + 1, dtype=numpy.float64)
    jvector = numpy.arange(1, self.P_glszm.shape[1] + 1, dtype=numpy.float64)

    self.coefficients['sumP_glszm'] = sumP_glszm
    self.coefficients['pr'] = pr
    self.coefficients['pg'] = pg
    self.coefficients['ivector'] = ivector
    self.coefficients['jvector'] = jvector

  def getSmallAreaEmphasisFeatureValue(self):
    """Calculate and return the Small Area Emphasis (SAE) value.

    A measure of the distribution of small size zones, with a greater value indicative
    of more smaller size zones and more fine textures.
    """
    try:
      sae = numpy.sum(self.coefficients['pr']/(self.coefficients['jvector']**2)) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      sae = numpy.core.numeric.NaN
    return (sae)

  def getLargeAreaEmphasisFeatureValue(self):
    """Calculate and return the Large Area Emphasis (LAE) value.

    A measure of the distribution of large area size zones, with a greater value indicative
    of more larger size zones and more coarse textures.
    """
    try:
      lae = numpy.sum(self.coefficients['pr']*(self.coefficients['jvector']**2)) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      lae = numpy.core.numeric.NaN
    return (lae)

  def getIntensityVariabilityFeatureValue(self):
    """Calculate and return the Intensity Variability (IV) value.

    Measures the variability of gray-level intensity values in the image, where a lower IV value
    correlates with more homogeneity in intensity values.
    """
    try:
      iv = numpy.sum(self.coefficients['pg']**2) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      iv = numpy.core.numeric.NaN
    return (iv)

  def getSizeZoneVariabilityFeatureValue(self):
    """Calculate and return the Size-Zone Variability (SZV) value.

    Measures the variability of size zone volumes in the image, where a lower SZV value
    correlates with more homogeneity in size zone volumes.
    """
    try:
      szv = numpy.sum(self.coefficients['pr']**2) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      szv = numpy.core.numeric.NaN
    return (szv)

  def getZonePercentageFeatureValue(self):
    """Calculate and return the Zone Percentage (ZP) value.

    Measures the homogeneity of the distribution of size zones in an image among the observed gray-levels.
    """
    try:
      zp = self.coefficients['sumP_glszm'] / numpy.sum(self.coefficients['pr']*self.coefficients['jvector'])
    except ZeroDivisionError:
      zp = numpy.core.numeric.NaN
    return (zp)

  def getLowIntensityEmphasisFeatureValue(self):
    """Calculate and return the Low Intensity Emphasis (LIE) value.

    Measures the distribution of lower gray-level size zones, with a higher value indicating a greater
    proportion of lower gray-level values and size zones in the image.
    """
    try:
      lie = numpy.sum( (self.coefficients['pg']/(self.coefficients['ivector']**2)) ) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      lie = numpy.core.numeric.NaN
    return (lie)

  def getHighIntensityEmphasisFeatureValue(self):
    """Calculate and return the High Intensity Emphasis (HIE) value.

    Measures the distribution of the higher gray-level values, with a higher value indicating
    a greater proportion of higher gray-level values and size zones in the image.
    """
    try:
      hie = numpy.sum( (self.coefficients['pg']*(self.coefficients['ivector']**2)) ) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      hie = numpy.core.numeric.NaN
    return (hie)

  def getLowIntensitySmallAreaEmphasisFeatureValue(self):
    """Calculate and return the Low Intensity Small Area Emphases (LISAE) value.

    Measures the proportion in the image of the joint distribution of smaller size zones with lower gray-level values.
    """
    try:
      lisae = numpy.sum( (self.P_glszm/((self.coefficients['ivector'][:,None]**2)*(self.coefficients['jvector'][None,:]**2))) , (0, 1) ) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      lisae = numpy.core.numeric.NaN
    return (lisae)

  def getHighIntensitySmallAreaEmphasisFeatureValue(self):
    """Calculate and return the High Intensity Small Area Emphases (HISAE) value.

    Measures the proportion in the image of the joint distribution of smaller size zones with higher gray-level values.
    """
    try:
      hisae = numpy.sum( (self.P_glszm*(self.coefficients['ivector'][:,None]**2)/(self.coefficients['jvector'][None,:]**2)) , (0, 1) ) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      hisae = numpy.core.numeric.NaN
    return (hisae)

  def getLowIntensityLargeAreaEmphasisFeatureValue(self):
    """Calculate and return the Low Intensity Large Area Emphases (LILAE) value.

    Measures the proportion in the image of the joint distribution of larger size zones with lower gray-level values.
    """
    try:
      lilae = numpy.sum( (self.P_glszm*(self.coefficients['jvector'][None,:]**2)/(self.coefficients['ivector'][:,None]**2)) , (0, 1) ) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      lilae = numpy.core.numeric.NaN
    return (lilae)

  def getHighIntensityLargeAreaEmphasisFeatureValue(self):
    """Calculate and return the High Intensity Large Area Emphases (HILAE) value.

    Measures the proportion in the image of the joint distribution of larger size zones with higher gray-level values.
    """
    try:
      hilae = numpy.sum( (self.P_glszm*((self.coefficients['jvector'][None,:]**2)*(self.coefficients['ivector'][:,None]**2))) , (0, 1) ) / self.coefficients['sumP_glszm']
    except ZeroDivisionError:
      hilae = numpy.core.numeric.NaN
    return (hilae)
