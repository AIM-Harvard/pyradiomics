import numpy
import collections
from radiomics import base, preprocessing
import SimpleITK as sitk

class RadiomicsGLCM(base.RadiomicsFeaturesBase):
    """
    GLCM feature calculation.
    """
    def __init__(self, inputImage, inputMask, binWidth):
        super(RadiomicsGLCM,self).__init__(inputImage, inputMask)

        self.imageArray = sitk.GetArrayFromImage(inputImage)
        self.maskArray = sitk.GetArrayFromImage(inputMask)
        self.binWidth = binWidth

        (self.matrix, self.matrixCoordinates) = \
          preprocessing.RadiomicsHelpers.padTumorMaskToCube(self.imageArray,self.maskArray)

        self.targetVoxelArray = self.matrix[self.matrixCoordinates]
        self.coefficients = {}
        self.P_glcm = {}

        # binning
        self.matrix, self.histogram = preprocessing.RadiomicsHelpers.binImage(self.binWidth, self.targetVoxelArray, self.matrix, self.matrixCoordinates)
        self.coefficients['Ng'] = len(self.histogram[0])

        self.createGLCM()

        self.calculateCoefficients()

    def calculateGLCM(self, distances, angles):
        """
        13 GLCM matrices for each image for every direction from the voxel
        (13 for each neighboring voxel from a reference voxel centered in a 3x3 cube)
        for GLCM matrices P(i,j;gamma, a), gamma = 1, a = 1...13
        """

        maxrows = self.matrix.shape[2]
        maxcols = self.matrix.shape[1]
        maxheight = self.matrix.shape[0]
        indices = zip(*self.matrixCoordinates)

        numIndices = len(indices)
        index = 0

        for h, c, r in indices:
            if index % 100 == 0:
                percentProgress = (float(index) / float(numIndices)) * 100.0
                print ('calculate GLCM: %.2f%s' % (percentProgress, '%'))
            index = index + 1
            for angles_idx, angle in enumerate(angles):
                for distances_idx, distance in enumerate(distances):
                    i = self.matrix[h, c, r]
                    i_idx = int(i-1)

                    row = r + angle[2]
                    col = c + angle[1]
                    height = h + angle[0]

                    if row >= 0 and row < maxrows and col >= 0 and col < maxcols:
                        if tuple((height, col, row)) in indices:
                            j = self.matrix[height, col, row]
                            j_idx = int(j-1)
                            self.P_glcm[i_idx, j_idx, distances_idx, angles_idx] += 1

    def createGLCM(self):
        """
        Generate container for GLCM Matrices: P_glcm
        """
        Ng = self.coefficients['Ng']

        distances=numpy.array([1])
        angles = numpy.array([ (0, 1, 0),
                               (-1, 1, 0),
                               (-1, 0, 0),
                               (-1, -1, 0),
                               (0, 1, -1),
                               (0, 0, -1),
                               (0, -1, -1),
                               (-1, 0, -1),
                               (1, 0, -1),
                               (-1, 1, -1),
                               (1, -1, -1),
                               (-1, -1, -1),
                               (1, 1, -1) ])
        self.P_glcm = numpy.zeros( (Ng, Ng, distances.size, int(angles.shape[0])), dtype='float32' )
        self.calculateGLCM(distances, angles)

        # Normalize each glcm
        for q in xrange(int(angles.shape[0])):
            self.P_glcm[:,:,0,q] = self.P_glcm[:,:,0,q]/(self.P_glcm[:,:,0,q].sum())

    # check if ivector and jvector can be replaced
    def calculateCoefficients(self):
        """
        Fill in the coefficients array
        """
        Ng = self.coefficients['Ng']

        eps = numpy.spacing(1)

        #
        # Calculate GLCM Coefficients
        #
        # shape = (Ng, distances.size, angles)
        ivector = numpy.arange(1,Ng+1)
        # shape = (Ng, distances.size, angles)
        jvector = numpy.arange(1,Ng+1)


        i,j = numpy.meshgrid(numpy.arange(1,self.P_glcm.shape[0]+1), numpy.arange(1,self.P_glcm.shape[1]+1))

        # shape = (Ng, Ng)
        prodMatrix = numpy.multiply.outer(ivector, jvector)
        # shape = (Ng, Ng)
        # sum-1 so 111 in corner
        sumMatrix = numpy.add.outer(ivector, jvector)
        # shape = (Ng, Ng)
        # absval(subtract)+1 so 1-1 and 56-56 in corners
        diffMatrix = numpy.absolute(numpy.subtract.outer(ivector, jvector))

        # shape = (2*Ng-1)
        # 1-111
        kValuesSum = numpy.arange(2, (Ng*2)+1) - 1
        # shape = (Ng-1)
        # 0-55
        kValuesDiff = numpy.arange(0,Ng) + 1

        # shape = (distances.size, angles)
        u = self.P_glcm.mean(0).mean(0)
        # marginal row probabilities #shape = (Ng, distances.size, angles)
        px = self.P_glcm.sum(1)
        # marginal column probabilities #shape = (Ng, distances.size, angles)
        py = self.P_glcm.sum(0)

        # shape = (distances.size, angles)
        ux = numpy.sum(numpy.sum( ivector[None,:,None,None]*self.P_glcm, 0 ), 0 )
        uy = numpy.sum(numpy.sum( jvector[:,None,None,None]*self.P_glcm, 0 ), 0 )

        # shape = (distances.size, angles)
        sigx = numpy.sum(numpy.sum(self.P_glcm*((i[:,:,None,None]-ux[None,None,:,:])**2), 0 ), 0 )**0.5
        # shape = (distances.size, angles)
        sigy = numpy.sum(numpy.sum(self.P_glcm*((j[:,:,None,None]-uy[None,None,:,:])**2), 0 ), 0 )**0.5

        # shape = (2*Ng-1, distances.size, angles)
        # 111 long
        pxAddy = numpy.array([ numpy.sum(self.P_glcm[i+j == k], 0) for k in kValuesSum ])
        # shape = (Ng, distances.size, angles)
        # 56 long
        pxSuby = numpy.array([ numpy.sum(self.P_glcm[numpy.abs(i-j) == k], 0) for k in kValuesDiff ])

        # entropy of px #shape = (distances.size, angles)
        HX = (-1) * numpy.sum( (px * numpy.log(px+eps)), 0)
        # entropy of py #shape = (distances.size, angles)
        HY = (-1) * numpy.sum( (py * numpy.log(py+eps)), 0)
        # shape = (distances.size, angles)
        HXY = (-1) * numpy.sum( numpy.sum( (self.P_glcm * numpy.log(self.P_glcm+eps)), 0 ), 0 )

        #shape = (distances.size, angles)
        HXY1 = (-1) * numpy.sum( numpy.sum( (self.P_glcm * numpy.log(self.P_glcm+eps)), 0 ), 0)
        # shape = (distances.size, angles)
        HXY2 = (-1) * numpy.sum( numpy.sum( ((px*py) * numpy.log(px*py+eps)), 0 ), 0)

        self.coefficients['Ng'] = Ng
        self.coefficients['eps'] = eps
        self.coefficients['i'] = i
        self.coefficients['j'] = j
        self.coefficients['kValuesSum'] = kValuesSum
        self.coefficients['kValuesDiff'] = kValuesDiff
        self.coefficients['u'] = u
        self.coefficients['px'] = px
        self.coefficients['py'] = py
        self.coefficients['ux'] = ux
        self.coefficients['uy'] = uy
        self.coefficients['sigx'] = sigx
        self.coefficients['sigy'] = sigy
        self.coefficients['pxAddy'] = pxAddy
        self.coefficients['pxSuby'] = pxSuby
        self.coefficients['HX'] = HX
        self.coefficients['HY'] = HY
        self.coefficients['HXY'] = HXY
        self.coefficients['HXY1'] = HXY1
        self.coefficients['HXY2'] = HXY2

    def getAutocorrelationFeatureValue(self):
        """
        Using the i and j arrays, calculate and return the
        autocorrelation feature value
        """
        i = self.coefficients['i']
        j = self.coefficients['j']
        ac = numpy.sum(numpy.sum(self.P_glcm * (i*j)[:,:,None,None], 0 ), 0 )
        return (ac.mean())

    def getClusterProminenceFeatureValue(self):
        """
        Using coefficients i, j, ux, uy, calculate the cluster prominence and return it.
        """
        i = self.coefficients['i']
        j = self.coefficients['j']
        ux = self.coefficients['ux']
        uy = self.coefficients['uy']
        cp = numpy.sum( numpy.sum( (self.P_glcm * (((i+j)[:,:,None,None] - ux[None,None,:,:] - uy[None,None,:,:])**4)), 0 ), 0 )
        return (cp.mean())

    def getClusterShadeFeatureValue(self):
        """
        Using coefficients i, j, ux, uy, calculate the cluster shade and return it.
        """
        i = self.coefficients['i']
        j = self.coefficients['j']
        ux = self.coefficients['ux']
        uy = self.coefficients['uy']
        cs = numpy.sum( numpy.sum( (self.P_glcm * (((i+j)[:,:,None,None] - ux[None,None,:,:] - uy[None,None,:,:])**3)), 0 ), 0 )
        return (cs.mean())

    def getClusterTendencyFeatureValue(self):
        """
        Using coefficients i, j, ux, uy, calculate the cluster tendency and return it.
        """
        i = self.coefficients['i']
        j = self.coefficients['j']
        ux = self.coefficients['ux']
        uy = self.coefficients['uy']
        ct = numpy.sum( numpy.sum( (self.P_glcm * (((i+j)[:,:,None,None] - ux[None,None,:,:] - uy[None,None,:,:])**2)), 0 ), 0 )
        return (ct.mean())

    def getContrastFeatureValue(self):
        """
        Using coefficients i, j, calculate the contrast and return it.
        """
        i = self.coefficients['i']
        j = self.coefficients['j']
        cont = numpy.sum( numpy.sum( (self.P_glcm * ((numpy.abs(i-j))[:,:,None,None]**2)), 0 ), 0 )
        return (cont.mean())

    def getCorrelationFeatureValue(self):
        """
        Using coefficients i, j, ux, uy, sigx, sigy, calculate the correlation and return it.
        """
        i = self.coefficients['i']
        j = self.coefficients['j']
        ux = self.coefficients['ux']
        uy = self.coefficients['uy']
        sigx = self.coefficients['sigx']
        sigy = self.coefficients['sigy']
        corm = numpy.sum(numpy.sum(self.P_glcm*(i[:,:,None,None]-ux[None,None,:,:])*(j[:,:,None,None]-uy[None,None,:,:]), 0 ), 0 )
        corr = corm/(sigx*sigy)
        return (corr.mean())

    def getDifferenceEntropyFeatureValue(self):
        """
        Using coefficients pxSuby, eps, calculate the difference entropy and return it.
        """
        pxSuby = self.coefficients['pxSuby']
        eps = self.coefficients['eps']
        difent = (-1) * numpy.sum( (pxSuby*numpy.log(pxSuby+eps)), 0 )
        return (difent.mean())

    def getDissimilarityFeatureValue(self):
        """
        Using coefficients i, j, calculate the dissimilarity and return it.
        """
        i = self.coefficients['i']
        j = self.coefficients['j']
        dis = numpy.sum( numpy.sum( (self.P_glcm * (numpy.abs(i-j))[:,:,None,None]), 0 ), 0 )
        return (dis.mean())

    def getEnergyFeatureValue(self):
        """
        Using P_glcm, calculate the energy and return it.
        """
        ene = numpy.sum( numpy.sum( (self.P_glcm**2), 0), 0 )
        return (ene.mean())

    def getEntropyFeatureValue(self):
        """
        Using coefficients eps, calculate the entropy and return it.
        """
        eps = self.coefficients['eps']
        ent = (-1) * numpy.sum( numpy.sum( (self.P_glcm * numpy.log(self.P_glcm+eps)), 0), 0)
        return (ent.mean())

    def getHomogeneity1FeatureValue(self):
        """
        Using coefficients i, j, calculate the homogeneity 1 and return it.
        """
        i = self.coefficients['i']
        j = self.coefficients['j']
        homo1 = numpy.sum( numpy.sum( (self.P_glcm / (1 + (numpy.abs(i-j))[:,:,None,None])), 0 ), 0 )
        return (homo1.mean())

    def getHomogeneity2FeatureValue(self):
        """
        Using coefficients i, j, calculate the homogeneity 2 and return it.
        """
        i = self.coefficients['i']
        j = self.coefficients['j']
        homo2 = numpy.sum( numpy.sum( (self.P_glcm / (1 + (numpy.abs(i-j))[:,:,None,None]**2)), 0 ), 0 )
        return (homo2.mean())

    def getImc1FeatureValue(self):
        """
        Using coefficients HX, HY, HXY, HXY1, HXY2, calculate the imc 1 and return it.
        """
        HX = self.coefficients['HX']
        HY = self.coefficients['HY']
        HXY = self.coefficients['HXY']
        HXY1 = self.coefficients['HXY1']
        HXY2 = self.coefficients['HXY2']
        imc1 = (HXY - HXY1)/numpy.max(([HX,HY]),0)
        return (imc1.mean())

    def getImc2FeatureValue(self):
        """
        Using coefficients HXY, HXY2, calculate the imc 2 and return it.
        NOT IMPLEMENTED
        """
        # Produces Error

        HXY = self.coefficients['HXY']
        HXY2 = self.coefficients['HXY2']

        # imc2[g,a] = ( 1-numpy.e**(-2*(HXY2[g,a]-HXY[g,a])) )**(0.5) #nan value too high
        # matlab:(1-exp(-2*(hxy2-hxy)))^0.5;

        #produces Nan(square root of a negative)
        #exponent = decimal.Decimal( -2*(HXY2[g,a]-HXY[g,a]) )
        #imc2.append( ( decimal.Decimal(1)-decimal.Decimal(numpy.e)**(exponent) )**(decimal.Decimal(0.5)) )

        #if meanFlag:
          #return (homo2.mean())
        #else:
          #return homo2

    def getIdmnFeatureValue(self):
        """
        Using coefficients i, j, Ng, calculate the idmn  and return it.
        """
        i = self.coefficients['i']
        j = self.coefficients['j']
        Ng = self.coefficients['Ng']
        idmn = numpy.sum( numpy.sum( (self.P_glcm / (1 + (((numpy.abs(i-j))[:,:,None,None]**2)/(Ng**2)))), 0 ), 0 )
        return (idmn.mean())

    def getIdnFeatureValue(self):
        """
        Using coefficients i, j, Ng, calculate the idn and return it.
        """
        i = self.coefficients['i']
        j = self.coefficients['j']
        Ng = self.coefficients['Ng']
        idn  = numpy.sum( numpy.sum( (self.P_glcm / (1 + ((numpy.abs(i-j))[:,:,None,None]/Ng))), 0 ), 0 )
        return (idn.mean())

    def getInverseVarianceFeatureValue(self):
        """
        Use the i, j, Ng coeffients, calculate the inverse variance and return it.
        """
        i = self.coefficients['i']
        j = self.coefficients['j']
        Ng = self.coefficients['Ng']
        maskDiags = numpy.ones((numpy.abs(i-j)).shape, dtype = bool)
        maskDiags[numpy.diag_indices(Ng)] = False
        inv = numpy.sum( (self.P_glcm[maskDiags] / ((numpy.abs(i-j))[:,:,None,None]**2)[maskDiags]), 0 )
        return (inv.mean())

    def getMaximumProbabilityFeatureValue(self):
        """
        Using P_glcm, calculate the maximum probability and return it.
        """
        maxprob = self.P_glcm.max(0).max(0)
        return (maxprob.mean())

    def getSumAverageFeatureValue(self):
        """
        Using coefficients pxAddy, kValuesSum, calculate the sum average and return it.
        """
        pxAddy = self.coefficients['pxAddy']
        kValuesSum = self.coefficients['kValuesSum']
        sumavg =  numpy.sum( (kValuesSum[:,None,None]*pxAddy), 0 )
        return (sumavg.mean())

    def getSumEntropyFeatureValue(self):
        """
        Using coefficients pxAddy, eps, calculate the sum entropy and return it.
        """
        pxAddy = self.coefficients['pxAddy']
        eps = self.coefficients['eps']
        sumentr = (-1) * numpy.sum( (pxAddy * numpy.log(pxAddy+eps)), 0 )
        return (sumentr.mean())

    def getSumVarianceFeatureValue(self):
        """
        Using coefficients pxAddy, kValuesSum, calculate the sum variance and return it.
        """
        pxAddy = self.coefficients['pxAddy']
        kValuesSum = self.coefficients['kValuesSum']
        sumvar = numpy.sum( (pxAddy*((kValuesSum[:,None,None] - kValuesSum[:,None,None]*pxAddy)**2)), 0 )
        return (sumvar.mean())

    def getSumSquaresFeatureValue(self):
        """
        Using coefficients j, u, calculate the sum of the squares and return it.
        """
        j = self.coefficients['j']
        u = self.coefficients['u']
        # Also known as Variance
        ss = numpy.sum( numpy.sum( (self.P_glcm * ((j[:,:,None,None]-u[None,None,:,:])**2)), 0 ), 0 )
        return (ss.mean())
