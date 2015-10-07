import numpy
import collections
import pdb

def smallAreaEmphasis(self, P_glszm, jvector, sumP_glszm):
    try:
        #sre = numpy.sum( numpy.sum( (P_glszm/(self.j[:,:,None]**2)) , 0 ), 0 ) / (sumP_glszm[None,None,:])
        sre = numpy.sum( (self.pr/(jvector[:,None]**2)), 0 ) / sumP_glszm
    except ZeroDivisionError:
        sre = 0    
    return (sre.mean())
  
def largeAreaEmphasis(self, P_glszm, jvector, sumP_glszm):
    try:
        #lre = numpy.sum( numpy.sum( (P_glszm*(self.j[:,:,None]**2)) , 0 ), 0 ) / (sumP_glszm[None,None,:])
        lre = numpy.sum( (self.pr*(jvector[:,None]**2)), 0 ) / sumP_glszm
    except ZeroDivisionError:
        lre = 0
    return (lre.mean())

def intensityVariability(self, P_glszm, sumP_glszm):
    try:
        gln = numpy.sum( (self.pg**2) , 0 ) / sumP_glszm   
    except ZeroDivisionError:
        gln = 0  
    return (gln.mean())
  
def sizeZoneVariability(self, P_glszm, sumP_glszm):
    try:
        rln = numpy.sum( (self.pr**2) , 0 ) / sumP_glszm 
    except ZeroDivisionError:
        rln = 0
    return (rln.mean())

def zonePercentage(self, P_glszm, Np):
    try:
        rp = numpy.sum( numpy.sum( (P_glszm/(Np)) , 0 ), 0 )
    except ZeroDivisionError:
        rp = 0
    return (rp.mean())

def lowIntensityEmphasis(self, P_glszm, ivector, sumP_glszm):
    try:
        #pdb.set_trace()
        #lglre = numpy.sum( numpy.sum( (P_glszm/(self.i[:,:,None]**2)) , 0 ), 0 )/ (sumP_glszm[None,None,:])
        lglre = numpy.sum( (self.pg/(ivector[:,None]**2)) , 0 ) / sumP_glszm
    except ZeroDivisionError:
        lglre = 0 
    return (lglre.mean())
    
def highIntensityEmphasis(self, P_glszm, ivector, sumP_glszm):
  try:
      #hglre = numpy.sum( numpy.sum( (P_glszm*(self.i[:,:,None]**2)) , 0 ), 0 ) / (sumP_glszm[None,None,:])
      hglre = numpy.sum( (self.pg*(ivector[:,None]**2)) , 0 ) / sumP_glszm
  except ZeroDivisionError:
      hglre = 0 
  return (hglre.mean())
  
def lowIntensitySmallAreaEmphasis(self, P_glszm, ivector, jvector, sumP_glszm):
  try:
      srlgle = numpy.sum( numpy.sum( (P_glszm/((ivector[:,None,None]**2)*(jvector[None,:,None]**2))) , 0 ), 0 ) / sumP_glszm
  except ZeroDivisionError:
      srlgle = 0
  return (srlgle.mean())
   
def highIntensitySmallAreaEmphasis(self, P_glszm, ivector, jvector, sumP_glszm):
  try:
      srhgle = numpy.sum( numpy.sum( (P_glszm*(ivector[:,None,None]**2)/(jvector[None,:,None]**2)) , 0 ), 0 ) / sumP_glszm
  except ZeroDivisionError:
      srhgle = 0   
  return (srhgle.mean())
      
def lowIntensityLargeAreaEmphasis(self, P_glszm, ivector, jvector, sumP_glszm):
  try:
      lrlgle = numpy.sum( numpy.sum( (P_glszm*(jvector[None,:,None]**2)/(ivector[:,None,None]**2)) , 0 ), 0 ) / sumP_glszm
  except ZeroDivisionError:
      lrlgle = 0
  return (lrlgle.mean())
          
def highIntensityLargeAreaEmphasis(self, P_glszm, ivector, jvector, sumP_glszm):
  try:
      lrhgle = numpy.sum( numpy.sum( (P_glszm*((jvector[None,:,None]**2)*(ivector[:,None,None]**2))) , 0 ), 0 ) / sumP_glszm
  except ZeroDivisionError:
      lrhgle = 0 
  return (lrhgle.mean())