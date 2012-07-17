import os
import numpy as np
import pyfits
from momentfunc import moment
import pymconvolve

class gini:
    """Calculate gini coefficient at different radii. This will also call 
       the function to compute M20. The algorithm is as follows

       1. Find the pixels in the image which belong to the galaxy, ie. make
          a segmentation map. This can be done by smoothing the image by a 
          boxcar of size r(\eta)/5 

       2. The surface brightness at r(\eta), I_(\eta) is measured and 
          pixels in the smoothed image with flux values greater than 
          I_(\eta) and less than 10(\sigma) is assigned to the galaxy. 
          \sigma is the sky deviation and which removes any remaining 
          cosmic rays or spurious noise pixels in the image and are not 
          included in the segmentation map.

       3. The Gini coefficient can be computed by the equation
          G = (1 / Avg(X) * n * (n-1)) * Sum over pixel[(2 * i - n -1 ) * X] 
    """
    def __init__(self, z, xcntr, ycntr, pa, eg, r20, r50, r80, ext_rad, \
                 sky, skysig):
        self.z            = z
        self.xcntr        = xcntr
        self.ycntr        = ycntr
        self.ext_rad      = ext_rad
        self.r20          = r20
        self.r50          = r50
        self.r80          = r80
        self.pa           = pa
        self.eg           = eg
        self.sky          = sky
        self.skysig       = skysig
        self.sigma        = 0.2 * self.ext_rad / 1.5#the size of the boxcar
        self.segmentation = segmentation(self.z, self.xcntr, self.ycntr, \
                            self.pa, self.eg, self.sky, self.r20, self.r50, \
                            self.r80, self.ext_rad, self.sigma, self.skysig, 1)

def gauss_kern(size, sizey=None):
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
    g = np.exp(-(x**2/np.float(size)+y**2/np.float(sizey)))
    return g / g.sum()


def segmentation(zextract, xcntr, ycntr, pa, eg, sky, r20, r50, r80, \
                 ext_rad, sigma, skysig, lower):
    """This function find the segmentation map"""
    def gini_coef(I):
        """This will calculate Gini Coefficient"""
        oneD_I = np.reshape(I, (1,-1))
        sorted_I = np.sort(oneD_I)
        sorted_I = sorted_I[sorted_I > 0]
        nn = sorted_I.size
        average = np.abs(sorted_I).mean()
        sumI = 0.0
        for i in range(nn):
            sumI += (2 * (i + 1) - nn - 1) * abs(sorted_I[i])
        G = sumI / (average * nn * (nn - 1))
        return G
    sigma = int(sigma)
    co= np.cos(pa*np.pi/180.0)
    si= np.sin(pa*np.pi/180.0)
    one_minus_eg_sq = (1.0 - eg)**2.0
    SizeY = zextract.shape[0] # Y size
    SizeX = zextract.shape[1] # X size
    x = np.reshape(np.arange(SizeX * SizeY), (SizeY, SizeX)) % SizeX
    x = x.astype(np.float32)
    y = np.reshape(np.arange(SizeX * SizeY), (SizeY, SizeX)) / SizeX
    y = y.astype(np.float32)
    rx = (x - xcntr)* co + (y - ycntr) * si
    ry = (xcntr - x) * si + (y - ycntr) * co
    R = np.sqrt(rx**2.0 + ry**2.0 / one_minus_eg_sq)
    boxcar = np.reshape(np.ones(sigma * sigma), (sigma, sigma))    
    I_sigma = pymconvolve.Convolve(zextract, boxcar) 
    lower_thre = (I_sigma[R <= ext_rad / 1.5].sum() - \
                  I_sigma[R <= ext_rad / 1.5 - 1].sum()) / \
                 (3.14 * ((ext_rad / 1.5)**2. - (ext_rad / 1.5 - 1)**2.))
    upper_thre = 10.0 * skysig # the upper threshold
    #The below will finds the segmentation image
    I0 = np.where(R > (2.5 / 1.5) * ext_rad, 0, zextract)
    I1 = np.where(I_sigma < lower_thre, 0, I0)
    I2 = np.where(R > (ext_rad / 1.5 - (0.02 * ext_rad)), 0, I1)
    I3 = np.where(R <= (ext_rad / 1.5 - (0.02 * ext_rad)), 0, I1)
    I4 = np.where(I_sigma > upper_thre, 0, I3)
    I = I2 + I4
    G = gini_coef(I)
    mo = moment(I, xcntr, ycntr)
    M = mo.moment_of_light[0]
    # FIX the following will calculate the Gini coef. & M20 of the residual
    # image of asymmetry. Now it is turning off. In the future it can be 
    # replaced with the Gini coef. of residual galfit image
    # f = pyfits.open(c.outdir +'AResidual.fits')
    # res =  f[0].data
    # f.close()
    # if os.access(c.outdir +'AResidual.fits', os.F_OK):
    #    os.remove(c.outdir +'AResidual.fits')
    # res1=np.where(I == 0, 0, res)
    # G_res = gini_coef(res1)
    # mo1=moment(res, xcntr, ycntr)
    # M_res = mo1.moment_of_light[0]
    G_res = 0.0
    M_res = 0.0
    # END
    I80 = np.where(R > r80, 0, zextract)
    G80 = gini_coef(I80)
    mo2 = moment(I80, xcntr, ycntr)
    M80 = mo2.moment_of_light[0]

    I50 = np.where(R > r50, 0, zextract)
    G50 = gini_coef(I50)
    mo3 = moment(I50, xcntr, ycntr)
    M_50 = mo3.moment_of_light[0]

    I20 = np.where(R > r20, 0, zextract)
    G20 = gini_coef(I20)
    mo4 = moment(I20, xcntr, ycntr)
    M_20 = mo4.moment_of_light[0]
    return G, G_res, G80, G50, G20, M, M_res, M80, M_50, M_20
