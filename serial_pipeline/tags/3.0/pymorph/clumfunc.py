import numpy as np
import pyfits
import pymconvolve
from rotate import ImSec

class clumpness:
    """The clumpness parameter and the algorithm used as follows
       1. The image is smoothed by a boxcar of width
          0.25 * r(Petrosian parameter = 0.2)

       2. The smoothness is computed with the radius 1.5 by using
          S = 10 * Sum(I_0 - I_S) / Sum(I_0)
          where I_0 is the galaxy pixels and I_S that of smoothed image

       3. Compute the average smoothness of the background and subtract 
          from S.

       4. The inner region of the galaxy is not considered in the
          computation of S as these are often unresolved.

       5. Use only the positive pixels for the computation."""
    def __init__(self, z, xcntr, ycntr, pa, eg, ext_rad, sigma, sky, flag_image):
        self.z                = z
        self.xcntr            = xcntr
        self.ycntr            = ycntr
        self.ext_rad          = ext_rad
        self.pa               = pa
        self.eg               = eg
        self.sky              = sky
        self.flag_image       = flag_image
        self.sigma            = np.int(sigma) #the size of the boxcar
        self.clumpness        = CLUMPNESS(self.z, self.xcntr, self.ycntr, \
                                self.pa, self.eg, self.sky, self.ext_rad, \
                                self.sigma, self.flag_image)

def CLUMPNESS(z, xcntr, ycntr, pa, eg, sky, ext_rad, sigma, flag_image):
    CutImDa, cut_xcntr, cut_ycntr, SizeY, SizeX, ymin, ymax, xmin, \
    xmax, flag_out = ImSec(z, xcntr, ycntr, ext_rad)
    co = np.cos(pa * np.pi / 180.0)
    si = np.sin(pa * np.pi / 180.0)  
    one_minus_eg_sq = (1.0 - eg)**2.0 
    x = np.reshape(np.arange(SizeX * SizeY), (SizeY, SizeX)) % SizeX
    x = x.astype(np.float32)
    y = np.reshape(np.arange(SizeX * SizeY), (SizeY, SizeX)) / SizeX
    y = y.astype(np.float32)
    tx = (x - cut_xcntr) * co + (y - cut_ycntr) * si
    ty = (cut_xcntr - x) * si + (y - cut_ycntr) * co
    R = np.sqrt(tx**2.0 + ty**2.0 / one_minus_eg_sq)
    boxcar = np.reshape(np.ones(sigma * sigma), (sigma, sigma)) 
    I_sigma = pymconvolve.Convolve(CutImDa, boxcar)
    res = CutImDa - I_sigma #the residual image
    if(flag_image):
        # the below will find the image portion which is an anulus of inner 
        # radius 0.3 * eta(.2) and outer radius 1.5 * eta(.2)
        # making the residual value equal to zero inside the 
        # ext_rad/20)
        res[R <= ext_rad * (1 / 20.0)] = 0
        res[R >= ext_rad] = 0
        # sum of positive values of residue
        res_inside_anulus_sum = res[res > 0].sum() 
        # Average inside ext_rad. res_inside_anulus_sum will be divided by 
        # z_inside_R_sum in casgm module
        z_inside_R_sum = CutImDa[R <= ext_rad].sum() / (3.14 * \
                                 ext_rad * ext_rad * np.sqrt(1 - eg**2.0))
        # FIX I dont know why 1/6 instead of 1/20. 
        area = 3.14 * (ext_rad * ext_rad * np.sqrt(1 - eg**2.0)) - \
               3.14 * (ext_rad * ext_rad * (1 / 6.0) * (1 / 6.0) * \
               np.sqrt(1 - eg**2.0))
        # END
        S = res_inside_anulus_sum / area 
        e1sq = CutImDa[res > 0].sum() + I_sigma[res > 0].sum() + \
               4. * CutImDa[res > 0].size * sky
    else:
        res[R >= ext_rad] = 0
        res_inside_anulus_sum = res[np.where(res > 0)].sum()
        area = 3.14 * ext_rad**2.0
        S = res_inside_anulus_sum / area
        z_inside_R_sum = 0 # just to return the value in the end
        e1sq = CutImDa[res > 0].sum() + I_sigma[res > 0].sum() + \
               2. * CutImDa[res > 0].size * sky
    e2sq = res_inside_anulus_sum**2.0
    e3sq = CutImDa[R <= ext_rad].sum() + 2 * CutImDa[R <= ext_rad].size * sky
    e4sq = (CutImDa[R <= ext_rad].sum())**2.0
    if(e2sq!=0):
        error = e1sq / e2sq
    else:
        print "Could not find error"
        error = 0.0
    return S, error, z_inside_R_sum, e3sq, e4sq    
