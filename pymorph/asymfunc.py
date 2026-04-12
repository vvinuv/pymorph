import os
import numpy as np    
import numpy.ma as ma
import pyfits
from rotate import rotate, ImSec
import config as c

class asymmetry:
    """Finding Asymmetry parameter.
       A = min(Sum(abs(I_0 - I_r) / Sum(abs(I_0)) - Noise_Correction
    """        
    def __init__(self, cutimage, maskimage, ini_xcntr, ini_ycntr, pa, eg, r50, extraction_radius, background, angle, flag_image, ABS_ZSUM):
        self.cutimage        = cutimage
        self.maskimage       = maskimage
        self.ini_xcntr       = ini_xcntr 
        self.ini_ycntr       = ini_ycntr
        self.extraction_radius = extraction_radius
        self.pa              = pa
        self.eg              = eg
        self.r50             = r50
        self.background      = background
        self.one_minus_eg_sq = (1.0 - eg)**2.0
        self.flag_image      = flag_image
        self.angle           = angle
        self.ABS_ZSUM        = ABS_ZSUM
        
        # Call the ASYM function and store results
        self.image_asymm = ASYM(self.cutimage, self.maskimage, 
                                self.ini_xcntr, self.ini_ycntr, self.pa, 
                                self.one_minus_eg_sq, r50, self.background, 
                                self.extraction_radius, self.angle, 
                                self.flag_image, self.ABS_ZSUM)
        return 

def ASYM(cutimage, maskimage, ini_xcntr, ini_ycntr, pa, one_minus_eg_sq, r50, background, extraction_radius, angle, flag_image, ABS_ZSUM):
    if flag_image == 0:
        maskimage = 'BackMask.fits'
    
    co = np.cos(pa * np.pi / 180.0)
    si = np.sin(pa * np.pi / 180.0)
    
    Aabs = np.zeros(9, dtype=np.float32)
    rot_sum = np.zeros(9, dtype=np.float32)
    abs_zsum = np.zeros(9, dtype=np.float32)
    sh_sum = np.zeros(9, dtype=np.float32)
    absres_sum = np.zeros(9, dtype=np.float32)
    
    # Load Image
    with pyfits.open(os.path.join(c.datadir, cutimage)) as f:
        z = f[0].data.astype(np.float32)
        z = z - background

    # Load Mask
    with pyfits.open(os.path.join(c.outdir, maskimage)) as fm:
        zm = fm[0].data
        
    # Rotate the mask
    rzm = rotate(zm, 180.0, ini_xcntr, ini_ycntr)
    mask = zm + rzm 
    
    center_rad = 0.01 * r50 
    flag_center = 0 
    nn = 0 
    
    while flag_center == 0:
        flag_out = 0
        xcntr_list = np.reshape(np.array([[ini_xcntr] * 3, 
                                          [ini_xcntr - center_rad] * 3, 
                                          [ini_xcntr + center_rad] * 3]), (9, ))
        ycntr_list = np.array([ini_ycntr, ini_ycntr - center_rad, ini_ycntr + center_rad] * 3)
        
        for it in range(9):
            # Extract Image Section
            CutImDa, cut_xcntr, cut_ycntr, SizeY, SizeX, ymin, ymax, xmin, xmax, f_out = \
                       ImSec(z, xcntr_list[it], ycntr_list[it], extraction_radius)
            
            if flag_image != 0:
                flag_out = f_out
            
            # Create coordinate grid
            y_grid, x_grid = np.mgrid[0:SizeY, 0:SizeX].astype(np.float32)
            
            tx = (x_grid - cut_xcntr) * co + (y_grid - cut_ycntr) * si
            ty = (cut_xcntr - x_grid) * si + (y_grid - cut_ycntr) * co
            R = np.sqrt(tx**2.0 + (ty**2.0 / one_minus_eg_sq))
            
            # Rotate Cutout
            rz = rotate(CutImDa, 180.0, cut_xcntr, cut_ycntr)
            res = CutImDa - rz
            masksub = mask[int(ymin):int(ymax), int(xmin):int(xmax)]

            # Masked calculations
            m_cut = ma.masked_array(CutImDa, masksub)
            m_rz = ma.masked_array(rz, masksub)
            m_res = ma.masked_array(res, masksub)
            
            radial_condition = (R <= extraction_radius)
            
            sh_sum[it] = m_cut[radial_condition].sum()
            rot_sum[it] = m_rz[radial_condition].sum()
            abs_zsum[it] = np.abs(m_cut[radial_condition]).sum()
            absres_sum[it] = np.abs(m_res[radial_condition]).sum()
            
            if flag_image:
                Aabs[it] = absres_sum[it] / abs_zsum[it] if abs_zsum[it] != 0 else 99.0
            else:
                Aabs[it] = absres_sum[it] / ABS_ZSUM if ABS_ZSUM != 0 else 99.0

        # Find minimum asymmetry in the 3x3 grid
        Aabso = Aabs.min()
        index = Aabs.argmin()
        
        if nn > 20:
            flag_center = 1
        elif xcntr_list[index] == ini_xcntr and ycntr_list[index] == ini_ycntr:
            flag_center = 1
            
            # Recalculate final error for the winning center
            final_absres = absres_sum[index]
            final_absz = abs_zsum[index] if flag_image else ABS_ZSUM
            final_sh = sh_sum[index]
            final_rot = rot_sum[index]
            
            # Pixel count for background noise scaling
            pix_count = np.sum(R <= extraction_radius)
            
            error_asym = np.sqrt(Aabso**2 * (
                ((final_sh + final_rot + 4 * background * pix_count) / (final_absres**2 if final_absres != 0 else 1)) +
                ((final_sh + 2.0 * background * pix_count) / (final_absz**2 if final_absz != 0 else 1))
            ))
        else:
            ini_xcntr = xcntr_list[index]
            ini_ycntr = ycntr_list[index]
            nn += 1

    return Aabso, error_asym, ini_xcntr, ini_ycntr, nn, flag_out, abs_zsum[index], sh_sum[index], extraction_radius