import os
import sys
import pyfits
from os.path import exists
import numpy as np
from flagfunc import GetFlag, isset, SetFlag
import config as c 
import traceback
from mask_or_fit import *

class ConfigFunc:
    """The class making configuration file for GALFIT."""
    def __init__(self, cutimage, whtimage, xcntr, ycntr, NXPTS, NYPTS, line_s, psffile, sex_cata):
        self.cutimage = cutimage
        self.line_s  = line_s
        self.whtimage = whtimage
        self.xcntr = xcntr
        self.ycntr = ycntr
        self.NXPTS = NXPTS
        self.NYPTS = NYPTS 
        self.psffile = psffile
        self.conff = self.write_conff(line_s, sex_cata, xcntr, ycntr)
        return

    def write_conff(self, line_s, sex_cata, xcntr, ycntr):
        imagefile = c.imagefile
        threshold = c.threshold
        thresh_area = c.thresh_area
        avoidme = c.avoidme

        try:
            ComP = c.components 
        except:
            ComP = ['bulge', 'disk']
        if len(ComP) == 0:
            ComP = ['bulge', 'disk']
            
        target = SEx_obj(self.NXPTS, self.NYPTS, line_s)
        target_imcenter = [target.xcntr, target.ycntr]
        target.set_center(xcntr, ycntr)
        
        outfile = 'O_' + c.fstring + '.fits'
        mask_file = 'M_' + c.fstring + '.fits'
        config_file = 'G_' + c.fstring + '.in'
        constrain_file = c.fstring + '.con'
        
        try:
            c_center_constrain = c.center_constrain
        except:
            c_center_constrain = 2.0
            
        MakeConstrain = 0 if exists(constrain_file) else 1
        confiles = con_G_writer(constrain_file, config_file, MakeConstrain, c.galfitv)

        mag_zero = c.mag_zero

        # Write header image parameters
        confiles.f_G.write('# IMAGE PARAMETERS\n')
        confiles.f_G.writelines(['A) ', c.datadir + str(self.cutimage), '    # Input data image (FITS file)\n'])
        confiles.f_G.writelines(['B) ', str(outfile), '        # Name for the output image\n'])
        confiles.f_G.writelines(['C) ', c.datadir + str(self.whtimage), '        # Noise image name\n'])
        confiles.f_G.writelines(['D) ', c.datadir + str(self.psffile), '            # Input PSF image\n'])
        confiles.f_G.writelines(['E) 1            # PSF oversampling factor\n'])
        confiles.f_G.writelines(['F) ', str(mask_file), '        # Bad pixel mask\n'])
        confiles.f_G.writelines(['G) ', str(constrain_file), '       # File with constraints\n'])
        confiles.f_G.writelines(['H) 1 ', str(self.NXPTS), ' 1 ', str(self.NYPTS), '        # Image region\n'])
        confiles.f_G.writelines(['I) 100 100        # Size of convolution box\n'])
        confiles.f_G.writelines(['J) ', str(mag_zero), '        # Photometric zeropoint\n'])
        confiles.f_G.writelines(['O) regular            # Display type\n'])
        confiles.f_G.writelines(['P) 0            # Create output image only?\n'])
        confiles.f_G.writelines(['S) 0            # Interactively?\n\n\n'])

        for Co in ComP:
            if Co == 'bulge':
                confiles.write_bulge(target)
                c.Flag = SetFlag(c.Flag, GetFlag('FIT_BULGE'))
            elif Co == 'disk':
                confiles.write_disk(target)
                c.Flag = SetFlag(c.Flag, GetFlag('FIT_DISK'))   
            elif Co == 'point':
                confiles.write_point(target)
                c.Flag = SetFlag(c.Flag, GetFlag('FIT_POINT'))
            elif Co == 'bar':
                confiles.write_bar(target)
                c.Flag = SetFlag(c.Flag, GetFlag('FIT_BULGE'))

        confiles.write_sky(target)
        
        if c.center_deviated:
            c.center_deviated = 0

        isneighbour = 0
        for line_j in open(sex_cata, 'r'):
            if line_j[0] != '#':
                neighbor = SEx_obj(self.NXPTS, self.NYPTS, line_j)
                if target.mask_or_fit(neighbor, threshold, thresh_area, avoidme) == 0:
                    isneighbour = 1
                    xn = xcntr - target_imcenter[0] + neighbor.xcntr
                    yn = ycntr - target_imcenter[1] + neighbor.ycntr
                    neighbor.set_center(xn, yn)
                    confiles.write_neighbor(neighbor)

        if isneighbour:
            c.Flag = SetFlag(c.Flag, GetFlag('NEIGHBOUR_FIT'))

        confiles.close_files()

class con_G_writer():
    """This class handles the writing of the configuration file and G file"""
    def __init__(self, confile, Gfile, MakeConstrain, galfitv):
        self.obj_counter = 0
        self.MakeConstrain = MakeConstrain
        self.set_galfitv(galfitv)
        self.set_files(confile, Gfile)
        return

    def set_galfitv(self, galfitv):
        # np.float is deprecated in newer numpy versions
        self.galfitv = float(galfitv.split('.')[0])
        return

    def set_files(self, confile, Gfile):
        self.confile = confile
        self.Gfile = Gfile
        if self.MakeConstrain == 1:
            self.f_constrain = open(confile, 'w')
        self.f_G = open(Gfile, 'w')
        return

    def close_files(self):
        if self.MakeConstrain == 1:
            self.f_constrain.close()
        self.f_G.close()

    def write_bulge(self, target):
        self.obj_counter += 1
        if self.MakeConstrain == 1:
            self.f_constrain.write(f"{self.obj_counter}      n      {c.LN} to {c.UN}\n")
            if c.center_deviated:
                dev = c.center_deviation - c.center_deviation / 4.0
                self.f_constrain.write(f"{self.obj_counter}      x      -{dev}     {dev}\n")
                self.f_constrain.write(f"{self.obj_counter}      y      -{dev}     {dev}\n")
            else:
                self.f_constrain.write(f"{self.obj_counter}      x      {-c.center_constrain}     {c.center_constrain}\n")
                self.f_constrain.write(f"{self.obj_counter}      y      {-c.center_constrain}     {c.center_constrain}\n")

            self.f_constrain.write(f"{self.obj_counter}     mag     {c.UMag} to {c.LMag}\n")
            self.f_constrain.write(f"{self.obj_counter}      re     {c.LRe} to {c.URe}\n")
            self.f_constrain.write(f"{self.obj_counter}      q       0.0 to 1.0\n")
            self.f_constrain.write(f"{self.obj_counter}      pa       -360.0 to 360.0\n")

        self.f_G.write('# Sersic function\n\n')
        self.f_G.writelines([' 0) sersic        # Object type\n'])
        self.f_G.writelines([' 1) ', str(target.xcntr), ' ', str(target.ycntr), ' ', str(c.fitting[0]), ' ', str(c.fitting[0]), '   # position x, y\n'])
        self.f_G.writelines([' 3) ', str(target.mag), ' 1        # total magnitude\n'])
        self.f_G.writelines([' 4) ', str(target.radius), ' 1        # R_e [Pixels]\n'])
        self.f_G.writelines([' 5) 4.0 ', str(int(not c.devauc)), '        # Sersic exponent\n'])
        
        param_idx = ' 9)' if self.galfitv >= 3.0 else ' 8)'
        self.f_G.writelines([param_idx, ' ', str(target.axis_rat), ' 1        # axis ratio\n'])
        param_idx = ' 10)' if self.galfitv >= 3.0 else ' 9)'
        self.f_G.writelines([param_idx, ' ', str(target.pos_ang_galfit), ' 1        # position angle\n'])
        
        if self.galfitv < 3.0:
            box_fit = '1' if (c.bdbox or c.bbox) else '0'
            self.f_G.writelines(['10) 0.0 ', box_fit, '        # diskiness/boxiness\n'])
        self.f_G.writelines([' Z) 0             # output image\n\n\n']) 

    def write_disk(self, target):
        self.obj_counter += 1
        if self.MakeConstrain == 1:
            if c.center_deviated:
                dev = c.center_deviation - c.center_deviation / 4.0
                self.f_constrain.write(f"{self.obj_counter}      x      -{dev}     {dev}\n")
                self.f_constrain.write(f"{self.obj_counter}      y      -{dev}     {dev}\n")
            else:
                self.f_constrain.write(f"{self.obj_counter}       x       {-c.center_constrain}     {c.center_constrain}\n")
                self.f_constrain.write(f"{self.obj_counter}       y       {-c.center_constrain}     {c.center_constrain}\n")
            self.f_constrain.write(f"{self.obj_counter}     mag     {c.UMag} to {c.LMag}\n")
            self.f_constrain.write(f"{self.obj_counter}      rs     {c.LRd} to {c.URd}\n")
            self.f_constrain.write(f"{self.obj_counter}      q       0.0 to 1.0\n")
            self.f_constrain.write(f"{self.obj_counter}      pa       -360.0 to 360.0\n")

        self.f_G.writelines(['# Exponential function\n\n'])
        self.f_G.writelines([' 0) expdisk         # Object type\n'])
        self.f_G.writelines([' 1) ', str(target.xcntr), ' ', str(target.ycntr), ' ', str(c.fitting[1]), ' ', str(c.fitting[1]), '    # position x, y\n'])
        self.f_G.writelines([' 3) ', str(target.mag), ' 1         # total magnitude\n'])
        self.f_G.writelines([' 4) ', str(target.radius), ' 1         # R_e [Pixels]\n'])
        
        if self.galfitv >= 3.0:
            self.f_G.writelines([' 9) ', str(target.axis_rat), ' 1        # axis ratio\n'])
            self.f_G.writelines([' 10) ', str(target.pos_ang_galfit), ' 1        # position angle\n'])
        else:
            self.f_G.writelines([' 8) ', str(target.axis_rat), ' 1        # axis ratio\n'])
            self.f_G.writelines([' 9) ', str(target.pos_ang_galfit), ' 1        # position angle\n'])
            box_fit = '1' if (c.bdbox or c.bbox) else '0'
            self.f_G.writelines(['10) 0.0 ', box_fit, '        # diskiness/boxiness\n'])
        self.f_G.writelines([' Z) 0             # output image\n\n\n']) 

    def write_point(self, target):
        self.obj_counter += 1
        if self.MakeConstrain == 1:
            self.f_constrain.write(f"{self.obj_counter}       x       -2.0      2.0\n")
            self.f_constrain.write(f"{self.obj_counter}       y       -2.0      2.0\n")
        pmag = target.mag + 2.5 * np.log10(6.0)
        self.f_G.writelines(['# point source\n\n'])
        self.f_G.writelines([' 0) psf              # Object type\n'])
        self.f_G.writelines([' 1) ', str(target.xcntr), ' ', str(target.ycntr), ' 1 1  # position x, y\n'])
        self.f_G.writelines([' 3) ', str(pmag), ' 1             # total magnitude\n'])
        self.f_G.writelines([' Z) 0                    # output image\n\n\n'])

    def write_bar(self, target):
        self.obj_counter += 1
        if self.MakeConstrain == 1:
            self.f_constrain.write(f"{self.obj_counter}      n      0.1 to 2.2\n")
            dev = c.center_deviation - c.center_deviation / 4.0
            self.f_constrain.write(f"{self.obj_counter}      x      -{dev}     {dev}\n")
            self.f_constrain.write(f"{self.obj_counter}      y      -{dev}     {dev}\n")
            self.f_constrain.write(f"{self.obj_counter}     mag     {c.UMag} to {c.LMag}\n")
            self.f_constrain.write(f"{self.obj_counter}      re     {c.LRe} to {c.URe}\n")
            self.f_constrain.write(f"{self.obj_counter}      q       0.0 to 0.5\n")
            self.f_constrain.write(f"{self.obj_counter}      pa       -360.0 to 360.0\n")

        barmag = target.mag + 2.5 * np.log10(3.0)
        self.f_G.write('# Sersic function for bar\n\n')
        self.f_G.writelines([' 0) sersic        # Object type\n'])
        self.f_G.writelines([' 1) ', str(target.xcntr), ' ', str(target.ycntr), ' ', str(c.fitting[3]), ' ', str(c.fitting[3]), '   # position x, y\n'])
        self.f_G.writelines([' 3) ', str(barmag), ' 1        # total magnitude\n'])
        self.f_G.writelines([' 4) ', str(target.radius), ' 1        # R_e\n'])
        self.f_G.writelines([' 5) 0.5 1        # Sersic exponent\n'])
        self.f_G.writelines([' 8) 0.3 1    # axis ratio\n'])
        self.f_G.writelines([' 9) ', str(target.pos_ang_galfit), ' 1        # position angle\n'])
        self.f_G.writelines(['10) 0.0 0        # diskiness/boxiness\n'])
        self.f_G.writelines([' Z) 0             # output image\n\n\n']) 

    def write_sky(self, target):
        self.obj_counter += 1
        self.f_G.writelines(['# sky\n\n']) 
        self.f_G.writelines([' 0) sky\n'])
        self.f_G.writelines([' 1) ', str(c.SexSky), ' ', str(c.fitting[2]), '    # sky background\n'])
        self.f_G.writelines([' 2) 0.000      0       # sky gradient in x\n',
                             ' 3) 0.000      0       # sky gradient in y\n',
                             ' Z) 0                  # output image\n\n\n'])

    def write_neighbor(self, target):
        self.obj_counter += 1
        if self.MakeConstrain == 1:
            self.f_constrain.write(f"{self.obj_counter}      n      0.02 to 20.0  \n")
            self.f_constrain.write(f"{self.obj_counter}     mag    -100.0 to 100.0\n")
            self.f_constrain.write(f"{self.obj_counter}      re      0.0 to 500.0\n")
            self.f_constrain.write(f"{self.obj_counter}      q       0.0 to 1.0\n")
            self.f_constrain.write(f"{self.obj_counter}      pa    -360.0 to 360.0\n")

        self.f_G.writelines(['# neighbor sersic function\n\n'])
        self.f_G.writelines([' 0) sersic               # Object type\n'])
        self.f_G.writelines([' 1) ', str(target.xcntr), ' ', str(target.ycntr), ' 1 1  # position x, y\n'])
        self.f_G.writelines([' 3) ', str(target.mag), ' 1         # total magnitude\n'])
        self.f_G.writelines([' 4) ', str(target.radius), ' 1          # R_e\n'])
        self.f_G.writelines([' 5) 4.0 1            # Sersic exponent\n'])
        
        if self.galfitv >= 3.0:
            self.f_G.writelines([' 9) ', str(target.axis_rat), ' 1        # axis ratio\n'])
            self.f_G.writelines([' 10) ', str(target.pos_ang_galfit), ' 1            # position angle\n'])
        else:
            self.f_G.writelines([' 8) ', str(target.axis_rat), ' 1        # axis ratio\n'])
            self.f_G.writelines([' 9) ', str(target.pos_ang_galfit), ' 1            # position angle\n'])
            self.f_G.writelines(['10) 0.0 0             # diskiness/boxiness\n'])
        self.f_G.writelines([' Z) 0                    # output image\n\n\n'])