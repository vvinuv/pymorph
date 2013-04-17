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
    """The class making configuration file for GALFIT. The configuration file 
       consists of bulge and disk component of the object and only Sersic 
       component for the neighbours, if any. The sky is always fixed and has
       the value of SExtractor. The disk/boxy parameter is also fixed to zero.
       The initial value for Sersic index 'n' is 4.The configuration file has 
       the name G_string(galid).in. The output image has the name 
       O_string(galid).fits"""
    def __init__(self, cutimage, whtimage, xcntr, ycntr, NXPTS, NYPTS, line_s, psffile, sex_cata):
        self.cutimage = cutimage
        self.line_s  = line_s
	self.whtimage = whtimage
        self.xcntr = xcntr
        self.ycntr = ycntr
        self.NXPTS = NXPTS
        self.NYPTS = NYPTS 
        self.psffile = psffile
        self.conff    = self.write_conff(line_s, sex_cata)
        return

    def write_conff(self, line_s, sex_cata):
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
        outfile   = 'O_' + c.fstring + '.fits'
        mask_file = 'M_' + c.fstring + '.fits'
        config_file = 'G_' + c.fstring + '.in' #Name of the GALFIT configuration file
        constrain_file = c.fstring + '.con'
        try:
            c.center_constrain = c.center_constrain
        except:
            c.center_constrain = 2.0
        if exists(constrain_file):
            MakeConstrain = 0
        else:
            MakeConstrain = 1

        confiles = con_G_writer(constrain_file, config_file, MakeConstrain,c.galfitv)


        mag_zero = c.mag_zero #magnitude zero point

        #Write configuration file
        confiles.f_G.write('# IMAGE PARAMETERS\n')
        confiles.f_G.writelines(['A) ', c.datadir + str(self.cutimage), '	# Input data image',\
                      ' (FITS file)\n'])
        confiles.f_G.writelines(['B) ', str(outfile), '		# Name for',\
                      ' the output image\n'])
        confiles.f_G.writelines(['C) ', c.datadir + str(self.whtimage), '		# Noise image name', \
                      ' (made from data if blank or "none")\n'])
        confiles.f_G.writelines(['D) ', c.datadir + str(self.psffile), '			# Input PSF', \
                      ' image for convolution (FITS file)\n'])
        confiles.f_G.writelines(['E) 1			# PSF oversampling factor relative',
                      ' to data\n'])
        confiles.f_G.writelines(['F) ', str(mask_file), '		# Bad pixel',
                      ' mask(FITS image or ASCII coord list)\n'])
        confiles.f_G.writelines(['G) ', str(constrain_file), '       # File with parameter',\
                      ' constraints (ASCII file)\n'])
        confiles.f_G.writelines(['H) 1 ', str(self.NXPTS), ' 1 ', str(self.NYPTS), '		#',\
                      ' Image region to fit (xmin xmax ymin ymax)\n'])
        #confiles.f_G.writelines(['I) ', str(self.NXPTS), ' ', str(self.NYPTS),	'		#',\
        #              ' Size of convolution box (x y)\n'])
        # This really shouldn't be hardcoded!!!!
        confiles.f_G.writelines(['I) ', str(100), ' ', str(100),	'		#',\
                      ' Size of convolution box (x y)\n'])
        confiles.f_G.writelines(['J) ', str(mag_zero), '		# Magnitude',\
                      ' photometric zeropoint\n'])
        confiles.f_G.writelines(['O) regular			# Display type',\
                      ' (regular, curses, both)\n'])
        confiles.f_G.writelines(['P) 0			# Create output image only?',\
                      ' (1=yes; 0=optimize)\n'])
        confiles.f_G.writelines(['S) 0			# Modify/create',\
                     ' objects interactively?\n\n\n'])


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
        for line_j in open(sex_cata,'r'):
            if line_j[0] != '#': #line is not a comment
                try:
                    neighbor = SEx_obj(self.NXPTS, self.NYPTS, line_j)
                    if mask_or_fit(target,neighbor,threshold,thresh_area,avoidme)==0:
                        isneighbour = 1
                        confiles.write_neighbor(neighbor)
                except:
                    pass

        if isneighbour:
            c.Flag  = SetFlag(c.Flag, GetFlag('NEIGHBOUR_FIT'))

        confiles.close_files()
        
    #    f_fit = open('fit2.log','a')
    #    if exists('fit.log'):
    #        os.system('rm fit.log')
    #Here the user should tell the location of the GALFIT excutable
    #    os.system('/Vstr/vstr/vvinuv/galfit/modified/galfit "' + config_file + '"')
    #    if exists('fit.log'):
    #        for line in open('fit.log','r'):
    #            f_fit.writelines([str(line)])
    #    f_fit.close()

class con_G_writer():
    """This class handles the writing of the configuration file and G file"""
    def __init__(self, confile, Gfile, MakeConstrain, galfitv):
        self.obj_counter = 0
        self.MakeConstrain = MakeConstrain
        self.set_files(confile, Gfile)
        self.set_galfitv(galfitv)
        return

    def set_galfitv(self, galfitv):
        self.galfitv = np.float(galfitv.split('.')[0])
        return

    def set_files(self, confile, Gfile):
        """Sets file names and opens those files"""
        try:
            self.f_constrain.close()
        except AttributeError, NameError:
            pass
        try:
            self.f_G.close()
        except AttributeError, NameError:
            pass

        self.confile = confile
        self.Gfile = Gfile

        if self.MakeConstrain == 1:
            self.f_constrain = open(confile, 'w')

        self.f_G = open(Gfile,'w')

        return

    def close_files(self):
        """Closes configuration file and G file"""
        if self.MakeConstrain == 1:
            self.f_constrain.close()
        self.f_G.close()
        return

    def write_bulge(self, target):
        self.obj_counter+=1

        if self.MakeConstrain == 1:
            # write constraints
            self.f_constrain.write(str(self.obj_counter) + '      n      ' + str(c.LN) + \
                                   ' to ' + str(c.UN) +  '\n')
            if c.center_deviated:
                self.f_constrain.write(str(self.obj_counter) + '      x      -' + \
                                       str(c.center_deviation - c.center_deviation / 4.0) + \
                                  '     ' + \
                                       str(c.center_deviation - c.center_deviation / 4.0) + '\n')
                self.f_constrain.write(str(self.obj_counter) + '      y      -' + \
                                       str(c.center_deviation - c.center_deviation / 4.0) + \
                                       '     ' + \
                                       str(c.center_deviation - c.center_deviation / 4.0) + '\n')
            else:
                self.f_constrain.write(str(self.obj_counter) + '      x      ' + str(-c.center_constrain) + '     ' + str(c.center_constrain) + '\n')
                self.f_constrain.write(str(self.obj_counter) + '      y      ' + str(-c.center_constrain) + '     ' + str(c.center_constrain) + '\n')

            self.f_constrain.write(str(self.obj_counter) + '     mag     ' + str(c.UMag) + \
                              ' to ' + str(c.LMag) + '\n')
            self.f_constrain.write(str(self.obj_counter) + '      re     ' + str(c.LRe) +\
                              ' to ' + str(c.URe) + '\n')
            self.f_constrain.write(str(self.obj_counter) + '      q       0.0 to 1.0\n')
            self.f_constrain.write(str(self.obj_counter) + '      pa       -360.0 to 360.0\n')

        # write config file
        self.f_G.write('# Sersic function\n\n')
        self.f_G.writelines([' 0) sersic		# Object type\n'])
        self.f_G.writelines([' 1) ', str(target.xcntr), ' ', str(target.ycntr),' ', \
                             str(c.fitting[0]), ' ', str(c.fitting[0]), '   #',\
                             ' position x, y [pixel]\n'])
        self.f_G.writelines([' 3) ', str(target.mag), ' 1		# total magnitude\n'])
        self.f_G.writelines([' 4) ', str(target.radius), ' 1		# R_e [Pixels]\n'])
        self.f_G.writelines([' 5) 4.0 ', str(int(not c.devauc)) ,\
                             '		#Sersic exponent (deVauc=4, expdisk=1)\n'])
        if self.galfitv >= 3.0:
            self.f_G.writelines([' 9) ', str(target.axis_rat), ' 1		', \
                                 '# axis ratio (b/a)\n'])
            self.f_G.writelines([' 10) ', str(target.pos_ang_galfit), ' 1		',\
                                 '# position angle (PA)',\
                                 '[Degrees: Up=0, Left=90]\n'])
        else:
            self.f_G.writelines([' 8) ', str(target.axis_rat), ' 1		',\
                                 '# axis ratio (b/a)\n'])
            self.f_G.writelines([' 9) ', str(target.pos_ang_galfit), ' 1		',\
                                 '# position angle (PA)',\
                                 '[Degrees: Up=0, Left=90]\n'])
            if c.bdbox or c.bbox:
                self.f_G.writelines(['10) 0.0 1		# diskiness (< 0) or ' \
                                     'boxiness (> 0)\n'])
            else:
                self.f_G.writelines(['10) 0.0 0            # diskiness (< 0) or ' \
                              'boxiness (> 0)\n'])
        self.f_G.writelines([' Z) 0 			# output image',\
                             ' (see above)\n\n\n']) 

        return
    

    def write_disk(self, target):
        self.obj_counter+=1

        if self.MakeConstrain == 1:
            # write constraints
            if c.center_deviated:
                self.f_constrain.write(str(self.obj_counter) + '      x      -' + \
                                  str(c.center_deviation - c.center_deviation / 4.0) + \
                                  '     ' + \
                                  str(c.center_deviation - c.center_deviation / 4.0) + '\n')
                self.f_constrain.write(str(self.obj_counter) + '      y      -' + \
                                  str(c.center_deviation - c.center_deviation / 4.0) + \
                                  '     ' + \
                                  str(c.center_deviation - c.center_deviation / 4.0) + '\n')
            else:
                self.f_constrain.write(str(self.obj_counter) + '       x       ' + str(-c.center_constrain) + '     ' + str(c.center_constrain) + '\n')
                self.f_constrain.write(str(self.obj_counter) + '       y       ' + str(-c.center_constrain) + '     ' + str(c.center_constrain) + '\n')
            self.f_constrain.write(str(self.obj_counter) + '     mag     ' + str(c.UMag) + \
                              ' to ' + str(c.LMag) + '\n')
            self.f_constrain.write(str(self.obj_counter) + '      rs     ' + str(c.LRd) + \
                              ' to ' + str(c.URd) + '\n')
            self.f_constrain.write(str(self.obj_counter) + '      q       0.0 to 1.0\n')
            self.f_constrain.write(str(self.obj_counter) + '      pa       -360.0 to 360.0\n')

        # write config file
        self.f_G.writelines(['# Exponential function\n\n'])
        self.f_G.writelines([' 0) expdisk 		# Object type\n'])
        self.f_G.writelines([' 1) ', str(target.xcntr), ' ', str(target.ycntr),' ', \
                      str(c.fitting[1]), ' ', str(c.fitting[1]), '    #',\
                      ' position x, y [pixel]\n'])
        self.f_G.writelines([' 3) ', str(target.mag), ' 1     	# total magnitude\n'])
        self.f_G.writelines([' 4) ', str(target.radius), ' 1 		# R_e [Pixels]\n'])
        if self.galfitv >= 3.0:
            self.f_G.writelines([' 9) ', str(target.axis_rat), ' 1		', \
                          '# axis ratio (b/a)\n'])
            self.f_G.writelines([' 10) ', str(target.pos_ang_galfit), ' 1		',\
                          '# position angle (PA)',\
                          '[Degrees: Up=0, Left=90]\n'])
            self.f_G.writelines([' Z) 0 			# output image',\
                              ' (see above)\n\n\n']) 
        else:
            self.f_G.writelines([' 8) ', str(target.axis_rat), ' 1		',\
                          '# axis ratio (b/a)\n'])
            self.f_G.writelines([' 9) ', str(target.pos_ang_galfit), ' 1		',\
                          '# position angle (PA)',\
                          '[Degrees: Up=0, Left=90]\n'])
            if c.bdbox or c.bbox:
                self.f_G.writelines(['10) 0.0 1		# diskiness (< 0) or ' \
                              'boxiness (> 0)\n'])
            else:
                self.f_G.writelines(['10) 0.0 0            # diskiness (< 0) or ' \
                                  'boxiness (> 0)\n'])
        self.f_G.writelines([' Z) 0 			# output image',\
                             ' (see above)\n\n\n']) 

        return
    
                         
    def write_point(self, target):
        self.obj_counter+=1

        if self.MakeConstrain == 1:
            # write constraints
            self.f_constrain.write(str(self.obj_counter) + '       x       -2.0      2.0\n')
            self.f_constrain.write(str(self.obj_counter) + '       y       -2.0      2.0\n')
            
        # write config
        pmag = target.mag + 2.5 * np.log10(6.0)
        self.f_G.writelines(['#point source\n\n'])
        self.f_G.writelines([' 0) psf              # Object type\n'])
        self.f_G.writelines([' 1) ', str(target.xcntr), ' ', str(target.ycntr),' 1 1  #',\
                             ' position x, y [pixel]\n'])
        self.f_G.writelines([' 3) ', str(pmag), ' 1             # total magnitude\n'])
        self.f_G.writelines([' Z) 0                    # output image '\
                             '(see above)\n\n\n'])

        return
    
    def write_gaussian(self, target):
        self.obj_counter+=1

        if self.MakeConstrain == 1:
            # write constraints
            self.f_constrain.write(str(self.obj_counter) + '       x       -2.0      2.0\n')
            self.f_constrain.write(str(self.obj_counter) + '       y       -2.0      2.0\n')
            

        # write config
        gmag = target.mag + 2.5 * np.log10(2.0)
        self.f_G.writelines(['# Gaussian function\n\n'])
        self.f_G.writelines([' 0) gaussian              # Object type\n'])
        self.f_G.writelines([' 1) ', str(target.xcntr), ' ', str(target.ycntr),' 1 1  #',\
                      ' position x, y [pixel]\n'])
        self.f_G.writelines([' 3) ', str(gmag), ' 1             # total magnitude\n'])
        self.f_G.writelines([' 4) 0.50 0             #FWHM\n'])
        self.f_G.writelines([' 8) 1 0     # axis ratio (b/a)\n'])
        self.f_G.writelines([' 9) 10.0 0                 # position '\
                      'angle(PA) [Degrees: Up=0, Left=90]\n'])
        self.f_G.writelines(['10) 0.0 0                # diskiness (< 0) or '\
                      'boxiness (> 0)\n'])
        self.f_G.writelines([' Z) 0                    # output image '\
                      '(see above)\n\n\n'])

        return

    def write_bar(self, target):
        self.obj_counter+=1
        
        if self.MakeConstrain == 1:
            # write constraints
            self.f_constrain.write(str(self.obj_counter) + '      n      ' + str('0.1') + \
                                   ' to ' + str('2.2') +  '\n')
            if c.center_deviated:
                self.f_constrain.write(str(self.obj_counter) + '      x      -' + \
                                       str(c.center_deviation - c.center_deviation / 4.0) + \
                                       '     ' + \
                                       str(c.center_deviation - c.center_deviation / 4.0) + '\n')
                self.f_constrain.write(str(self.obj_counter) + '      y      -' + \
                                       str(c.center_deviation - c.center_deviation / 4.0) + \
                                       '     ' + \
                                       str(c.center_deviation - c.center_deviation / 4.0) + '\n')
            else:
                self.f_constrain.write(str(self.obj_counter) + '      x      ' + str(-c.center_constrain) + '     ' + str(c.center_constrain) + '\n')
                self.f_constrain.write(str(self.obj_counter) + '      y      ' + str(-c.center_constrain) + '     ' + str(c.center_constrain) + '\n')
            self.f_constrain.write(str(self.obj_counter) + '     mag     ' + str(c.UMag) + \
                                   ' to ' + str(c.LMag) + '\n')
            self.f_constrain.write(str(self.obj_counter) + '      re     ' + str(c.LRe) +\
                                   ' to ' + str(c.URe) + '\n')
            self.f_constrain.write(str(self.obj_counter) + '      q       0.0 to 0.5\n')
            self.f_constrain.write(str(self.obj_counter) + '      pa       -360.0 to 360.0\n')

        barmag = target.mag + 2.5 * np.log10(3.0)
        f.write('# Sersic function for bar\n\n')
        f.writelines([' 0) sersic		# Object type\n'])
        f.writelines([' 1) ', str(target.xcntr), ' ', str(target.ycntr),' ', \
                      str(c.fitting[3]), ' ', str(c.fitting[3]), '   #',\
                      ' position x, y [pixel]\n'])
        f.writelines([' 3) ', str(barmag), ' 1		# total magnitude\n'])
        f.writelines([' 4) ', str(target.radius), ' 1		# R_e [Pixels]\n'])
        f.writelines([' 5) 0.5 1		#Sersic exponent',\
                      ' (deVauc=4, expdisk=1)\n'])
        f.writelines([' 8) ', str('0.3'), ' 1	# axis ratio (b/a)\n'])
        f.writelines([' 9) ', str(target.pos_ang_galfit), ' 1		# position angle (PA)',\
                          '[Degrees: Up=0, Left=90]\n'])
        f.writelines(['10) 0.0 0		# diskiness (< 0) or ' \
                      'boxiness (> 0)\n'])
        f.writelines([' Z) 0 			# output image',\
                      ' (see above)\n\n\n']) 

        return
    
    def write_sky(self, target):
        self.obj_counter+=1
            
        self.f_G.writelines(['# sky\n\n']) 
        self.f_G.writelines([' 0) sky\n'])#str(2.22604*(1.0+0.01)), ' ', str(0.0)
    #    self.f_G.writelines([' 1) ', str(2.41154*(1.0+0.01)),' ', str(c.fitting[2]), \
    #                  '	# sky background [ADU counts\n'])
        self.f_G.writelines([' 1) ', str(c.SexSky),' ', str(c.fitting[2]), \
                      '	# sky background [ADU counts\n'])
    #    self.f_G.writelines([' 1) ', str(0.0), ' ', str(c.fitting[2]), \
    #                  '	# sky background [ADU counts\n'])
    #    self.f_G.writelines([' 1) ', str(c.SkyMin), '      ', str(c.fitting[2]), \
    #                  '	# sky background [ADU counts\n'])
        self.f_G.writelines([' 2) 0.000      0       # dsky/dx (sky gradient in x)\n',\
                      ' 3) 0.000      0       # dsky/dy (sky gradient in y)\n',\
                      ' Z) 0                  # output image\n\n\n'])

        return


    def write_neighbor(self, target):
        self.obj_counter+=1
        
        if self.MakeConstrain == 1:
            # write constraints
            self.f_constrain.write(str(self.obj_counter) + '      n      0.02 to 20.0  \n')
            self.f_constrain.write(str(self.obj_counter) + '     mag    -100.0 to 100.0\n')
            self.f_constrain.write(str(self.obj_counter) + '      re      0.0 to 500.0\n')
            self.f_constrain.write(str(self.obj_counter) + '      q       0.0 to 1.0\n')
            self.f_constrain.write(str(self.obj_counter) + '      pa    -360.0 to 360.0\n')

        # write config
        self.f_G.writelines(['# neighbor sersic function\n\n'])
        self.f_G.writelines([' 0) sersic               # Object type\n'])
        self.f_G.writelines([' 1) ', str(target.xcntr), ' ', str(target.ycntr), \
                      ' 1 1  # position x, y [pixel]\n'])
        self.f_G.writelines([' 3) ', str(target.mag), ' 1     	#',\
                      ' total magnitude\n'])
        self.f_G.writelines([' 4) ', str(target.radius), ' 1  		',\
                      '# R_e [Pixels]\n'])
        self.f_G.writelines([' 5) 4.0 1        	#Sersic exponent', \
                      ' (deVauc=4, expdisk=1)\n'])
        if self.galfitv >= 3.0:
            self.f_G.writelines([' 9) ', str(target.axis_rat), ' 1        # axis',\
                          ' ratio (b/a)\n'])
            self.f_G.writelines([' 10) ', str(target.pos_ang_galfit), ' 1 	       ',\
                          ' # position angle (PA)  [Degrees: Up=0,'\
                          ' Left=90]\n'])
        else:
            self.f_G.writelines([' 8) ', str(target.axis_rat), ' 1        # axis',\
                          ' ratio (b/a)\n'])
            self.f_G.writelines([' 9) ', str(target.pos_ang_galfit), ' 1 	       ',\
                          ' # position angle (PA)  [Degrees: Up=0,'\
                          ' Left=90]\n'])
            self.f_G.writelines(['10) 0.0 0         	# diskiness',\
                          ' (< 0) or boxiness (> 0)\n'])
        self.f_G.writelines([' Z) 0 	           	# output',\
                      ' image (see above)\n\n\n'])

        return
    
