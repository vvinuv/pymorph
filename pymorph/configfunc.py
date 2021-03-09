import os
import sys
import fitsio
import numpy as np
from flagfunc import GetFlag, isset, SetFlag
import traceback
import mask_or_fit as mf

class GalfitConfigFunc:
    
    """

    The class making configuration file for GALFIT. The configuration file 
    consists of bulge and disk component of the object and only Sersic 
    component for the neighbours, if any. The sky is always fixed and has
    the value of SExtractor. The disk/boxy parameter is also fixed to zero.
    The initial value for Sersic index 'n' is 4.The configuration file has 
    the name G_string(galid).in. The output image has the name 
    O_string(galid).fits

    """

    def __init__(self, datadir, cutimage, whtimage,  
                 xcntr, ycntr, NXPTS, NYPTS, 
                 components, fitting, psffile, line_s, sex_cata, skyval, 
                 mag_zero, flag):

        self.datadir = datadir
        self.cutimage = cutimage
        self.whtimage = whtimage
        self.xcntr = xcntr
        self.ycntr = ycntr
        self.NXPTS = NXPTS
        self.NYPTS = NYPTS 
        self.components = components
        self.fitting = fitting
        self.psffile = psffile
        self.line_s  = line_s
        self.sex_cata = sex_cata
        self.skyval = skyval
        self.mag_zero = mag_zero
        self.flag = flag

    def write_config(self, fstring, threshold, thresh_area, 
                    center_deviated, galfitv, center_constrain=2.0,
                    avoidme=0, 
                    LMag=500., UMag=-500.,
                    LN=0.1, UN=20.,
                    LRe=0., URe=500.,
                    LRd=0., URd=500.,
                    bdbox=False, bbox=False, dbox=False, devauc=False):

        target = mf.GetSExObj(self.NXPTS, self.NYPTS, self.line_s)
        target_imcenter =[target.xcntr, target.ycntr]
        target.set_center(self.xcntr, self.ycntr)
        
        config_file = 'G_{}.in'.format(fstring) #Name of the GALFIT configuration file
        constrain_file = '{}.con'.format(fstring)
        self.mask_file = 'M_{}.fits'.format(fstring) 
        outfile   = 'O_{}.fits'.format(fstring)

        if os.path.exists(constrain_file):
            make_constrain = 0
        else:
            make_constrain = 1

        gconfig = GConfigWrite(config_file, constrain_file, self.fitting,
                                make_constrain, center_deviated, 
                                center_constrain, self.skyval,  galfitv,
                                LMag=500., UMag=-500.,
                                LN=0.1, UN=20.,
                                LRe=0., URe=500.,
                                LRd=0., URd=500.,
                                bdbox=False, bbox=False, 
                                dbox=False, devauc=False)


        #Write configuration file
        gconfig.f_G.write('# IMAGE PARAMETERS\n')
        gconfig.f_G.writelines(['A) ', os.path.join(self.datadir, self.cutimage), '	# Input data image',\
                      ' (FITS file)\n'])
        gconfig.f_G.writelines(['B) ', str(outfile), '		# Name for',\
                      ' the output image\n'])
        gconfig.f_G.writelines(['C) ', os.path.join(self.datadir, self.whtimage), '		# Noise image name', \
                      ' (made from data if blank or "none")\n'])
        gconfig.f_G.writelines(['D) ', os.path.join(self.datadir, self.psffile), '			# Input PSF', \
                      ' image for convolution (FITS file)\n'])
        gconfig.f_G.writelines(['E) 1			# PSF oversampling factor relative',
                      ' to data\n'])
        gconfig.f_G.writelines(['F) ', str(self.mask_file), '		# Bad pixel',
                      ' mask(FITS image or ASCII coord list)\n'])
        gconfig.f_G.writelines(['G) ', str(constrain_file), '       # File with parameter',\
                      ' constraints (ASCII file)\n'])
        gconfig.f_G.writelines(['H) 1 ', str(self.NXPTS), ' 1 ', str(self.NYPTS), '		#',\
                      ' Image region to fit (xmin xmax ymin ymax)\n'])
        #gconfig.f_G.writelines(['I) ', str(self.NXPTS), ' ', str(self.NYPTS),	'		#',\
        #              ' Size of convolution box (x y)\n'])
        # This really shouldn't be hardcoded!!!!
        gconfig.f_G.writelines(['I) ', str(100), ' ', str(100),	'		#',\
                      ' Size of convolution box (x y)\n'])
        gconfig.f_G.writelines(['J) ', str(self.mag_zero), '		# Magnitude',\
                      ' photometric zeropoint\n'])
        gconfig.f_G.writelines(['O) regular			# Display type',\
                      ' (regular, curses, both)\n'])
        gconfig.f_G.writelines(['P) 0			# Create output image only?',\
                      ' (1=yes; 0=optimize)\n'])
        gconfig.f_G.writelines(['S) 0			# Modify/create',\
                     ' objects interactively?\n\n\n'])


        for comp in self.components:
            if comp == 'bulge':
                gconfig.write_bulge(target)
                self.flag = SetFlag(self.flag, GetFlag('FIT_BULGE'))
            elif comp == 'disk':
                gconfig.write_disk(target)
                self.flag = SetFlag(self.flag, GetFlag('FIT_DISK'))   
            elif comp == 'point':
                gconfig.write_point(target)
                self.flag = SetFlag(self.flag, GetFlag('FIT_POINT'))
            elif comp == 'bar':
                gconfig.write_bar(target)
                self.flag = SetFlag(self.flag, GetFlag('FIT_BAR'))

        gconfig.write_sky(target)
        
        if center_deviated:
            center_deviated = 0

        isneighbour = 0
        for line_j in open(self.sex_cata,'r'):
            neighbor = mf.GetSExObj(self.NXPTS, self.NYPTS, line_j)
            if target.get_mask(neighbor, threshold, thresh_area,
                                  avoidme) == 0:
                isneighbour = 1
                # recenter in chip coordinates
                xn = self.xcntr - target_imcenter[0] + neighbor.xcntr
                yn = self.xcntr - target_imcenter[1] + neighbor.ycntr
                neighbor.set_center(xn, yn)

                gconfig.write_neighbor(neighbor)

        if isneighbour:
            self.flag  = SetFlag(self.flag, GetFlag('NEIGHBOUR_FIT'))

        gconfig.close_files()
        
    #    f_fit = open('fit2.log','a')
    #    if exists('fit.log'):
    #        os.system('rm fit.log')
    #Here the user should tell the location of the GALFIT excutable
    #    os.system('/Vstr/vstr/vvinuv/galfit/modified/galfit "' + config_file + '"')
    #    if exists('fit.log'):
    #        for line in open('fit.log','r'):
    #            f_fit.writelines([str(line)])
    #    f_fit.close()

class GConfigWrite():
    """

    This class handles the writing of the configuration file and G file

    """

    def __init__(self, config_file, constrain_file, fitting,
                 make_constrain, center_deviated,
                 center_constrain, skyval,  galfitv,
                 LMag=500., UMag=-500., 
                 LN=0.1, UN=20., 
                 LRe=0., URe=500., 
                 LRd=0., URd=500.,
                 bdbox=False, bbox=False, dbox=False, devauc=False
                 ):

        self.obj_counter = 0
        self.make_constrain = make_constrain
        self.set_files(config_file, constrain_file)
        self.fitting = fitting
        self.make_constrain = center_constrain
        self.center_deviated = center_deviated
        self.center_constrain = center_constrain
        self.skyval = skyval
        self.set_galfitv(galfitv)
        
        self.LMag = LMag
        self.UMag = UMag
        self.LN = LN
        self.UN = UN
        self.LRe = LRe
        self.URe = URe
        self.LRd = LRd
        self.URd = URd
        self.bdbox = bdbox
        self.bbox = bbox
        self.dbox = dbox
        self.devauc = devauc


    def set_galfitv(self, galfitv):
        self.galfitv = np.float(galfitv.split('.')[0])

    def set_files(self, config_file, constrain_file):
        """Sets file names and opens those files"""
        try:
            self.f_constrain.close()
        except AttributeError as e:
            print('f_constrain', e)

        try:
            self.f_g.close()
        except AttributeError as e:
            print('f_G', e)

        if self.make_constrain == 1:
            self.f_constrain = open(constrain_file, 'w')

        self.f_G = open(config_file, 'w')

    print(1)
    def close_files(self):
        """Closes configuration file and G file"""
        if self.make_constrain == 1:
            self.f_constrain.close()
        self.f_G.close()

    def write_bulge(self, target):
        self.obj_counter+=1

        if self.make_constrain == 1:
            # write constraints
            sline = '{}     n       {}      to      {}\n'.format(
                    self.obj_counter, self.LN, self.UN)
            self.f_constrain.write(sline)

            if self.center_deviated:
                self.f_constrain.write(str(self.obj_counter) + '      x      -' + \
                                       str(self.center_deviation - self.center_deviation / 4.0) + \
                                  '     ' + \
                                       str(self.center_deviation - self.center_deviation / 4.0) + '\n')
                self.f_constrain.write(str(self.obj_counter) + '      y      -' + \
                                       str(self.center_deviation - self.center_deviation / 4.0) + \
                                       '     ' + \
                                       str(self.center_deviation - self.center_deviation / 4.0) + '\n')
            else:
                self.f_constrain.write(str(self.obj_counter) + '      x      ' + str(-self.center_constrain) + '     ' + str(self.center_constrain) + '\n')
                self.f_constrain.write(str(self.obj_counter) + '      y      ' + str(-self.center_constrain) + '     ' + str(self.center_constrain) + '\n')

            self.f_constrain.write(str(self.obj_counter) + '     mag     ' + str(self.UMag) + \
                              ' to ' + str(self.LMag) + '\n')
            self.f_constrain.write(str(self.obj_counter) + '      re     ' + str(self.LRe) +\
                              ' to ' + str(self.URe) + '\n')
            self.f_constrain.write(str(self.obj_counter) + '      q       0.0 to 1.0\n')
            self.f_constrain.write(str(self.obj_counter) + '      pa       -360.0 to 360.0\n')

        # write config file
        self.f_G.write('# Sersic function\n\n')
        self.f_G.writelines([' 0) sersic		# Object type\n'])
        self.f_G.writelines([' 1) ', str(target.xcntr), ' ', str(target.ycntr),' ', \
                             str(self.fitting[0]), ' ', str(self.fitting[0]), '   #',\
                             ' position x, y [pixel]\n'])
        self.f_G.writelines([' 3) ', str(target.mag), ' 1		# total magnitude\n'])
        self.f_G.writelines([' 4) ', str(target.radius), ' 1		# R_e [Pixels]\n'])
        self.f_G.writelines([' 5) 4.0 ', str(int(not self.devauc)) ,\
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
            if self.bdbox or self.bbox:
                self.f_G.writelines(['10) 0.0 1		# diskiness (< 0) or ' \
                                     'boxiness (> 0)\n'])
            else:
                self.f_G.writelines(['10) 0.0 0            # diskiness (< 0) or ' \
                              'boxiness (> 0)\n'])
        self.f_G.writelines([' Z) 0 			# output image',\
                             ' (see above)\n\n\n']) 

    

    def write_disk(self, target):
        self.obj_counter+=1

        if self.make_constrain == 1:
            # write constraints
            if self.center_deviated:
                self.f_constrain.write(str(self.obj_counter) + '      x      -' + \
                                  str(self.center_deviation - self.center_deviation / 4.0) + \
                                  '     ' + \
                                  str(self.center_deviation - self.center_deviation / 4.0) + '\n')
                self.f_constrain.write(str(self.obj_counter) + '      y      -' + \
                                  str(self.center_deviation - self.center_deviation / 4.0) + \
                                  '     ' + \
                                  str(self.center_deviation - self.center_deviation / 4.0) + '\n')
            else:
                self.f_constrain.write(str(self.obj_counter) + '       x       ' + str(-self.center_constrain) + '     ' + str(self.center_constrain) + '\n')
                self.f_constrain.write(str(self.obj_counter) + '       y       ' + str(-self.center_constrain) + '     ' + str(self.center_constrain) + '\n')
            self.f_constrain.write(str(self.obj_counter) + '     mag     ' + str(self.UMag) + \
                              ' to ' + str(self.LMag) + '\n')
            self.f_constrain.write(str(self.obj_counter) + '      rs     ' + str(self.LRd) + \
                              ' to ' + str(self.URd) + '\n')
            self.f_constrain.write(str(self.obj_counter) + '      q       0.0 to 1.0\n')
            self.f_constrain.write(str(self.obj_counter) + '      pa       -360.0 to 360.0\n')

        # write config file
        self.f_G.writelines(['# Exponential function\n\n'])
        self.f_G.writelines([' 0) expdisk 		# Object type\n'])
        self.f_G.writelines([' 1) ', str(target.xcntr), ' ', str(target.ycntr),' ', \
                      str(self.fitting[1]), ' ', str(self.fitting[1]), '    #',\
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
            if self.bdbox or self.bbox:
                self.f_G.writelines(['10) 0.0 1		# diskiness (< 0) or ' \
                              'boxiness (> 0)\n'])
            else:
                self.f_G.writelines(['10) 0.0 0            # diskiness (< 0) or ' \
                                  'boxiness (> 0)\n'])
        self.f_G.writelines([' Z) 0 			# output image',\
                             ' (see above)\n\n\n']) 

    
                         
    def write_point(self, target):
        self.obj_counter+=1

        if self.make_constrain == 1:
            # write constraints
            self.f_constrain.write(str(self.obj_counter) + '       x       -2.0      2.0\n')
            self.f_constrain.write(str(self.obj_counter) + '       y       -2.0      2.0\n')
            
        # write config
        pmag = target.mag + 2.5 * np.log10(6.0)
        self.f_G.writelines(['#point source\n\n'])
        self.f_G.writelines([' 0) psf              # Object type\n'])
        self.f_G.writelines([' 1) ', str(target.xcntr), ' ', str(target.ycntr),' ', \
                               str(self.fitting[4]), ' ', str(self.fitting[4]), '    #',\
                      ' position x, y [pixel]\n'])
        self.f_G.writelines([' 3) ', str(pmag), ' 1             # total magnitude\n'])
        self.f_G.writelines([' Z) 0                    # output image '\
                             '(see above)\n\n\n'])

    
    def write_gaussian(self, target):
        self.obj_counter+=1

        if self.make_constrain == 1:
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


    def write_bar(self, target):
        self.obj_counter+=1
        
        if self.make_constrain == 1:
            # write constraints
            self.f_constrain.write(str(self.obj_counter) + '      n      ' + str('0.1') + \
                                   ' to ' + str('2.2') +  '\n')
            if self.center_deviated:
                self.f_constrain.write(str(self.obj_counter) + '      x      -' + \
                                       str(self.center_deviation - self.center_deviation / 4.0) + \
                                       '     ' + \
                                       str(self.center_deviation - self.center_deviation / 4.0) + '\n')
                self.f_constrain.write(str(self.obj_counter) + '      y      -' + \
                                       str(self.center_deviation - self.center_deviation / 4.0) + \
                                       '     ' + \
                                       str(self.center_deviation - self.center_deviation / 4.0) + '\n')
            else:
                self.f_constrain.write(str(self.obj_counter) + '      x      ' + str(-self.center_constrain) + '     ' + str(self.center_constrain) + '\n')
                self.f_constrain.write(str(self.obj_counter) + '      y      ' + str(-self.center_constrain) + '     ' + str(self.center_constrain) + '\n')
            self.f_constrain.write(str(self.obj_counter) + '     mag     ' + str(self.UMag) + \
                                   ' to ' + str(self.LMag) + '\n')
            self.f_constrain.write(str(self.obj_counter) + '      re     ' + str(self.LRe) +\
                                   ' to ' + str(self.URe) + '\n')
            self.f_constrain.write(str(self.obj_counter) + '      q       0.0 to 0.5\n')
            self.f_constrain.write(str(self.obj_counter) + '      pa       -360.0 to 360.0\n')

        barmag = target.mag + 2.5 * np.log10(3.0)
        self.f_G.write('# Sersic function for bar\n\n')
        self.f_G.writelines([' 0) sersic		# Object type\n'])
        self.f_G.writelines([' 1) ', str(target.xcntr), ' ', str(target.ycntr),' ', \
                      str(self.fitting[3]), ' ', str(self.fitting[3]), '   #',\
                      ' position x, y [pixel]\n'])
        self.f_G.writelines([' 3) ', str(barmag), ' 1		# total magnitude\n'])
        self.f_G.writelines([' 4) ', str(target.radius), ' 1		# R_e [Pixels]\n'])
        self.f_G.writelines([' 5) 0.5 1		#Sersic exponent',\
                      ' (deVauc=4, expdisk=1)\n'])
        self.f_G.writelines([' 8) ', str('0.3'), ' 1	# axis ratio (b/a)\n'])
        self.f_G.writelines([' 9) ', str(target.pos_ang_galfit), ' 1		# position angle (PA)',\
                          '[Degrees: Up=0, Left=90]\n'])
        self.f_G.writelines(['10) 0.0 0		# diskiness (< 0) or ' \
                      'boxiness (> 0)\n'])
        self.f_G.writelines([' Z) 0 			# output image',\
                      ' (see above)\n\n\n']) 

    
    def write_sky(self, target):
        self.obj_counter+=1
            
        self.f_G.writelines(['# sky\n\n']) 
        self.f_G.writelines([' 0) sky\n'])
        self.f_G.writelines([' 1) ', str(self.skyval),' ', str(self.fitting[2]), \
                      '	# sky background [ADU counts\n'])
        self.f_G.writelines([' 2) 0.000      0       # dsky/dx (sky gradient in x)\n',\
                      ' 3) 0.000      0       # dsky/dy (sky gradient in y)\n',\
                      ' Z) 0                  # output image\n\n\n'])



    def write_neighbor(self, target):
        self.obj_counter+=1
        
        if self.make_constrain == 1:
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

    
