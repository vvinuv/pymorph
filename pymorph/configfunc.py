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

    def __init__(self, datadir, 
                 cutimage, whtimage,  
                 xcntr_img, ycntr_img,
                 good_object, 
                 half_size, 
                 components, fitting, 
                 psffile, 
                 sex_cata_neighbor, 
                 SexSky,
                 fstring,
                 threshold, thresh_area,
                 center_deviated, center_constrain,
                 avoidme,
                 LMag, UMag,
                 LN, UN,
                 LRadius, URadius,
                 bdbox, bbox, dbox, devauc,
                 galfitv,
                 mag_zero, flag):

        self.datadir = datadir
        self.cutimage = cutimage
        self.whtimage = whtimage
        self.xcntr_img = xcntr_img
        self.ycntr_img = ycntr_img
        self.half_size = half_size
        self.components = components
        self.fitting = fitting
        self.psffile = psffile
        self.good_object  = good_object
        #print(self.good_object.shape)
        self.sex_cata_neighbor = sex_cata_neighbor
        self.SexSky = SexSky
        self.fstring = fstring
        self.threshold = threshold
        self.thresh_area = thresh_area
        self.center_deviated = center_deviated
        self.center_constrain = center_constrain
        self.avoidme = avoidme
        self.LMag = LMag
        self.UMag = UMag
        self.LN = LN
        self.UN = UN
        self.LRe = LRadius
        self.URe = URadius
        self.LRd = LRadius
        self.URd = URadius
        self.bdbox = bdbox
        self.bbox = bbox
        self.dbox = dbox
        self.devauc = devauc
        self.galfitv = galfitv
        self.mag_zero = mag_zero
        self.flag = flag

        self.obj_counter = 0

        self.make_constrain = 1

    def _write_sersic_bar_constrain(self, component, cO):

        f_constrain = open(self.constrain_file, 'a')
        print(cO, self.LN, self.UN)
        if component == 'sersic_main':
            f_constrain.write('{}  n  {} to {} \n'.format(cO, self.LN, self.UN))
        elif component == 'bar':
            f_constrain.write('{}  n  0.1 to 2.2 \n'.format(cO))
        elif component == 'sersic_neighbor':
            f_constrain.write('{}  n  0.5 to 15.0 \n'.format(cO))
            center_constrain = center_contrain * 10
        else:
            print('No sersic is provided')

        if self.center_deviated:
            xy_lim = self.center_deviation - self.center_deviation / 4.0
            xlim_l = self.xcntr_img - xy_lim
            xlim_u = self.xcntr_img + xy_lim
            ylim_l = self.ycntr_img - xy_lim
            ylim_u = self.ycntr_img + xy_lim
            f_constrain.write('{}  x  {} to {} \n'.format(cO, xlim_l,
                                                          xlim_u))
            f_constrain.write('{}  y  {} to {} \n'.format(cO, ylim_l,
                                                          ylim_u))
        else:
            xlim_l = self.xcntr_img - self.center_constrain
            xlim_u = self.xcntr_img + self.center_constrain
            ylim_l = self.ycntr_img - self.center_constrain
            ylim_u = self.ycntr_img + self.center_constrain

            f_constrain.write('{}  x  {} to {} \n'.format(cO, 
                                                          xlim_l, xlim_u))
            f_constrain.write('{}  y  {} to {} \n'.format(cO,
                                                          ylim_l, ylim_u))

        f_constrain.write('{}  mag  {} to {} \n'.format(cO, self.UMag, 
                                                        self.LMag))
        f_constrain.write('{}  re  {} to {} \n'.format(cO, self.LRe, self.URe))
        f_constrain.write('{}  q  {} to {} \n'.format(cO, 0.0, 1.0))
        f_constrain.write('{}  pa  {} to {} \n'.format(cO, -360.0, 360.0))
        f_constrain.close()

    def _sersic_main_constrain(self, cO):
        '''
        Use _write_sersic_bar_constrain
        '''
        f_constrain = open(self.constrain_file, 'a')
        f_constrain.write('{}  n  {} to {} \n'.format(cO, self.LN, self.UN))
        f_constrain.write('{}  x  {} to {} \n'.format(cO, -center_constrain,
                                                      center_constrain))
        f_constrain.write('{}  y  {} to {} \n'.format(cO, -center_constrain,
                                                      center_constrain))
        f_constrain.write('{}  mag  {} to {} \n'.format(cO, self.UMag, self.LMag))
        f_constrain.write('{}  re  {} to {} \n'.format(cO, self.LRe, self.URe))
        f_constrain.write('{}  q  {} to {} \n'.format(cO, Lq, Uq))
        f_constrain.write('{}  pa  {} to {} \n'.format(cO, Lpa, Upa))
        f_constrain.close()

    def _bar_constrain(self, cO):
        '''
        Use _write_sersic_bar_constrain
        '''
        f_constrain = open(self.constrain_file, 'a')
        f_constrain.write('{}  n  0.1 to 2.2 \n'.format(cO))
        f_constrain.write('{}  x  {} to {} \n'.format(cO, -center_constrain,
                                                      center_constrain))
        f_constrain.write('{}  y  {} to {} \n'.format(cO, -center_constrain,
                                                      center_constrain))
        f_constrain.write('{}  mag  {} to {} \n'.format(cO, self.UMag, self.LMag))
        f_constrain.write('{}  re  {} to {} \n'.format(cO, self.LRe, self.URe))
        f_constrain.write('{}  q  {} to {} \n'.format(cO, Lq, Uq))
        f_constrain.write('{}  pa  {} to {} \n'.format(cO, Lpa, Upa))
        f_constrain.close()

    def _expdisk_constrain(self, cO):
        f_constrain = open(self.constrain_file, 'a')
        if self.center_deviated:
            xy_lim = self.center_deviation - self.center_deviation / 4.0
            xlim_l = self.xcntr_img - xy_lim
            xlim_u = self.xcntr_img + xy_lim
            ylim_l = self.ycntr_img - xy_lim
            ylim_u = self.ycntr_img + xy_lim
            f_constrain.write('{}  x  {} to {} \n'.format(cO, xlim_l,
                                                          xlim_u))
            f_constrain.write('{}  y  {} to {} \n'.format(cO, ylim_l,
                                                          ylim_u))
        else:
            #The below line is used to make the self.center_constrain shorter
            #to ccen
            xlim_l = self.xcntr_img - self.center_constrain
            xlim_u = self.xcntr_img + self.center_constrain
            ylim_l = self.ycntr_img - self.center_constrain
            ylim_u = self.ycntr_img + self.center_constrain

            f_constrain.write('{}  x  {} to {} \n'.format(cO,
                                                          xlim_l, xlim_u))
            f_constrain.write('{}  y  {} to {} \n'.format(cO,
                                                          ylim_l, ylim_u))

        f_constrain.write('{}  mag  {} to {} \n'.format(cO, self.UMag, 
                                                        self.LMag))
        f_constrain.write('{}  rs  {} to {} \n'.format(cO, self.LRd, self.URd))
        f_constrain.write('{}  q  {} to {} \n'.format(cO, 0.0, 1.0))
        f_constrain.write('{}  pa  {} to {} \n'.format(cO, -360.0, 360.0))
        f_constrain.close()

    def _sersic_neighbor_constrain(self, cO):
        '''
        Use _write_sersic_bar_constrain
        '''
        f_constrain = open(self.constrain_file, 'a')
        f_constrain.write('{} n 0.2 to 20. \n'.format(cO))
        f_constrain.write('{} mag 100. to -100 \n'.format(cO))
        f_constrain.write('{} re 0.1 to 500.0 \n'.format(cO))
        f_constrain.write('{} q 0.01 to 1.0 \n'.format(cO))
        f_constrain.write('{} pa -360. to 360. \n'.format(cO))
        f_constrain.close()

       
    def _write_bulge(self, target, fcon):
        self.obj_counter += 1

        if self.make_constrain == 1:
            self._write_sersic_bar_constrain('sersic_main', 1)            

        # write config file
        fcon.write('# Sersic function\n\n')
        comment = 'Object type\n'
        fcon.writelines([' 0) sersic # {}'.format(comment)])

        comment = ' position x, y [pixel]\n'
        fcon.writelines([' 1) {:.2f} {:.2f} {} {} # {}'.format(target.xcntr, 
                                                          target.ycntr,
                                                          self.fitting[0],
                                                          self.fitting[0],
                                                          comment)])
        comment = 'total magnitude\n'
        fcon.writelines([' 3) {:.2f} 1 # {}'.format(target.mag, comment)])

        comment = 'R_e [Pixels]\n'
        fcon.writelines([' 4) {:.2f} 1 # {}'.format(target.radius, comment)])

        comment = 'Sersic exponent (deVauc=4, expdisk=1)\n'
        fcon.writelines([' 5) 4.0 {} # {}'.format(int(not self.devauc),
                                                     comment)])

        if self.galfitv >= 3.0:
            comment1 = 'axis ratio (b/a)\n'
            fcon.writelines([' 9) {:.2f} 1 # {}'.format(target.axis_rat, 
                                                        comment1)])
            
            comment2 = 'position angle (PA) [Degrees: Up=0, Left=90]\n'
            fcon.writelines([' 10) {:.2f} 1 # {}'.format(target.pos_ang_galfit,
                                                         comment2)])
        else:
            fcon.writelines([' 8) {:.2f} 1 # {}'.format(target.axis_rat, 
                                                        comment1)])

            comment2 = 'position angle (PA) [Degrees: Up=0, Left=90]\n'
            fcon.writelines([' 9) {:.2f} 1 # {}'.format(target.pos_ang_galfit,                                                          comment2)])

            if self.bdbox or self.bbox:
                comment = 'diskiness (< 0) or boxiness (> 0)\n'
                fcon.writelines(['10) 0.0 1	# {}'.format(comment)])
            else:
                fcon.writelines(['10) 0.0 0	# {}'.format(comment)])

        comment = 'output image (see above)\n\n\n'
        fcon.writelines([' Z) 0 # {}'.format(comment)])
    

    def _write_disk(self, target, fcon):
        self.obj_counter += 1

        if self.make_constrain == 1:
            self._expdisk_constrain(2)

        # write config file
        fcon.writelines(['# Exponential function\n\n'])

        comment = 'Object type\n'
        fcon.writelines([' 0) expdisk # {}'.format(comment)])

        comment = ' position x, y [pixel]\n'
        fcon.writelines([' 1) {:.2f} {:.2f} {} {} # {}'.format(target.xcntr,
                                                          target.ycntr,
                                                          self.fitting[0],
                                                          self.fitting[0],
                                                          comment)])
        comment = 'total magnitude\n'
        fcon.writelines([' 3) {:.2f} 1 # {}'.format(target.mag, comment)])

        comment = 'R_e [Pixels]\n'
        fcon.writelines([' 4) {:.2f} 1 # {}'.format(target.radius, comment)])

        comment1 = 'axis ratio (b/a)\n'
        comment2 = 'position angle (PA) [Degrees: Up=0, Left=90]\n'
        if self.galfitv >= 3.0:
            fcon.writelines([' 9) {:.2f} 1 # {}'.format(target.axis_rat,
                                                        comment1)])

            fcon.writelines([' 10) {:.2f} 1 # {}'.format(target.pos_ang_galfit,
                                                         comment2)])        
        else:
            fcon.writelines([' 8) {:.2f} 1 # {}'.format(target.axis_rat,
                                                        comment1)])

            comment2 = 'position angle (PA) [Degrees: Up=0, Left=90]\n'
            fcon.writelines([' 9) {:.2f} 1 # {}'.format(target.pos_ang_galfit,                                                          comment2)])

            if self.bdbox or self.bbox:
                comment = 'diskiness (< 0) or boxiness (> 0)\n'
                fcon.writelines(['10) 0.0 1 # {}'.format(comment)])
            else:
                fcon.writelines(['10) 0.0 0 # {}'.format(comment)])

        comment = 'output image (see above)\n\n\n'
        fcon.writelines([' Z) 0 # {}'.format(comment)])


    
                         
    def _write_point(self, target, fcon):
        self.obj_counter += 1

        if self.make_constrain == 1:
            f_constrain = open(self.constrain_file, 'a')
            # write constraints
            f_constrain.write('{} x -2.0 -2.0'.format(self.obj_counter)) 
            f_constrain.write('{} y -2.0 -2.0'.format(self.obj_counter))
            f_constrain.close()
 
        # write config
        pmag = target.mag + 2.5 * np.log10(6.0)
        fcon.writelines(['#point source\n\n'])
        comment = 'Object type\n'
        fcon.writelines([' 0) psf # {}'.format(comment)])

        comment = 'position x, y [pixel]\n'
        #print('Point', target.xcntr, target.ycntr, self.fitting[4], self.fitting[4])
        fcon.writelines([' 1) {:.2f} {:.2f} {} {} # {}'.format(target.xcntr, 
                                                           target.ycntr,
                                                           self.fitting[4],
                                                           self.fitting[4],
                                                           comment)])

        comment = 'total magnitude\n'
        fcon.writelines([' 3) {:.2f} # {}'.format(pmag, comment)])

        comment = 'output image (see above)\n\n\n'
        fcon.writelines([' Z) 0 {}'.format(comment)])

    
    def _write_gaussian(self, target, fcon):
        self.obj_counter += 1

        if self.make_constrain == 1:
            # write constraints
            f_constrain = open(self.constrain_file, 'a')
            f_constrain.write('{} x -2.0 2.0\n'.format(self.obj_counter))
            f_constrain.write('{} y -2.0 2.0\n'.format(self.obj_counter))
            f_constrain.close()            

        # write config
        gmag = target.mag + 2.5 * np.log10(2.0)
        fcon.writelines(['# Gaussian function\n\n'])
        fcon.writelines([' 0) gaussian              # Object type\n'])

        comment = ' position x, y [pixel]\n'
        fcon.writelines([' 1) {:.2f} {:.2f} 1 1 # {}'.format(target.xcntr,
                                                          target.ycntr,
                                                          comment)])
        comment = 'total magnitude\n'
        fcon.writelines([' 3) {:.2f} 1 # {}'.format(gmag, comment)])

        comment = 'FWHM\n'
        fcon.writelines([' 4) 0.50 1 # {}'.format(comment)])


        if self.galfitv >= 3.0:
            comment1 = 'axis ratio (b/a)\n'
            fcon.writelines([' 9) {:.2f} 1 # {}'.format(target.axis_rat,
                                                            comment1)])

            comment2 = 'position angle (PA) [Degrees: Up=0, Left=90]\n'
            fcon.writelines([' 10) {:.2f} 1 # {}'.format(target.pos_ang_galfit,
                                                             comment2)])         
        else:
            fcon.writelines([' 8) {:.2f} 1 # {}'.format(target.axis_rat,
                                                        comment1)])

            fcon.writelines([' 9) {:.2f} 1 # {}'.format(target.pos_ang_galfit,                                                          comment2)])

            comment = 'diskiness (< 0) or boxiness (> 0)\n'
            fcon.writelines(['10) 0.0 0 # {}'.format(comment)])

        comment = 'output image (see above)\n\n\n'
        fcon.writelines([' Z) 0 {}'.format(comment)])


    def _write_bar(self, target, fcon):
        self.obj_counter += 1
        
        if self.make_constrain == 1:
            self._write_sersic_bar_constrain('bar', self.obj_counter)
        
        barmag = target.mag + 2.5 * np.log10(3.0)
        fcon.write('# Sersic function for bar\n\n')

        comment = 'Object type\n'
        fcon.writelines([' 0) sersic # {}'.format(comment)])

        comment = ' position x, y [pixel]\n'
        fcon.writelines([' 1) {:.2f} {:.2f} {} {} # {}'.format(target.xcntr,
                                                          target.ycntr,
                                                          self.fitting[3],
                                                          self.fitting[3],
                                                          comment)])
        comment = 'total magnitude\n'
        fcon.writelines([' 3) {:.2f} 1 # {}'.format(barmag, comment)])

        comment = 'R_e [Pixels]\n'
        fcon.writelines([' 4) {:.2f} 1 # {}'.format(target.radius, comment)])

        comment = 'Sersic exponent using a lower value of n=0.05\n'
        fcon.writelines([' 5) 0.5 1 # {}'.format(comment)])

        if self.galfitv >= 3.0:
            comment1 = 'smaller axis ratio (b/a) = 0.3\n'
            fcon.writelines([' 9) 0.3 1 # {}'.format(comment1)])

            comment2 = 'position angle (PA) [Degrees: Up=0, Left=90]\n'
            fcon.writelines([' 10) {:.2f} 1 # {}'.format(target.pos_ang_galfit,
                                                         comment2)])
        else:
            fcon.writelines([' 8) 0.3 1 # {}'.format(comment1)])

            comment2 = 'position angle (PA) [Degrees: Up=0, Left=90]\n'
            fcon.writelines([' 9) {:.2f} 1 # {}'.format(target.pos_ang_galfit,                                                          comment2)])

            if self.bdbox or self.bbox:
                comment = 'diskiness (< 0) or boxiness (> 0)\n'
                fcon.writelines(['10) 0.0 1 # {}'.format(comment)])
            else:
                fcon.writelines(['10) 0.0 0 # {}'.format(comment)])

        comment = 'output image (see above)\n\n\n'
        fcon.writelines([' Z) 0 # {}'.format(comment)])


    
    def _write_sky(self, target, fcon):
        self.obj_counter += 1
            
        fcon.writelines(['# sky\n\n']) 
        fcon.writelines([' 0) sky\n'])
        comment = 'sky background [ADU counts\n'
        fcon.writelines([' 1) {:.5f} {} # {}'.format(self.SexSky, 
                                                     self.fitting[2], 
                                                     comment)])
        comment = 'dsky/dy (sky gradient in x)\n'
        fcon.writelines([' 2) 0.0 0 # {}'.format(comment)])

        comment = 'dsky/dy (sky gradient in y)\n'
        fcon.writelines([' 3) 0.0 0 # {}'.format(comment)])   

        comment = 'output image (see above)\n\n\n'
        fcon.writelines([' Z) 0 # {}'.format(comment)])



    def _write_neighbor(self, neighbor, fcon):
        self.obj_counter += 1
        
        if self.make_constrain == 1:
            # write constraints
            self._sersic_neighbor_constrain(self.obj_counter)

        # write config file
        fcon.write('# Neighbor Sersic function\n\n')
        comment = 'Object type\n'
        fcon.writelines([' 0) sersic # {}'.format(comment)])

        comment = ' position x, y [pixel]\n'
        fcon.writelines([' 1) {:.2f} {:.2f} 1 1 # {}'.format(neighbor.xcntr, 
                                                     neighbor.ycntr,
                                                     comment)])
        comment = 'total magnitude\n'
        fcon.writelines([' 3) {:.2f} 1 # {}'.format(neighbor.mag, comment)])

        comment = 'R_e [Pixels]\n'
        fcon.writelines([' 4) {:.2f} 1 # {}'.format(neighbor.radius, comment)])

        comment = 'Sersic exponent (deVauc=4, expdisk=1)\n'
        fcon.writelines([' 5) 4.0 1 # {}'.format(comment)])

        if self.galfitv >= 3.0:
            comment1 = 'axis ratio (b/a)\n'
            fcon.writelines([' 9) {:.2f} 1 # {}'.format(neighbor.axis_rat, 
                                                    comment1)])
            
            comment2 = 'position angle (PA) [Degrees: Up=0, Left=90]\n'
            fcon.writelines([' 10) {:.2f} 1 # {}'.format(neighbor.pos_ang_galfit,
                                                     comment2)])
        else:
            fcon.writelines([' 8) {:.2f} 1 # {}'.format(neighbor.axis_rat, 
                                                    comment1)])

            comment2 = 'position angle (PA) [Degrees: Up=0, Left=90]\n'
            fcon.writelines([' 9) {:.2f} 1 # {}'.format(neighbor.pos_ang_galfit,                                                          comment2)])

            if self.bdbox or self.bbox:
                comment = 'diskiness (< 0) or boxiness (> 0)\n'
                fcon.writelines(['10) 0.0 1	# {}'.format(comment)])
            else:
                fcon.writelines(['10) 0.0 0	# {}'.format(comment)])

        comment = 'output image (see above)\n\n\n'
        fcon.writelines([' Z) 0 # {}'.format(comment)])

    def write_config(self):

        #print(self.good_object.shape)
        target = mf.GetSExObj(NXPTS=2 * self.half_size, 
                              NYPTS=2 * self.half_size, 
                              values=self.good_object)
        #print(self.good_object.shape)
        target_imcenter = [target.xcntr, target.ycntr]
        #print('target_imcenter', target_imcenter)
        target.set_center(self.xcntr_img, self.ycntr_img)
        #print(target.xcntr, target.ycntr)
        target.pos_ang
        
        self.config_file = 'G_{}.in'.format(self.fstring) #GALFIT configuration file
        self.constrain_file = '{}.con'.format(self.fstring)
        self.mask_file = 'M_{}.fits'.format(self.fstring) 
        self.oimg   = 'O_{}.fits'.format(self.fstring)

        if os.path.exists(self.constrain_file):
            make_constrain = 0
        else:
            make_constrain = 1

        if os.path.exists(self.config_file):
            pass
        else:
            fcon = open(self.config_file, 'w')

            #Write configuration file
            fcon.write('# IMAGE PARAMETERS\n')

            self.cutimage = os.path.split(self.cutimage)[-1]
            comment = 'Input data image(FITS file)\n'
            fcon.writelines(['A) {} # {}'.format(self.cutimage, comment)])

            comment = 'Name of the output image\n'
            fcon.writelines(['B) {} # {}'.format(self.oimg, comment)])

            self.whtimage = os.path.split(self.whtimage)[-1]
            comment = 'Noise image name (made from data if blank or "none")\n'
            fcon.writelines(['C) {} # {}'.format(self.whtimage, comment)])

            cpsf = os.path.join(self.datadir, self.psffile)
            comment = 'Input PSF image for convolution (FITS file)\n'
            fcon.writelines(['D) {} # {}'.format(self.psffile, comment)])

            comment = 'PSF oversampling factor relative to data\n'
            fcon.writelines(['E) 1 # {}'.format(comment)])

            comment = 'Bad pixel mask(FITS image or ASCII coord list)\n'
            fcon.writelines(['F) {} # {}'.format(self.mask_file, comment)])

            comment = 'File with parameter constraints (ASCII file)\n'
            fcon.writelines(['G) {} # {}'.format(self.constrain_file, comment)])


            comment = 'Image region to fit (xmin xmax ymin ymax)\n'
            fcon.writelines(['H) 1 {} 1 {} # {}'.format(self.half_size * 2, 
                                                        self.half_size * 2,
                                                               comment)])
            #comment = 'Size of convolution box (x y)\n'
            #fcon.writelines(['I) {} {} # {}'.format(self.NXPTS, 
                                                               #self.NYPTS,
                                                               #comment)])
            # This below line really shouldn't be hardcoded!!!!
            comment = 'Size of convolution box (x y)\n'
            fcon.writelines(['I) 100 100 # {}'.format(comment)])

            comment = 'Magnitude photometric zeropoint\n'
            fcon.writelines(['J) {} # {}'.format(self.mag_zero, comment)])

            comment = 'Display type (regular, curses, both)\n'
            fcon.writelines(['O) regular # {}'.format(comment)])

            comment = 'Create output image only? (1=yes; 0=optimize)\n'
            fcon.writelines(['P) 0 # {}'.format(comment)])

            comment = 'Modify/create objects interactively?\n\n\n'
            fcon.writelines(['S) 0 # {}'.format(comment)]) 


            #print('self.components', self.components)
            for comp in self.components:
                if comp == 'bulge':
                    self._write_bulge(target, fcon)
                    self.flag = SetFlag(self.flag, GetFlag('FIT_BULGE'))
                elif comp == 'disk':
                    self._write_disk(target, fcon)
                    self.flag = SetFlag(self.flag, GetFlag('FIT_DISK'))   
                elif comp == 'point':
                    self._write_point(target, fcon)
                    self.flag = SetFlag(self.flag, GetFlag('FIT_POINT'))
                elif comp == 'bar':
                    self._write_bar(target, fcon)
                    self.flag = SetFlag(self.flag, GetFlag('FIT_BAR'))

            self._write_sky(target, fcon)
            
            if self.center_deviated:
                self.center_deviated = 0

            isneighbour = 0
            #print('CCCC')
            #Neighbour in the cutimage needs to be fitted
            fit_neighbor_cutimage = [[self.xcntr_img, self.ycntr_img]]
            target.set_center(target_imcenter[0], target_imcenter[1])
            for line_neigh in open(self.sex_cata_neighbor, 'r'):
                #print(line_neigh) 
                values_neigh = np.fromstring(line_neigh, sep=' ')
                neighbor = mf.GetSExObj(NXPTS=self.half_size * 2, 
                                        NYPTS=self.half_size * 2, 
                        values=values_neigh)
                #print('target area', target.area)
                #print('target xcntr', target.xcntr, target.ycntr)
                #print('neigh xcntr', neighbor.xcntr, neighbor.ycntr)
                #print('self.xcntr_img', self.xcntr_img, self.ycntr_img)
                xn = self.xcntr_img - target_imcenter[0] + neighbor.xcntr
                yn = self.ycntr_img - target_imcenter[1] + neighbor.ycntr
                #print(xn, yn)
                if target.get_mask(neighbor, self.threshold, self.thresh_area,
                                   self.avoidme) == 0:
                    isneighbour = 1
                    # recenter in chip coordinates
                    xn = self.xcntr_img - target_imcenter[0] + neighbor.xcntr
                    yn = self.ycntr_img - target_imcenter[1] + neighbor.ycntr
                    neighbor.set_center(xn, yn)
                    
                    self._write_neighbor(neighbor, fcon)
                    fit_neighbor_cutimage.append([xn, yn])
            
            #print('CCCC')
            if isneighbour:
                self.flag  = SetFlag(self.flag, GetFlag('NEIGHBOUR_FIT'))
            
            fcon.close()

            self.fit_neighbor_cutimage = np.array(fit_neighbor_cutimage).astype(int)
            #print('neighbor_cutimage', self.fit_neighbor_cutimage)
            
   
