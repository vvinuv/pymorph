import os
import configparser
import sys
import subprocess
import numpy as np
#from .flagfunc import GetFlag, isset, SetFlag
import traceback
#from .mask_or_fit import GetSExObj

class GalfitConfigRunFunc:
    
    """

    The class making configuration file for GALFIT. The configuration file 
    consists of bulge and disk component of the object and only Sersic 
    component for the neighbours, if any. The sky is always fixed and has
    the value of SExtractor. The disk/boxy parameter is also fixed to zero.
    The initial value for Sersic index 'n' is 4.The configuration file has 
    the name G_string(galid).in. The output image has the name 
    O_string(galid).fits

    """

    def __init__(self, config_file):

        config = configparser.ConfigParser()
        config.read(config_file)

        self.NMag = config.getfloat('params_limit', 'NMag')
        self.NRadius = config.getfloat('params_limit', 'NRadius')

        self.LN = config.getfloat('params_limit', 'LN')
        self.UN = config.getfloat('params_limit', 'UN')

        self.center_constrain = config.getfloat('params_limit', 'center_constrain')

        self.bbox = config.getfloat('params_limit', 'bbox')
        self.barbox = config.getfloat('params_limit', 'barbox')

        components = config.get('galfit', 'components').split(',')
        self.components = [cm.strip() for cm in components]

        fitting = config.get('galfit', 'fitting')
        self.fitting = [int(tf) for tf in fitting.split(',')]

        self.devauc = config.getboolean('galfit', 'devauc')

        self.GALFIT_PATH = config.get('external', 'GALFIT_PATH')

        self.obj_counter = 0

        self.make_constrain = 1

    def _write_sersic_bar_constrain(self, component, cnum):

        new = '\n'

        if component == 'sersic_main':
            constraints = f'{cnum} n {self.LN} to {self.UN}' 
        elif component == 'bar':
            constraints = f'{cnum} n 0.1 to 2.2' 
        lx = self.xcntr_img - 5
        ux = self.xcntr_img + 5
        ly = self.ycntr_img - 5
        uy = self.ycntr_img + 5

        ure = self.NRadius * self.flux_radius

        constraints += f'{new} {cnum} x {lx:.2f} to {ux:.2f}'
        constraints += f'{new} {cnum} y {ly:.2f} to {uy:.2f}'
        constraints += f'{new} {cnum} mag {self.mag_auto - self.NMag:.2f} to {self.mag_auto + self.NMag:.2f}'
        constraints += f'{new} {cnum} re 0.2 to {ure:.2f}'
        constraints += f'{new} {cnum} q 0.1 to 0.9'
        constraints += f'{new} {cnum} pa {0} to {360}'

        f_constrain = open(self.constrain_file, 'a')
        f_constrain.write(constraints)
        f_constrain.close()

    def _expdisk_constrain(self, cnum):

        new = '\n'

        lx = self.xcntr_img - 5
        ux = self.xcntr_img + 5
        ly = self.ycntr_img - 5
        uy = self.ycntr_img + 5

        ure = self.NRadius * self.flux_radius

        constraints = f'{new} {cnum} x {lx:.2f} to {ux:.2f}'
        constraints += f'{new} {cnum} y {ly:.2f} to {uy:.2f}'
        constraints += f'{new} {cnum} mag {self.mag_auto - self.NMag:.2f} to {self.mag_auto + self.NMag:.2f}'
        constraints += f'{new} {cnum} re 0.2 to {ure:.2f}'
        constraints += f'{new} {cnum} q 0.1 to 0.9'
        constraints += f'{new} {cnum} pa {0} to {360}'
        
        f_constrain = open(self.constrain_file, 'a')
        f_constrain.write(constraints)
        f_constrain.close()

    def _sersic_neighbor_constrain(self, cnum):

        new = '\n'

        constraints = f'{new} {cnum} n 0.2 to 10'
        constraints += f'{new} {cnum} mag {self.mag_auto - self.NMag:.2f} to {self.mag_auto + self.NMag:.2f}'
        constraints += f'{new} {cnum} re 0.1 to 500'
        constraints += f'{new} {cnum} q 0.1 to 0.9'
        constraints += f'{new} {cnum} pa {0} to {360}'

        f_constrain = open(self.constrain_file, 'a')
        f_constrain.write(constraints)
        f_constrain.close()

       
    def _write_bulge(self):

        new = '\n'

        self.obj_counter += 1

        if 1: #self.make_constrain == 1:
            self._write_sersic_bar_constrain('sersic_main', self.obj_counter)            

        # write config file
        bulge_conf = f'# Sersic function {new}{new}'

        bulge_conf += f' 0) sersic # Object type {new}'

        bulge_conf += f' 1) {self.xcntr_img:.2f} {self.ycntr_img:.2f} '
        bulge_conf += f'{self.fitting[0]} {self.fitting[0]}'
        bulge_conf += f' # position x, y [pixel] {new}'


        bulge_conf += f' 3) {self.mag_auto:.2f} 1 # total magnitud {new}'

        bulge_conf += f' 4) {self.flux_radius:.2f} 1 # R_e [Pixels] {new}'
        bulge_conf += f' 5) 4.0 {int(not self.devauc)} '
        bulge_conf += f'#Sersic exponent (deVauc=4, expdisk=1) {new}'

        bulge_conf += f' 9) {self.axis_ratio:.2f} 1 # axis ratio (b/a) {new}'
        bulge_conf += f' 10) {self.pos_ang:.2f} 1 ' 
        bulge_conf += f'# position angle (PA) [Degrees: Up=0, Left=90] {new}'

        if self.bbox == -99:
            pass
        else:
            bulge_conf += f'C0) {self.bbox} 1 '
            bulge_conf += f'# diskiness (< 0) or boxiness (> 0){new}'

        bulge_conf += f' Z) 0 # output image (see above) {new}{new}'
    
        return bulge_conf

    def _write_disk(self):

        new = '\n'

        self.obj_counter += 1

        if 1:#self.make_constrain == 1:
            self._expdisk_constrain(2)

        # write config file
        disk_conf = f'# Exponential function {new}{new}'

        disk_conf += f' 0) expdisk # Object type {new}'

        disk_conf += f' 1) {self.xcntr_img:.2f} {self.ycntr_img:.2f} '
        disk_conf += f'{self.fitting[1]} {self.fitting[1]}'
        disk_conf += f' # position x, y [pixel] {new}'


        disk_conf += f' 3) {self.mag_auto:.2f} 1 # total magnitud {new}'

        disk_conf += f' 4) {self.flux_radius:.2f} 1 # R_e [Pixels] {new}'

        disk_conf += f' 9) {self.axis_ratio:.2f} 1 # axis ratio (b/a) {new}'
        disk_conf += f' 10) {self.pos_ang:.2f} 1 '
        disk_conf += f'# position angle (PA) [Degrees: Up=0, Left=90] {new}'

        disk_conf += f' Z) 0 # output image (see above) {new}{new}'

        return disk_conf
                         
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


    def _write_bar(self):

        new = '\n'

        self.obj_counter += 1
        
        if 1:#self.make_constrain == 1:
            self._write_sersic_bar_constrain('bar', self.obj_counter)

        
        # write config file
        bar_conf = f'# Sersic function {new}{new}'

        bar_conf += f' 0) sersic # Object type {new}'

        bar_conf += f' 1) {self.xcntr_img:.2f} {self.ycntr_img:.2f} '
        bar_conf += f'{self.fitting[3]} {self.fitting[3]}'
        bar_conf += f' # position x, y [pixel] {new}'


        bar_conf += f' 3) {self.barmag:.2f} 1 # total magnitud {new}'

        bar_conf += f' 4) {self.bar_radius:.2f} 1 # R_e [Pixels] {new}'
        bar_conf += f' 5) {self.barn} 1 #Sersic exponent {new}'

        bar_conf += f' 9) {self.axis_ratio:.2f} 1 # axis ratio (b/a) {new}'
        bar_conf += f' 10) {self.pos_ang:.2f} 1 '
        bar_conf += f'# position angle (PA) [Degrees: Up=0, Left=90] {new}'

        if self.barbox == -99:
            pass
        else:
            bar_conf += f'C0) {self.barbox} 1 '
            bar_conf += f'# diskiness (< 0) or boxiness (> 0){new}'

        bar_conf += f' Z) 0 # output image (see above) {new}'

        return bar_conf
    
    def _write_sky(self):

        new = '\n'

        self.obj_counter += 1
            
        # write config file
        sky_conf = f'# Sky {new}{new}'

        sky_conf += f' 0) sky {new}'
        sky_conf += f' 1) {self.SexSky:.5f} {self.fitting[2]} '
        sky_conf += f'# sky background [ADU counts] {new}'


        sky_conf += f' 2) 0.0 0 # dsky/dx (sky gradient in x) {new}'
        sky_conf += f' 3) 0.0 0 # dsky/dy (sky gradient in y) {new}'

        sky_conf += f' Z) 0 # output image (see above) {new}{new}{new}'
        
        return sky_conf

    def _write_neighbor(self, neigh):

        new = '\n'

        self.obj_counter += 1
        self._sersic_neighbor_constrain(self.obj_counter)

        # write config file
        neigh_conf = f'# Neighbor Sersic function {new}{new}'

        neigh_conf += f' 0) sersic # Object type {new}'

        neigh_conf += f' 1) {neigh["X_IMAGE"]:.2f} {neigh["Y_IMAGE"]:.2f} '
        neigh_conf += f' 1 1  # position x, y [pixel] {new}'


        neigh_conf += f' 3) {neigh["MAG_AUTO"]:.2f} 1 # total magnitud {new}'

        neigh_conf += f' 4) {neigh["FLUX_RADIUS"]:.2f} 1 # R_e [Pixels] {new}'
        neigh_conf += f' 5) 4.0 1 #Sersic exponent {new}'

        neigh_conf += f' 9) {1/neigh["ELONGATION"]:.2f} 1 # axis ratio (b/a) {new}'
        neigh_conf += f' 10) {neigh["THETA_IMAGE"] - 90:.2f} 1 '
        neigh_conf += f'# position angle (PA) [Degrees: Up=0, Left=90] {new}'

        neigh_conf += f' Z) 0 # output image (see above) {new}{new}'
        
        return neigh_conf

    def write_config(self, galaxies):
        
        new = '\n'


        self.xcntr_img = galaxies.get('target')["X_IMAGE"]
        self.ycntr_img = galaxies.get("target")["Y_IMAGE"]
        self.mag_auto = galaxies.get("target")["MAG_AUTO"]
        self.flux_radius = galaxies.get("target")["FLUX_RADIUS"]
        self.axis_ratio = 1 / galaxies.get("target")["ELONGATION"]
        self.pos_ang = galaxies.get("target")["GALFIT_ANGLE"]
        self.SexSky = galaxies.get("target")["BACKGROUND"]

        self.barmag = self.mag_auto + 0.5
        self.barn = 1.0
        self.bar_radius = self.flux_radius 

        self.size = galaxies.get("target")["IMG_SIZE"]
        self.psffile = galaxies.get("target")["psf_file"]
        self.mag_zero = galaxies.get("target")["mag_zero"]

        rootname = galaxies.get("target")["rootname"]
        gal_id = galaxies.get("target")["gal_id"]
        self.fstring = f"{rootname}_{gal_id}"

        self.galfit_file = f'G_{self.fstring}.in' #GALFIT configuration file
        self.constrain_file = f'{self.fstring}.con'
        self.mask_file = f'M_{self.fstring}.fits'
        self.oimg   = f'O_{self.fstring}.fits'
        self.cutimage   = f'I_{self.fstring}.fits'
        self.whtimage   = f'W_{self.fstring}.fits'
        

        #if os.path.exists(self.constrain_file):
        #    make_constrain = 0
        #else:
        #    make_constrain = 1

        if 0:#os.path.exists(self.galfit_file):
            self.fit_neighbor_cutimage = np.array([[self.xcntr_img, 
                                              self.ycntr_img]]).astype(int)
        else:
            
            galfit_config = f'# IMAGE PARAMETERS {new}'
            # self.cutimage
            comment = 'Input data image(FITS file)'
            galfit_config += f"A) {self.cutimage} # {comment}{new}"

            # self.oimg
            comment = 'Name of the output image'
            galfit_config += f"B) {self.oimg} # {comment}{new}"

            # self.whtimage
            comment = 'Noise image name (made from data if blank or \"none\")'
            galfit_config += f"C) {self.whtimage} # {comment}{new}"

            # self.psffile
            comment = 'Input PSF image for convolution (FITS file)'
            galfit_config += f"D) {self.psffile} # {comment}{new}"

            # constant
            comment = 'PSF oversampling factor relative to data'
            galfit_config += f"E) 1 # {comment}{new}"

            # mask file
            comment = 'Bad pixel mask(FITS image or ASCII coord list)'
            galfit_config += f"F) {self.mask_file} # {comment}{new}"

            # constraint file
            comment = 'File with parameter constraints (ASCII file)'
            galfit_config += f"G) {self.constrain_file} # {comment}{new}"

            # image region
            comment = 'Image region to fit (xmin xmax ymin ymax)'
            galfit_config += f"H) 1 {self.size} 1 {self.size} # {comment}{new}"

            # convolution box
            comment = 'Size of convolution box (x y)'
            galfit_config += f"I) 100 100 # {comment}{new}"

            # zeropoint
            comment = 'Magnitude photometric zeropoint'
            galfit_config += f"J) {self.mag_zero} # {comment}{new}"

            # display type
            comment = 'Display type (regular, curses, both)'
            galfit_config += f"O) regular # {comment}{new}"

            # output mode
            comment = '# Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps'
            galfit_config += f"P) 0 # {comment}{new}"

            # interactive mode
            comment = 'Modify/create objects interactively?'
            galfit_config += f"S) 0 # {comment}{new}{new}{new}"

            #Creating a constraint file and contrain function will append
            open(self.constrain_file, "w").close()

            #self.flag  = SetFlag(self.flag, GetFlag('NEIGHBOUR_FIT'))
            #print('self.components', self.components)
            yes_bar = False
            for comp in self.components:
                if comp == 'bulge':
                    galfit_config += self._write_bulge()
                    #self.flag = SetFlag(self.flag, GetFlag('FIT_BULGE'))
                elif comp == 'disk':
                    galfit_config += self._write_disk()
                    #self.flag = SetFlag(self.flag, GetFlag('FIT_DISK'))   
                elif comp == 'point':
                    self._write_point(target, fcon)
                    #self.flag = SetFlag(self.flag, GetFlag('FIT_POINT'))
                elif comp == 'bar':
                    yes_bar = True
                    galfit_config += self._write_bar()
                    #self.flag = SetFlag(self.flag, GetFlag('FIT_BAR'))
            galfit_config += self._write_sky()
            
            galaxies.get('target')['YES_BAR'] = yes_bar

            isneighbour = 0

            neighbours = galaxies.get('neighbours')
            for _, neigh in neighbours.iterrows():
                galfit_config += self._write_neighbor(neigh)
                isneighbour += 1

            #if isneighbour > 0:
                #self.flag  = SetFlag(self.flag, GetFlag('NEIGHBOUR_FIT'))
            
            g = open(self.galfit_file, 'w')
            g.write(galfit_config) 
            g.close()

    def GalfitRun(self):    

        try:
            subprocess.run([self.GALFIT_PATH, self.galfit_file], check=True)
        except subprocess.CalledProcessError as e:
            print("GALFIT failed")
            print("Return code:", e.returncode)
            print("Command:", e.cmd)
 
