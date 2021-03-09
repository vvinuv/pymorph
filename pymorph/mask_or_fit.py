import numpy as np

class GetSExObj():

    def __init__(self, NXPTS=None, NYPTS=None, line_s=None):

        if NXPTS is None and NYPTS is None and line_s is None:
            print('No SExtractor lines (GetSExObject)')

        if NXPTS is not None and NYPTS is not None:
            self.set_IM_dim(NXPTS, NYPTS)

        if line_s is None:
            values = [-999]

        values = line_s.split()
        if len(values) == 19:
            values = [float(a) for a in values]
        else:
            print(values)
            print("Non-standard SEx Cat line!!\nAll values set to -999!!")
            values = 20*[-999.0]
        
        self.set_sex_num(int(values[0]))
        self.set_center(values[1], values[2])
        self.set_center_world(values[3], values[4])
        self.set_mag(values[7]) 
        self.set_mag_err(values[8]) 
        self.set_radius(values[9]) 
        self.set_sky(values[10]) 
        self.set_pos_ang(values[11]) 
        self.set_bbya(values[12])
        self.set_axis_rat(1.0/values[12]) 
        self.set_area(values[13]) 
        self.set_maj_axis(values[14])
        self.set_star_prob(values[16])

    def set_IM_dim(self, NXPTS, NYPTS):
        """Stores the image dimensions"""
        self.NXPTS = NXPTS
        self.NYPTS = NYPTS

    def set_sex_num(self, obj_num):
        """Stores the object number"""
        self.sex_num = obj_num
    
    def set_center(self, xcntr, ycntr):
        """set the object center in pixels"""
        self.xcntr   = xcntr 
        self.ycntr   = ycntr

    def set_center_world(self, ra_cntr, dec_cntr):
        """set the object center in world coordinate"""
        self.ra_cntr   = ra_cntr
        self.dec_cntr   = dec_cntr

    def set_mag(self, mag):
        """set the object magnitude"""
        self.mag = mag
    
    def set_mag_err(self, mag_err):
        """set the object magnitude"""
        self.mag_err = mag_err

    def set_radius(self, rad):
        """set the object radius in pixels"""
        self.radius = rad

    def set_sky(self, sky):
        """set the object sky"""
        self.sky = sky
    
    def set_pos_ang(self, pos_ang):
        """set the pos ang"""
        self.pos_ang = pos_ang
        self.pos_ang_galfit = pos_ang -90.0
        self.si = np.sin(self.pos_ang * np.pi / 180.0)
        self.co = np.cos(self.pos_ang * np.pi / 180.0)

    def set_bbya(self, bbya):
        '''Set minor/major ratio'''
        self.bbya = bbya

    def set_axis_rat(self, b_by_a):
        """set the object center axis ratio (b/a)"""
        self.axis_rat = b_by_a
        self.eg = 1.0 - self.axis_rat
        self.one_minus_eg_sq    = (1.0-self.eg)**2.0

    def set_area(self, area):
        """set the object area"""
        self.area = area

    def set_maj_axis(self, maj_axis):
        """set the object major axis"""
        self.maj_axis = maj_axis

    def set_star_prob(self, star_prob):
        """set the object major axis"""
        self.star_prob = star_prob

    def calc_rad(self, xloc, yloc):
        """calculate the radius of the given pixels"""
        tx =(xloc-self.xcntr+1.0)*self.co + (yloc-self.ycntr+1.0)*self.si
        ty = (self.xcntr-1.0-xloc)*self.si + (yloc-self.ycntr+1.0)*self.co
        R = np.sqrt(tx**2.0 + ty**2.0 / self.one_minus_eg_sq)
        return R

    
    def get_mask(self, neighbor, threshold, thresh_area, avoidme):

        """

        Returns 1 if the neighbor should be masked, 0 if it should be fit,
        -1 if the object is beyond the bounds of the image

        """

        mask_it = -999

        if(abs(neighbor.xcntr - self.xcntr) < self.NXPTS / 2.0 + avoidme and \
           abs(neighbor.ycntr - self.ycntr) < self.NYPTS / 2.0 + avoidme and \
           np.sqrt((neighbor.xcntr - self.xcntr)**2.0 + \
           (neighbor.ycntr - self.ycntr)**2.0) > 5.0):
            
        
            if(abs(neighbor.xcntr - self.xcntr) > threshold * (neighbor.maj_axis + \
                                                               self.maj_axis) or \
               abs(neighbor.ycntr - self.ycntr) > threshold * (neighbor.maj_axis + \
                                                               self.maj_axis) or \
               neighbor.area < thresh_area * self.area):

                mask_it = 1 #Mask it!!!!

            else:
                mask_it = 0 #It is a neighbor that needs to be fit!!!

        else:
            mask_it = -1 # It is an object beyond the bounds that should be ignored

        return mask_it
