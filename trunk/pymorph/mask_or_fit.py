import os, sys, pyfits
import numpy as np
import config as c
import ndimage as im

class SEx_obj():
    def __init__(self, NXPTS, NYPTS, line_s = ''):

        self.set_IM_dim(NXPTS, NYPTS)

        values = line_s.split()
        if len(values) == 19:
            values = [float(a) for a in values]
        else:
            print values
            print "Non-standard SEx Cat line!!\nAll values set to -999!!"
            values = 20*[-999.0]
        
        self.set_sex_num(int(values[0]))
        self.set_center(values[1], values[2])
        self.set_mag(values[7]) 
        self.set_radius(values[9]) 
        self.set_sky(values[10]) 
        self.set_pos_ang(values[11]) 
        self.set_axis_rat(1.0/values[12]) 
        self.set_area(values[13]) 
        self.set_maj_axis(values[14])
        self.set_star_prob(values[16])
        return

    def set_IM_dim(self, NXPTS, NYPTS):
        """Stores the image dimensions"""
        self.NXPTS = NXPTS
        self.NYPTS = NYPTS
        return

    def set_sex_num(self, obj_num):
        """Stores the object number"""
        self.sex_num = obj_num
        return
    
    def set_center(self, xcntr, ycntr):
        """set the object center in pixels"""
        self.xcntr   = xcntr 
        self.ycntr   = ycntr
        return

    def set_mag(self, mag):
        """set the object magnitude"""
        self.mag = mag
        return
    
    def set_radius(self, rad):
        """set the object radius in pixels"""
        self.radius = rad
        return

    def set_sky(self, sky):
        """set the object sky"""
        self.sky = sky
        return
    
    def set_pos_ang(self, pos_ang):
        """set the pos ang"""
        self.pos_ang = pos_ang
        self.pos_ang_galfit = pos_ang -90.0
        self.si = np.sin(self.pos_ang * np.pi / 180.0)
        self.co = np.cos(self.pos_ang * np.pi / 180.0)
        return

    def set_axis_rat(self, b_by_a):
        """set the object center axis ratio (b/a)"""
        self.axis_rat = b_by_a
        self.eg = 1.0 - self.axis_rat
        self.one_minus_eg_sq    = (1.0-self.eg)**2.0
        return

    def set_area(self, area):
        """set the object area"""
        self.area = area
        return

    def set_maj_axis(self, maj_axis):
        """set the object major axis"""
        self.maj_axis = maj_axis
        return

    def set_star_prob(self, star_prob):
        """set the object major axis"""
        self.star_prob = star_prob
        return

    def calc_rad(self, xloc, yloc):
        """calculate the radius of the given pixels"""
        tx =(xloc-self.xcntr+1.0)*self.co + (yloc-self.ycntr+1.0)*self.si
        ty = (self.xcntr-1.0-xloc)*self.si + (yloc-self.ycntr+1.0)*self.co
        R = np.sqrt(tx**2.0 + ty**2.0 / self.one_minus_eg_sq)
        return R

    
    def mask_or_fit(self, neighbor, threshold, thresh_area, avoidme):
        """Returns 1 if the neighbor should be masked, 0 if it should be fit,
        -1 if the object is beyond the bounds of the image"""
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
