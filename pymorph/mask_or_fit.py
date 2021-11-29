import numpy as np

class GetSExObj():

    def __init__(self, NXPTS=None, NYPTS=None, values=None):

        
        self.values = values
        #self.values = np.fromstring(values, sep=' ')
        if NXPTS is None and self.values is None:
            print('No SExtractor lines (GetSExObject)')

        if NXPTS is not None:
            self.NXPTS = NXPTS
            self.NYPTS = NYPTS
            self.IM_dim()

        if values is None:
            self.values = [-999]
        
        print(self.values.shape)
        ##values = line_s.split()
        if len(self.values) != 19:
            ##values = [float(a) for a in values]
            print(self.values)
            print("Non-standard SEx Cat line!!\nAll values set to -999!!")
            self.values = 20*[-999.0]
        #It will be nice it we can have @property decorator        
#         self.sex_num(int(values[0]))
        #self.sex_center#values[1], values[2])
        #self.center_world#values[3], values[4])
        self.set_center(self.values[1], self.values[2])
#         self.mag(values[7]) 
#         self.mag_err(values[8]) 
#         self.radius(values[9]) 
#         self.sky(values[10]) 
#        self.pos_ang() 
#         self.bbya(values[12])
        self.axis_rat() 
#         self.area(values[13]) 
#         self.maj_axis(values[14])
#         self.star_prob(values[16])

    #@property
    def IM_dim(self):
        """Stores the image dimensions"""
        return self.NXPTS, self.NXPTS

    @property
    def sex_num(self):
        """Stores the object number"""
        return int(self.values[0])
    
    @property
    def sex_center(self):#, xcntr, ycntr):
        """set the object center in pixels"""  
        self.xcntr = self.values[1]
        self.ycntr = self.values[2]
        return self.xcntr, self.ycntr


    @property
    def center_world(self):#, ra_cntr, dec_cntr):
        """set the object center in world coordinate"""      
        self.ra_cntr = self.values[3]
        self.dec_cntr = self.values[4]
        return self.ra_cntr, self.dec_cntr
        
    def set_center(self, xcntr, ycntr):
        """set the object center in pixels"""   
        self.xcntr   = xcntr
        self.ycntr   = ycntr


    def set_center_world(self, ra_cntr, dec_cntr):
        """set the object center in world coordinate"""  
        self.ra_cntr   = ra_cntr
        self.dec_cntr   = dec_cntr
    
    @property
    def mag(self):#, mag):
        """set the object magnitude"""
        return self.values[7]
    
    @property
    def mag_err(self):#, mag_err):
        """set the object magnitude"""
        return self.values[8]

    @property
    def radius(self):#, rad):
        """set the object radius in pixels"""
        return self.values[9]

    @property
    def sky(self):#, sky):
        """set the object sky"""
        return self.values[10]
    
    @property
    def pos_ang_galfit(self):
        return self.values[11] - 90.0
        
    @property
    def pos_ang(self):
        """set the pos ang"""
        self.si = np.sin(self.values[11] * np.pi / 180.0)
        self.co = np.cos(self.values[11] * np.pi / 180.0)
        return self.values[11]
        
    @property
    def bbya(self):#, bbya):
        '''Set minor/major ratio'''
        return 1 / self.values[12]

    def axis_rat(self):#, b_by_a):
        """set the object center axis ratio (b/a)"""        
        self.axis_rat = self.bbya
        self.eg = 1.0 - self.axis_rat
        self.one_minus_eg_sq = (1.0 - self.eg)**2.0

    @property
    def area(self):#, area):
        """set the object area"""
        return self.values[13]

    @property
    def maj_axis(self):#, maj_axis):
        """set the object major axis"""
        return self.values[14]

    @property
    def star_prob(self):#, star_prob):
        """set the object major axis"""
        return self.values[16]

    def calc_rad(self, xloc, yloc):
        """calculate the radius of the given pixels"""
        tx = (xloc - self.xcntr + 1.0) * self.co + (yloc - self.ycntr + 1.0) * self.si
        ty = (self.xcntr - 1.0 - xloc) * self.si + (yloc - self.ycntr + 1.0) * self.co
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
