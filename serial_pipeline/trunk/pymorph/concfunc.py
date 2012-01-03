import numpy as np
import pyfits
import os
import numpy.ma as ma
class concentration:
    """The class for finding concentration parameter. The algorithm used
       here is as follows
       1. It calculates the average light at different radii. i.e., 
          the average light in an annular ring at different radii. 

       2. It calculates the average light inside the apertures of 
          different radii.  

       3. From the above two it calculates the petrosian eta(r) value. 
          Petrosian ratio at a radius r from the center of an object 
          to be the ratio of the local surface brightness in an annulus 
          at r to the mean surface brightness within r
          eta(r) =  I(r) / <I(r)>
       4. Find the radius at which the Petrosian equal to  0.2 
 
       5. Compute the light inside the aperture of radius 1.5 times the   
          Petrosian radius, that contains more than 90% of the galaxy's 
          total light. 

       6. Find the 20%, 50% and 80% light radii. Linear interpolation 
          is used for this.

       7. Compute concentration parameter as 5*log(r(80%)/r(20%))         
    """        
    def __init__(self, z, mask, xcntr, ycntr, nxpts, nypts, pa, eg, background):
        self.z              = z
        self.mask           = mask
        self.xcntr          = xcntr
        self.ycntr          = ycntr
        self.nxpts          = nxpts
        self.nypts          = nypts
        self.eta_radius     = 2.0
        self.pa             = pa
        self.eg             = eg
        self.oversamp       = 10
        self.incr           = 1.0
        self.incrl          = 1.0
        self.background     = background
        self.oversamp       = self.oversamp * 1.0
        
        xmin = np.floor(xcntr - 4)
        ymin = np.floor(ycntr - 4)
        xmax = np.floor(xcntr + 4)
        ymax = np.floor(ycntr + 4)
        CntrZ = self.z[ymin:ymax, xmin:xmax].copy()
        # zooming the image
        z = z / (self.oversamp * self.oversamp)
        ZoomZ = np.repeat(z, self.oversamp, axis=1)
        ZoomZ = np.repeat(ZoomZ, self.oversamp, axis=0)
        # zooming the mask
        ZoomM = np.repeat(mask, self.oversamp, axis=1)
        ZoomM = np.repeat(ZoomM, self.oversamp, axis=0)
        ZoomZM = ma.masked_array(ZoomZ, ZoomM)
        pa = self.pa / 180 * angle
        co = np.cos(pa)
        si = np.sin(pa)
        SizeY = ZoomZ.shape[0] # Y size
        SizeX = ZoomZ.shape[1] # X size
        axis_rat = (1 - eg) #eg is defind in this way in pymorph.py
        one_minus_eg_sq = axis_rat * axis_rat 
        x = np.reshape(np.arange(SizeX * SizeY), (SizeY, SizeX)) % SizeX
        x = x.astype(np.float32)
        x /= self.oversamp 
        y = np.reshape(np.arange(SizeX * SizeY), (SizeY, SizeX)) / SizeX
        y = y.astype(np.float32)
        y /= self.oversamp
        rx = (x - xcntr)* co + (y - ycntr) * si
        ry = (xcntr - x) * si + (y - ycntr) * co
        r = np.sqrt(rx**2.0 + ry**2.0 / one_minus_eg_sq)
        AvgIAtR = []
        AvgIInR = []
        IInRArr = []
        RArr = []
        Rbin1 = np.linspace(0.2, 10, num= 9.8 * self.oversamp) # finner bins 
                                                               # at small R
        Rbin2 = np.arange(10, np.max(ZoomZ.shape))
        Rbin = np.concatenate((Rbin1, Rbin2))      
        for ri in Rbin:
            con = (r > ri - 1/self.oversamp) & (r < ri + 1/self.oversamp)
            IAtR = ZoomZM[con].sum()
            NAtR = ma.count(ZoomZM[con]) * 1.0
            AvgIAtR.append(IAtR / NAtR)
            con = (r < ri)
            IInR = ZoomZM[con].sum()
            NInR = ma.count(ZoomZM[con]) * 1.0
            AvgIInR.append(IInR / NInR)
            IInRArr.append(IInR)
            RArr.append(ri) 
        AvgIAtR = np.asarray(AvgIAtR)
        AvgIInR = np.asarray(AvgIInR)
        IInR = np.asarray(IInR)
        RArr = np.asarray(RArr) 
        Eta = AvgIAtR / AvgIInR
        EtaRAd = RArr[np.argmin(np.abs(Eta - 0.2))]
        TotRad = EtaRad * 1.5
        con = (r < TotRad)
        TotI = ZoomZM[con].sum() 
        FracI = IInR / TotI
        # IInRe = np.sqrt(IInR)
        # TotI = np.sqrt(TotI)
        FracIe = FracI * np.sqrt((1 / IInR) + (1 / TotI))
        self.r20 = RArr[np.argmin(np.abs(FracI - 0.2))]
        self.r50 = RArr[np.argmin(np.abs(FracI - 0.5))]
        self.r80 = RArr[np.argmin(np.abs(FracI - 0.8))]
        self.r90 = RArr[np.argmin(np.abs(FracI - 0.9))]
        self.r20e = (self.r20 / 0.2) * FracIe[np.argmin(np.abs(FracI - 0.2))]
        print self.r20, self.r20e, self.r50

f = pyfits.open(   
concentration(z, xcntr, ycntr, nxpts, nypts, pa, eg, background)
def none():
        Flag_check = 1
#    calculating the eta parameter
#        print self.divide_r,self.r,self.nxpts,self.nypts,self.background,self.incr,self.eta_radius
        try:
             self.total_rad = eta_radius_fnc(self.divide_r,self.r,self.mask,self.nxpts,self.nypts,self.background,self.incr,self.eta_radius)
        except:
            self.total_rad = 9999
            print 'Cannot Find PetroRadius. May be problem with mask'
#        print self.total_rad
        if(self.total_rad == 9999):
            self.concen=9999
            self.error_con=9999
        else:
            self.total_I = self.divide_r[np.where(self.r<=self.total_rad)].sum()

        #calculating the 20% (r20) light radius and 80% (r80) light radius
            tempI90 = 0.0
            tempI80 = 0.0
            tempI50 = 0.0
            tempI20 = 0.0
            tempr20 = 0.0
            tempr50 = 0.0
            tempr80 = 0.0
            tempr90 = 0.0
            flag90 = 0
            flag80 = 0
            flag50 = 0
            flag20 = 0 #this flags equal 1 when r80 and r20 find
            orad=self.total_rad*1.0 #orad is the initial radius to find the r80 and hopes the r80 is inside the eta_radius
            Inner = 1
            I20 = 0.2 * self.total_I
            I50 = 0.5 * self.total_I
            I80 = 0.8 * self.total_I
            I90 = 0.9 * self.total_I
#            ttempI20 = self.divide_r[np.where(self.r<=self.eta_radius-1.0)].sum()
            ttempI20 = DivideCntr[np.where(DivideCntrR<=self.eta_radius-1.0)].sum()
            while flag90==0 or flag80==0 or flag50==0 or flag20==0:
                if(flag90==0):
                    tempI90 = self.divide_r[np.where(self.r<=orad)].sum()
                if(flag80==0):
                    tempI80 = self.divide_r[np.where(self.r<=orad)].sum()
                if(flag50==0):
                    tempI50 = self.divide_r[np.where(self.r<=orad)].sum()
                if(flag20==0):
                    if irad > DivideCntr.shape[0] - 1.0 and Inner:
                        Inner = 0
                        irad = self.eta_radius
                        ttempI20 = self.divide_r[np.where(self.r<=self.eta_radius-1.0)].sum()
                    if Inner:
                        tempI20 = DivideCntr[np.where(DivideCntrR<=irad)].sum()
                    else:
                        tempI20 = self.divide_r[np.where(self.r<=irad)].sum()
                if(flag90==0 and tempI90<I90):
                    flag90=1
                    alpha90=error(I90,tempI90,ttempI90,orad,self.total_I,self.total_rad,self.background)
                    self.r90=orad+alpha90[0]
#                    print 'r90', self.r90, orad
                if(flag80==0 and tempI80<I80):
                    flag80=1
                    alpha80=error(I80,tempI80,ttempI80,orad,self.total_I,self.total_rad,self.background)
                    self.r80=orad+alpha80[0]
#                    print 'r80', self.r80, orad
                if(flag50==0 and tempI50<I50):
                    flag50=1
                    alpha50=error(I50,tempI50,ttempI50,orad,self.total_I,self.total_rad,self.background)
                    self.r50=orad + alpha50[0]
#                    print 'r50', self.r50, orad
                if(flag20==0 and tempI20>I20):
                    flag20=1
                    alpha20=error(I20,ttempI20,tempI20,irad-self.incr,self.total_I,self.total_rad,self.background)
                    if Inner:
                        self.r20=(irad-self.incrl+alpha20[0])/ 10.0
                    else:
                        self.r20=irad-self.incrl+alpha20[0]
                orad-=self.incr
                irad+=self.incrl
                ttempI90=tempI90
                ttempI80=tempI80
                ttempI50=tempI50
                ttempI20=tempI20

            self.r20_error=alpha20[1]
            self.r50_error=alpha50[1]
            self.r80_error=alpha80[1]
            self.r90_error=alpha90[1]
            #calculate concentration parameter "concen"
            self.concen=5*np.log10(self.r80/self.r20)
            self.error_ratio=np.sqrt((self.r80/self.r20)**2*((alpha80[1]/self.r80)**2+(alpha20[1]/self.r20)**2))
            self.error_con=5*self.error_ratio/(self.r80/self.r20)
        #econcen = 5*log10 (r80/r20) - 5*log10 ((r80-0.5)/(r20+0.5)) This is the way Conselice finds error in concentration
#            print self.concen, self.error_con

#------This function returns r------#

def return_r(nxpts,nypts,xcntr,ycntr,pa,eg,pixel_division):
    divide_x=(np.reshape(np.arange(nxpts*pixel_division*pixel_division*nypts),(nxpts*pixel_division,nypts*pixel_division))/(nypts*pixel_division))*(1/(pixel_division*1.0))

    divide_y=(np.reshape(np.arange(nxpts*pixel_division*pixel_division*nypts),(nxpts*pixel_division,nypts*pixel_division))%(nypts*pixel_division))*(1/(pixel_division*1.0))
    axis_rat = (1 - eg) #eg is defind in this way in pymorph.py
    one_minus_eg_sq = axis_rat * axis_rat # eg is axis ratio. check fitgal code(1.0-eg)**2.0
    co = np.cos(pa * np.pi / 180.0)
    si = np.sin(pa * np.pi / 180.0)
    
    #radius corresponds to the pixel index
    tx = (divide_x-xcntr) * co + (divide_y-ycntr) * si
    ty = (xcntr-divide_x) * si + (divide_y-ycntr) * co
    tx = tx.astype(np.float32)
    tx = tx.astype(np.float32)
    r = np.sqrt(tx**2.0 + ty**2.0 / one_minus_eg_sq)
    return r



#-----Function which divide image-----#



#------The function which returns the total radius ie. 1.5*r(eta=0.2)-----#
#I_ave_radius is the intensity at different anulus, I_ave - intensity inside different radius, no_I_ave_radius - the no of pixels in the anular region, no_I_ave - the no. of pixels inside the different radius

def eta_radius_fnc(divide_r, r, mask, nxpts, nypts, background, incr, eta_radius):
    divide_r = ma.masked_array(divide_r, mask)
    FLAG_ETA, nn = 0, 0 # The n is using because the program won't calculate eta radius if the eta<0.2 in the first run itself. FLAG_ETA will tell whether the condition eta<0.2 is reached.
    eta, I_ave_radius, I_ave, no_I_ave, no_I_ave_radius, r50, temp_eta = \
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
#Here it try to find the eta radius. It will run until FLAG_ETA=0 and eta_radius< min(NXPTS,NYPTS)
    while FLAG_ETA == 0:
        temp_eta = eta
        I_ave_radius = divide_r[np.where(r <= eta_radius + 1.0)].sum()\
                               - divide_r[np.where(r < eta_radius - 1.0)].sum()
        no_I_ave_radius = ma.count(divide_r[np.where(r <= eta_radius + \
                                  1.0)]) - ma.count(divide_r[np.where(r <= \
                                  eta_radius - 1.0)])
        I_ave = divide_r[np.where(r <= eta_radius)].sum()
        no_I_ave = ma.count(divide_r[np.where(r<=eta_radius)])    
        z = ma.filled(divide_r, 0.0)
#        os.system('rm -f WriteTestZ.fits WriteTestR.fits WriteTestM.fits')
#        HdU = pyfits.PrimaryHDU(z.astype(np.float32))
#        HdU.writeto('WriteTestZ.fits')
#        HdU = pyfits.PrimaryHDU(r.astype(np.float32))
#        HdU.writeto('WriteTestR.fits')
#        HdU = pyfits.PrimaryHDU(mask.astype(np.float32))
#        HdU.writeto('WriteTestM.fits')
        if(I_ave==0 or no_I_ave==0 or I_ave_radius == 0 or no_I_ave_radius==0):
            FLAG_ETA=0
        else :
            I_ave=I_ave/no_I_ave #average intensity inside a radius
            I_ave_radius=I_ave_radius/no_I_ave_radius 
                        #average intensity inside the annulus
            #print I_ave,I_ave_radius,eta_radius,eta
#            I_ave = I_ave / (3.14 * eta_radius * eta_radius)
#            I_ave_radius = I_ave_radius / (3.14 * ((eta_radius \
#                                       + 1.0) * (eta_radius + 1.0) - \
#                                       (eta_radius - 1.0) * (eta_radius - 1.0)))
            eta = (I_ave_radius) / (I_ave) #eta parameter
            if(eta < 0.2 and nn == 0):
                FLAG_ETA = 1
                print "Finding eta radius failed. \
                                       Exiting CASGM function"
                eta_radius_corre = 9999
#                os._exit(0)
            elif(eta < 0.2 and nn > 0):
                #check = check_eta(eta_radius,eta,divide_r,r)
                #eta = check[0]
                #eta_radius = check[1]
                #Flag_check = check[2]
            #if(Flag_check==0):
                FLAG_ETA = 1
                ERROR1 = np.sqrt(eta**2 * (((I_ave_radius * \
                                         no_I_ave_radius + background * \
                                         no_I_ave_radius) / (I_ave_radius * \
                                         no_I_ave_radius)**2) + ((I_ave * \
                                         no_I_ave + background * no_I_ave) / \
                                         (I_ave * no_I_ave)**2)))        
                ERROR2 = np.sqrt(temp_eta**2 * \
                                         (((temp_I_ave_radius * \
                                         temp_no_I_ave_radius + background * \
                                         temp_no_I_ave_radius) / \
                                         (temp_I_ave_radius * \
                                         temp_no_I_ave_radius)**2) + \
                                         ((temp_I_ave * temp_no_I_ave + \
                                         background * temp_no_I_ave) / \
                                         (temp_I_ave * temp_no_I_ave)**2)))
                alpha = (0.2 - temp_eta) / (eta - temp_eta)
                eta_radius_corre = eta_radius - 1 + alpha
                error_eta_radius = np.sqrt(alpha * alpha * \
                                                   (((ERROR2)**2 / (0.2 - \
                                                   temp_eta)**2) + ((ERROR1 + \
                                                   ERROR2)**2 / (eta - \
                                                   temp_eta)**2)))
                TEMP_RAD = eta_radius_corre
        #print eta,eta_radius,I_ave,I_ave_radius
        temp_I_ave_radius = I_ave_radius
        temp_I_ave = I_ave
        temp_no_I_ave = no_I_ave
        temp_no_I_ave_radius = no_I_ave_radius
        eta_radius += incr
        if(eta_radius > np.min([nxpts, nypts]) / 2.0):
            print 'The radius measurment exceeded the size '\
                               'of the frame. Exiting CASGM module'
            eta_radius_corre = 9999
            FLAG_ETA = 1
#            os._exit(0)
        I_ave_radius = 0.0
        I_ave = 0.0
        no_I_ave = 0.0
        no_I_ave_radius = 0.0
        nn += 1
    #Calculating the total light which the light inside 1.5*r(eta=.2). conselice 2003
    if(eta_radius_corre == 9999):
        total_rad = eta_radius_corre
    else:
        total_rad = 1.5*(eta_radius_corre)    
    return total_rad



#-----Function for error estimation for concentration-----#
def error(CONST, x0, x1, rad, total_I, total_rad, background):
    """Function to find error"""
    alpha = (CONST - x0) / (x1 - x0)
    error_alpha = np.sqrt(alpha * alpha * (((x0 + total_I + (total_rad * \
                      background) + (rad * background)) / (CONST - x0)**2) + \
                      ((x1 + x0 + (rad * background) + ((rad + 1) * \
                      background)) / (x1 - x0)**2)))
    return alpha, error_alpha



#----Function which check whether the eta is .2 locally or it is global
def check_eta(eta_radius, eta, divide_r, r):
    ETA = eta
    ETA_RADIUS = eta_radius
    eta_radius = eta_radius + 5
    Flag_check = 0
    while Flag_check == 0 and eta_radius > ETA_RADIUS:
        temp_eta = eta
        I_ave_radius = divide_r[np.where(r <= eta_radius + 0.5)].sum()\
                                - divide_r[np.where(r <= eta_radius - 0.5)].sum()
        no_I_ave_radius = divide_r[np.where(r <= eta_radius + 0.1)].size\
                                  - divide_r[np.where(r <= eta_radius - 0.1)].\
                                  size
        I_ave = divide_r[np.where(r <= eta_radius)].sum()
        no_I_ave = divide_r[np.where(r <= eta_radius)].size            
        if(I_ave == 0 or no_I_ave == 0 or I_ave_radius == 0 or \
                   no_I_ave_radius == 0):
            eta_radius -= 1
        else:
            I_ave = I_ave / no_I_ave #average intensity inside a radius
            I_ave_radius = I_ave_radius / no_I_ave_radius #average intensity inside the annulus
            eta = (I_ave_radius) / (I_ave)
            eta_radius -= 1
        if(abs(eta - 0.2) < abs(ETA - 0.2)):
            Flag_check = 1
            ETA = eta
            ETA_RADIUS = eta_radius + 1
    return ETA, ETA_RADIUS, Flag_check

#f=pyfits.open('n5585_lR.fits')
#z=f[0].data
#header = f[0].header
#if (header.has_key('sky')):
#    sky = header['sky']
#f.close()
#xcntr=192.03
#ycntr=157.42
#pa=0.0
#eg=0.0
#z=z-sky
#background=1390.377
#nxpts=z.shape[0]
#nypts=z.shape[1]
#concentration(z, xcntr, ycntr, nxpts, nypts, pa, eg, background)
