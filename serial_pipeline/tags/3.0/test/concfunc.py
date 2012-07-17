import numpy as np
import pyfits
import os
import numpy.ma as ma
import pylab
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
    def __init__(self, z, mask, xcntr, ycntr, pa, eg, sky):
        self.z              = z
        self.mask           = mask
        self.xcntr          = xcntr
        self.ycntr          = ycntr
        self.pa             = pa
        self.eg             = eg
        self.oversamp       = 10.
        self.sky     = sky
        xmin = np.floor(self.xcntr - 10)
        ymin = np.floor(self.ycntr - 10)
        xmax = np.floor(self.xcntr + 10)
        ymax = np.floor(self.ycntr + 10) 
        # Center region of image and mask
        CntrZ = self.z[ymin:ymax, xmin:xmax].copy()
        CntrM = self.mask[ymin:ymax, xmin:xmax].copy()
        #hdu = pyfits.PrimaryHDU(CntrZ)
        #os.system('rm -f test.fits')
        #hdu.writeto('test.fits')
        # zooming the image
        CntrZ = CntrZ / (self.oversamp * self.oversamp)
        ZoomZ = np.repeat(CntrZ, self.oversamp, axis=1)
        ZoomZ = np.repeat(ZoomZ, self.oversamp, axis=0)
        # zooming the mask
        ZoomM = np.repeat(CntrM, self.oversamp, axis=1)
        ZoomM = np.repeat(ZoomM, self.oversamp, axis=0)
        ZoomZM = ma.masked_array(ZoomZ, ZoomM)
        ZM = ma.masked_array(self.z, self.mask)
        pa = self.pa * np.pi / 180.0 
        self.co = np.cos(pa)
        self.si = np.sin(pa)
        axis_rat = (1 - eg) #eg is defind in this way in pymorph.py
        self.one_minus_eg_sq = axis_rat * axis_rat
        def ReturnIs(marray, xc, yc, rbins, oversamp, total=0):
            """Returns the average quantities at different radius of a
              masked array. total=1 just return a total count within rbins""" 
            SizeY = marray.shape[0] # Y size
            SizeX = marray.shape[1] # X size
            x = np.reshape(np.arange(SizeX * SizeY), (SizeY, SizeX)) % SizeX
            x = x.astype(np.float32)
            x /= oversamp 
            y = np.reshape(np.arange(SizeX * SizeY), (SizeY, SizeX)) / SizeX
            y = y.astype(np.float32)
            y /= oversamp
            rx = (x - xc)* self.co + (y - yc) * self.si
            ry = (xc - x) * self.si + (y - yc) * self.co
            r = np.sqrt(rx**2.0 + ry**2.0 / self.one_minus_eg_sq)
            if total:
                con = (r < rbins)
                TotI = marray[con].sum()
                TotN = ma.count(marray[con]) / (oversamp * oversamp * 1.0)
                return TotI, TotN
            else:
                AvgIAtR = []
                AvgIInR = []
                IInRArr = []
                RArr = []
                NInRArr = []
                letbreak = 0 # this will be used to break the loop if eta is 
                             # less than 0.2 for 20 ri's
                for ri in rbins:
                    con = (r > ri - 1/oversamp) & (r < ri + 1/oversamp)
                    IAtR = marray[con].sum()
                    NAtR = ma.count(marray[con]) * 1.0 
                    con = (r < ri)
                    IInR = marray[con].sum()
                    NInR = ma.count(marray[con]) * 1.0
                    if NAtR == 0 or NInR == 0 or ri > 20 and NAtR < 30 or \
                       ri > 20 and NInR < 30:
                        pass
                    else:
                        AvgIAtR.append(IAtR / NAtR)
                        AvgIInR.append(IInR / NInR)
                        IInRArr.append(IInR)
                        RArr.append(ri)
                        NInRArr.append(NInR)
                        if IAtR * NInR / (NAtR * IInR) < 0.2:
                            letbreak += 1
                    if letbreak > 20:
                        break
                AvgIAtR = np.asarray(AvgIAtR)
                AvgIInR = np.asarray(AvgIInR)
                IInRArr = np.asarray(IInRArr)
                RArr = np.asarray(RArr) 
                NInRArr = np.asarray(NInRArr) / (oversamp * oversamp * 1.0)
                return AvgIAtR, AvgIInR, IInRArr, RArr, NInRArr
        Rbin1 = np.linspace(0.2, 10, num= 9.8 * self.oversamp) # finner bins 
                                                               # at small R
        AvgIAtR, AvgIInR, IInRArr, RArr, NInRArr = ReturnIs(ZoomZM, 9., 9., \
                                                Rbin1, self.oversamp)
        Rbin2 = np.arange(10, np.max(ZM.shape))
        AvgIAtR1, AvgIInR1, IInRArr1, RArr1, NInRArr1 = ReturnIs(ZM, \
                                           self.xcntr, self.ycntr, Rbin2, 1)
        AvgIAtR = np.concatenate((AvgIAtR, AvgIAtR1))
        AvgIInR = np.concatenate((AvgIInR, AvgIInR1))
        IInRArr = np.concatenate((IInRArr, IInRArr1))
        RArr = np.concatenate((RArr, RArr1)) 
        NInRArr = np.concatenate((NInRArr, NInRArr1))
        Eta = AvgIAtR / AvgIInR
        polz = np.poly1d(np.polyfit(RArr, Eta - 0.2, 3))
        print polz.r
        EtaRad = RArr[np.argmin(np.abs(Eta - 0.2))]
        self.TotRad = EtaRad * 1.5
        if self.TotRad < 10:
            TotI, TotN = ReturnIs(ZoomZM, 9.0, 9.0, self.TotRad, 1, total=1)
        else:
            TotI, TotN = ReturnIs(ZM, self.xcntr, self.ycntr, \
                                  self.TotRad, 1, total=1)
        print RArr[-1], IInRArr[-1], TotI, EtaRad, self.TotRad

        FracI = IInRArr / TotI
        print RArr.shape, FracI.shape
        # IInRe = np.sqrt(IInR)
        # TotI = np.sqrt(TotI)
        FracIe = FracI * np.sqrt(((IInRArr + self.sky * NInRArr) / \
                 IInRArr**2.) + ((TotI + self.sky * TotN) / TotI**2.))
        pylab.plot(RArr, FracI)
        pylab.show()

        self.r20 = RArr[np.argmin(np.abs(FracI - 0.2))]
        self.r50 = RArr[np.argmin(np.abs(FracI - 0.5))]
        self.r80 = RArr[np.argmin(np.abs(FracI - 0.8))]
        self.r90 = RArr[np.argmin(np.abs(FracI - 0.9))]
        # Errors are calculated as follows
        # The radius = a factor * fraction of light, for eg.
        # the radius at which the fraction of light becomes 20% 
        # r20 = alpha * f20, r20 = alpha * 0.2, alpha = r20 / 0.2
        # so the error in r20, alpha times error in f20
        self.r20e = (self.r20 / 0.2) * FracIe[np.argmin(np.abs(FracI - 0.2))]
        self.r50e = (self.r50 / 0.5) * FracIe[np.argmin(np.abs(FracI - 0.5))]
        self.r80e = (self.r80 / 0.8) * FracIe[np.argmin(np.abs(FracI - 0.8))]
        self.r90e = (self.r90 / 0.9) * FracIe[np.argmin(np.abs(FracI - 0.9))]
        print self.r20, self.r20e, self.r50, self.r50e, self.r80, \
         self.r80e, self.r90, self.r90e
        rat2080 = self.r80 / self.r20
        self.concen = 5 * np.log10(rat2080)
        rat2080e = rat2080 * np.sqrt((self.r20e / self.r20)**2. + \
                   (self.r80e / self.r80)**2.)
        self.error_con = 5 * 0.434 * rat2080e / rat2080
        print self.concen, self.error_con
f = pyfits.open('n5585_lR.fits')
z = f[0].data
f.close() 
f = pyfits.open('BMask.fits')
mask = f[0].data
f.close() 
z = z - 1390.0
concentration(z, mask, 192.03,157.42, 313, 313, 0, 0, 1390.)

