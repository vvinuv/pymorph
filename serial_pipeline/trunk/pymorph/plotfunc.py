import sys, pyfits
from pylab import *
import numpy as n
import numpy.core.ma as ma

class PlotFunc:
    """The class for plotting. It will plot the galaxy image, the model galaxy
       and the residual image. It also plot the histogram of residual image 
       and compute the Goodness parameter out of it. It is defind as the ratio 
       of number pixels within n times sky sigma around zero value to the total
       number of pixels. The writehtmlfunc.py will check whether the fit is
       good or bad according to some conditions. This will also plot the 1D
       profile of the galaxy and model galaxy from ellipse fitting. The other
       plot it gives is the mask image. After plotting it saves the image
       as png file with name P_string(galid).png"""
    def __init__(self, cutimage, outimage, maskimage, xcntr, ycntr, skysig):
        self.cutimage  = cutimage
        self.outimage  = outimage
        self.maskimage = maskimage
        self.xcntr     = xcntr
        self.ycntr     = ycntr
        self.skysig    = skysig
        self.plot_profile = plot_profile(cutimage, outimage, maskimage, xcntr, ycntr, skysig)
        return
def get_data(ticker):
    """ Returns the values from the ellipse output table"""
    class C: pass
    def get_ticker(ticker):
        vals = []
        lines = file( '%s' % ticker ).readlines()
        for line in lines[1:]:
            try:
                vals.append([float(val) for val in line.split()[0:]])
            except:
                pass
        M = array(vals)
        c = C()
        c.sma = M[:,0]
        c.flux = M[:,1]
        c.flux_err = M[:,2]
        c.mag = M[:,3]
        c.mag_uerr = M[:,4]
        c.mag_lerr = M[:,5]
        return c
    c1 = get_ticker(ticker)
    return c1

def plot_profile(cutimage, outimage, maskimage, xcntr, ycntr, skysig):
    try:
        #Read the GALFIT output
        f=pyfits.open(outimage)
        galaxy = f[1].data 
        model = f[2].data
        residual = f[3].data
        f.close()
        anorm = normalize(galaxy.min(), 3.0*skysig) #Color normalization
        residual0 = residual
        #Read mask
        f_mask = pyfits.open(maskimage)
        mask = f_mask[0].data 
        f_mask.close()
        size = galaxy.shape[0]
        maskedModel = ma.masked_array(model, mask)
        model = ma.filled(maskedModel, value=9999)
        #The calculations for Goodness is starting here
        maskedresidual = ma.masked_array(residual, mask)
        residual = ma.filled(maskedresidual, value=9999)
        #colorbar(shrink=0.90)
        valid_pixels = ma.count(maskedresidual)
        print 'valid_pixels', valid_pixels
        pixels_in_skysig = residual[n.where(abs(residual) <= skysig / 2.0)].size
        print 'pixels_in_skysig', pixels_in_skysig
        goodness = pixels_in_skysig / float(valid_pixels)
        hist_mask = n.zeros((size, size))
        hist_mask[n.where(abs(residual) > 2.0 * skysig)] = 1
        hist_res = ma.masked_array(residual, hist_mask)
        #The procedure for chi2nu wrt radius is starting here
        x = n.reshape(n.arange(size*size),(size,size)) % size
        x = x.astype(n.float32)
        y = n.reshape(n.arange(size*size),(size,size)) / size
        y = y.astype(n.float32)
        tx = x - xcntr + 0.5
        ty = y - ycntr + 0.5
        R = n.sqrt(tx**2.0 + ty**2.0)
        Chi2Nu = []
        Chi2NuRad = []
        StartRad = 2.0
        TempRad = 0.0
        while StartRad <= size:
#            Chi2NuEle = (ma.sum(abs(maskedresidual[n.where(R <= \
#                         StartRad)])) - ma.sum(abs(maskedresidual[n.where(R <= \
#                         TempRad)])))**2.0 / \
#                         ((galaxy[n.where(R <= StartRad)].sum() - \
#                         galaxy[n.where(R <= TempRad)].sum()) * \
#                         (ma.count(maskedresidual[n.where(R <= StartRad)]) -\
#                          ma.count(maskedresidual[n.where(R <= TempRad)])))
            Chi2NuEle = (ma.sum(abs(maskedresidual[n.where(R <= \
                         StartRad)])) - ma.sum(abs(maskedresidual[n.where(R <= \
                         TempRad)])))**2.0 / \
                         (galaxy[n.where(R <= StartRad)].sum() - \
                         galaxy[n.where(R <= TempRad)].sum())
            try:
                Chi2Nu.append(float(Chi2NuEle))
                Chi2NuRad.append(StartRad)
            except:
                pass
            TempRad = StartRad
            StartRad += 1.0
    except:
        pass
    try:
        data = get_data('E_' + str(cutimage)[:-4] + 'txt')
        sma = data.sma		#sma from ellise fitting
        flux = data.flux	#Flux at various sma
        flux_err =data.flux_err	#Error in Flux
        mag = data.mag		#Magnitude at various sma
        mag_uerr = data.mag_uerr	#Upper error in magnitude
        mag_lerr = data.mag_lerr	#lower error in Magnitude
        GalEll = 1
    except:
        GalEll = 0
        pass
    try:
        data1 = get_data('OE_' + str(cutimage)[:-4] + 'txt')
        sma1 = data1.sma		#sma from ellise fitting
        flux1 = data1.flux	#Flux at various sma
        flux_err1 =data1.flux_err	#Error in Flux
        mag1 = data1.mag		#Magnitude at various sma
        mag_uerr1 = data1.mag_uerr	#Upper error in magnitude
        mag_lerr1 = data1.mag_lerr	#lower error in Magnitude
        ModelEll = 1
    except:
        ModelEll = 0
        pass
    try:
        SmaCommon = []
        MagDev = []
        MagLErr = []
        MagUErr = []
        for i in range(len(sma)):
            for j in range(len(sma1)):
                if sma[i] == sma1[j]:
                    try:
                        SmaCommon.append(sma[i])
                        MagDev.append(mag[i] - mag1[j])
                        FluxErr = n.sqrt((flux_err[i] / flux[i])**2.0 + \
                                  (flux_err1[j]/flux1[j])**2.0)
                        MagLErr.append((n.log10(flux[i]/flux1[j]) - \
                                        n.log10((flux[i]/flux1[j]) - FluxErr)) \
                                        * -2.5)
                        MagUErr.append((n.log10((flux[i]/flux1[j]) + FluxErr) -\
                                       (n.log10(flux[i]/flux1[j]))) * -2.5) 
                    except:
                        pass
    except:
        pass
    #Plotting Starts
    figure(1, figsize=(8.0, 4.6))
    rect1 = [0.125, 0.5, 0.225, 1.5*0.225]
    rect2 = [3.25*0.125, 0.5, 0.225, 1.5*0.225]
    rect3 = [5.5*0.125, 0.5, 0.225, 1.5*0.225]
    rect4 = [5.5*0.125, 0.075, 0.225, 1.5*0.225]
    rect5 = [3.25*0.125, 0.075, 0.225, 1.5*0.225]
#    rect6 = [0.125, 0.75*0.225 + 0.075, 0.225, 1.5*0.225]
    rect7 = [0.125, 0.075, 0.225, 1.5*0.225]
    try:
        axUL = axes(rect1)
        image1 = imshow(n.flipud(galaxy), \
                        extent=[0, size, 0, size], norm=anorm)
        title('Original Galaxy')
        axUM = axes(rect2)
        image1 = imshow(n.flipud(model), cmap=cm.jet, \
                        extent=[0, size, 0, size], norm=anorm)
        title('Model Galaxy + Mask')
        axUR = axes(rect3)
        image1 = imshow(n.flipud(residual0), cmap=cm.jet, \
                        extent=[0, size, 0, size], norm=anorm)
        title('Residual')
        axLR = axes(rect4)
        nn, bins, patches = hist(hist_res, 50, normed=0)
        nMaxArg = nn.argmax()
        if(nMaxArg < 16):
            ArgInc = nMaxArg
        else:
            ArgInc = 16
#    print trapz(bins, nn)
        nMax = max(nn) 
        binmin = bins[nMaxArg - ArgInc]
        binmax = bins[nMaxArg + ArgInc]
        axis([binmin, binmax, 0.0, nMax])
        setp(patches, 'facecolor', 'g', 'alpha', 0.75)
        xticks((binmin, binmin /2.0, 0.0, binmax/2.0, binmax),\
                (str(binmin)[:5], str(binmin /2.0)[:5], str(0.0)[:3], \
                 str(binmax/2.0)[:5], str(binmax)[:5]))
        grid(True)
        title('Difference Histogram')
        Dx = abs(axLR.get_xlim()[0]-axLR.get_xlim()[1])
        Dy = abs(axLR.get_ylim()[0]-axLR.get_ylim()[1])
        axLR.set_aspect(Dx/Dy)
    except:
        pass
    try:
        if GalEll and ModelEll:
            ymin = min(min(mag), min(mag1))
            ymax = max(max(mag), max(mag1))
            xmax = max(max(sma), max(sma1)) 
            axLL = axes(rect7)
            errorbar(sma, mag, [mag_uerr,mag_lerr], fmt='o',ecolor='r', ms=3)
            plot(sma1, mag1,color='g',lw=2)
            axLL.set_ylim(ymax, ymin)
            axLL.set_xlim(0, xmax)
            xlabel(r'Radius', size='medium')
            ylabel(r'Surface Brightness', size='medium')
            title('1-D Profile Comparison')
            grid(True)
            Dx = abs(axLL.get_xlim()[0]-axLL.get_xlim()[1])
            Dy = abs(axLL.get_ylim()[0]-axLL.get_ylim()[1])
            axLL.set_aspect(Dx/Dy)
            axLM = axes(rect5)
            try:
#                xmax = max(max(SmaCommon), max(Chi2NuRad))
                errorbar(SmaCommon, MagDev, [MagUErr, MagLErr], \
                         fmt='o',ecolor='r', ms=3)
                ylabel('Magnitude Deviation', size='medium')
                xlabel('Radius', size='medium')
                grid(True)
                Dx = abs(axLM.get_xlim()[0]-axLM.get_xlim()[1])
                Dy = abs(axLM.get_ylim()[0]-axLM.get_ylim()[1])
                axLM.set_aspect(Dx/Dy)
            except:
                pass
        elif GalEll:
            ymin = min(mag)
            ymax = max(mag)
            xmax = max(sma)
            axLL = axes(rect7)
            errorbar(sma, mag, [mag_uerr,mag_lerr], fmt='o',ecolor='r', ms=3)
            axLL.set_ylim(ymax, ymin)
            axLL.set_xlim(0, xmax)
            xlabel(r'Radius', size='medium')
            ylabel(r'Surface Brightness', size='medium')
            title('1-D Profile Comparison')
            grid(True)
            Dx = abs(axLL.get_xlim()[0]-axLL.get_xlim()[1])
            Dy = abs(axLL.get_ylim()[0]-axLL.get_ylim()[1])
            axLL.set_aspect(Dx/Dy)
        elif ModelEll:
            ymin = min(mag1)
            ymax = max(mag1)
            xmax = max(sma1)
            axLL = axes(rect7)
            plot(sma1, mag1,color='g',lw=2)
            axLL.set_ylim(ymax, ymin)
            axLL.set_xlim(0, xmax)
            xlabel(r'Radius', size='medium')
            ylabel(r'Surface Brightness', size='medium')
            title('1-D Profile Comparison')
            grid(True)
            Dx = abs(axLL.get_xlim()[0]-axLL.get_xlim()[1])
            Dy = abs(axLL.get_ylim()[0]-axLL.get_ylim()[1])
            axLL.set_aspect(Dx/Dy)
    except:
        pass     
#    show()
    savefig('P_' + str(cutimage)[:-4] + 'png')
    close()
    return goodness
#PlotFunc('LFC1208I_1038.fits', 'O_LFC1208I_1038.fits', 'M_LFC1208I_1038.fits', 40, 40, 0.02)
