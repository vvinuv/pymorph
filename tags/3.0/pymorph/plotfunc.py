import sys, pyfits
import numpy as np
import config as c
from pylab import *
import numpy.ma as ma

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
    def __init__(self, outimage, maskimage, xcntr, ycntr, sky, skysig, save_name = '999'):
        self.outimage  = outimage
        self.maskimage = maskimage
        self.xcntr     = xcntr
        self.ycntr     = ycntr
        self.sky       = sky
        self.skysig    = skysig
        if save_name == '999':
            save_name = 'P_' + c.fstring + '.png'
        self.plot_profile = plot_profile(outimage, maskimage, xcntr, ycntr, sky, skysig, save_name)
        return

def plot_profile(outimage, maskimage, xcntr, ycntr, sky, skysig, save_name):
    try:
        # Read the GALFIT output
        f=pyfits.open(outimage)
        galaxy = f[1].data 
        model = f[2].data
        residual = f[3].data
        f.close()
        anorm = normalize(sky - 2*skysig , sky + 12.0*skysig) # Color 
                                                              # normalization
        # Read mask
        f_mask = pyfits.open(maskimage)
        mask = f_mask[0].data 
        f_mask.close()
        maskedModel = ma.masked_array(model, mask)
        model = ma.filled(maskedModel, 9999)
        #The calculations for Goodness is starting here
        maskedresidual = ma.masked_array(residual, mask)
        anormRes = normalize(-2 * skysig, 3 * skysig) 
        residual = ma.filled(maskedresidual, 9999)
        #colorbar(shrink=0.90)
        valid_pixels = ma.count(maskedresidual)
        print 'No of valid pixels >>> ', valid_pixels
        pixels_in_skysig = residual[where(abs(residual) <= skysig)].size
        print 'No of pixels within sky sigma >>> ', pixels_in_skysig
        try:
            goodness = pixels_in_skysig / float(valid_pixels)
        except:
            goodness = 9999
        NYPTS, NXPTS = galaxy.shape
        hist_mask = zeros((NXPTS, NYPTS))
        hist_mask[where(abs(residual) > 12.0 * skysig)] = 1
        hist_res = ma.masked_array(residual, hist_mask)
    except:
        pass
    try:
        data = np.genfromtxt('E_' + c.fstring + '.txt', delimiter=' ', \
               names=True)
        sma = data['sma']		#sma from ellise fitting
        flux = data['inte']	#Flux at various sma
        flux_err =data['intee']	#Error in Flux
        mag = data['mag'] + float(c.mag_zero) #Magnitude at various sma
        mag_uerr = data['magu']	#Upper error in magnitude
        mag_lerr = data['magl']	#lower error in Magnitude
        GalEll = 1
    except:
        GalEll = 0
    try:
        data1 = np.genfromtxt('OE_' + c.fstring + '.txt', delimiter=' ', \
               names=True)
        sma1 = data1['sma']		#sma from ellise fitting
        flux1 = data1['inte']	#Flux at various sma
        flux_err1 =data1['intee']	#Error in Flux
        mag1 = data1['mag'] + float(c.mag_zero)#Magnitude at various sma
        mag_uerr1 = data1['magu']	#Upper error in magnitude
        mag_lerr1 = data1['magl']	#lower error in Magnitude
        ModelEll = 1
    except:
        ModelEll = 0
    try:
        try:
            MaxRad = sma.max()
        except:
            MaxRad = sma1.max()
        NoOfPoints = int(30 * MaxRad / 50.)
        SmaCommon = np.logspace(0, np.log10(MaxRad), \
                                NoOfPoints, endpoint=True)
        MagI = np.interp(SmaCommon, sma, mag)
        MagI1 = np.interp(SmaCommon, sma1, mag1)
        FluxI = np.interp(SmaCommon, sma, flux)
        FluxI1 = np.interp(SmaCommon, sma1, flux1)
        FluxErrI = np.interp(SmaCommon, sma, flux_err)
        FluxErrI1 = np.interp(SmaCommon, sma1, flux_err1)
        MagDev = MagI - MagI1
        FluxErr = np.sqrt((FluxErrI / FluxI)**2.0 + (FluxErrI1 / FluxI1)**2.0)
        MagLErr = (np.log10(FluxI / FluxI1) - \
                        np.log10((FluxI / FluxI1) - FluxErr)) * -2.5
        MagUErr = (np.log10((FluxI / FluxI1) + FluxErr) -\
                       (np.log10(FluxI / FluxI1))) * -2.5
    except:
        pass
    #Plotting Starts
    FigSize = [8.0, 4.6]
    MatPlotParams = {'axes.titlesize': 10,
                     'axes.labelsize': 10,
                     'xtick.labelsize': 5,
                     'ytick.labelsize': 5,
                     'figure.figsize': FigSize}
    rcParams.update(MatPlotParams)
    rect1 = [0.125, 0.5, 0.225, 1.5*0.225]
    rect2 = [3.25*0.125, 0.5, 0.225, 1.5*0.225]
    rect3 = [5.5*0.125, 0.5, 0.225, 1.5*0.225]
    rect4 = [5.5*0.125, 0.075, 0.225, 1.5*0.225]
    rect5 = [3.25*0.125, 0.075, 0.225, 1.5*0.225]
#    rect6 = [0.125, 0.75*0.225 + 0.075, 0.225, 1.5*0.225]
    rect7 = [0.125, 0.075, 0.225, 1.5*0.225]
    try:
        axUL = axes(rect1)
        image1 = imshow(galaxy, extent=[0, NXPTS, 0, NYPTS], norm=anorm)
        Dx = abs(axUL.get_xlim()[0]-axUL.get_xlim()[1])
        Dy = abs(axUL.get_ylim()[0]-axUL.get_ylim()[1])
        colorbar(shrink=0.9, format='%.2f')
        axUL.set_aspect(Dx/Dy)

        title('Original Galaxy')
        axUM = axes(rect2)
        image1 = imshow(model, cmap=cm.jet, \
                        extent=[0, NXPTS, 0, NYPTS], norm=anorm)
        colorbar(shrink=0.9, format='%.2f')
        Dx = abs(axUM.get_xlim()[0]-axUM.get_xlim()[1])
        Dy = abs(axUM.get_ylim()[0]-axUM.get_ylim()[1])
        axUM.set_aspect(Dx/Dy)

        title('Model Galaxy + Mask')
        axUR = axes(rect3)
        image1 = imshow(residual, cmap=cm.jet, \
                        extent=[0, NXPTS, 0, NYPTS], norm=anormRes)
        colorbar(shrink=0.9)
        Dx = abs(axUR.get_xlim()[0]-axUR.get_xlim()[1])
        Dy = abs(axUR.get_ylim()[0]-axUR.get_ylim()[1])
        axUR.set_aspect(Dx/Dy)

        title('Residual')
        axLR = axes(rect4)
	hist_res1d = hist_res.compressed()
        nn, bins, patches = hist(hist_res1d, 50, normed=0)
        nMaxArg = nn.argmax() 
        if(nMaxArg < 16):
            ArgInc = nMaxArg
        else:
            ArgInc = 16
        # print trapz(bins, nn)
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
            xlabel(r'Radius')
            ylabel(r'Surface Brightness')
            title('1-D Profile Comparison')
            grid(True)
            Dx = abs(axLL.get_xlim()[0]-axLL.get_xlim()[1])
            Dy = abs(axLL.get_ylim()[0]-axLL.get_ylim()[1])
            axLL.set_aspect(Dx/Dy)
            axLM = axes(rect5)
            try:
                errorbar(SmaCommon, MagDev, [MagUErr, MagLErr], \
                         fmt='o',ecolor='r', ms=3)
                ylabel('Magnitude Deviation')
                xlabel('Radius')
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
            xlabel(r'Radius')
            ylabel(r'Surface Brightness')
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
            xlabel(r'Radius')
            ylabel(r'Surface Brightness')
            title('1-D Profile Comparison')
            grid(True)
            Dx = abs(axLL.get_xlim()[0]-axLL.get_xlim()[1])
            Dy = abs(axLL.get_ylim()[0]-axLL.get_ylim()[1])
            axLL.set_aspect(Dx/Dy)
    except:
        pass     
#    show()
    savefig(save_name)
    close()
    return goodness
#c.fstring = 'test_n5585_lR'
#PlotFunc('O_test_n5585_lR.fits', 'M_test_n5585_lR.fits', 82, 82, 1434., 40.)
