import sys
import pyfits
import numpy as np
import config as c
from pylab import *
import numpy.ma as ma

class PlotFunc:
    """The class for plotting galaxy images, models, and residuals."""
    def __init__(self, outimage, maskimage, xcntr, ycntr, sky, skysig, save_name='999'):
        self.outimage = outimage
        self.maskimage = maskimage
        self.xcntr = xcntr
        self.ycntr = ycntr
        self.sky = sky
        self.skysig = skysig
        if save_name == '999':
            save_name = 'P_' + c.fstring + '.png'
        self.plot_profile = plot_profile(outimage, maskimage, xcntr, ycntr, sky, skysig, save_name)
        return

def plot_profile(outimage, maskimage, xcntr, ycntr, sky, skysig, save_name):
    goodness = 9999
    try:
        # Read the GALFIT output
        f = pyfits.open(outimage)
        galaxy = f[1].data 
        model = f[2].data
        residual = f[3].data
        f.close()
        
        anorm = normalize(sky - 2*skysig, sky + 12.0*skysig)
        
        # Read mask
        f_mask = pyfits.open(maskimage)
        mask = f_mask[0].data 
        f_mask.close()
        
        maskedModel = ma.masked_array(model, mask)
        model = ma.filled(maskedModel, 9999)
        
        # Calculations for Goodness
        maskedresidual = ma.masked_array(residual, mask)
        anormRes = normalize(-2 * skysig, 3 * skysig) 
        residual_filled = ma.filled(maskedresidual, 9999)
        
        valid_pixels = ma.count(maskedresidual)
        pixels_in_skysig = residual_filled[where(abs(residual_filled) <= skysig)].size
        
        try:
            goodness = pixels_in_skysig / float(valid_pixels)
        except:
            goodness = 9999
            
        NYPTS, NXPTS = galaxy.shape
        hist_mask = zeros((NYPTS, NXPTS)) # Note: galaxy.shape is (Y, X)
        hist_mask[where(abs(residual_filled) > 12.0 * skysig)] = 1
        hist_res = ma.masked_array(residual_filled, hist_mask)
    except Exception as e:
        print(f"Error in data reading/processing: {e}")
        pass

    # Profile Data Loading
    GalEll, ModelEll = 0, 0
    try:
        data = np.genfromtxt('E_' + c.fstring + '.txt', delimiter=' ', names=True)
        sma, flux, flux_err = data['sma'], data['inte'], data['intee']
        mag = data['mag'] + float(c.mag_zero)
        mag_uerr, mag_lerr = data['magu'], data['magl']
        GalEll = 1
    except: pass

    try:
        data1 = np.genfromtxt('OE_' + c.fstring + '.txt', delimiter=' ', names=True)
        sma1, flux1, flux_err1 = data1['sma'], data1['inte'], data1['intee']
        mag1 = data1['mag1'] + float(c.mag_zero) if 'mag1' in data1.dtype.names else data1['mag'] + float(c.mag_zero)
        mag_uerr1, mag_lerr1 = data1['magu'], data1['magl']
        ModelEll = 1
    except: pass

    # Interp logic
    try:
        MaxRad = max(sma.max() if GalEll else 0, sma1.max() if ModelEll else 0)
        NoOfPoints = int(30 * MaxRad / 50.)
        SmaCommon = np.logspace(0, np.log10(MaxRad), NoOfPoints, endpoint=True)
        MagI = np.interp(SmaCommon, sma, mag)
        MagI1 = np.interp(SmaCommon, sma1, mag1)
        FluxI = np.interp(SmaCommon, sma, flux)
        FluxI1 = np.interp(SmaCommon, sma1, flux1)
        FluxErrI = np.interp(SmaCommon, sma, flux_err)
        FluxErrI1 = np.interp(SmaCommon, sma1, flux_err1)
        MagDev = MagI - MagI1
        FluxErr = np.sqrt((FluxErrI / FluxI)**2.0 + (FluxErrI1 / FluxI1)**2.0)
        MagLErr = (np.log10(FluxI / FluxI1) - np.log10((FluxI / FluxI1) - FluxErr)) * -2.5
        MagUErr = (np.log10((FluxI / FluxI1) + FluxErr) - (np.log10(FluxI / FluxI1))) * -2.5
    except: pass

    # Plotting Starts
    FigSize = [8.0, 4.6]
    MatPlotParams = {'axes.titlesize': 10, 'axes.labelsize': 10, 'xtick.labelsize': 5, 'ytick.labelsize': 5, 'figure.figsize': FigSize}
    rcParams.update(MatPlotParams)
    
    rect1, rect2, rect3 = [0.125, 0.5, 0.225, 0.337], [0.406, 0.5, 0.225, 0.337], [0.687, 0.5, 0.225, 0.337]
    rect4, rect5, rect7 = [0.687, 0.075, 0.225, 0.337], [0.406, 0.075, 0.225, 0.337], [0.125, 0.075, 0.225, 0.337]

    try:
        # Original Galaxy
        axUL = axes(rect1)
        imshow(galaxy, extent=[0, NXPTS, 0, NYPTS], norm=anorm)
        colorbar(shrink=0.9, format='%.2f')
        title('Original Galaxy')

        # Model
        axUM = axes(rect2)
        imshow(model, cmap=cm.jet, extent=[0, NXPTS, 0, NYPTS], norm=anorm)
        colorbar(shrink=0.9, format='%.2f')
        title('Model Galaxy + Mask')

        # Residual
        axUR = axes(rect3)
        imshow(residual_filled, cmap=cm.jet, extent=[0, NXPTS, 0, NYPTS], norm=anormRes)
        colorbar(shrink=0.9)
        title('Residual')

        # Histogram
        axLR = axes(rect4)
        # FIXED: Moved inside the try block
        hist_res1d = hist_res.compressed()
        nn, bins, patches = hist(hist_res1d, 50, density=False) 
        nMaxArg = nn.argmax() 
        ArgInc = min(nMaxArg, 16)
        nMax = max(nn) 
        axis([bins[max(0, nMaxArg - ArgInc)], bins[min(len(bins)-1, nMaxArg + ArgInc)], 0.0, nMax])
        setp(patches, 'facecolor', 'g', 'alpha', 0.75)
        title('Difference Histogram')
        grid(True)
    except Exception as e:
        print(f"Error in plotting top row: {e}")

    try:
        if GalEll or ModelEll:
            axLL = axes(rect7)
            if GalEll:
                errorbar(sma, mag, [mag_uerr, mag_lerr], fmt='o', ecolor='r', ms=3)
            if ModelEll:
                plot(sma1, mag1, color='g', lw=2)
            
            # Surface Brightness convention: y-axis is inverted
            all_mag = (mag if GalEll else []) + (mag1 if ModelEll else [])
            axLL.set_ylim(max(all_mag), min(all_mag))
            xlabel(r'Radius')
            ylabel(r'Surface Brightness')
            title('1-D Profile Comparison')
            grid(True)

            if GalEll and ModelEll:
                axLM = axes(rect5)
                errorbar(SmaCommon, MagDev, [MagUErr, MagLErr], fmt='o', ecolor='r', ms=3)
                ylabel('Mag Deviation')
                xlabel('Radius')
                grid(True)
    except Exception as e:
        print(f"Error in plotting profiles: {e}")

    savefig(save_name)
    close()
    return goodness