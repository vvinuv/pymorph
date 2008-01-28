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
    def __init__(self, cutimage, outimage, maskimage, skysig):
        self.cutimage  = cutimage
        self.outimage  = outimage
        self.maskimage = maskimage
        self.skysig = skysig
        self.plot_profile = plot_profile(cutimage, outimage, maskimage, skysig)
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

def plot_profile(cutimage, outimage, maskimage, skysig):
    data = get_data('E_' + str(cutimage)[:-4] + 'txt')
    data1 = get_data('OE_' + str(cutimage)[:-4] + 'txt')
    sma = data.sma		#sma from ellise fitting
    flux = data.flux	#Flux at various sma
    flux_err =data.flux_err	#Error in Flux
    mag = data.mag		#Magnitude at various sma
    mag_uerr = data.mag_uerr	#Upper error in magnitude
    mag_lerr = data.mag_lerr	#lower error in Magnitude
    sma1 = data1.sma		#sma from ellise fitting
    flux1 = data1.flux	#Flux at various sma
    flux_err1 =data1.flux_err	#Error in Flux
    mag1 = data1.mag		#Magnitude at various sma
    mag_uerr1 = data1.mag_uerr	#Upper error in magnitude
    mag_lerr1 = data1.mag_lerr	#lower error in Magnitude
    sc1=subplot(234)
    #sc1.scatter(sma, mag, s=10, alpha=0.75, c='r')
    sc1.errorbar(sma, mag, [mag_uerr,mag_lerr], fmt='o',ecolor='r', ms=3) 
    #fmt point type, ecolor point color, ms, point size
    #sc1.errorbar(sma1, mag1, [mag_uerr1,mag_lerr1], fmt='o',ecolor='g', ms=3)
    sc1.plot(sma1, mag1,color='g',lw=2)
    ymin = min(min(mag), min(mag1))
    ymax = max(max(mag), max(mag1))
    sc1.set_ylim(ymax, ymin)
    Dx = abs(sc1.get_xlim()[0]-sc1.get_xlim()[1])
    Dy = abs(sc1.get_ylim()[0]-sc1.get_ylim()[1])
    sc1.set_aspect(Dx/Dy)
    xlabel(r'Radius', size='medium')
    ylabel(r'Surface Brightness', size='medium')
    title('1-D Profile Comparison')
    grid(True)
    #savefig('plot_' + str(cutimage)[6:-4] + 'png')
    #colorbar()
    #show()
    f=pyfits.open(outimage)
    galaxy = f[1].data 
    model = f[2].data
    residual = f[3].data
    f.close()
    f_mask = pyfits.open(maskimage)
    mask = f_mask[0].data 
    f_mask.close()
    size = galaxy.shape[0]
    subplot(231)
    title('Original Galaxy')
#    anorm = normalize(galaxy.min() + (galaxy.max()-galaxy.min())/ 60.0, \
#            galaxy.max() - (galaxy.max() - galaxy.min()) / 1.1)
    anorm = normalize(0.01,.2)
#    print skysig, galaxy.min(), galaxy.max()
#    anorm = normalize(0.0, skysig / 2.0)
    image1 = imshow(n.flipud(galaxy), \
                  extent=[0, size, 0, size])
    image1.autoscale()
    subplot(232)
    title('Model Galaxy')
    image1 = imshow(n.flipud(model), cmap=cm.jet, \
                  extent=[0, size, 0, size])
    image1.autoscale()
    subplot(233)
    title('Residual')
    image1 = imshow(n.flipud(residual), cmap=cm.jet, \
                  extent=[0, size, 0, size])
    image1.autoscale()
    maskedresidual = ma.masked_array(residual, mask)
    residual = ma.filled(maskedresidual, value=0.0)
#    a = axes([0, 0, 150, 150], axisbg='y')
    subplot(235)
    nn, bins, patches = hist(residual, 100, normed=0)
    goodness = float(nn[n.where(abs(bins) <= skysig)].sum()) / \
               float(nn.sum())
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
    grid(True)
    title('Difference Histogram')
#    setp(a, xticks=[], yticks=[])
    subplot(236)
    title('Mask')
    image1 = imshow(n.flipud(mask), norm=anorm, cmap=cm.jet, \
                  extent=[0, size, 0, size])
    image1.autoscale()
    savefig('P_' + str(cutimage)[:-4] + 'png')
#    figure()
    close()
    return goodness
