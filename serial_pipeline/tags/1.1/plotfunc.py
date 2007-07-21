import sys, pyfits
from pylab import *
from matplotlib.numerix import fromstring, argsort, take, array, resize
import random, numarray.random_array,numarray.mlab


class PlotFunc:
    """The class for plotting"""
    def __init__(self, files):
        self.files        = files
        self.plot_profile = plot_profile(files)

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

def plot_profile(files):
    data = get_data('elli_' + str(files)[6:-4] + 'txt')
    data1 = get_data('out_elli_' + str(files)[6:-4] + 'txt')
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
    sc1=subplot(224)
    #sc1.scatter(sma, mag, s=10, alpha=0.75, c='r')
    sc1.errorbar(sma, mag, [mag_uerr,mag_lerr], fmt='o',ecolor='r', ms=3) 
    #fmt point type, ecolor point color, ms, point size
    #sc1.errorbar(sma1, mag1, [mag_uerr1,mag_lerr1], fmt='o',ecolor='g', ms=3)
    sc1.plot(sma1, mag1,color='g',lw=2)
    ymin = min(min(mag), min(mag1))
    ymax = max(max(mag), max(mag1))
    sc1.set_ylim(ymax, ymin)
    Dx=abs(sc1.get_xlim()[0]-sc1.get_xlim()[1])
    Dy=abs(sc1.get_ylim()[0]-sc1.get_ylim()[1])
    sc1.set_aspect(Dx/Dy)
    xlabel(r'Radius', size='medium')
    ylabel(r'Surface Brightness', size='medium')
    title('1-D Profile Comparison')
    grid(True)
    #savefig('plot_' + str(files)[6:-4] + 'png')
    #colorbar()
    #show()
    f=pyfits.open('out_' + str(files)[6:-4] + 'fits')
    galaxy 		= f[1].data 
    model  		= f[2].data
    residual	= f[3].data
    f.close()
    subplot(221)
    anorm = normalize(0.0,.04)
    image1=imshow(numarray.mlab.rot90(numarray.swapaxes(galaxy,0,1)),norm=anorm)
    image1.autoscale()
    subplot(222)
    image1=imshow(numarray.mlab.rot90(numarray.swapaxes(model,0,1)),norm=anorm)
    image1.autoscale()
    subplot(223)
    image1=imshow(numarray.mlab.rot90(numarray.swapaxes(residual,0,1)),norm=anorm)
    image1.autoscale()
    savefig('plot_' + str(files)[6:-4] + 'png')
    figure()
