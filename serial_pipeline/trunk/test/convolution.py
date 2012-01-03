from numpy.fft import fft2, ifft2
import pyfits
from numpy import log, array
import os
import convolve as conv
def Convolve(image1, image2, MinPad=True, pad=True):
    """ Not so simple convolution """

    #Just for comfort:
    FFt = fft2
    iFFt = ifft2

    #The size of the images:
    r1,c1 = image1.shape
    r2,c2 = image2.shape

    #MinPad results simpler padding,smaller images:
    if MinPad:
        r = r1+r2
        c = c1+c2
    else:
        #if the Numerical Recipies says so:
        r = 2*max(r1,r2)
        c = 2*max(c1,c2)
    
    #For nice FFT, we need the power of 2:
    if pad:
        pr2 = int(log(r)/log(2.0) + 1.0 )
        pc2 = int(log(c)/log(2.0) + 1.0 )
        rOrig = r
        cOrig = c
        r = 2**pr2
        c = 2**pc2
    #end of if pad
    
    #numpy fft has the padding built in, which can save us some steps
    #here. The thing is the s(hape) parameter:
    fftimage = FFt(image1, s=(r,c))*FFt(image2[::-1,::-1],s=(r,c))
#    fftimage = FFt(image1,s=(r,c)) * FFt(image2,s=(r,c))

    if pad:
        return (iFFt(fftimage))[:rOrig,:cOrig].real
    else:
        return (iFFt(fftimage)).real

f = pyfits.open('n5585_lR.fits')
z = f[0].data
f.close()
ker = array([[1, 1], [1, 1]])
zc = Convolve(z, ker) / 4.0
hdu = pyfits.PrimaryHDU(zc)
os.system('rm -f conv.fits conv1.fits')
hdu.writeto('conv.fits')
zc = conv.boxcar(z, (2,2))
hdu = pyfits.PrimaryHDU(zc)
hdu.writeto('conv1.fits')

