import pyfits
import numpy as np
from numpy import fft
def Convolve(image, kernal, zeropad=True):
    """ Convolution using numpy fft2. zeropad=True always in this function"""
    kernal = np.asarray(kernal)
    kernal = kernal.astype(np.float32) 
    kernal = kernal / kernal.sum() # Normalizing the kernal
    # The shape of the input including zero pading
    s0 = image.shape[0] + kernal.shape[0] - 1
    s1 = image.shape[1] + kernal.shape[1] - 1
    # The extra size
    es0 = (kernal.shape[0] - 1) / 2. 
    es1 = (kernal.shape[1] - 1) / 2.
    # The product of fourier transform of image and kernal. [s0, s1] is the 
    # dimension of the image with zero pading
    FouTra = np.multiply(fft.fft2(image, [s0, s1]), fft.fft2(kernal,[s0, s1]))
    # The inverse fourier transform of the products of trasforms gives the
    # convolved image
    ConvIm = fft.ifft2(FouTra).real
    # Removing the zero padded region
    ymin = np.floor(es0)
    ymax = np.floor(s0 - es0)
    xmin = np.floor(es1)
    xmax = np.floor(s1 - es1)
    ConvIm = ConvIm[ymin:ymax, xmin:xmax]
    return ConvIm
