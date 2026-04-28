import numpy as np
from astropy.io import fits
from astropy.convolution import convolve, Gaussian2DKernel


def return_sigma(image_data, outname='W.fits',  gain=4.71, ncombine=1):
    '''
    This is exactly return sigma image as GALFIT in sigma.c  
    '''
    # To get the 20 and 80 percentage of counts 
    ql = 0.2
    qu = 0.8
    frac = qu - ql
    gaus_radius = 5
    gaus_sigma = 2
    detect_sigma = 5

    #gain = 4.71
    #ncombine = 3
    effgain = gain * ncombine

    id_sort = np.quantile(image_data, (0.2, 0.8))

    image_data_quant = image_data[(image_data >= id_sort[0]) & (image_data <= id_sort[1])]

    med = np.median(image_data)


    #image_data.size * 0.6 means that the shape of image_data_quant which have only 60% of all points because of the quantile
    std_sky = np.sqrt(np.sum((image_data_quant - med)**2) / (frac * image_data.size - 1))

    avg = np.mean(image_data[abs(image_data - med) <= detect_sigma * std_sky])

    std_sky2 = np.std(image_data[abs(image_data - med) <= detect_sigma * std_sky])


    kernel = Gaussian2DKernel(x_stddev=gaus_sigma, x_size=gaus_radius, 
                              y_size=gaus_radius)
    image_data_smooth = convolve(image_data, kernel, nan_treatment='interpolate')

    #image_data_smooth = gaussian_filter(image_data, sigma=gaus_sigma, radius=gaus_radius)

    #Here varp is the variance as it is the squre of sqrt(count * effgain).
    #i.e. count * gain
    varp = (image_data_smooth - avg) * effgain
    #Below varp is rms
    varp[varp >= 0] = np.sqrt(varp[varp >= 0] / effgain / effgain + 1 * std_sky2 * std_sky2)
    varp[varp < 0] = std_sky2
    std = varp.copy()
    fits.writeto(outname, varp, overwrite=True)
    #return std

if __name__  == '__main__':
    ffits = "/Users/vinu/github/pymorph/examples/small_image/data/Icl1358_9.fits"
    hdul = fits.open(ffits)
    img_data = hdul[0].data
    header = hdul[0].header


    gain = header['GAIN']
    ncombine = header['NCOMBINE']
    varp = return_sigma(img_data, gain=gain, ncombine=ncombine)

    fits.writeto("Isigma.fits", varp, overwrite=True)

