import numpy as np
import fitsio

def return_sigma(image_data, gain=4.71, ncombine=1):
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


    image_data_smooth = gaussian_filter(image_data, sigma=gaus_sigma, radius=gaus_radius)

    #Here varp is the variance as it is the squre of sqrt(count * effgain).
    #i.e. count * gain
    varp = (image_data_smooth - avg) * effgain
    varp[varp >= 0] = np.sqrt(varp[varp >= 0] / effgain / effgain + 1 * std_sky2 * std_sky2)
    varp[varp < 0] = std_sky2

    return varp

if __name__  == '__main__':
    root = '../examples/small_image/data'
    header = fitsio.read_header(f'{root}/I.fits')
    gain = header['GAIN']
    ncombine = header['NCOMBINE']
    image_data = fitsio.FITS(f'{root}/I.fits')[0].read()
    varp = return_sigma(image_data, gain=gain, ncombine=ncombine)

    os.system(f'/bin/rm -f {root}/Isigma.fits')
    f = fitsio.FITS(f'{rmroot}/Isigma.fits', 'rw')
    f.write(varp)
    f.close()

