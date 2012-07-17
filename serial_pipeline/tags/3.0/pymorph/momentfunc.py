import pyfits
import numpy as np
import numpy.ma as ma

class moment:
    """This clas is for moment calculation. The algorithm used here is as 
       follows

       1. The total second-order moment M_tot is the flux in each pixel 
          f_i multiplied by the squared distance to the center of the 
          galaxy, summed over all the galaxy pixels assigned by the 
          segmentation map.

          M_tot = Sum(f_i * [(x_i - x_c)^2 + (y_i - y_c)^2])

          Where xc, yc is the galaxy's center.

       2. The center is computed by finding xc, yc such that M_tot is 
          minimized.

       3. Define M20 as the brightest 20% of the galaxy's flux.

       4. To compute M20, sort the pixels by flux, sum M_i over the 
          brightest pixels until the sum of the brightest pixels equals 20%
          of the total galaxy flux, and then normalize by M_tot.
    """
    def __init__(self, z, ixcntr, iycntr):
        self.z               = z
        self.ixcntr          = ixcntr
        self.iycntr          = iycntr
        self.moment_of_light = moment_of_light(self.z, self.ixcntr, self.iycntr)


def moment_of_light(zextract, ixcntr, iycntr):
    flag_center = 0
    center_rad = 1.0
    nn = 0
    SizeY = zextract.shape[0] # Y size
    SizeX = zextract.shape[1] # X size
    x = np.reshape(np.arange(SizeX * SizeY), (SizeY, SizeX)) % SizeX
    x = x.astype(np.float32)
    y = np.reshape(np.arange(SizeX * SizeY), (SizeY, SizeX)) / SizeX
    y = y.astype(np.float32)
    M_tot = np.zeros([9])
    M_tot = M_tot.astype(np.float32)
    while flag_center == 0:
        xcntr = np.reshape(np.array([[ixcntr] * 3, [ixcntr - center_rad] * 3, \
                [ixcntr + center_rad] * 3]), (9, ))
        xcntr = xcntr.astype(np.float32) 
        ycntr = np.array([iycntr, iycntr - center_rad, iycntr + center_rad] * 3)
        ycntr = ycntr.astype(np.float32)
        for iter in range(9):
            tx = x - xcntr[iter]
            ty = y - ycntr[iter]
            R = tx**2.0 + ty**2.0
            RZ = R * zextract
            M_tot[iter] = RZ.sum()
        M_tot_min = M_tot.min()
        index = M_tot.argmin()
        if(nn > 100):
            xcntr[index] = ixcntr
            ycntr[index] = iycntr
        if(xcntr[index] == ixcntr and ycntr[index] == iycntr):
            flag_center = 1
        else:
            ixcntr = xcntr[index]
            iycntr = ycntr[index]
            nn += 1
    oneD_I = zextract.flat
    sorted_I = np.sort(oneD_I)
    argument = np.argsort(oneD_I)
    tx = x - xcntr[index]
    ty = y - ycntr[index]
    R = tx**2.0 + ty**2.0
    RZ = R * zextract
    oneD_RZ = RZ.flat
    oneD_RZ = oneD_RZ[argument]
    N = sorted_I.size
    total = sorted_I.sum()
    sum = 0.0
    sumM = 0.0
    i = N - 1
    while sum < 0.2 * total:
        sum += sorted_I[i]
        sumM += oneD_RZ[i]
        i -= 1
    M20 = np.log10(sumM / M_tot_min)
    return M20, M_tot[index], sumM, total, sum, nn, xcntr[index], ycntr[index]
