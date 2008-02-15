import pyfits
import numpy as n
import numpy.core.ma as ma

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
	def __init__(self, z, ini_xcntr, ini_ycntr):
		self.z			= z
		self.ini_xcntr		= ini_xcntr
		self.ini_ycntr		= ini_ycntr
		self.moment_of_light	= moment_of_light(self.z, self.ini_xcntr, self.ini_ycntr)


def moment_of_light(zextract, ini_xcntr, ini_ycntr):
	flag_center = 0
	center_rad = 1.0
	nn = 0
	NXPTS = zextract.shape[0]
	NYPTS = zextract.shape[1]
	x = n.reshape(n.arange(NXPTS * NYPTS),(NXPTS, NYPTS)) / NYPTS
	x = x.astype(n.float32)
	y = n.reshape(n.arange(NXPTS * NYPTS),(NXPTS, NYPTS)) % NYPTS
	y = y.astype(n.float32)
	M_tot = n.zeros([9])
        M_tot = M_tot.astype(n.float32)
	while flag_center == 0:
		xcntr = n.reshape(n.array([[ini_xcntr] * 3, [ini_xcntr - \
                        center_rad] * 3, [ini_xcntr + center_rad] * 3]), (9, ))
                xcntr = xcntr.astype(n.float32) 
		ycntr = n.array([ini_ycntr, ini_ycntr - center_rad, ini_ycntr + \
                             center_rad] * 3)
                ycntr = ycntr.astype(n.float32)
		for iter in range(9):
			tx = x - xcntr[iter]
			ty = y - ycntr[iter]
			R = tx**2.0 + ty**2.0
			RZ = R * zextract
			M_tot[iter] = RZ.sum()
		M_tot_min = M_tot.min()
		index = M_tot.argmin()
		if(nn > 100):
			xcntr[index] = ini_xcntr
			ycntr[index] = ini_ycntr
		if(xcntr[index] == ini_xcntr and ycntr[index] == ini_ycntr):
			flag_center = 1
		else:
			ini_xcntr = xcntr[index]
			ini_ycntr = ycntr[index]
			nn += 1
			#print n
	oneD_I = zextract.flat
	sorted_I = n.sort(oneD_I)
	argument = n.argsort(oneD_I)
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
	M20 = n.log10(sumM / M_tot_min)
	return M20, M_tot[index], sumM, total, sum, nn, xcntr[index],ycntr[index]

#f=pyfits.open('n5585_lR.fits')
#z=f[0].data
#header = f[0].header
#if (header.has_key('sky')):
#    sky = header['sky']
#f.close()
#xcntr=192.03
#ycntr=157.42
#pa=0.0
#eg=0.0
#z=z-sky
#moment(z, xcntr, ycntr)
