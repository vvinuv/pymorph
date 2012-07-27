import numpy as n
import pyfits
from momentfunc import *
#import scipy.signal
import convolve as conv

class gini:
	"""Calculate gini coefficient at different radii. This will also call 
           the function to compute M20. The algorithm is as follows

           1. Find the pixels in the image which belong to the galaxy, ie. make
              a segmentation map. This can be done by smoothing the image by a 
               boxcar of size r(\eta)/5 

           2. The surface brightness at r(\eta), I_(\eta) is measured and 
              pixels in the smoothed image with flux values greater than 
              I_(\eta) and less than 10(\sigma) is assigned to the galaxy. 
              \sigma is the sky deviation and which removes any remaining 
              cosmic rays or spurious noise pixels in the image and are not 
              included in the segmentation map.

           3. The Gini coefficient can be computed by the equation
             G = (1 / Avg(X) * n * (n-1)) * Sum over pixel[(2 * i - n -1 ) * X] 
           """
	def __init__(self,z,ini_xcntr,ini_ycntr,pa,eg,r20,r50,r80,extraction_radius,background,skysig):
		self.z			= z
		self.ini_xcntr		= ini_xcntr
		self.ini_ycntr		= ini_ycntr
		self.extraction_radius	= extraction_radius
		self.r20		= r20
                self.r50                = r50
		self.r80		= r80
		self.pa			= pa
		self.eg			= eg
		self.background		= background
		self.skysig		= skysig
		self.sigma		= 0.2*self.extraction_radius/1.5#the size of the boxcar
		self.segmentation	= segmentation(self.z,self.ini_xcntr,self.ini_ycntr,self.pa,self.eg,self.background,self.r20,self.r50,self.r80,self.extraction_radius,self.sigma,self.skysig,1)
#		self.gini_coef		= gini_coef(self.segmentation)


def gauss_kern(size, sizey=None):
	size = int(size)
	if not sizey:
		sizey = size
	else:
		sizey = int(sizey)
	x, y = n.mgrid[-size:size+1, -sizey:sizey+1]
	g = n.exp(-(x**2/n.float(size)+y**2/n.float(sizey)))
	return g / g.sum()


def segmentation(zextract,ini_xcntr,ini_ycntr,pa,eg,background,r20,r50,r80,extraction_radius,sigma,skysig,lower):
	"""This function is resposible for find the segmentation map"""
	def gini_coef(I):
		"""This will calculate Gini Coefficient"""
		oneD_I=n.reshape(I,(1,-1))
		sorted_I=n.sort(oneD_I)
		sorted_I=sorted_I[n.where(sorted_I>0)]
		nn=sorted_I.size
		average=n.abs(sorted_I).mean()
		sum=0.0
		for i in range(nn):
			sum+=(2*(i+1)-nn-1)*abs(sorted_I[i])
		G=sum/(average*nn*(nn-1))
		return G
	NXPTS=zextract.shape[0]
	NYPTS=zextract.shape[1]
	co= n.cos(pa*n.pi/180.0)
	si= n.sin(pa*n.pi/180.0)	
	x=n.reshape(n.arange(NXPTS*NYPTS),(NXPTS,NYPTS))%NYPTS
	x = x.astype(n.float32)
	y=n.reshape(n.arange(NXPTS*NYPTS),(NXPTS,NYPTS))/NYPTS
	y = y.astype(n.float32)
	tx = (x-ini_xcntr)*co + (y-ini_ycntr)*si
	ty = (ini_xcntr-x)*si + (y-ini_ycntr)*co
	R=n.sqrt(tx**2.0+ty**2.0/(1.0-eg)**2.0)

#	I_sigma=image.gaussian_filter(zextract, int(sigma), order=0, mode='nearest')#the convolve image with the boxcar
#	g = gauss_kern(sigma, sizey=None)
#	I_sigma=scipy.signal.convolve(zextract, g, mode='same')	
	I_sigma = conv.boxcar(zextract,(int(sigma),int(sigma)),mode='nearest') 
	lower_thre=(I_sigma[n.where(R<=extraction_radius/1.5)].sum()-I_sigma[n.where(R<=extraction_radius/1.5-1)].sum())/(3.14*((extraction_radius*extraction_radius/1.5**2)-(extraction_radius/1.5-1)*(extraction_radius/1.5-1)))
	upper_thre= 10*skysig#the upper threshold
	#print lower_thre,upper_thre,extraction_radius
	#print zextract,ini_xcntr,ini_ycntr,pa,eg,background,extraction_radius,sigma
	#The below will finds the segmentation image
	I0=n.where(R>(2.5/1.5)*extraction_radius,0,zextract)
	I1=n.where(I_sigma<lower_thre,0,I0)
	I2=n.where(R>(extraction_radius/1.5-(0.02*extraction_radius)),0,I1)
	I3=n.where(R<=(extraction_radius/1.5-(0.02*extraction_radius)),0,I1)
	I4=n.where(I_sigma>upper_thre,0,I3)
	I=I2+I4
	G = gini_coef(I)
	mo=moment(I, ini_xcntr, ini_ycntr)
	M = mo.moment_of_light[0]
	f=pyfits.open('AResidual.fits')
	res =  f[0].data
	f.close()
	res1=n.where(I == 0, 0, res)
	G_res = gini_coef(res1)
	mo1=moment(res, ini_xcntr, ini_ycntr)
	M_res = mo1.moment_of_light[0]

        I80=n.where(R>r80,0,zextract)
        G80 = gini_coef(I80)
	mo2=moment(I80, ini_xcntr, ini_ycntr)
	M80 = mo2.moment_of_light[0]

        I50=n.where(R>r50,0,zextract)
        G50 = gini_coef(I50)
        mo3=moment(I50, ini_xcntr, ini_ycntr)
        M_50 = mo3.moment_of_light[0]

	I20=n.where(R>r20,0,zextract)
        G20 = gini_coef(I20)
	mo4=moment(I20, ini_xcntr, ini_ycntr)
	M_20 = mo4.moment_of_light[0]
#	print G, G_res, G80, G20, M, M_res, M80, M_20
	return G, G_res, G80, G50, G20, M, M_res, M80, M_50, M_20

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
#background=1390.377
#nxpts=z.shape[0]
#nypts=z.shape[1]
#extraction_radius=100
#sigma=20.0
#gini(z,xcntr,ycntr,pa,eg,20,80,extraction_radius,background,3)
