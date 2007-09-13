from numarray import *
import numarray.convolve as conv
import pyfits
from pyraf import iraf
class clumpness:
	"""clumpness parameter"""
	def __init__(self,z,ini_xcntr,ini_ycntr,pa,eg,extraction_radius,sigma,background,flag_image):
		self.z			= z
		self.ini_xcntr		= ini_xcntr
		self.ini_ycntr		= ini_ycntr
		self.extraction_radius	= extraction_radius
		self.pa			= pa
		self.eg			= eg
		self.background		= background
		self.flag_image		= flag_image
		self.sigma		= sigma#the size of the boxcar
		self.image_clumpness	= CLUMPNESS(self.z,self.ini_xcntr,self.ini_ycntr,self.pa,self.eg,self.background,self.extraction_radius,self.sigma,self.flag_image)

def CLUMPNESS(z,ini_xcntr,ini_ycntr,pa,eg,background,extraction_radius,sigma,flag_image):
	zextract=z[int(ini_ycntr-extraction_radius):int(ini_ycntr+extraction_radius),int(ini_xcntr-extraction_radius):int(ini_xcntr+extraction_radius)]
	
	NXPTS=zextract.shape[0]
	NYPTS=zextract.shape[1]
	#print NXPTS, NYPTS,ini_xcntr,ini_ycntr
	co= math.cos(pa*pi/180.0)
	si= math.sin(pa*pi/180.0)	
	x=reshape(arange(NXPTS*NYPTS),(NXPTS,NYPTS))/NYPTS
	x = x.astype(Float32)
	y=reshape(arange(NXPTS*NYPTS),(NXPTS,NYPTS))%NYPTS
	y = y.astype(Float32)
	tx = (x-NXPTS/2-1)*co + (y-NYPTS/2-1)*si
	ty = (NXPTS/2-1-x)*si + (y-NYPTS/2-1)*co
	R=sqrt(tx**2.0+ty**2.0/(1.0-eg)**2.0)
#	hdu = pyfits.PrimaryHDU(zextract.astype(Float32))
#	hdu.writeto('sImage.fits')
#	input = 'sImage.fits'
#	iraf.images(_doprint=0)
#	iraf.images.imfilter(_doprint=0)		
#	iraf.boxcar("".join(input), output = "SImage.fits", xwindow = sigma, ywindow = sigma, boundar = "nearest", constant = 0.0)
	I_sigma=conv.boxcar(zextract, (int(sigma),int(sigma)),mode='nearest')#the convolve image with the boxcar
#	f=pyfits.open("SImage.fits")
#	I_sigma = f[0].data
#	f.close()
#	for myfile in ['SImage.fits','sImage.fits']:
#		if os.access(myfile,os.F_OK):
#			os.remove(myfile)
	res=zextract-I_sigma #the residual image

	if(flag_image):
	#the below will find the image portion which is an anulus of inner radius .3*eta(.2) and outer radius 1.5*eta(.2)
#		res[where(R<=extraction_radius*(0.25/1.5))]=0 #making the residual value equal to zero inside the extraction_radius*(.25/1.5)
		res[where(R<=extraction_radius*(1/20.0))]=0
		res[where(R>=extraction_radius)]=0
		res_inside_anulus_sum=res[where(res>0)].sum() #the sum of residue inside the anulus
#                print res_inside_anulus_sum, 3.14*extraction_radius*extraction_radius
		z_inside_R_sum=zextract[where(R<=extraction_radius)].sum()/(3.14*extraction_radius*extraction_radius*sqrt(1-eg**2.0))
		area=3.14*(extraction_radius*extraction_radius*sqrt(1-eg**2.0))-3.14*(extraction_radius*extraction_radius*(1/6.0)*(1/6.0)*sqrt(1-eg**2.0))
		S=res_inside_anulus_sum/area#-(0.25*extraction_radius/1.5)*(0.25*extraction_radius/1.5)*sqrt(1-eg**2.0)))
		e1sq=zextract[where(res>0)].sum()+I_sigma[where(res>0)].sum()+4*zextract[where(res>0)].nelements()*background

	#no_res_inside_anulus=res[where(res>0)].nelements()
	else:
		res[where(R>=extraction_radius)]=0
		res_inside_anulus_sum=res[where(res>0)].sum()
		area=3.14*extraction_radius**2.0
		S=res_inside_anulus_sum/area
		z_inside_R_sum=0#just to return the value in the end
		e1sq=zextract[where(res>0)].sum()+I_sigma[where(res>0)].sum()+2*zextract[where(res>0)].nelements()*background
	e2sq=res_inside_anulus_sum**2.0
	e3sq=zextract[where(R<=extraction_radius)].sum()+2*zextract[where(R<=extraction_radius)].nelements()*background
 	e4sq=(zextract[where(R<=extraction_radius)].sum())**2.0
	if(e2sq!=0):
		error=e1sq/e2sq
	else:
		print "Could not find error"
		error=0.0

	return S,error,z_inside_R_sum,e3sq,e4sq	
