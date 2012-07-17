import numpy as n

class gini:
	"""gini coefficient"""
	def __init__(self,z,ini_xcntr,ini_ycntr,pa,eg,extraction_radius,background,skysig):
		self.z			= z
		self.ini_xcntr		= ini_xcntr
		self.ini_ycntr		= ini_ycntr
		self.extraction_radius	= extraction_radius
		self.pa			= pa
		self.eg			= eg
		self.background		= background
		self.skysig		= skysig
		self.sigma		= 0.2*self.extraction_radius/1.5#the size of the boxcar
		self.segmentation	= segmentation(self.z,self.ini_xcntr,self.ini_ycntr,self.pa,self.eg,self.background,self.extraction_radius,self.sigma,self.skysig,1)
		self.gini_coef		= gini_coef(self.segmentation)


def segmentation(zextract,ini_xcntr,ini_ycntr,pa,eg,background,extraction_radius,sigma,skysig,lower):
	NXPTS=zextract.shape[0]
	NYPTS=zextract.shape[1]
	co= math.cos(pa*pi/180.0)
	si= math.sin(pa*pi/180.0)	
	x=reshape(arange(NXPTS*NYPTS),(NXPTS,NYPTS))%NYPTS
	x = x.astype(Float32)
	y=reshape(arange(NXPTS*NYPTS),(NXPTS,NYPTS))/NYPTS
	y = y.astype(Float32)
	tx = (x-ini_xcntr)*co + (y-ini_ycntr)*si
	ty = (ini_xcntr-x)*si + (y-ini_ycntr)*co
	R=sqrt(tx**2.0+ty**2.0/(1.0-eg)**2.0)
	I_sigma=image.gaussian_filter(zextract, int(sigma), order=0, mode='nearest')#the convolve image with the boxcar
	lower_thre=(I_sigma[where(R<=extraction_radius/1.5)].sum()-I_sigma[where(R<=extraction_radius/1.5-1)].sum())/(3.14*((extraction_radius*extraction_radius/1.5**2)-(extraction_radius/1.5-1)*(extraction_radius/1.5-1)))
	upper_thre= 10*skysig#the upper threshold
	#print lower_thre,upper_thre,extraction_radius
	#print zextract,ini_xcntr,ini_ycntr,pa,eg,background,extraction_radius,sigma
	#The below will finds the segmentation image
	I0=where(R>(2.5/1.5)*extraction_radius,0,zextract)
	I1=where(I_sigma<lower_thre,0,I0)
	I2=where(R>(extraction_radius/1.5-(0.02*extraction_radius)),0,I1)
	I3=where(R<=(extraction_radius/1.5-(0.02*extraction_radius)),0,I1)
	I4=where(I_sigma>upper_thre,0,I3)
	I=I2+I4
	return I

def gini_coef(I):
	oneD_I=reshape(I,(1,-1))
	sorted_I=sort(oneD_I)
	sorted_I=sorted_I[where(sorted_I>0)]
	n=sorted_I.nelements()
	average=sorted_I.mean()
	sum=0.0
	for i in range(n):
		sum+=(2*(i+1)-n-1)*sorted_I[i]
	G=sum/(average*n*(n-1))
	return G
