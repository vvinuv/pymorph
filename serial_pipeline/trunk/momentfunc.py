from numarray import *
import numarray.ma as ma
import random, numarray.random_array,numarray.mlab
import numarray.nd_image as image
import numarray.convolve as conv
class moment:
	"""moment calculation"""
	def __init__(self,z,ini_xcntr,ini_ycntr):
		self.z			= z
		self.ini_xcntr		= ini_xcntr
		self.ini_ycntr		= ini_ycntr
		self.moment_of_light	= moment_of_light(self.z,self.ini_xcntr,self.ini_ycntr)


def moment_of_light(zextract,ini_xcntr,ini_ycntr):
	flag_center=0
	center_rad=1.0
	n=0
	NXPTS=zextract.shape[0]
	NYPTS=zextract.shape[1]
	x=reshape(arange(NXPTS*NYPTS),(NXPTS,NYPTS))/NYPTS
	x = x.astype(Float32)
	y=reshape(arange(NXPTS*NYPTS),(NXPTS,NYPTS))%NYPTS
	y = y.astype(Float32)
	M_tot=zeros([9],type=Float32)
	while flag_center==0:
		xcntr=reshape(array([[ini_xcntr]*3,[ini_xcntr-center_rad]*3,[ini_xcntr+center_rad]*3],type=Float32),(9,)) 
		ycntr=array([ini_ycntr,ini_ycntr-center_rad,ini_ycntr+center_rad]*3,type=Float32)
		for iter in range(9):
			tx = x-xcntr[iter]
			ty = y-ycntr[iter]
			R=tx**2.0+ty**2.0
			RZ=R*zextract
			M_tot[iter]=RZ.sum()
		M_tot_min=M_tot.min()
		index=M_tot.argmin()
		if(n>100):
			xcntr[index] = ini_xcntr
			ycntr[index] = ini_ycntr
		if(xcntr[index] == ini_xcntr and ycntr[index] == ini_ycntr):
			flag_center=1
		else:
			ini_xcntr=xcntr[index]
			ini_ycntr=ycntr[index]
			n+=1
			#print n
	oneD_I=zextract.flat
	sorted_I=sort(oneD_I)
	argument=argsort(oneD_I)
	tx = x-xcntr[index]
	ty = y-ycntr[index]
	R=tx**2.0+ty**2.0
	RZ=R*zextract
	oneD_RZ=RZ.flat
	oneD_RZ=oneD_RZ[argument]
	N=sorted_I.nelements()
	total=sorted_I.sum()
	sum=0.0
	sumM=0.0
	i=N-1
	while sum<0.2*total:
		sum+=sorted_I[i]
		sumM+=oneD_RZ[i]
		i-=1
	
	M20=log10(sumM/M_tot_min)
	return M20,M_tot[index],sumM,total,sum,n,xcntr[index],ycntr[index]
