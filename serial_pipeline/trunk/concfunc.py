from numarray import *
import pyfits
class concentration:
	"""concentration parameter"""		
	def __init__(self, z, xcntr, ycntr, nxpts, nypts, pa, eg, background):
		self.z			= z
		self.xcntr		= xcntr
		self.ycntr		= ycntr
		self.nxpts		= nxpts
		self.nypts		= nypts
		self.eta_radius		= 2.0
		self.pa			= pa
		self.eg			= eg
		self.pixel_division	= 1
		self.incr		= 1.0
		self.background		= background
		
		#Find the radius to each pixel	
		self.r=return_r(self.nxpts,self.nypts,self.xcntr,self.ycntr,self.pa,self.eg,self.pixel_division)
		#Initiallizing the divided image
		#self.divide_r=divide_image(self.z,self.nxpts,self.nypts,self.pixel_division)
		self.divide_r=self.z
		Flag_check = 1
		
#	calculating the eta parameter
		self.total_rad = eta_radius_fnc(self.divide_r,self.r,self.nxpts,self.nypts,self.background,self.incr,self.eta_radius)
		if(self.total_rad == 9999):
			self.concen=9999
			self.error_con=9999
		else:
			self.total_I = self.divide_r[where(self.r<=self.total_rad)].sum()

		#calculating the 20% (r20) light radius and 80% (r80) light radius
			tempI80=tempI20=0.0
			tempr20=tempr80=0.0
			flag80=flag20=0 #this flags equal 1 when r80 and r20 find
			orad=self.total_rad*1.0 #orad is the initial radius to find the r80 and hopes the r80 is inside the eta_radius
			irad=self.eta_radius#irad is the initial radius to find the r20
			I20=0.2*self.total_I
			I80=0.8*self.total_I
			ttempI20=self.divide_r[where(self.r<=self.eta_radius-1.0)].sum()

			while flag80==0 or flag20==0:
				if(flag80==0):
					tempI80 = self.divide_r[where(self.r<=orad)].sum()
				if(flag20==0):
					tempI20 = self.divide_r[where(self.r<=irad)].sum()
				if(flag80==0 and tempI80<I80):
					flag80=1
					alpha80=error(I80,tempI80,ttempI80,orad,self.total_I,self.total_rad,self.background)
					self.r80=orad+alpha80[0]
				if(flag20==0 and tempI20>I20):
					flag20=1
					alpha20=error(I20,ttempI20,tempI20,irad-self.incr,self.total_I,self.total_rad,self.background)
					self.r20=irad-self.incr+alpha20[0]
				orad-=self.incr
				irad+=self.incr
				ttempI80=tempI80
				ttempI20=tempI20

			self.r20_error=alpha20[1]
			self.r80_error=alpha80[1]
			#calculate concentration parameter "concen"
			self.concen=5*log10(self.r80/self.r20)
			self.error_ratio=sqrt((self.r80/self.r20)**2*((alpha80[1]/self.r80)**2+(alpha20[1]/self.r20)**2))
			self.error_con=5*self.error_ratio/(self.r80/self.r20)
		#econcen = 5*log10 (r80/r20) - 5*log10 ((r80-0.5)/(r20+0.5)) This is the way Conselice finds error in concentration


#------This function returns r------#

def return_r(nxpts,nypts,xcntr,ycntr,pa,eg,pixel_division):
	divide_x=(reshape(arange(nxpts*pixel_division*pixel_division*nypts),(nxpts*pixel_division,nypts*pixel_division))%(nypts*pixel_division))*(1/(pixel_division*1.0))

	divide_y=(reshape(arange(nxpts*pixel_division*pixel_division*nypts),(nxpts*pixel_division,nypts*pixel_division))/(nypts*pixel_division))*(1/(pixel_division*1.0))
	
	one_minus_eg_sq		= (1.0-eg)**2.0
	co			= math.cos(pa*pi/180.0)
	si			= math.sin(pa*pi/180.0)
	
	#radius corresponds to the pixel index
	tx = (divide_x-xcntr)*co + (divide_y-ycntr)*si
	ty = (xcntr-divide_x)*si + (divide_y-ycntr)*co
	tx = tx.astype(Float32)
	tx = tx.astype(Float32)
	r=sqrt(tx**2.0+ty**2.0/one_minus_eg_sq)
	return r



#-----Function which make divided image-----#
def divide_image(z,nxpts,nypts,pixel_division):
	divide_r=zeros([nxpts*pixel_division,nypts*pixel_division],type=Float32)
	I_at_diff_radii=zeros([nypts*nxpts],type=Float32)
	radii=zeros([nypts*nxpts],type=Float32)
		#filling values in the divided images
	a=0
	for i in range(nypts):
		for m in range(pixel_division):
			b=0
			for j in range(nxpts):
				for n in range(pixel_division):
					divide_r[a,b]=z[i,j]/(pixel_division*pixel_division*1.0)
					b+=1
			a+=1
	return divide_r



#------The function which returns the total radius ie. 1.5*r(eta=0.2)-----#
#I_ave_radius is the intensity at different anulus, I_ave - intensity inside different radius, no_I_ave_radius - the no of pixels in the anular region, no_I_ave - the no. of pixels inside the different radius

def eta_radius_fnc(divide_r,r,nxpts,nypts,background,incr,eta_radius):
	FLAG_ETA=n=0 # The n is using because the program won't calculate eta radius if the eta<0.2 in the first run itself. FLAG_ETA will tell whether the condition eta<0.2 is reached.
	eta=I_ave_radius=I_ave=no_I_ave=no_I_ave_radius=r50=temp_eta=0.0
#Here it try to find the eta radius. It will run until FLAG_ETA=0 and eta_radius< min(NXPTS,NYPTS)
		
	while FLAG_ETA==0:
		temp_eta = eta
		I_ave_radius = 	divide_r[where(r<=eta_radius+1.0)].sum()-divide_r[where(r<eta_radius-1.0)].sum()
		no_I_ave_radius = divide_r[where(r<=eta_radius+1.0)].nelements()-divide_r[where(r<=eta_radius-1.0)].nelements()
		I_ave = divide_r[where(r<=eta_radius)].sum()
		no_I_ave = divide_r[where(r<=eta_radius)].nelements()			
		if(I_ave==0 or no_I_ave==0 or I_ave_radius==0 or no_I_ave_radius==0):
			FLAG_ETA=0
		else :
			#I_ave=I_ave/no_I_ave #average intensity inside a radius
			#I_ave_radius=I_ave_radius/no_I_ave_radius #average intensity inside the annulus
			#print I_ave,I_ave_radius,eta_radius,eta
			I_ave=I_ave/(3.14*eta_radius*eta_radius)
			I_ave_radius=I_ave_radius/(3.14*((eta_radius+1.0)*(eta_radius+1.0)-(eta_radius-1.0)*(eta_radius-1.0)))
			eta=(I_ave_radius)/(I_ave) #eta parameter
			if(eta<0.2 and n==0):
				FLAG_ETA=1
				print "Finding eta radius failed. Exiting program"
				eta_radius_corre = 9999
#				os._exit(0)
			elif(eta<0.2 and n>0):
				#check = check_eta(eta_radius,eta,divide_r,r)
				#eta = check[0]
				#eta_radius = check[1]
				#Flag_check = check[2]
			#if(Flag_check==0):
				FLAG_ETA=1
				ERROR1=sqrt(eta**2*(((I_ave_radius*no_I_ave_radius+background*no_I_ave_radius)/(I_ave_radius*no_I_ave_radius)**2)+((I_ave*no_I_ave+background*no_I_ave)/(I_ave*no_I_ave)**2)))		
				ERROR2=sqrt(temp_eta**2*(((temp_I_ave_radius*temp_no_I_ave_radius+background*temp_no_I_ave_radius)/(temp_I_ave_radius*temp_no_I_ave_radius)**2)+((temp_I_ave*temp_no_I_ave+background*temp_no_I_ave)/(temp_I_ave*temp_no_I_ave)**2)))
				alpha=(0.2-temp_eta)/(eta-temp_eta)
				eta_radius_corre=eta_radius-1+alpha
				error_eta_radius=sqrt(alpha*alpha*(((ERROR2)**2/(0.2-temp_eta)**2)+((ERROR1+ERROR2)**2/(eta-temp_eta)**2)))
				TEMP_RAD=eta_radius_corre
		#print eta,eta_radius,I_ave,I_ave_radius
		temp_I_ave_radius=I_ave_radius
		temp_I_ave=I_ave
		temp_no_I_ave=no_I_ave
		temp_no_I_ave_radius=no_I_ave_radius
		eta_radius += incr
		if(eta_radius>min(nxpts,nypts)/2.0):
			print "The radius measurment exceeded the size of the frame. Exiting program"
			eta_radius_corre = 9999
#			os._exit(0)
		I_ave_radius=I_ave=no_I_ave=no_I_ave_radius=0.0
		n+=1
	#Calculating the total light which the light inside 1.5*r(eta=.2). conselice 2003
	if(eta_radius_corre == 9999):
		total_rad = eta_radius_corre
	else:
		total_rad = 1.5*(eta_radius_corre)	
	return total_rad



#-----Function for error estimation for concentration-----#
def error(CONST,x0,x1,rad,total_I,total_rad,background):
	"""Function to find error"""
	alpha=(CONST-x0)/(x1-x0)
	error_alpha=sqrt(alpha*alpha*(((x0+total_I+(total_rad*background)+(rad*background))/(CONST-x0)**2)+((x1+x0+(rad*background)+((rad+1)*background))/(x1-x0)**2)))
	return alpha,error_alpha



#----Function which check whether the eta is .2 locally or it is global
def check_eta(eta_radius,eta,divide_r,r):
	ETA=eta
	ETA_RADIUS=eta_radius
	eta_radius=eta_radius+5
	Flag_check=0
	while Flag_check == 0 and eta_radius>ETA_RADIUS:
		temp_eta = eta
		I_ave_radius = 	divide_r[where(r<=eta_radius+0.5)].sum()-divide_r[where(r<=eta_radius-0.5)].sum()
		no_I_ave_radius = divide_r[where(r<=eta_radius+0.1)].nelements()-divide_r[where(r<=eta_radius-0.1)].nelements()
		I_ave = divide_r[where(r<=eta_radius)].sum()
		no_I_ave = divide_r[where(r<=eta_radius)].nelements()			
		if(I_ave==0 or no_I_ave==0 or I_ave_radius==0 or no_I_ave_radius==0):
			eta_radius-=1
		else:
			I_ave=I_ave/no_I_ave #average intensity inside a radius
			I_ave_radius=I_ave_radius/no_I_ave_radius #average intensity inside the annulus
			eta=(I_ave_radius)/(I_ave)
			eta_radius-=1
		if(abs(eta-0.2)<abs(ETA-0.2)):
			Flag_check=1
			ETA=eta
			ETA_RADIUS=eta_radius+1
	return ETA, ETA_RADIUS, Flag_check

