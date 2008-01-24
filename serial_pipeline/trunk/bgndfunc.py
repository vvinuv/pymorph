import numpy as n
import ndimage as im
class bgnd:
	"""background parameters"""
	def __init__(self,z,back_ini_xcntr,back_ini_ycntr,back_extraction_radius):
		self.z=z
		self.back_ini_xcntr=back_ini_xcntr
		self.back_ini_ycntr=back_ini_ycntr
		self.back_extraction_radius=back_extraction_radius
		self.zextract=self.z[int(self.back_ini_ycntr-self.back_extraction_radius):int(self.back_ini_ycntr+self.back_extraction_radius),int(self.back_ini_xcntr-self.back_extraction_radius):int(self.back_ini_xcntr+self.back_extraction_radius)]
		self.mean=im.mean(self.zextract)
		self.standard_deviation=im.standard_deviation(self.zextract)
