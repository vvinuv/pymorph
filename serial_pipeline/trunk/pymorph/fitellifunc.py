import sys,os
from pyraf import iraf
import config as c
import pyfits
import numpy as n

class FitElliFunc:
    """The class which will run ellipse task automatically. This will not
       keep any parameter fixed during the fit. The output of the ellipse
       fit will be E_+string(galid).txt for the galaxy and 
       OE_+string(galid).txt for the model galaxy"""
    def __init__(self, clus_id, line_s):
        self.clus_id   = clus_id
        self.line_s    = line_s
        self.fit_elli  = fit_elli(clus_id, line_s)
       

def fit_elli(clus_id, line_s):
    #try:
    size = c.size
    values = line_s.split()
    mask_file = 'ell_mask_' + str(c.rootname) + '_'  + str(clus_id) + '.fits'
    image_file = 'image_' + str(c.rootname) + '_'  + str(clus_id) + '.fits' 
    print image_file
    xcntr_o  = float(values[1]) #x center of the object
    ycntr_o  = float(values[2]) #y center of the object
    xcntr = (size/2.0) + 1.0 + xcntr_o - int(xcntr_o)
    ycntr = (size/2.0) + 1.0 + ycntr_o - int(ycntr_o)
    mag    = float(values[7]) #Magnitude
    radius = float(values[9]) #Half light radius
    mag_zero = 25.256 #magnitude zero point
    sky	 = float(values[10]) #sky
    if(float(values[11])>=0 and float(values[11])<=180.0): 
        pos_ang = float(values[11]) - 90.0 #position angle
    if(float(values[11])<0 and float(values[11])>=-180.0):
        pos_ang = 90.0 - abs(float(values[11]))  #position angle
    if(float(values[11])>180 and float(values[11])<=360.0):
        pos_ang = float(values[11]) - 360.0 + 90.0 #position angle
    if(float(values[11])>=-360 and float(values[11])<-180.0):
        pos_ang = float(values[11]) + 360.0 - 90.0 #position angle	
    axis_rat = 1.0 / float(values[12]) #axis ration b/a
    eg = 1 - axis_rat
    if(eg<=0.05):
        eg = 0.07
    major_axis = float(values[14])#major axis of the object
    iraf.imcopy(mask_file, 'image'+str(mask_file)[8:]+'.pl')
    run_elli(image_file, xcntr, ycntr, eg, pos_ang, major_axis)


def run_elli(input, output, xcntr, ycntr, eg, pa, sma, sky):#,radd,background):
    """The function responsible for fit ellipse"""
    EllGal = 'GalEllFit.fits'
    fEl = pyfits.open(input)
    GaLaXy = fEl[0].data
    fEl.close()
    GaLaXy = GaLaXy - sky
    hdu = pyfits.PrimaryHDU(GaLaXy.astype(n.float32))
    hdu.writeto(EllGal)
    try:
        iraf.stsdas(_doprint=0)
        iraf.tables(_doprint=0)
        iraf.stsdas.analysis(_doprint=0)
        iraf.stsdas.analysis.isophote(_doprint=0)
    except:
        iraf.stsdas()
	iraf.tables()
	iraf.stsdas.analysis()
	iraf.stsdas.analysis.isophote()
    image_exist = 1
    iraf.unlearn('geompar')
    iraf.geompar.x0=xcntr
    iraf.geompar.y0=ycntr
    iraf.geompar.ellip0=eg
    iraf.geompar.pa0=pa
    iraf.geompar.sma0=6.0
    iraf.geompar.minsma=0.1
    iraf.geompar.maxsma=sma*5.0
    iraf.geompar.step=0.1
    iraf.geompar.recente="no"
    iraf.geompar.xylearn="no"
    iraf.unlearn('controlpar')
    iraf.controlpar.conver=0.05
    iraf.controlpar.minit=10
    iraf.controlpar.maxit=50
    iraf.controlpar.hcenter="no"
    iraf.controlpar.hellip="no"
    iraf.controlpar.hpa="no"
    iraf.controlpar.wander="INDEF"
    iraf.controlpar.maxgerr=0.5
    iraf.controlpar.olthres=1
    iraf.controlpar.soft="no"
    iraf.samplepar.integrm="bi-linear"
    iraf.samplepar.usclip=3
    iraf.samplepar.lsclip=3
    iraf.samplepar.nclip=0
    iraf.samplepar.fflag=0.5
    iraf.unlearn('ellipse')
    iraf.ellipse("".join(EllGal), output="test", interac="no",Stdout="ellip", \
                 Stderr="err")
    iraf.tprint("test.tab", prparam="no", prdata="yes", pwidth=80, plength=0, \
                showrow="no", orig_row="no", showhdr="no", showunits="no", \
                columns="SMA, INTENS, INT_ERR, MAG, MAG_LERR, MAG_UERR, \
                TFLUX_E", rows="-", \
                option="plain", align="yes", sp_col="", lgroup=0, Stdout=output)
    for myfile in ['ellip','err','test.tab', EllGal]:
        if os.access(myfile,os.F_OK):
            os.remove(myfile)

def run_elli_full(input, output, xcntr, ycntr, eg, pa, sma, sky):#,radd,background):
    """The function responsible for fit ellipse"""
    EllGal = 'GalEllFit.fits'
    fEl = pyfits.open(input)
    GaLaXy = fEl[0].data
    fEl.close()
    GaLaXy = GaLaXy - sky
    hdu = pyfits.PrimaryHDU(GaLaXy.astype(n.float32))
    hdu.writeto(EllGal)
    try:
        iraf.stsdas(_doprint=0)
        iraf.tables(_doprint=0)
        iraf.stsdas.analysis(_doprint=0)
        iraf.stsdas.analysis.isophote(_doprint=0)
    except:
        iraf.stsdas()
	iraf.tables()
	iraf.stsdas.analysis()
	iraf.stsdas.analysis.isophote()
    image_exist = 1
    iraf.unlearn('geompar')
    iraf.geompar.x0=xcntr
    iraf.geompar.y0=ycntr
    iraf.geompar.ellip0=eg
    iraf.geompar.pa0=pa
    iraf.geompar.sma0=6.0
    iraf.geompar.minsma=0.1
    iraf.geompar.maxsma=sma*5.0
    iraf.geompar.step=0.1
    iraf.geompar.recente="no"
    iraf.geompar.xylearn="no"
    iraf.unlearn('controlpar')
    iraf.controlpar.conver=0.05
    iraf.controlpar.minit=10
    iraf.controlpar.maxit=50
    iraf.controlpar.hcenter="no"
    iraf.controlpar.hellip="no"
    iraf.controlpar.hpa="no"
    iraf.controlpar.wander="INDEF"
    iraf.controlpar.maxgerr=0.5
    iraf.controlpar.olthres=1
    iraf.controlpar.soft="no"
    iraf.samplepar.integrm="bi-linear"
    iraf.samplepar.usclip=3
    iraf.samplepar.lsclip=3
    iraf.samplepar.nclip=0
    iraf.samplepar.fflag=0.5
    iraf.unlearn('ellipse')
    iraf.ellipse("".join(EllGal), output="test", interac="no",Stdout="ellip", \
                 Stderr="err")
    iraf.tprint("test.tab", prparam="no", prdata="yes", pwidth=80, plength=0, \
                showrow="no", orig_row="no", showhdr="no", showunits="no", \
                columns="SMA, INTENS, INT_ERR, MAG, MAG_LERR, MAG_UERR, \
                TFLUX_E", rows="-", \
                option="plain", align="yes", sp_col="", lgroup=0, Stdout=output)
    for myfile in ['ellip','err','test.tab', EllGal]:
        if os.access(myfile,os.F_OK):
            os.remove(myfile)
