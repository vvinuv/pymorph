import os
import sys
import pyfits
import config as c
from os.path import exists
from numpy import log10
from readlog import ReadLog
from runsexfunc import *
import copy

class ConfigIter:
    """The class making configuration file for GALFIT. The configuration file 
       consists of bulge and disk component of the object and only Sersic 
       component for the neighbours, if any. The sky is always fixed and has
       the value of SExtractor. The disk/boxy parameter is also fixed to zero.
       The initial value for Sersic index 'n' is 4.The configuration file has 
       the name G_string(galid).in. The output image has the name 
       O_string(galid).fits"""
    def __init__(self, cutimage, whtimage, xcntr, ycntr, NXPTS, NYPTS, line_s, psffile):
        self.cutimage = cutimage
        self.line_s  = line_s
	self.whtimage = whtimage
        self.xcntr = xcntr
        self.ycntr = ycntr
        self.NXPTS = NXPTS
        self.NYPTS = NYPTS 
        self.psffile = psffile
        self.confiter    = confiter(cutimage, whtimage, xcntr, ycntr, NXPTS, NYPTS, line_s, psffile)
		

def confiter(cutimage, whtimage, xcntr, ycntr, NXPTS, NYPTS, line_s, psffile):
    RunSex(cutimage, whtimage, 'TEMP.SEX.cat', 9999, 9999, 0)
    imagefile = c.imagefile
    sex_cata = c.sex_cata
    threshold = c.threshold
    thresh_area = c.thresh_area
    mask_reg = c.mask_reg
    try:
        ComP = c.components 
    except:
        ComP = ['bulge', 'disk']
    if len(ComP) == 0:
        ComP = ['bulge', 'disk']
    values = line_s.split()
    outfile   = 'O_' + str(cutimage)[:-5] + '.fits'
    mask_file = 'M_' + str(cutimage)[:-5] + '.fits'
    config_file = 'G_' + str(cutimage)[:-5] + '.in' #Name of the GALFIT configuration file
    constrain_file = str(cutimage)[:-5] + '.con'
    try:
	c.center_constrain = c.center_constrain
    except:
	c.center_constrain = 2.0
    def SersicMainConstrain(constrain_file, cO, cen_con, re_con):
        f_constrain = open(constrain_file, 'ab')
        f_constrain.write(str(cO) + '      n      ' + str(c.LN) + \
                          ' to ' + str(c.UN) +  '\n')
        f_constrain.write(str(cO) + '      x      ' + \
                          str(-cen_con) + '     ' + \
                          str(cen_con) + '\n')
        f_constrain.write(str(cO) + '      y      ' + \
                          str(-cen_con) + '     ' + \
                          str(cen_con) + '\n')
        f_constrain.write(str(cO) + '     mag     ' + str(c.UMag) + \
                          ' to ' + str(c.LMag) + '\n')
        if re_con == 0:
            f_constrain.write(str(cO) + '      re     ' + str(0.1) +\
                          ' to ' + str(c.URe) + '\n')
        else:
            f_constrain.write(str(cO) + '      re     ' + str(0.1) +\
                          ' to ' + str(re_con) + '\n')

        f_constrain.write(str(cO) + '      q       0.05 to 1.0\n')
        f_constrain.write(str(cO) + '      pa       -360.0 to 360.0\n')
        f_constrain.close()

    def ExpdiskConstrain(constrain_file, cO):
        f_constrain = open(constrain_file, 'ab')
        f_constrain.write(str(cO) + '       x       ' + \
                          str(-c.center_constrain) + '     ' + \
                          str(c.center_constrain) + '\n')
        f_constrain.write(str(cO) + '       y       ' + \
                          str(-c.center_constrain) + '     ' + \
                          str(c.center_constrain) + '\n')
        f_constrain.write(str(cO) + '     mag     ' + str(c.UMag) + \
                          ' to ' + str(c.LMag) + '\n')
        f_constrain.write(str(cO) + '      rs     ' + str(c.LRd) + \
                          ' to ' + str(c.URd) + '\n')
        f_constrain.write(str(cO) + '      q       0.0 to 1.0\n')
        f_constrain.write(str(cO) + '      pa       -360.0 to 360.0\n')
        f_constrain.close()

    def BarConstrain(constrain_file, cO):
        f_constrain = open(constrain_file, 'ab')
        f_constrain.write(str(cO) + '      n      ' + str(0.01) + \
                          ' to ' + str(2.2) +  '\n')
        f_constrain.write(str(cO) + '      x      ' + \
                          str(-(c.center_constrain - 1.0)) + '     ' + \
                          str(c.center_constrain - 1.0) + '\n')
        f_constrain.write(str(cO) + '      y      ' + \
                          str(-(c.center_constrain - 1.0)) + '     ' + \
                          str(c.center_constrain + 1.0) + '\n')
        f_constrain.write(str(cO) + '     mag     ' + str(c.UMag) + \
                          ' to ' + str(c.LMag) + '\n')
        f_constrain.write(str(cO) + '      re     ' + str(c.LRe) +\
                          ' to ' + str(c.URe) + '\n')
        f_constrain.write(str(cO) + '      q       0.0 to 0.6\n')
        f_constrain.write(str(cO) + '      pa       -360.0 to 360.0\n')
        f_constrain.close()

    def PsfConstrain(constrain_file, cO):
        f_constrain = open(constrain_file, 'ab')
        f_constrain.write(str(cO) + '       x       ' + \
                          str(-c.center_constrain) + '     ' + \
                          str(c.center_constrain) + '\n')
        f_constrain.write(str(cO) + '       y       ' + \
                          str(-c.center_constrain) + '     ' + \
                          str(c.center_constrain) + '\n')
        f_constrain.write(str(cO) + '     mag     ' + str(c.UMag) + \
                          ' to ' + str(c.LMag) + '\n')
        f_constrain.close()

    def SersicConstrain(constrain_file, cO):
        f_constrain = open(constrain_file, 'ab')
        f_constrain.write(str(cO) + '      n      0.02 to 20.0  \n')
        f_constrain.write(str(cO) + '     mag    -100.0 to 100.0\n')
        f_constrain.write(str(cO) + '      re      0.0 to 500.0\n')
        f_constrain.write(str(cO) + '      q       0.0 to 1.0\n')
        f_constrain.write(str(cO) + '      pa    -360.0 to 360.0\n')
        f_constrain.close()
    
    def SkyConstrain(constrain_file, cO):
        f_constrain = open(constrain_file, 'ab')
        f_constrain.write(str(cO) + '      sky      ' + str(c.SexSky - c.SexSky * 0.001) + ' to ' + str(c.SexSky + c.SexSky * 0.001) + '  \n')
        f_constrain.close()
    xcntr_o  = xcntr #float(values[1]) #x center of the object
    ycntr_o  = ycntr #float(values[2]) #y center of the object
    mag    = float(values[7]) #Magnitude
    radius = float(values[9]) #Half light radius
    mag_zero = c.mag_zero #magnitude zero point
    sky	 = float(values[10]) #sky 
    pos_ang = float(values[11]) - 90.0 #position angle
    axis_rat = 1.0/float(values[12]) #axis ration b/a
    area_o = float(values[13])   # object's area
    major_axis = float(values[14])	#major axis of the object
    ParamDict = {}
    ParamDict[0] = {}
    #Add components
    AdComp = 1
    if 'bulge' in ComP:
        c.Flag += 512
        ParamDict[0][AdComp] = {}
        #Bulge Parameters
        ParamDict[0][AdComp][1] = 'sersic'
        ParamDict[0][AdComp][2] = [xcntr_o, ycntr_o]
        ParamDict[0][AdComp][3] = mag 
        ParamDict[0][AdComp][4] = radius 
        ParamDict[0][AdComp][5] = 4.0
        ParamDict[0][AdComp][6] = axis_rat
        ParamDict[0][AdComp][7] = pos_ang
        ParamDict[0][AdComp][8] = 0
        ParamDict[0][AdComp][9] = 0
        ParamDict[0][AdComp][11] = 'Main'
        AdComp += 1
    if 'disk' in ComP:
        c.Flag += 1024
        #Disk parameters
        ParamDict[0][AdComp] = {}
        ParamDict[0][AdComp][1] = 'expdisk'
        ParamDict[0][AdComp][2] = [xcntr_o, ycntr_o]
        ParamDict[0][AdComp][3] = mag 
        ParamDict[0][AdComp][4] = radius
        ParamDict[0][AdComp][5] = axis_rat
        ParamDict[0][AdComp][6] = pos_ang
        ParamDict[0][AdComp][7] = 0
        ParamDict[0][AdComp][8] = 0
        ParamDict[0][AdComp][11] = 'Main'
        AdComp += 1
    if 'bar' in ComP:
#        c.Flag += 512
        ParamDict[0][AdComp] = {}
        #Bulge Parameters
        ParamDict[0][AdComp][1] = 'bar'
        ParamDict[0][AdComp][2] = [xcntr_o, ycntr_o]
        ParamDict[0][AdComp][3] = mag + 2.5 * log10(2.0)
        ParamDict[0][AdComp][4] = radius
        ParamDict[0][AdComp][5] = 0.5
        ParamDict[0][AdComp][6] = 0.3
        ParamDict[0][AdComp][7] = pos_ang
        ParamDict[0][AdComp][8] = 0
        ParamDict[0][AdComp][9] = 0
        ParamDict[0][AdComp][11] = 'Main'
        AdComp += 1
    if 'point' in ComP:
#        c.Flag += 512
        ParamDict[0][AdComp] = {}
        #Point Parameters
        ParamDict[0][AdComp][1] = 'psf'
        ParamDict[0][AdComp][2] = [xcntr_o, ycntr_o]
        ParamDict[0][AdComp][3] = mag + 2.5 * log10(6.0)
        ParamDict[0][AdComp][4] = 0
        ParamDict[0][AdComp][11] = 'Main'
        AdComp += 1

    isneighbour = 0
    f_constrain = open(constrain_file, 'ab')
    for line_j in open('TEMP.SEX.cat','r'):
        try:
            values = line_j.split()
            xcntr_n  = float(values[1]) #x center of the neighbour
            ycntr_n  = float(values[2]) #y center of the neighbour
            mag    = float(values[7]) #Magnitude
            radius = float(values[9]) #Half light radius
            sky      = float(values[10]) #sky
            pos_ang = float(values[11]) - 90.0 #position angle
            axis_rat = 1.0/float(values[12]) #axis ration b/a
            area_n = float(values[13]) # neighbour area
            maj_axis = float(values[14])#major axis of neighbour
            NotFitNeigh = 0
            if abs(xcntr_n - xcntr_o) > NXPTS / 2.0 + c.avoideme or \
               abs(ycntr_n - ycntr_o) > NYPTS / 2.0 + c.avoideme or \
               abs(xcntr_n - xcntr_o) < 5.0 and abs(ycntr_n - ycntr_o) < 5.0:
                NotFitNeigh = 1
            if(abs(xcntr_n - xcntr_o) <= (major_axis + maj_axis) * \
               threshold and \
               abs(ycntr_n - ycntr_o) <= (major_axis  + maj_axis) * \
               threshold and area_n >= thresh_area * area_o and \
               xcntr_n != xcntr_o and ycntr_n != ycntr_o and NotFitNeigh == 0):
                if((xcntr_o - xcntr_n) < 0):
                    xn = xcntr + abs(xcntr_n - xcntr_o)
                if((ycntr_o - ycntr_n) < 0):
                    yn = ycntr + abs(ycntr_n - ycntr_o)
                if((xcntr_o - xcntr_n) > 0):
                    xn = xcntr - (xcntr_o - xcntr_n)
                if((ycntr_o - ycntr_n) > 0):
                    yn = ycntr - (ycntr_o - ycntr_n)
                ParamDict[0][AdComp] = {}
                ParamDict[0][AdComp][1] = 'sersic'
                ParamDict[0][AdComp][2] = [xn, yn]
                ParamDict[0][AdComp][3] = mag
                ParamDict[0][AdComp][4] = radius
                ParamDict[0][AdComp][5] = 4.0
                ParamDict[0][AdComp][6] = axis_rat
                ParamDict[0][AdComp][7] = pos_ang
                ParamDict[0][AdComp][8] = 0
                ParamDict[0][AdComp][9] = 0
                ParamDict[0][AdComp][11] = 'Other'

                isneighbour = 1
                AdComp += 1
        except:
            pass
    f_constrain.close()
    if isneighbour:
        c.Flag  += 4096
    #Sky component
    ParamDict[0][AdComp] = {}
    ParamDict[0][AdComp][1] = 'sky'
    ParamDict[0][AdComp][2] = c.SexSky
    ParamDict[0][AdComp][3] = 0
    ParamDict[0][AdComp][4] = 0
    ParamDict[0][AdComp][5] = 0
    ParamDict[0][AdComp][11] = 'Other'

    #Write Sersic function
    def SersicFunc(conffile, ParamDict, FitDict, No, RunNo):
        f=open(config_file, 'ab')
        f.write('# Sersic function\n\n')
        f.writelines([' 0) sersic \n'])
        f.writelines([' 1) ', str(ParamDict[RunNo][No][2][0]), ' ', \
                              str(ParamDict[RunNo][No][2][1]), ' ', \
                              str(FitDict[No][1][0]),   ' ', \
                              str(FitDict[No][1][1]),   '\n'])
        f.writelines([' 3) ', str(ParamDict[RunNo][No][3]), ' ', \
                              str(FitDict[No][2]),  '\n'])
        f.writelines([' 4) ', str(ParamDict[RunNo][No][4]), ' ', \
                              str(FitDict[No][3]),  '\n'])
        f.writelines([' 5) ', str(ParamDict[RunNo][No][5]), ' ',\
                              str(FitDict[No][4]),  '\n'])
        f.writelines([' 8) ', str(ParamDict[RunNo][No][6]), ' ', \
                              str(FitDict[No][5]),  '\n'])
        f.writelines([' 9) ', str(ParamDict[RunNo][No][7]), ' ', \
                              str(FitDict[No][6]),  '\n'])
        if c.bdbox or c.bbox:
            f.writelines(['10) 0.0 1		\n'])
        else:
            f.writelines(['10) 0.0 0            \n'])
        f.writelines([' Z) 0 			\n\n\n'])
        f.close()
    def ExpFunc(conffile, ParamDict, FitDict, No, RunNo):
        f=open(config_file, 'ab')
        f.writelines(['# Exponential function\n\n'])
        f.writelines([' 0) expdisk \n'])
        f.writelines([' 1) ', str(ParamDict[RunNo][No][2][0]), ' ', \
                              str(ParamDict[RunNo][No][2][1]),' ', \
                              str(FitDict[No][1][0]), ' ', \
                              str(FitDict[No][1][1]), '\n'])
        f.writelines([' 3) ', str(ParamDict[RunNo][No][3]),  ' ', \
                              str(FitDict[No][2]),    '\n'])
        f.writelines([' 4) ', str(ParamDict[RunNo][No][4]),  ' ', \
                              str(FitDict[No][3]),    '\n'])
        f.writelines([' 8) ', str(ParamDict[RunNo][No][5]),  ' ', \
                              str(FitDict[No][4]),    '\n'])
        f.writelines([' 9) ', str(ParamDict[RunNo][No][6]),  ' ', \
                              str(FitDict[No][5]),    '\n'])
        if c.bdbox or c.dbox:
            f.writelines(['10) 0.0 1   \n']) 
        else:
            f.writelines(['10) 0.0 0   \n'])
        f.writelines([' Z) 0           \n\n\n'])
        f.close()
    def BarFunc(conffile, ParamDict, FitDict, No, RunNo):
        f=open(config_file, 'ab')
        f.write('# Sersic function\n\n')
        f.writelines([' 0) sersic \n'])
        f.writelines([' 1) ', str(ParamDict[RunNo][No][2][0]), ' ', \
                              str(ParamDict[RunNo][No][2][1]), ' ', \
                              str(FitDict[No][1][0]),   ' ', \
                              str(FitDict[No][1][1]),   '\n'])
        f.writelines([' 3) ', str(ParamDict[RunNo][No][3]), ' ', \
                              str(FitDict[No][2]),  '\n'])
        f.writelines([' 4) ', str(ParamDict[RunNo][No][4]), ' ', \
                              str(FitDict[No][3]),  '\n'])
        f.writelines([' 5) ', str(ParamDict[RunNo][No][5]), ' ',\
                              str(FitDict[No][4]),  '\n'])
        f.writelines([' 8) ', str(ParamDict[RunNo][No][6]), ' ', \
                              str(FitDict[No][5]),  '\n'])
        f.writelines([' 9) ', str(ParamDict[RunNo][No][7]), ' ', \
                              str(FitDict[No][6]),  '\n'])
        if c.bdbox or c.bbox:
            f.writelines(['10) 0.0 1		\n'])
        else:
            f.writelines(['10) 0.0 0            \n'])
        f.writelines([' Z) 0 			\n\n\n'])
        f.close()
    def PsfFunc(conffile, ParamDict, FitDict, No, RunNo):
        f=open(config_file, 'ab')
        f.writelines([' 0) psf\n'])
        f.writelines([' 1) ', str(ParamDict[RunNo][No][2][0]), ' ', \
                              str(ParamDict[RunNo][No][2][1]), ' ', \
                              str(FitDict[No][1][0]),   ' ', \
                              str(FitDict[No][1][1]),   '\n'])
        f.writelines([' 3) ', str(ParamDict[RunNo][No][3]), ' ', \
                              str(FitDict[No][2]),  '\n'])
        f.writelines([' Z) 0                    \n\n\n'])
        f.close()
    def SkyFunc(conffile, ParamDict, FitDict, No, RunNo):
        f=open(config_file, 'ab')
        f.writelines([' 0) sky\n'])
        f.writelines([' 1) ', str(ParamDict[RunNo][No][2]), \
                    '      ', str(FitDict[No][1]), '\n'])
        f.writelines([' 2) 0.000      0       \n',\
                      ' 3) 0.000      0       \n',\
                      ' Z) 0                  \n\n\n'])
        f.writelines(['# Neighbour sersic function\n\n'])
        f.close()
    
    def DecideFitting(ParamDict, RunNo):
        No = RunNo + 1
        FitDict = {}
        #print ParamDict 
        if No == 1:
            for j in range(len(ParamDict[RunNo])):
                i = j + 1
                FitDict[i] = {} 
                if ParamDict[RunNo][i][1] == 'sersic' and ParamDict[RunNo][i][11] == 'Main':
                    FitDict[i][1] = [1, 1]
                    FitDict[i][2] = 1 
                    FitDict[i][3] = 1 
                    FitDict[i][4] = 0
                    FitDict[i][5] = 1       
                    FitDict[i][6] = 1    
                if ParamDict[RunNo][i][1] == 'expdisk' and ParamDict[RunNo][i][11] == 'Main':
                    FitDict[i][1] = [1, 1]
                    FitDict[i][2] = 1 
                    FitDict[i][3] = 1 
                    FitDict[i][4] = 1       
                    FitDict[i][5] = 1  
                #In the first run the bar and point will not be fitted.
                #The following is to keep the order of ParamDict and FitDict 
                if ParamDict[RunNo][i][1] == 'bar' and ParamDict[RunNo][i][11] == 'Main':
                    FitDict[i][1] = [1, 1]
                    FitDict[i][2] = 1 
                    FitDict[i][3] = 1 
                    FitDict[i][4] = 1
                    FitDict[i][5] = 1       
                    FitDict[i][6] = 1    
                if ParamDict[RunNo][i][1] == 'psf' and ParamDict[RunNo][i][11] == 'Main':
                    FitDict[i][1] = [1, 1]
                    FitDict[i][2] = 1 
                if ParamDict[RunNo][i][1] == 'sky':
                    FitDict[i][1] = 1
                    FitDict[i][2] = 0 
                    FitDict[i][3] = 0 
                if ParamDict[RunNo][i][1] == 'sersic' and ParamDict[RunNo][i][11] == 'Other':
                    FitDict[i][1] = [1, 1]
                    FitDict[i][2] = 1 
                    FitDict[i][3] = 1 
                    FitDict[i][4] = 1
                    FitDict[i][5] = 1       
                    FitDict[i][6] = 1    
        if No == 2:
            for j in range(len(ParamDict[RunNo])):
                i = j + 1
                FitDict[i] = {}  
                if ParamDict[RunNo][i][1] == 'sersic' and ParamDict[RunNo][i][11] == 'Main':
                    FitDict[i][1] = [1, 1]
                    FitDict[i][2] = 1 
                    FitDict[i][3] = 1 
                    FitDict[i][4] = 0
                    FitDict[i][5] = 1       
                    FitDict[i][6] = 1    
                if ParamDict[RunNo][i][1] == 'expdisk' and ParamDict[RunNo][i][11] == 'Main':
                    FitDict[i][1] = [1, 1]
                    FitDict[i][2] = 1 
                    FitDict[i][3] = 1 
                    FitDict[i][4] = 1       
                    FitDict[i][5] = 1  
                if ParamDict[RunNo][i][1] == 'bar' and ParamDict[RunNo][i][11] == 'Main':
                    FitDict[i][1] = [1, 1]
                    FitDict[i][2] = 1 
                    FitDict[i][3] = 1 
                    FitDict[i][4] = 1
                    FitDict[i][5] = 1       
                    FitDict[i][6] = 1    
                if ParamDict[RunNo][i][1] == 'psf' and ParamDict[RunNo][i][11] == 'Main':
                    FitDict[i][1] = [1, 1]
                    FitDict[i][2] = 1 
                if ParamDict[RunNo][i][1] == 'sky':
                    FitDict[i][1] = 1
                    FitDict[i][2] = 0 
                    FitDict[i][3] = 0 
                if ParamDict[RunNo][i][1] == 'sersic' and ParamDict[RunNo][i][11] == 'Other':
                    FitDict[i][1] = [1, 1]
                    FitDict[i][2] = 1 
                    FitDict[i][3] = 1 
                    FitDict[i][4] = 1
                    FitDict[i][5] = 1       
                    FitDict[i][6] = 1    
        if No == 3:
            for j in range(len(ParamDict[RunNo])):
                i = j + 1
                FitDict[i] = {}  
                if ParamDict[RunNo][i][1] == 'sersic' and ParamDict[RunNo][i][11] == 'Main':
                    FitDict[i][1] = [1, 1]
                    FitDict[i][2] = 1 
                    FitDict[i][3] = 1 
                    FitDict[i][4] = 1
                    FitDict[i][5] = 1       
                    FitDict[i][6] = 1    
                if ParamDict[RunNo][i][1] == 'expdisk' and ParamDict[RunNo][i][11] == 'Main':
                    FitDict[i][1] = [1, 1]
                    FitDict[i][2] = 1 
                    FitDict[i][3] = 1 
                    FitDict[i][4] = 1       
                    FitDict[i][5] = 1    
                if ParamDict[RunNo][i][1] == 'bar' and ParamDict[RunNo][i][11] == 'Main':
                    FitDict[i][1] = [1, 1]
                    FitDict[i][2] = 1 
                    FitDict[i][3] = 1 
                    FitDict[i][4] = 1
                    FitDict[i][5] = 1       
                    FitDict[i][6] = 1    
                if ParamDict[RunNo][i][1] == 'psf' and ParamDict[RunNo][i][11] == 'Main':
                    FitDict[i][1] = [1, 1]
                    FitDict[i][2] = 1 
                if ParamDict[RunNo][i][1] == 'sky':
                    FitDict[i][1] = 1
                    FitDict[i][2] = 0 
                    FitDict[i][3] = 0 
                if ParamDict[RunNo][i][1] == 'sersic' and ParamDict[RunNo][i][11] == 'Other':
                    FitDict[i][1] = [0, 0]
                    FitDict[i][2] = 0 
                    FitDict[i][3] = 0 
                    FitDict[i][4] = 0
                    FitDict[i][5] = 0       
                    FitDict[i][6] = 0    
        if No == 4:
            for j in range(len(ParamDict[RunNo])):
                i = j + 1
                FitDict[i] = {}  
                if ParamDict[RunNo][i][1] == 'sersic' and ParamDict[RunNo][i][11] == 'Main':
                    FitDict[i][1] = [1, 1]
                    FitDict[i][2] = 1 
                    FitDict[i][3] = 1 
                    FitDict[i][4] = 1
                    FitDict[i][5] = 1       
                    FitDict[i][6] = 1    
                if ParamDict[RunNo][i][1] == 'expdisk' and ParamDict[RunNo][i][11] == 'Main':
                    FitDict[i][1] = [1, 1]
                    FitDict[i][2] = 1 
                    FitDict[i][3] = 1 
                    FitDict[i][4] = 1       
                    FitDict[i][5] = 1    
                if ParamDict[RunNo][i][1] == 'sky':
                    FitDict[i][1] = 1
                    FitDict[i][2] = 0 
                    FitDict[i][3] = 0 
                if ParamDict[RunNo][i][1] == 'sersic' and ParamDict[RunNo][i][11] == 'Other':
                    FitDict[i][1] = [0, 0]
                    FitDict[i][2] = 0 
                    FitDict[i][3] = 0 
                    FitDict[i][4] = 0
                    FitDict[i][5] = 0       
                    FitDict[i][6] = 0    
        return FitDict

    def DecideHowToMove(ParamDict, RunNo):
        ContinueLoop = 0
        if ParamDict[RunNo][1][5] > 7: #If n > 7, re = 1 and n = 1
            ParamDict[RunNo][1][4] = 1.0
            ParamDict[RunNo][1][5] = 1.0
            ContinueLoop = 1
        if ParamDict[RunNo][1][5] > 7 and ParamDict[RunNo][3][5] > 2.15:
            ParamDict[RunNo][1][4] = 1.0
            ParamDict[RunNo][1][5] = 1.0
            ParamDict[RunNo][3][4] = copy.deepcopy(ParamDict[RunNo - 1][3][4])
            ParamDict[RunNo][3][5] = 0.5
            ContinueLoop = 1
        if ParamDict[RunNo][1][6] < ParamDict[RunNo][3][6] or ParamDict[RunNo][3][5] > 2.15:
            ParamDict[RunNo][3][7] = copy.deepcopy(ParamDict[RunNo][1][7])
            ParamDict[RunNo][3][6] = copy.deepcopy(ParamDict[RunNo][1][6])
            ParamDict[RunNo][3][5] = 0.5
            ParamDict[RunNo][3][4] = copy.deepcopy(ParamDict[RunNo - 1][3][4]) 
            ParamDict[RunNo][3][3] = copy.deepcopy(ParamDict[RunNo - 1][3][3])
            ParamDict[RunNo][3][2][0] = copy.deepcopy(ParamDict[RunNo][1][2][0])
            ParamDict[RunNo][3][2][1] = copy.deepcopy(ParamDict[RunNo][1][2][1])
            ParamDict[RunNo][1][4] = 1.0
            ParamDict[RunNo][1][5] = 1.0
            ParamDict[RunNo][1][6] = copy.deepcopy(ParamDict[RunNo - 1][1][6])
            ContinueLoop = 1
        if ParamDict[RunNo][3][6] > 0.58:
            ParamDict[RunNo][3][7] = copy.deepcopy(ParamDict[RunNo][1][7])
            ParamDict[RunNo][3][6] = 0.3
            ParamDict[RunNo][3][5] = 0.5
            ContinueLoop = 1
        return ContinueLoop
    def DecideHowToMove1(ParamDict, RunNo):
        ContinueLoop = 0
        if ParamDict[RunNo][1][6] < ParamDict[RunNo][3][6] and ParamDict[RunNo][1][6] < 0.45:
            ParamDict[RunNo][3][7] = copy.deepcopy(ParamDict[RunNo][1][7])

            ParamDict[RunNo][1][2][0] = copy.deepcopy(ParamDict[RunNo - 1][1][2][0])
            ParamDict[RunNo][1][2][1] = copy.deepcopy(ParamDict[RunNo - 1][1][2][1])
            ParamDict[RunNo][1][3] = copy.deepcopy(ParamDict[RunNo - 1][1][3])
            ParamDict[RunNo][1][3] = ParamDict[RunNo][1][3] + 1.0 
            ParamDict[RunNo][1][4] = 1.0
            ParamDict[RunNo][1][5] = 1.0
            ParamDict[RunNo][1][6] = copy.deepcopy(ParamDict[RunNo - 1][1][6])
            ParamDict[RunNo][1][7] = copy.deepcopy(ParamDict[RunNo - 1][1][7])
 
#        ParamDict[RunNo][2][2][0] = copy.deepcopy(ParamDict[RunNo - 1][2][2][0])
#        ParamDict[RunNo][2][2][1] = copy.deepcopy(ParamDict[RunNo - 1][2][2][1])
#        ParamDict[RunNo][2][3] = copy.deepcopy(ParamDict[RunNo - 1][2][3])
#        ParamDict[RunNo][2][4] = copy.deepcopy(ParamDict[RunNo - 1][2][4])
#        ParamDict[RunNo][2][6] = copy.deepcopy(ParamDict[RunNo - 1][2][6])
#        ParamDict[RunNo][2][7] = copy.deepcopy(ParamDict[RunNo - 1][2][7])

            ParamDict[RunNo][3][2][0] = copy.deepcopy(ParamDict[RunNo - 1][3][2][0])
            ParamDict[RunNo][3][2][1] = copy.deepcopy(ParamDict[RunNo - 1][3][2][1])
            ParamDict[RunNo][3][3] = copy.deepcopy(ParamDict[RunNo - 1][3][3])
            ParamDict[RunNo][3][3] =  ParamDict[RunNo][3][3] 
            ParamDict[RunNo][3][4] = copy.deepcopy(ParamDict[RunNo - 1][3][4]) 
            ParamDict[RunNo][3][5] = 0.5
            ParamDict[RunNo][3][6] = copy.deepcopy(ParamDict[RunNo - 1][3][6])
            ContinueLoop = 1
        elif ParamDict[RunNo][1][5] > 7 :
            ParamDict[RunNo][1][2][0] = copy.deepcopy(ParamDict[RunNo - 1][1][2][0])
            ParamDict[RunNo][1][2][1] = copy.deepcopy(ParamDict[RunNo - 1][1][2][1])
            ParamDict[RunNo][1][3] = copy.deepcopy(ParamDict[RunNo - 1][1][3])
            ParamDict[RunNo][1][3] = ParamDict[RunNo][1][3] + 2.0 
            ParamDict[RunNo][1][4] = 1.0
            ParamDict[RunNo][1][5] = 1.0
            ParamDict[RunNo][1][6] = copy.deepcopy(ParamDict[RunNo - 1][1][6])
            ParamDict[RunNo][1][7] = copy.deepcopy(ParamDict[RunNo - 1][1][7])
 
#        ParamDict[RunNo][2][2][0] = copy.deepcopy(ParamDict[RunNo - 1][2][2][0])
#        ParamDict[RunNo][2][2][1] = copy.deepcopy(ParamDict[RunNo - 1][2][2][1])
#        ParamDict[RunNo][2][3] = copy.deepcopy(ParamDict[RunNo - 1][2][3])
#        ParamDict[RunNo][2][4] = copy.deepcopy(ParamDict[RunNo - 1][2][4])
#        ParamDict[RunNo][2][6] = copy.deepcopy(ParamDict[RunNo - 1][2][6])
#        ParamDict[RunNo][2][7] = copy.deepcopy(ParamDict[RunNo - 1][2][7])
            ContinueLoop = 1

        return ContinueLoop

#The following function is to check whether the large bulge radii is real for
#deva + disk fitting
    def DecideHowToMove2(ParamDict, RunNo):
        ContinueLoop = 0
        if ParamDict[RunNo][1][4] > ParamDict[RunNo][2][4] * 3.0 and ParamDict[RunNo][1][3] > ParamDict[RunNo][2][3]:
            ParamDict[RunNo][1][2][0] = copy.deepcopy(ParamDict[RunNo-1][1][2][0])
            ParamDict[RunNo][1][2][1] = copy.deepcopy(ParamDict[RunNo-1][1][2][1])
            ParamDict[RunNo][1][3] = copy.deepcopy(ParamDict[RunNo - 1][1][3] + 0.5)
            ParamDict[RunNo][1][4] = copy.deepcopy(ParamDict[RunNo - 1][1][4] / 3.0)
            ParamDict[RunNo][1][6] = copy.deepcopy(ParamDict[RunNo - 1][1][6])
            ParamDict[RunNo][1][7] = copy.deepcopy(ParamDict[RunNo - 1][1][7])
            ParamDict[RunNo][2][2][0] = copy.deepcopy(ParamDict[RunNo-1][2][2][0])
            ParamDict[RunNo][2][2][1] = copy.deepcopy(ParamDict[RunNo-1][2][2][1])
            ParamDict[RunNo][2][3] = copy.deepcopy(ParamDict[RunNo - 1][2][3])
            ParamDict[RunNo][2][4] = copy.deepcopy(ParamDict[RunNo - 1][2][4])
            ParamDict[RunNo][2][5] = copy.deepcopy(ParamDict[RunNo - 1][2][5])
            ParamDict[RunNo][2][6] = copy.deepcopy(ParamDict[RunNo - 1][2][6])
            ParamDict[RunNo][3][2] = copy.deepcopy(ParamDict[RunNo - 1][3][2])
            ContinueLoop = 1
        return ContinueLoop

    #Write configuration file. RunNo is the number of iteration
    for RunNo in range(2):
        f_constrain = open(constrain_file, 'w')
        f_constrain.close()
        f=open(config_file,'w')
        f.write('# IMAGE PARAMETERS\n')
        f.writelines(['A) ', str(cutimage), '	# Input data image',\
                      ' (FITS file)\n'])
        f.writelines(['B) ', str(outfile), '		# Name for',\
                      ' the output image\n'])
        f.writelines(['C) ', str(whtimage), '		# Noise image name', \
                      ' (made from data if blank or "none")\n'])
        f.writelines(['D) ', str(psffile), '			# Input PSF', \
                      ' image for convolution (FITS file)\n'])
        f.writelines(['E) 1			# PSF oversampling factor '\
                      'relative to data\n'])
        f.writelines(['F) ', str(mask_file), '		# Bad pixel',
                      ' mask(FITS image or ASCII coord list)\n'])
        f.writelines(['G) ', str(constrain_file), '       # File with'\
                      ' parameter constraints (ASCII file)\n'])
        f.writelines(['H) 1 ', str(NXPTS), ' 1 ', str(NYPTS), '		#',\
                      ' Image region to fit (xmin xmax ymin ymax)\n'])
#        f.writelines(['I) ', str(NXPTS), ' ', str(NYPTS),	'	#',\
#                      ' Size of convolution box (x y)\n'])
        f.writelines(['I) ', str(120), ' ', str(120),      '       #',\
                      ' Size of convolution box (x y)\n'])
        f.writelines(['J) ', str(mag_zero), '		# Magnitude',\
                      ' photometric zeropoint\n'])
        f.writelines(['O) regular			# Display type',\
                      ' (regular, curses, both)\n'])
        f.writelines(['P) 0			# Create output image only?',\
                      ' (1=yes; 0=optimize)\n'])
        f.writelines(['S) 0			# Modify/create',\
                     ' objects interactively?\n\n\n'])
        f.close()
        FitDict = DecideFitting(ParamDict, RunNo)
        for i in range(len(ParamDict[RunNo])):
            if ParamDict[RunNo][i + 1][1] == 'sersic':
                SersicFunc(config_file, ParamDict, FitDict, i+1, RunNo)
                if ParamDict[RunNo][i + 1][11] == 'Main':
                    if RunNo == 1: #Second run
                        SersicMainConstrain(constrain_file, i + 1, 1.0, 0)#2. * ParamDict[RunNo][2][4])
                    else: #Second run
                        SersicMainConstrain(constrain_file, i + 1, c.center_constrain, 0)
                else:
                    SersicConstrain(constrain_file, i + 1)
 
            if ParamDict[RunNo][i + 1][1] == 'expdisk':
                if RunNo + 1 == 0:
                    pass
                else:
                    ExpFunc(config_file, ParamDict, FitDict, i + 1, RunNo)
                    ExpdiskConstrain(constrain_file, i + 1)
            if ParamDict[RunNo][i + 1][1] == 'bar':
                if RunNo + 1 == 0:
                    pass
                else:
                    BarFunc(config_file, ParamDict, FitDict, i+1, RunNo)
                    BarConstrain(constrain_file, i + 1)

            if ParamDict[RunNo][i + 1][1] == 'psf':
                if RunNo + 1 < 3:
                    pass
                else:
                    PsfFunc(config_file, ParamDict, FitDict, i+1, RunNo)
                    PsfConstrain(constrain_file, i + 1)

            if  ParamDict[RunNo][i + 1][1] == 'sky':
                if RunNo == 9: 
                    SkyConstrain(constrain_file, i + 1)
                SkyFunc(config_file, ParamDict, FitDict, i+1, RunNo) 

#        print ParamDict
#        raw_input('Waiting >>> ')
#        print 'Waiting'
        if exists('fit.log'):
            os.remove('fit.log') 
        try:
            cmd = str(c.GALFIT_PATH) + ' ' + config_file
            os.system(cmd)
        #ReadLog(ParamDict, 2) => this function reads fig.log  depends on 
        #whether, for example, the expdisk function is the only one function
        #for fitting. It can happends, for example, the first fitting where
        #we can find the best  initial values for the rest of it.
        #ReadLog(ParamDict, 1) reads the fit.log in the order, ie. sersic, 
        #expdisk, other sersic etc
            ParamDict = ReadLog(ParamDict, 1, RunNo)
            if DecideHowToMove2(ParamDict, RunNo + 1):
                pass
            else:
                break
            if RunNo == -10:
                for i in range(3):
                    try:
                        ParamDict[RunNo][i + 2][2][0]  = ParamDict[RunNo][1][2][0]
                        ParamDict[RunNo][i + 2][2][1]  = ParamDict[RunNo][1][2][1]
                    except:
                        pass
                ParamDict[RunNo][2][3] = ParamDict[RunNo][1][3] 
        except:
#            print 'Hello'
            ParamDict[RunNo + 1] = copy.deepcopy(ParamDict[RunNo])
#            ParamDict[RunNo + 1][3][5] = 0.5
#            ParamDict[RunNo + 1][1][4] = 0.5
#            ParamDict[RunNo + 1][1][5] = 1.5

#            for i in range(3):
#                if abs(ParamDict[RunNo][i + 1][2][0] - ParamDict[RunNo + 1][i + 1][2][0]) > 2.5 or abs(ParamDict[RunNo][i + 1][2][1] - ParamDict[RunNo + 1][i + 1][2][1]) > 2.5: 
                    
                    
#        print 'Updated \n', ParamDict, '\n\n'
                        

