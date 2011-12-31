#####################################
#
# NAME: config_twostep.py
#
# DESCRIPTION: This script will perform an
#  iterative fitting process when invoked using
#  "detail = True" command in the config file.
#  as it is currently configured, the program does
#  a serexp fit first to estimate the sky. Then fixes
#  the sky at that value. Then it does a 'ser' fit.
#  this tells it the total sersic index for the galaxy.
#  Then it fits a 'dev' profile. Finally it does a 'devexp'
#  and 'serexp' profile choosing the starting parameters
#  based on the 'ser' fit sersic index.
#
#########################################

# import relevant python module
import os
import sys
import pyfits
from os.path import exists
from numpy import log10
import numpy as n
import copy
import numpy.ma as ma
import time
import traceback

# import pymorph-specific modules 
import config as c
from runsexfunc import *
from readlog import ReadLog
from cosmocal import cal
from flagfunc import *
try:
    from utilities import WriteDbDetail
except:
    print 'No database'


# declare the class

class ConfigIter:
    """The class making configuration file for GALFIT. The configuration file 
       consists of bulge and disk component of the object and only Sersic 
       component for the neighbours, if any. The sky is always fixed and has
       the value of SExtractor. The disk/boxy parameter is also fixed to zero.
       The initial value for Sersic index 'n' is 4.The configuration file has 
       the name G_string(galid).in. The output image has the name 
       O_string(galid).fits"""
    def __init__(self, cutimage, whtimage, xcntr, ycntr, NXPTS, NYPTS, line_s, psffile, z):
        self.cutimage = cutimage
        self.line_s  = line_s
	self.whtimage = whtimage
        self.xcntr = xcntr
        self.ycntr = ycntr
        self.NXPTS = NXPTS
        self.NYPTS = NYPTS 
        self.psffile = psffile

        # check whether the variable is defined, if not, define it 
        try: 
            c.center_constrain = c.center_constrain
        except:
            c.center_constrain = 2.0
        # run the fitting iteration
        self.confiter    = confiter(cutimage, whtimage, xcntr, ycntr, NXPTS, NYPTS, line_s, psffile, z)
        
#end of class definitions

# function definitions 



def confiter(cutimage, whtimage, xcntr, ycntr,
             NXPTS, NYPTS, line_s, psffile, z):
    # Run SExtractor 
    RunSex(c.datadir+cutimage, c.datadir+whtimage,
           'TEMP.SEX.cat', 9999, 9999, 0)
    # define improtant variables in the namespace
    imagefile = c.imagefile
    sex_cata = c.sex_cata
    threshold = c.threshold
    thresh_area = c.thresh_area
    mask_reg = c.mask_reg
    c.sersic_loc  = 2
    c.sky_loc = 1
    try:
        ComP = c.components 
    except:
        ComP = ['bulge', 'disk']
    if len(ComP) == 0:
        ComP = ['bulge', 'disk']
    values = line_s.split()
    # create filenames
    outfile   = 'O_' + c.fstring + '.fits' # galfit output
    mask_file = 'M_' + c.fstring + '.fits' # mask file
    config_file = 'G_' + c.fstring + '.in' #GALFIT configuration file
    constrain_file = c.fstring + '.con' #galfit constraint file


    xcntr_o  = xcntr  #x center of the object
    ycntr_o  = ycntr  #y center of the object
    mag    = float(values[7]) #Magnitude
    radius = float(values[9]) #Half light radius
    mag_zero = c.mag_zero #magnitude zero point
    sky	 = float(values[10]) #sky 
    pos_ang = float(values[11]) - 90.0 #position angle
    axis_rat = 1.0/float(values[12]) #axis ration b/a
    area_o = float(values[13])   # object's area
    major_axis = float(values[14])  #major axis of the object

    # calculate a 40kpc bulge radius to use for maximum bulge size.
    # This will only work if z is supplied, otherwise, a default maximum
    # is used in the constraint file
    if z > 0 and z < 10:
        try:
            KpCArc = cal(z, c.H0, c.WM, c.WV, c.pixelscale)[3]
            max_rad = 40.0/KpCArc
        except:
            KpCArc = 9999.0
            max_rad = c.URe
            
    else:
        max_rad = c.URe
        KpCArc = 9999.0

    # Modify flags for the detailed fitting process
    if c.Flag & 2**GetFlag('FIT_BULGE_CNTR'):
        c.Flag -= 2**GetFlag('FIT_BULGE_CNTR')
    if c.Flag & 2**GetFlag('FIT_DISK_CNTR'):
        c.Flag -= 2**GetFlag('FIT_DISK_CNTR')
    if c.Flag & 2**GetFlag('FIT_SKY'):
        c.Flag -= 2**GetFlag('FIT_SKY')
    if c.Flag & 2**GetFlag('FIT_BULGE'):
        c.Flag -= 2**GetFlag('FIT_BULGE')
    if c.Flag & 2**GetFlag('FIT_DISK'):
        c.Flag -= 2**GetFlag('FIT_DISK')
    if c.Flag & 2**GetFlag('FIT_POINT'):
        c.Flag -= 2**GetFlag('FIT_POINT')
        

    # define initial fitting parameters for single sersic fit
    ParamDict = {}
    ErrDict = {}
    ParamDict[0] = {}
    ErrDict[0] = {}
    BlankDict = {} # to be used in cases where the program crashes
    BlankDict[0] = {}
    
    # Add components
    AdComp = 1
    if 'bulge' in ComP:
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
        
        BlankDict[0][AdComp] = {}
        #Bulge Parameters
        BlankDict[0][AdComp][1] = 'sersic'
        BlankDict[0][AdComp][2] = [xcntr_o, ycntr_o]
        BlankDict[0][AdComp][3] = 9999 
        BlankDict[0][AdComp][4] = 9999 
        BlankDict[0][AdComp][5] = 9999
        BlankDict[0][AdComp][6] = 9999
        BlankDict[0][AdComp][7] = 9999
        BlankDict[0][AdComp][8] = 9999
        BlankDict[0][AdComp][9] = 9999
        BlankDict[0][AdComp][11] = 'Main'
        AdComp += 1
    if 'disk' in ComP:
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

        BlankDict[0][AdComp] = {}
        BlankDict[0][AdComp][1] = 'expdisk'
        BlankDict[0][AdComp][2] = [xcntr_o, ycntr_o]
        BlankDict[0][AdComp][3] = 9999#mag 
        BlankDict[0][AdComp][4] = 9999#radius
        BlankDict[0][AdComp][5] = 9999#axis_rat
        BlankDict[0][AdComp][6] = 9999#pos_ang
        BlankDict[0][AdComp][7] = 9999
        BlankDict[0][AdComp][8] = 9999
        BlankDict[0][AdComp][11] = 'Main'
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
            if abs(xcntr_n - xcntr_o) > NXPTS / 2.0 + c.avoidme or \
               abs(ycntr_n - ycntr_o) > NYPTS / 2.0 + c.avoidme or \
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

                BlankDict[0][AdComp] = {}
                BlankDict[0][AdComp][1] = 'sersic'
                BlankDict[0][AdComp][2] = [xn, yn]
                BlankDict[0][AdComp][3] = mag
                BlankDict[0][AdComp][4] = radius
                BlankDict[0][AdComp][5] = 4.0
                BlankDict[0][AdComp][6] = axis_rat
                BlankDict[0][AdComp][7] = pos_ang
                BlankDict[0][AdComp][8] = 0
                BlankDict[0][AdComp][9] = 0
                BlankDict[0][AdComp][11] = 'Other'

                isneighbour = 1
                AdComp += 1
                c.Flag = c.Flag & 2**GetFlag('NEIGHBOUR_FIT')
        except:
            pass
    f_constrain.close()
    #Sky component
    ParamDict[0][AdComp] = {}
    ParamDict[0][AdComp][1] = 'sky'
    ParamDict[0][AdComp][2] = c.SexSky
    ParamDict[0][AdComp][3] = 0
    ParamDict[0][AdComp][4] = 0
    ParamDict[0][AdComp][5] = 0
    ParamDict[0][AdComp][11] = 'Other'

    BlankDict[0][AdComp] = {}
    BlankDict[0][AdComp][1] = 'sky'
    BlankDict[0][AdComp][2] = c.SexSky
    BlankDict[0][AdComp][3] = 0
    BlankDict[0][AdComp][4] = 0
    BlankDict[0][AdComp][5] = 0
    BlankDict[0][AdComp][11] = 'Other'
    #Write Sersic function
    ErrDict[0][1] = {}
    ErrDict[0][1][1] = [9999, 9999]
    ErrDict[0][1][2] = 9999
    ErrDict[0][1][3] = 9999
    ErrDict[0][1][4] = 9999
    ErrDict[0][2] = {}
    ErrDict[0][2][1] = [9999, 9999]
    ErrDict[0][2][2] = 9999
    ErrDict[0][2][3] = 9999

    c.SkyNo = AdComp
    c.Chi2DOFArr = []
    c.FitArr = []
    c.RadArr = []
    c.CntrDevArr = []
    c.ErrArr = []
    c.ParamDictBook = copy.deepcopy(ParamDict)
    bad_fit = 0
    failed_ser = 0
    bt_fit = -1
    failed_sky = 0
    bad_sky = 0
    
   #Write configuration file. RunNo is the number of iteration
    for RunNo, fit_type in zip(range(5), ['sky','ser','dev','devexp','serexp']):
        # note that the 'sky' fit is currently serexp
        #force batch fitting using standard operations
        
        run_flag = c.Flag
        print 'run_flag', run_flag
        if exists('fit.log'):
            print "removing fit log"
            os.system('rm fit.log')

        print "RunNo ", RunNo, "fit type ", fit_type
        f_constrain = open(constrain_file, 'w')
        f_constrain.close()
        f=open(config_file,'w')
        f.write('# IMAGE PARAMETERS\n')
        f.writelines(['A) ', c.datadir+str(cutimage), '	# Input data image',\
                      ' (FITS file)\n'])
        f.writelines(['B) ', str(outfile), '		# Name for',\
                      ' the output image\n'])
        f.writelines(['C) ', c.datadir + str(whtimage), '		# Noise image name', \
                      ' (made from data if blank or "none")\n'])
        f.writelines(['D) ', c.datadir+str(psffile), '			# Input PSF', \
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

        FitDict, run_flag = DecideFitting(ParamDict, RunNo,fit_type, bad_fit, failed_ser, run_flag, failed_sky, bad_sky)            

        for i in range(len(ParamDict[RunNo])):
            if ParamDict[RunNo][i + 1][1] == 'sersic':
                SersicFunc(config_file, ParamDict, FitDict, i+1, RunNo)
                if ParamDict[RunNo][i + 1][11] == 'Main':
                    run_flag += 2**GetFlag('FIT_BULGE')
                    print "ADDING FLAG FIT_BULGE"
                    if RunNo > 1 and not bad_fit and ParamDict[c.sersic_loc][1][4]*2 < max_rad: # later fits
                        SersicMainConstrain(constrain_file, i + 1, 2.0 , ParamDict[c.sersic_loc][1][4]*2)
                    else:
                        SersicMainConstrain(constrain_file, i + 1, c.center_constrain, max_rad)
                else:
                    SersicConstrain(constrain_file, i + 1)
 
            if ParamDict[RunNo][i + 1][1] == 'expdisk':
                if RunNo in [1 ,2]:
                    pass #dont fit a disk in the cases of ser or dev fits 
                else:
                    run_flag += 2**GetFlag('FIT_DISK')
                    print "ADDING FLAG FIT_DISK"
                    ExpFunc(config_file, ParamDict, FitDict, i + 1, RunNo)
                    if RunNo > 0 and not bad_fit: # later runs 
                        ExpdiskConstrain(constrain_file, i + 1, 2.0)
                    else:
                        ExpdiskConstrain(constrain_file, i + 1, c.center_constrain)

            if  ParamDict[RunNo][i + 1][1] == 'sky':
                if RunNo > 0 and not bad_sky: 
                    SkyConstrain(constrain_file, i + 1, ParamDict[c.sky_loc][i+1][2])
                else:
                    SkyConstrain(constrain_file, i + 1, c.SexSky)
                SkyFunc(config_file, ParamDict, FitDict, i+1, RunNo) 

        # only if forcing normal fit
        if RunNo > 2 or RunNo == 0: #later fits after single ser fit
            bt_range = add_constrain(constrain_file, bt_fit, fit_type)
        cmd = str(c.GALFIT_PATH) + ' ' + config_file



        os.system(cmd)
        # ReadLog(ParamDict, 2) => this function reads fig.log
        # depends on whether, for example, the expdisk function
        # is the only one function for fitting. It can happends,
        # for example, the first fitting where we can find the best
        # initial values for the rest of it. ReadLog(ParamDict, 1)
        # reads the fit.log in the order, ie. sersic,
        # expdisk, other sersic etc
        try:
            print 'readlog'
            ParamDict, ErrDict, Chi2DOF = ReadLog(ParamDict, ErrDict, 1, RunNo, detail = True)
            print 'log read'
            print 'chi2dof ',Chi2DOF
            
            c.ParamDictBook[RunNo+1] = copy.deepcopy(ParamDict[RunNo+1])

            # Set limit flags
            mag_b = ParamDict[RunNo + 1][1][3]
            re = ParamDict[RunNo+1][1][4]
            SersicIndex = ParamDict[RunNo+1][1][5]
            if abs(mag_b - c.UMag) < 0.2 or abs(mag_b - c.LMag) < 0.2 or \
                   abs(re - max_rad) < 1.0 or abs(re - c.LRe) < 0.1 or \
                   abs(SersicIndex - c.LN) < 0.03 or abs(SersicIndex - c.UN) < 0.5:
                run_flag += 2**GetFlag('BULGE_AT_LIMIT')
                print "ADDING FLAG BULGE_AT_LIMIT"

                if abs(re - max_rad) < 1.0 or abs(re - c.LRe) < 0.1:
                    run_flag += 2**GetFlag('RE_AT_LIMIT')
                    print "ADDING FLAG RE_AT_LIMIT"
                if abs(SersicIndex - c.LN) < 0.03 or abs(SersicIndex - c.UN) < 0.5:
                    run_flag += 2**GetFlag('N_AT_LIMIT')
                    print "ADDING FLAG N_AT_LIMIT"
            if fit_type in ['devexp', 'serexp']:
                mag_d = ParamDict[RunNo + 1][2][3]
                rd = ParamDict[RunNo+1][2][4]

                fb = 10**(-0.4*mag_b)
                fd = 10**(-0.4*mag_d)

                if abs(mag_d - c.UMag) < 0.2 or abs(mag_d - c.LMag) < 0.2 or \
                       abs(rd - c.LRd) < 0.1 or abs(rd - c.URd) < 1.0:
                    run_flag += 2**GetFlag('DISK_AT_LIMIT')
                    print "ADDING FLAG DISK_AT_LIMIT"

                if abs((re/rd) - 1.0) < 0.02 or abs((re/rd) - 0.1) < 0.02:
                    run_flag += 2**GetFlag('RERD_AT_LIMIT')
                    print "ADDING FLAG RERD_AT_LIMIT"
                try:
                    BT = fb / (fb + fd)
                except:
                    BT = 9999.0

                if abs(BT - bt_range[0]) < .02 or abs(BT - bt_range[1]) < .02:
                    run_flag += 2**GetFlag('BT_AT_LIMIT')
                    print "ADDING FLAG BT_AT_LIMIT"
        except:
            print "failure at readlog!!!"
            c.ParamDictBook[RunNo+1] = copy.deepcopy(BlankDict[0])
            ParamDict[RunNo + 1] = copy.deepcopy(BlankDict[0])
            ErrDict[RunNo + 1] = copy.deepcopy(ErrDict[0])
            if fit_type == 'ser':
                failed_ser = 1 #track the failure of the single sersic fit
                run_flag += 2**GetFlag('DETAIL_FAILED')
                c.Flag += 2**GetFlag('DETAIL_FAILED')
                print "ADDING FLAG DETAIL_FAILED"
            elif fit_type == 'sky':
                failed_sky = 1 #track that we did not fix the sky
                run_flag += 2**GetFlag('DETAIL_FAILED')
                c.Flag += 2**GetFlag('DETAIL_FAILED')
                print "ADDING FLAG DETAIL_FAILED"
        try:
            c.ErrArr.append(FractionalError(ParamDict, ErrDict, RunNo + 1))
        except:
            c.ErrArr.append(9999.0)
        c.Chi2DOFArr.append(Chi2DOF)

        if fit_type != 'serexp': #set flags for other fit types
            bad_sky, bad_fit, bt_fit = DecideHowToMove2(ParamDict, RunNo + 1,fit_type,KpCArc, failed_ser, failed_sky )

            try:
                plot_name = 'P_' + str(cutimage)[6:-5] \
                          +'_'+fit_type+ '.png'
                if exists(plot_name):	
                    os.system('rm ' + plot_name)
                GoodNess = PlotFunc(cutimage, outfile, mask_file, 
                                    ParamDict[RunNo + 1][1][2][0], #xctr
                                    ParamDict[RunNo + 1][1][2][1], #yctr
                                    ParamDict[RunNo + 1][1][c.SkyNo][2], #sky
                                    c.skysig, save_name = plot_name)
                Goodness = GoodNess.plot_profile
            except:
                Goodness = 9999
                run_flag += GetFlag('PLOT_FAIL')
                print "ADDING FLAG PLOT_FAIL"
                                
            if Goodness < c.Goodness:
                run_flag += 2**GetFlag('SMALL_GOODNESS')
                print "ADDING FLAG SMALL_GOODNESS"
            if Chi2DOF > c.chi2sq:
                if chi2nu != 9999:
                    run_flag += 2**GetFlag('LARGE_CHISQ')
                    print "ADDING FLAG LARGE_CHISQ"

            print 'write flags ', run_flag 
            print 'writing db'
            try:
                WriteDbDetail(cutimage.split('.')[0], c.ParamDictBook[RunNo+1], ErrDict[RunNo + 1], c.SexSky, c.ParamDictBook[RunNo+1][c.SkyNo][2], RunNo, run_flag, Chi2DOF, model_type = fit_type, goodness = Goodness)
            except:
                print 'No database'
                traceback.print_exc()
                time.sleep(10)
        fit_log = 'fit_%s.log' %fit_type

        print "saving fit log to fit log %s" %(fit_log)

        f_fit = open(fit_log,'a')
        if exists('fit.log'):
            for line in open('fit.log','r'):
                f_fit.writelines([str(line)])
        f_fit.close()

        
        for mv_file_nm in [outfile, config_file, constrain_file]:
            os.system('cp %s %s_%s.%s' %(mv_file_nm, mv_file_nm.split('.')[0],
                                          fit_type, mv_file_nm.split('.')[1]))

        if fit_type == 'serexp':
            c.Flag = run_flag

    return

def add_constrain(constrain_file, bt, fit_type):
    f_constrain = open(constrain_file, 'ab')
    f_constrain.write('1-2     x      -0.001 to 0.001\n')
    f_constrain.write('1-2     y      -0.001 to 0.001\n')

    if fit_type == 'sky':
        rerdmin = 0.05
        rerdmax = 1.0
        bt_min = 0
        bt_max = 1
    else:
        rerdmin = 0.05
        rerdmax = 1.0
        absmin_bt = .0000001
        absmax_bt = .9999999

        min_down = .3
        min_up = .3
        d_trans1 = .4
        d_trans2 = .49
        u_trans1 = .3
        u_trans2 = .49 

        down_slope = (d_trans1 - min_down)/(d_trans2 - d_trans1)
        up_slope = (u_trans2 - min_up)/(u_trans2 - u_trans1)

        if bt <= d_trans1:
            bt_min = absmin_bt
        elif bt <= d_trans2:
            bt_min = bt - min_down - down_slope*(d_trans2 - bt) 
        else:
            bt_min = bt - min_down

        if bt <= u_trans1:
            bt_max = min_up + bt
        elif bt <= u_trans2:
            bt_max = bt + min_up + up_slope*(bt - u_trans1)
        else:
            bt_max = absmax_bt

        
        bd_min = (1.0/bt_min - 1.0)**(-1)
        bd_max = (1.0/bt_max - 1.0)**(-1)
        mag_min = -2.5*n.log10(bd_max)
        mag_max = -2.5*n.log10(bd_min)
    
        f_constrain.write('1-2      mag     ' + str(mag_min) +\
                          ' to ' + str(mag_max) + '\n')

    f_constrain.write('1/2     re     %s   %s\n' %(rerdmin, rerdmax))
    f_constrain.close()

    return bt_min, bt_max

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
                      ' to ' + str(c.LMag - 2.0) + '\n')
    if re_con == 0:
        f_constrain.write(str(cO) + '      re     ' + str(0.2) +\
                          ' to ' + str(c.URe) + '\n')
    else:
        f_constrain.write(str(cO) + '      re     ' + str(0.5) +\
                          ' to ' + str(re_con) + '\n')
        
    f_constrain.write(str(cO) + '      q       0.05 to 0.95\n')
    f_constrain.write(str(cO) + '      pa       -360.0 to 360.0\n')
    f_constrain.close()

def ExpdiskConstrain(constrain_file, cO, cen_con, rs_con = 0):
    f_constrain = open(constrain_file, 'ab')
    f_constrain.write(str(cO) + '       x       ' + \
                      str(-cen_con) + '     ' + \
                      str(cen_con) + '\n')
    f_constrain.write(str(cO) + '       y       ' + \
                      str(-cen_con) + '     ' + \
                      str(cen_con) + '\n')
    f_constrain.write(str(cO) + '     mag     ' + str(c.UMag) + \
                      ' to ' + str(c.LMag) + '\n')
    if rs_con == 0:
        f_constrain.write(str(cO) + '      rs     ' + str(c.LRd) + \
                          ' to ' + str(c.URd) + '\n')
    else:
        f_constrain.write(str(cO) + '      rs     ' + str(c.LRd) + \
                          ' to ' + str(rs_con) + '\n')
    f_constrain.write(str(cO) + '      q       0.0 to 1.0\n')
    f_constrain.write(str(cO) + '      pa       -360.0 to 360.0\n')
    f_constrain.close()


def SersicConstrain(constrain_file, cO):
    f_constrain = open(constrain_file, 'ab')
    f_constrain.write(str(cO) + '      n      0.02 to 20.0  \n')
    f_constrain.write(str(cO) + '     mag    -100.0 to 100.0\n')
    f_constrain.write(str(cO) + '      re      0.0 to 500.0\n')
    f_constrain.write(str(cO) + '      q       0.0 to 1.0\n')
    f_constrain.write(str(cO) + '      pa    -360.0 to 360.0\n')
    f_constrain.close()
    
def SkyConstrain(constrain_file, cO, SkyValToCon):
    f_constrain = open(constrain_file, 'ab')
    f_constrain.write(str(cO) + '      sky      ' +
                      str(SkyValToCon * 1.0-0.1) + ' to   ' +
                      str(SkyValToCon * 1.0+0.1) + '  \n')
    f_constrain.close()


def SersicFunc(conffile, ParamDict, FitDict, No, RunNo):
    f=open(conffile, 'ab')
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
    f=open(conffile, 'ab')
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
def SkyFunc(conffile, ParamDict, FitDict, No, RunNo):
    f=open(conffile, 'ab')
    f.writelines([' 0) sky\n'])
    f.writelines([' 1) ', str(ParamDict[RunNo][No][2]), \
                '      ', str(FitDict[No][1]), '\n'])
    f.writelines([' 2) 0.000      0       \n',\
                  ' 3) 0.000      0       \n',\
                  ' Z) 0                  \n\n\n'])
    f.writelines(['# Neighbour sersic function\n\n'])
    f.close()

    
def DecideFitting(ParamDict, RunNo, fit_type, bad_fit, failed_ser, run_flag, failed_sky, bad_sky):
    FitDict = {}
    #print ParamDict 
    if fit_type =='sky': # Perform a serexp fit for sky value 
        for j in range(len(ParamDict[RunNo])):
            i = j + 1
            FitDict[i] = {}  
            if ParamDict[RunNo][i][1] == 'sersic' and ParamDict[RunNo][i][11] == 'Main':
                FitDict[i][1] = [1, 1] 
                FitDict[i][2] = 1 
                FitDict[i][3] = 1
                run_flag += 2**GetFlag('FIT_BULGE_CNTR')
                print "ADDING FLAG FIT_BULGE_CNTR"
                FitDict[i][4] = 1 # free n
                FitDict[i][5] = 1       
                FitDict[i][6] = 1    
            if ParamDict[RunNo][i][1] == 'expdisk' and ParamDict[RunNo][i][11] == 'Main':
                FitDict[i][1] = [1, 1]
                FitDict[i][2] = 1 
                FitDict[i][3] = 1 
                FitDict[i][4] = 1       
                FitDict[i][5] = 1
                run_flag += 2**GetFlag('FIT_DISK_CNTR')
                print "ADDING FLAG FIT_DISK_CNTR"
            if ParamDict[RunNo][i][1] == 'sky':
                FitDict[i][1] = 1
                FitDict[i][2] = 0 
                FitDict[i][3] = 0
                run_flag += 2**GetFlag('FIT_SKY')
                print "ADDING FLAG FIT_SKY"
            if ParamDict[RunNo][i][1] == 'sersic' and ParamDict[RunNo][i][11] == 'Other':
                #should this be all 1s or 0s?
                FitDict[i][1] = [1, 1]
                FitDict[i][2] = 1 
                FitDict[i][3] = 1 
                FitDict[i][4] = 1
                FitDict[i][5] = 1       
                FitDict[i][6] = 1    

    if fit_type == 'ser': # perform a single-sersic fit
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

                run_flag += 2**GetFlag('FIT_BULGE_CNTR')
                print "ADDING FLAG FIT_BULGE_CNTR"
            #In the first run the disk will not be fitted.
            #This is to keep the order of ParamDict and FitDict 
            if ParamDict[RunNo][i][1] == 'expdisk' and ParamDict[RunNo][i][11] == 'Main':
                FitDict[i][1] = [1, 1]
                FitDict[i][2] = 1 
                FitDict[i][3] = 1 
                FitDict[i][4] = 1       
                FitDict[i][5] = 1
                
            if ParamDict[RunNo][i][1] == 'sky':
                FitDict[i][1] = bad_sky
                FitDict[i][2] = 0
                FitDict[i][3] = 0
                if bad_fit:
                    run_flag += GetFlag('FIT_SKY')
                    print "ADDING FLAG FIT_SKY"
            if ParamDict[RunNo][i][1] == 'sersic' and ParamDict[RunNo][i][11] == 'Other':
                FitDict[i][1] = [1, 1]
                FitDict[i][2] = 1 
                FitDict[i][3] = 1 
                FitDict[i][4] = 1
                FitDict[i][5] = 1       
                FitDict[i][6] = 1    

    elif fit_type == 'dev': # Perform a single dev fit 
        for j in range(len(ParamDict[RunNo])):
            i = j + 1
            FitDict[i] = {}  
            if ParamDict[RunNo][i][1] == 'sersic' and ParamDict[RunNo][i][11] == 'Main':
                FitDict[i][1] = [failed_ser, failed_ser] 
                FitDict[i][2] = 1 
                FitDict[i][3] = 1 
                FitDict[i][4] = 0 # fix n
                FitDict[i][5] = 1       
                FitDict[i][6] = 1
                if failed_ser:
                    run_flag += 2**GetFlag('FIT_BULGE_CNTR')
                    print "ADDING FLAG FIT_BULGE_CNTR"
            if ParamDict[RunNo][i][1] == 'expdisk' and ParamDict[RunNo][i][11] == 'Main':
                FitDict[i][1] = [failed_ser, failed_ser]
                FitDict[i][2] = 1 
                FitDict[i][3] = 1 
                FitDict[i][4] = 1       
                FitDict[i][5] = 1  
            if ParamDict[RunNo][i][1] == 'sky':
                FitDict[i][1] = bad_sky
                FitDict[i][2] = 0 
                FitDict[i][3] = 0
                if bad_fit:
                    run_flag += 2**GetFlag('FIT_SKY')
                    print "ADDING FLAG FIT_SKY"
            if ParamDict[RunNo][i][1] == 'sersic' and ParamDict[RunNo][i][11] == 'Other':
                #should this be all 1s or 0s?
                FitDict[i][1] = [1, 1]
                FitDict[i][2] = 1 
                FitDict[i][3] = 1 
                FitDict[i][4] = 1
                FitDict[i][5] = 1       
                FitDict[i][6] = 1    

    elif fit_type in ['devexp', 'serexp']: # Perform a devexp or serexp fit 
        for j in range(len(ParamDict[RunNo])):
            i = j + 1
            FitDict[i] = {}  
            if ParamDict[RunNo][i][1] == 'sersic' and ParamDict[RunNo][i][11] == 'Main':
                FitDict[i][1] = [failed_ser, failed_ser] 
                FitDict[i][2] = 1 
                FitDict[i][3] = 1
                if failed_ser:
                    run_flag += 2**GetFlag('FIT_BULGE_CNTR')
                    print "ADDING FLAG FIT_BULGE_CNTR"
                if fit_type == 'devexp':
                    FitDict[i][4] = 0 # fix n
                else:
                    FitDict[i][4] = 1 # free n
                FitDict[i][5] = 1       
                FitDict[i][6] = 1    
            if ParamDict[RunNo][i][1] == 'expdisk' and ParamDict[RunNo][i][11] == 'Main':
                FitDict[i][1] = [failed_ser, failed_ser]
                FitDict[i][2] = 1 
                FitDict[i][3] = 1 
                FitDict[i][4] = 1       
                FitDict[i][5] = 1
                if failed_ser:
                    run_flag += 2**GetFlag('FIT_DISK_CNTR')
                    print "ADDING FLAG FIT_DISK_CNTR"
            if ParamDict[RunNo][i][1] == 'sky':
                FitDict[i][1] = bad_sky
                FitDict[i][2] = 0 
                FitDict[i][3] = 0
                if bad_fit:
                    run_flag += 2**GetFlag('FIT_SKY')
                    print "ADDING FLAG FIT_SKY"
            if ParamDict[RunNo][i][1] == 'sersic' and ParamDict[RunNo][i][11] == 'Other':
                #should this be all 1s or 0s?
                FitDict[i][1] = [1, 1]
                FitDict[i][2] = 1 
                FitDict[i][3] = 1 
                FitDict[i][4] = 1
                FitDict[i][5] = 1       
                FitDict[i][6] = 1    

    return FitDict, run_flag


def FractionalError(ParamDict, ErrDict, RunNo):
    Xberr = ErrDict[RunNo][1][1][0] / ParamDict[RunNo][1][2][0] 
    Yberr = ErrDict[RunNo][1][1][1] / ParamDict[RunNo][1][2][1] 
    Fb = 10**((c.mag_zero -  ParamDict[RunNo][1][3]) / 2.5)
    Fbe = abs(10**((c.mag_zero -  ParamDict[RunNo][1][3] + ErrDict[RunNo][1][2]) / 2.5) - 10**((c.mag_zero -  ParamDict[RunNo][1][3] - ErrDict[RunNo][1][2]) / 2.5))
    ffbe = Fbe / Fb
    reerr = ErrDict[RunNo][1][3] / ParamDict[RunNo][1][4]
    nerr = ErrDict[RunNo][1][4] / ParamDict[RunNo][1][5]
    Xderr = ErrDict[RunNo][2][1][0] / ParamDict[RunNo][2][2][0] 
    Yderr = ErrDict[RunNo][2][1][1] / ParamDict[RunNo][2][2][1] 
    Fd = 10**((c.mag_zero -  ParamDict[RunNo][2][3]) / 2.5)
    Fde = n.abs(10**((c.mag_zero -  ParamDict[RunNo][2][3] + ErrDict[RunNo][2][2]) / 2.5) - 10**((c.mag_zero -  ParamDict[RunNo][2][3] - ErrDict[RunNo][2][2]) / 2.5))
    ffde = Fde / Fd
    rderr = ErrDict[RunNo][2][3] / ParamDict[RunNo][2][4]
    #        print Xberr, Yberr, Fb, Fbe, reerr, nerr
    #        print Xderr, Yderr, Fd, Fde, rderr
    toterr = n.sqrt(ffbe**2.0 + reerr**2.0 + nerr**2.0 + ffde**2.0 + rderr**2.0)
    if n.median([ffbe, reerr, nerr, ffde, rderr]) < 0.1:
        toterr = 0.01
    else:
        toterr = n.median([ffbe, reerr, nerr, ffde, rderr])
    print 'median ' ,n.median([ffbe, reerr, nerr, ffde, rderr])
    return toterr


def DecideHowToMove2(ParamDict, RunNo,fit_type, KpCArc, failed_ser, failed_sky):
    sersic_loc = c.sersic_loc # index location of the single sersic fit in the ParamDict used for decision-making in the dev+exp and Ser+exp fits
    sky_loc = c.sky_loc
    bad_fit = 0
    bad_sky = 0
    HitLimitCheck = 0

    n_points=n.array([0.5, 3,4,5,6, 7, 8])
    b_points = n.array([0.2, 0.5, 0.5,0.5, 0.5,0.5,0.5])
    fit_out=n.polyfit(n_points, b_points, 3)
    p=n.poly1d(fit_out)

    if c.ParamDictBook[sky_loc][1][5] > 8.0 or c.ParamDictBook[sky_loc][1][4] *KpCArc > 40.0 or failed_sky:
        print "bad sky fit"
        bad_sky = 1
        sky_loc = 0
    else:
        bad_sky = 0
        
    if fit_type != 'sky':
        if abs(c.ParamDictBook[sersic_loc][1][3] - (c.LMag - 1.0)) < 0.05 or abs(c.ParamDictBook[sersic_loc][1][3] - c.UMag) < 0.05 or c.ParamDictBook[sersic_loc][1][4] < 0.21 or  c.ParamDictBook[sersic_loc][1][5] > 8.0 or c.ParamDictBook[sersic_loc][1][4] *KpCArc > 40.0 or failed_ser:
            print "bad sersic fit"
            bad_fit = 1
            
        # note that bad_fit will always be true when failed_ser is true, but bad_fit does not imply failed_ser
        
        if bad_fit:
            look_loc = 0 # look only at SExtractor values  
        else:
            look_loc = sersic_loc #look at Sersic values where reasonable then look at SExtractor

    if fit_type == 'sky': #prepares for the ser fit
        bt_fit = -1
                
        ParamDict[RunNo][1][2][0] = copy.deepcopy(c.ParamDictBook[0][1][2][0])
        ParamDict[RunNo][1][2][1] = copy.deepcopy(c.ParamDictBook[0][1][2][1])
        
        ParamDict[RunNo][1][3] = copy.deepcopy(c.ParamDictBook[0][1][3])
        ParamDict[RunNo][1][4] = copy.deepcopy(c.ParamDictBook[0][1][4])
        ParamDict[RunNo][1][5] = 4.0
        ParamDict[RunNo][1][6] = copy.deepcopy(c.ParamDictBook[0][1][6])
        ParamDict[RunNo][1][7] = copy.deepcopy(c.ParamDictBook[0][1][7])    
        ParamDict[RunNo][2][2][0] = 9999
        ParamDict[RunNo][2][2][1] = 9999
        ParamDict[RunNo][2][3] = 9999
        ParamDict[RunNo][2][4] = 9999
        ParamDict[RunNo][2][5] = 9999
        ParamDict[RunNo][2][6] = 9999


    elif fit_type == 'ser': #prepares for the dev fit
        bt_fit = -1
                
        ParamDict[RunNo][1][2][0] = copy.deepcopy(c.ParamDictBook[look_loc][1][2][0])
        ParamDict[RunNo][1][2][1] = copy.deepcopy(c.ParamDictBook[look_loc][1][2][1])
        
        ParamDict[RunNo][1][3] = copy.deepcopy(c.ParamDictBook[0][1][3])
        ParamDict[RunNo][1][4] = copy.deepcopy(c.ParamDictBook[0][1][4])
        ParamDict[RunNo][1][5] = 4.0
        ParamDict[RunNo][1][6] = copy.deepcopy(c.ParamDictBook[look_loc][1][6])
        ParamDict[RunNo][1][7] = copy.deepcopy(c.ParamDictBook[look_loc][1][7])    
        ParamDict[RunNo][2][2][0] = 9999
        ParamDict[RunNo][2][2][1] = 9999
        ParamDict[RunNo][2][3] = 9999
        ParamDict[RunNo][2][4] = 9999
        ParamDict[RunNo][2][5] = 9999
        ParamDict[RunNo][2][6] = 9999
            
    elif bad_fit:
        bt_fit = .5
            
        ParamDict[RunNo][1][2][0] = copy.deepcopy(c.ParamDictBook[look_loc][1][2][0])
        ParamDict[RunNo][1][2][1] = copy.deepcopy(c.ParamDictBook[look_loc][1][2][1])

        ParamDict[RunNo][1][3] = copy.deepcopy(c.ParamDictBook[0][1][3]) - 2.5 * n.log10(bt_fit)
        ParamDict[RunNo][1][4] = copy.deepcopy(c.ParamDictBook[0][1][4])
        ParamDict[RunNo][1][5] = 4.0
        ParamDict[RunNo][1][6] = copy.deepcopy(c.ParamDictBook[0][1][6])
        ParamDict[RunNo][1][7] = copy.deepcopy(c.ParamDictBook[0][1][7])
        ParamDict[RunNo][2][2][0] = copy.deepcopy(c.ParamDictBook[look_loc][1][2][0])
        ParamDict[RunNo][2][2][1] = copy.deepcopy(c.ParamDictBook[look_loc][1][2][1])
        ParamDict[RunNo][2][3] = copy.deepcopy(c.ParamDictBook[0][1][3])-2.5 * n.log10(1.0 - bt_fit)
        ParamDict[RunNo][2][4] = copy.deepcopy(c.ParamDictBook[0][1][4])
        ParamDict[RunNo][2][5] = copy.deepcopy(c.ParamDictBook[0][1][6])
        ParamDict[RunNo][2][6] = copy.deepcopy(c.ParamDictBook[0][1][7])
    
    elif c.ParamDictBook[sersic_loc][1][5] <= 2.0: # disk-like galaxy -> so set disk parameters to bulge 
        print 'disky'
        bt_fit = p(c.ParamDictBook[sersic_loc][1][5])
        
        ParamDict[RunNo][1][2][0] = copy.deepcopy(c.ParamDictBook[sersic_loc][1][2][0])
        ParamDict[RunNo][1][2][1] = copy.deepcopy(c.ParamDictBook[sersic_loc][1][2][1])
        ParamDict[RunNo][1][3] = copy.deepcopy(c.ParamDictBook[sersic_loc][1][3]) - 2.5 * n.log10(bt_fit)
        ParamDict[RunNo][1][4] = copy.deepcopy(c.ParamDictBook[sersic_loc][1][4])/4.0 #shrink bulge radius
        ParamDict[RunNo][1][5] = 4.0
        ParamDict[RunNo][1][6] = copy.deepcopy(c.ParamDictBook[sersic_loc][1][6])
        ParamDict[RunNo][1][7] = copy.deepcopy(c.ParamDictBook[sersic_loc][1][7])
        ParamDict[RunNo][2][2][0] = copy.deepcopy(c.ParamDictBook[sersic_loc][1][2][0])
        ParamDict[RunNo][2][2][1] = copy.deepcopy(c.ParamDictBook[sersic_loc][1][2][1])
        ParamDict[RunNo][2][3] = copy.deepcopy(c.ParamDictBook[sersic_loc][1][3])-2.5 * n.log10(1.0 - bt_fit)
        ParamDict[RunNo][2][4] = copy.deepcopy(c.ParamDictBook[sersic_loc][1][4])
        ParamDict[RunNo][2][5] = copy.deepcopy(c.ParamDictBook[sersic_loc][1][6])
        ParamDict[RunNo][2][6] = copy.deepcopy(c.ParamDictBook[sersic_loc][1][7])

    else:  
        print 'mix'
        bt_fit = p(c.ParamDictBook[sersic_loc][1][5])
        
        ParamDict[RunNo][1][2][0] = copy.deepcopy(c.ParamDictBook[sersic_loc][1][2][0])
        ParamDict[RunNo][1][2][1] = copy.deepcopy(c.ParamDictBook[sersic_loc][1][2][1])
        ParamDict[RunNo][1][3] = copy.deepcopy(c.ParamDictBook[sersic_loc][1][3]) - 2.5 * n.log10(bt_fit)
        ParamDict[RunNo][1][4] = copy.deepcopy(c.ParamDictBook[sersic_loc][1][4])/2.0 
        ParamDict[RunNo][1][5] = 4.0
        ParamDict[RunNo][1][6] = copy.deepcopy(c.ParamDictBook[sersic_loc][1][6])
        ParamDict[RunNo][1][7] = copy.deepcopy(c.ParamDictBook[sersic_loc][1][7])
        ParamDict[RunNo][2][2][0] = copy.deepcopy(c.ParamDictBook[sersic_loc][1][2][0])
        ParamDict[RunNo][2][2][1] = copy.deepcopy(c.ParamDictBook[sersic_loc][1][2][1])
        ParamDict[RunNo][2][3] = copy.deepcopy(c.ParamDictBook[sersic_loc][1][3])-2.5 * n.log10(1.0 - bt_fit)
        ParamDict[RunNo][2][4] = copy.deepcopy(c.ParamDictBook[sersic_loc][1][4])*2.0
        ParamDict[RunNo][2][5] = copy.deepcopy(c.ParamDictBook[sersic_loc][1][6])
        ParamDict[RunNo][2][6] = copy.deepcopy(c.ParamDictBook[sersic_loc][1][7])


    try:
#            ParamDict[RunNo][c.SkyNo][2] = n.mean([copy.deepcopy(c.ParamDictBook[sky_loc][c.SkyNo][2]), c.SexSky, c.SkyMin])
            ParamDict[RunNo][c.SkyNo][2] = copy.deepcopy(c.ParamDictBook[sky_loc][c.SkyNo][2])
    except:
        print 'EXCEPTION!!!!!!!'
        traceback.print_exc()
        try:
            ParamDict[RunNo][c.SkyNo][2] = copy.deepcopy(c.ParamDictBook[0][c.SkyNo][2])
        except:
            pass


    bad_sky = 0 # this will fix sky which we will calculate as the average of three skys        
    return bad_sky, bad_fit,bt_fit
