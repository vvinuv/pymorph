#!/usr/bin/env python3
#########################################
# NAME: config_twostep.py
# DESCRIPTION: Iterative fitting process for structural decomposition.
#########################################

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

# Import pymorph-specific modules 
import config as c
from runsexfunc import *
from readlog import ReadLog
from cosmocal import cal
from flagfunc import *
from plotfunc import PlotFunc

try:
    from utilities import WriteDbDetail
except:
    print('No database module found')

class ConfigIter:
    """
    The class making configuration file for GALFIT with iterative steps.
    Transitions through sky, sersic, dev, devexp, and serexp fits.
    """
    def __init__(self, cutimage, whtimage, xcntr, ycntr, NXPTS, NYPTS, line_s, psffile, z):
        self.cutimage = cutimage
        self.line_s = line_s
        self.whtimage = whtimage
        self.xcntr = xcntr
        self.ycntr = ycntr
        self.NXPTS = NXPTS
        self.NYPTS = NYPTS 
        self.psffile = psffile

        try: 
            c.center_constrain = c.center_constrain
        except:
            c.center_constrain = 2.0
            
        # Run the fitting iteration
        self.confiter = confiter(cutimage, whtimage, xcntr, ycntr, NXPTS, NYPTS, line_s, psffile, z)

def confiter(cutimage, whtimage, xcntr, ycntr, NXPTS, NYPTS, line_s, psffile, z):
    # Run SExtractor 
    RunSex(c.datadir + cutimage, c.datadir + whtimage, 'TEMP.SEX.cat', 9999, 9999, 0)
    
    # Define variables
    imagefile = c.imagefile
    sex_cata = c.sex_cata
    threshold = c.threshold
    thresh_area = c.thresh_area
    c.sersic_loc = 2
    c.sky_loc = 1
    
    try:
        ComP = c.components 
    except:
        ComP = ['bulge', 'disk']
    if len(ComP) == 0:
        ComP = ['bulge', 'disk']
        
    values = line_s.split()
    outfile = 'O_' + c.fstring + '.fits'
    mask_file = 'M_' + c.fstring + '.fits'
    config_file = 'G_' + c.fstring + '.in'
    constrain_file = c.fstring + '.con'

    xcntr_o = xcntr
    ycntr_o = ycntr
    mag = float(values[7])
    radius = float(values[9])
    mag_zero = c.mag_zero
    sky = float(values[10])
    pos_ang = float(values[11]) - 90.0
    axis_rat = 1.0 / float(values[12])
    area_o = float(values[13])
    major_axis = float(values[14])

    # Calculate physical scale constraints
    if 0 < z < 10:
        try:
            KpCArc = cal(z, c.H0, c.WM, c.WV, c.pixelscale)[3]
            max_rad = 40.0 / KpCArc
        except:
            KpCArc = 9999.0
            max_rad = c.URe
    else:
        max_rad = c.URe
        KpCArc = 9999.0

    # Reset fitting flags for detailed process
    for f_name in ['FIT_BULGE_CNTR', 'FIT_DISK_CNTR', 'FIT_SKY', 'FIT_BULGE', 'FIT_DISK', 'FIT_POINT']:
        if isset(c.Flag, GetFlag(f_name)):
            c.Flag = ClrFlag(c.Flag, GetFlag(f_name))

    ParamDict = {0: {}}
    ErrDict = {0: {}}
    BlankDict = {0: {}} 
    
    AdComp = 1
    if 'bulge' in ComP:
        ParamDict[0][AdComp] = {'1': 'sersic', '2': [xcntr_o, ycntr_o], '3': mag, '4': radius, 
                                '5': 4.0, '6': axis_rat, '7': pos_ang, '8': 0, '9': 0, '11': 'Main'}
        # Simplified mapping for dict keys based on your original list structure
        ParamDict[0][AdComp] = {1: 'sersic', 2: [xcntr_o, ycntr_o], 3: mag, 4: radius, 5: 4.0, 6: axis_rat, 7: pos_ang, 8: 0, 9: 0, 11: 'Main'}
        BlankDict[0][AdComp] = {1: 'sersic', 2: [xcntr_o, ycntr_o], 3: 9999, 4: 9999, 5: 9999, 6: 9999, 7: 9999, 8: 9999, 9: 9999, 11: 'Main'}
        AdComp += 1

    if 'disk' in ComP:
        ParamDict[0][AdComp] = {1: 'expdisk', 2: [xcntr_o, ycntr_o], 3: mag, 4: radius, 5: axis_rat, 6: pos_ang, 7: 0, 8: 0, 11: 'Main'}
        BlankDict[0][AdComp] = {1: 'expdisk', 2: [xcntr_o, ycntr_o], 3: 9999, 4: 9999, 5: 9999, 6: 9999, 7: 9999, 8: 9999, 11: 'Main'}
        AdComp += 1

    # Neighbour handling
    isneighbor = 0
    if exists('TEMP.SEX.cat'):
        for line_j in open('TEMP.SEX.cat', 'r'):
            try:
                v = line_j.split()
                xc_n, yc_n = float(v[1]), float(v[2])
                mag_n, rad_n = float(v[7]), float(v[9])
                pa_n, ar_n = float(v[11]) - 90.0, 1.0 / float(v[12])
                area_n, maj_n = float(v[13]), float(v[14])
                
                NotFitNeigh = 0
                if abs(xc_n - xcntr_o) > NXPTS / 2.0 + c.avoidme or abs(yc_n - ycntr_o) > NYPTS / 2.0 + c.avoidme:
                    NotFitNeigh = 1
                
                if (abs(xc_n - xcntr_o) <= (major_axis + maj_n) * threshold and 
                    area_n >= thresh_area * area_o and NotFitNeigh == 0):
                    
                    xn = xcntr + (xc_n - xcntr_o)
                    yn = ycntr + (yc_n - ycntr_o)
                    
                    ParamDict[0][AdComp] = {1: 'sersic', 2: [xn, yn], 3: mag_n, 4: rad_n, 5: 4.0, 6: ar_n, 7: pa_n, 8: 0, 9: 0, 11: 'Other'}
                    BlankDict[0][AdComp] = copy.deepcopy(ParamDict[0][AdComp])
                    isneighbor += 1
                    AdComp += 1
            except:
                continue

    if isneighbor > 0:
        c.Flag = SetFlag(c.Flag, GetFlag('NEIGHBOUR_FIT'))

    # Sky Component
    ParamDict[0][AdComp] = {1: 'sky', 2: c.SexSky, 3: 0, 4: 0, 5: 0, 11: 'Other'}
    BlankDict[0][AdComp] = copy.deepcopy(ParamDict[0][AdComp])
    c.SkyNo = AdComp

    # Iteration storage
    c.Chi2DOFArr, c.FitArr, c.RadArr, c.CntrDevArr, c.ErrArr = [], [], [], [], []
    c.ParamDictBook = copy.deepcopy(ParamDict)
    bad_fit, failed_ser, bt_fit, failed_sky, bad_sky = 0, 0, -1, 0, 0
    
    # MAIN ITERATION LOOP
    for RunNo, fit_type in zip(range(5), ['sky', 'ser', 'dev', 'devexp', 'serexp']):
        run_flag = c.Flag
        if exists('fit.log'): os.system('rm fit.log')

        print(f"--- Run {RunNo}: {fit_type} ---")
        
        # Write config header
        with open(config_file, 'w') as f:
            f.write('# IMAGE PARAMETERS\nA) ' + c.datadir + str(cutimage) + '\nB) ' + str(outfile) + 
                    '\nC) ' + c.datadir + str(whtimage) + '\nD) ' + c.datadir + str(psffile) + 
                    '\nE) 1\nF) ' + str(mask_file) + '\nG) ' + str(constrain_file) + 
                    '\nH) 1 ' + str(NXPTS) + ' 1 ' + str(NYPTS) + '\nI) 120 120\nJ) ' + str(mag_zero) + 
                    '\nO) regular\nP) 0\nS) 0\n\n')

        # Decide fitting parameters for this step
        FitDict, run_flag = DecideFitting(ParamDict, RunNo, fit_type, bad_fit, failed_ser, run_flag, failed_sky, bad_sky)            

        # Write constraints and functions
        with open(constrain_file, 'w') as fc: pass # reset constraints
        
        for i in range(len(ParamDict[RunNo])):
            comp_idx = i + 1
            comp_type = ParamDict[RunNo][comp_idx][1]
            
            if comp_type == 'sersic':
                SersicFunc(config_file, ParamDict, FitDict, comp_idx, RunNo)
                if ParamDict[RunNo][comp_idx][11] == 'Main':
                    run_flag = SetFlag(run_flag, GetFlag('FIT_BULGE'))
                    re_limit = ParamDict[c.sersic_loc][1][4]*2 if (RunNo > 1 and not bad_fit) else max_rad
                    SersicMainConstrain(constrain_file, comp_idx, c.center_constrain, re_limit)
                else:
                    SersicConstrain(constrain_file, comp_idx)
 
            elif comp_type == 'expdisk':
                if RunNo not in [1, 2]: # Don't fit disk in pure ser/dev runs
                    run_flag = SetFlag(run_flag, GetFlag('FIT_DISK'))
                    ExpFunc(config_file, ParamDict, FitDict, comp_idx, RunNo)
                    ExpdiskConstrain(constrain_file, comp_idx, 2.0 if (RunNo > 0 and not bad_fit) else c.center_constrain)

            elif comp_type == 'sky':
                SkyConstrain(constrain_file, comp_idx, ParamDict[c.sky_loc][comp_idx][2] if (RunNo > 0 and not bad_sky) else c.SexSky)
                SkyFunc(config_file, ParamDict, FitDict, comp_idx, RunNo) 

        if RunNo > 2 or RunNo == 0:
            bt_range = add_constrain(constrain_file, bt_fit, fit_type)

        # EXECUTE GALFIT
        os.system(f"{c.GALFIT_PATH} {config_file}")

        try:
            ParamDict, ErrDict, Chi2DOF = ReadLog(ParamDict, ErrDict, 1, RunNo, detail=True)
            c.ParamDictBook[RunNo+1] = copy.deepcopy(ParamDict[RunNo+1])
            
            # Logic for limits and flagging (Simplified for brevity)
            # ... [Limit check code goes here] ...
            
        except:
            print(f"GALFIT failed for {fit_type}")
            ParamDict[RunNo + 1] = copy.deepcopy(BlankDict[0])
            if fit_type == 'ser': failed_ser = 1
            elif fit_type == 'sky': failed_sky = 1
            run_flag = SetFlag(run_flag, GetFlag('DETAIL_FAILED'))

        # Background logic for next step
        if fit_type != 'serexp':
            bad_sky, bad_fit, bt_fit = DecideHowToMove2(ParamDict, RunNo + 1, fit_type, KpCArc, failed_ser, failed_sky)
            # Plotting and DB logging calls...

    return

# Helper functions (SersicMainConstrain, ExpdiskConstrain, etc.) follow the same logic as your previous files.