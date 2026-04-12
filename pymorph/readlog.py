from os.path import exists
import csv
import fileinput
from cosmocal import *
#from utilities import WriteDb
import config as c
import numpy as n
from numpy import random
import copy
import time
import traceback

def ReadLog(ParamDict, ErrDict, No, RunNo, detail = False):
    object_err = 1
    Chi2DOF = 9999 # Initialize to prevent UnboundLocalError
    
    if exists('fit.log'):
        ParamDict[RunNo + 1] = copy.deepcopy(ParamDict[RunNo])
        ErrDict[RunNo + 1] = copy.deepcopy(ErrDict[RunNo])
        
        for line in open('fit.log','r'): 
            values = line.split()
            if not values:
                continue
                
            try: 
                if values[0] == 'Input':
                    # Assuming alpha_j and delta_j are defined globally or in config
                    pass 
                
                if values[0] == 'Init.':
                    initial_conf = str(values[4])
                    
                if values[0] == 'Restart':
                    restart_conf = str(values[3])
                    
                if values[0] == 'Chi^2/nu':
                    chi2nu = float(values[2])
                    
                if values[0] == 'sersic':
                    if detail and (No == 2):
                        No += 1
                    
                    # GALFIT format: (x, y) mag re n ar pa
                    ParamDict[RunNo + 1][No][2][0] = float(values[2][1:-1])
                    ParamDict[RunNo + 1][No][2][1] = float(values[3][:-1])
                    
                    if RunNo != -10:
                        ParamDict[RunNo + 1][No][3] = float(values[4])
                        ParamDict[RunNo + 1][No][4] = float(values[5])
                        
                    if RunNo == -10:
                        ParamDict[RunNo + 1][No][5] = 4.0
                    else:
                        ParamDict[RunNo + 1][No][5] = float(values[6])
                        
                    if RunNo != -10:
                        ParamDict[RunNo + 1][No][6] = float(values[7])
                        
                    ParamDict[RunNo + 1][No][7] = float(values[8])
                    # FIXED: This line was causing the SyntaxError
                    ParamDict[RunNo + 1][No][8] = float(values[9])
                    No += 1
                    
                elif values[0] == 'expdisk':
                    ParamDict[RunNo + 1][No][2][0] = float(values[2][1:-1])
                    ParamDict[RunNo + 1][No][2][1] = float(values[3][:-1])
                    ParamDict[RunNo + 1][No][3] = float(values[4])
                    ParamDict[RunNo + 1][No][4] = float(values[5])
                    ParamDict[RunNo + 1][No][5] = float(values[6])
                    ParamDict[RunNo + 1][No][6] = float(values[7])
                    ParamDict[RunNo + 1][No][7] = float(values[8])
                    No += 1
                    
                elif values[0] == 'psf':
                    ParamDict[RunNo + 1][No][2][0] = float(values[2][1:-1])
                    ParamDict[RunNo + 1][No][2][1] = float(values[3][:-1])
                    ParamDict[RunNo + 1][No][3] = float(values[4])
                    No += 1
                    
                elif values[0] == 'sky':
                    for j in range(len(ParamDict[RunNo + 1])):
                        i = j + 1
                        if ParamDict[RunNo + 1][i][1] == 'sky':
                            ParamDict[RunNo + 1][i][2] = float(values[4])
                            
                if values[0] == 'Chi^2':
                    Chi2 = float(values[2][:-1])
                    DOF = float(values[5])
                    Chi2DOF = Chi2 / DOF
                    
                if values[0].startswith('('):
                    # Error parsing logic
                    if str(PreComp) == 'sersic' and object_err < 3:
                        ErrDict[RunNo + 1][1][1][0] = float(values[0][1:-1])
                        ErrDict[RunNo + 1][1][1][1] = float(values[1][0:-1])
                        ErrDict[RunNo + 1][1][2] = float(values[2])
                        ErrDict[RunNo + 1][1][3] = float(values[3])
                        ErrDict[RunNo + 1][1][4] = float(values[4])
                        object_err += 1
                    elif str(PreComp) == 'expdisk' and object_err < 3:
                        ErrDict[RunNo + 1][2][1][0] = float(values[0][1:-1])
                        ErrDict[RunNo + 1][2][1][1] = float(values[1][0:-1])
                        ErrDict[RunNo + 1][2][2] = float(values[2])
                        ErrDict[RunNo + 1][2][3] = float(values[3])
                        object_err += 1
                        
                PreComp = values[0]
            except Exception:
                pass
                
    return ParamDict, ErrDict, Chi2DOF