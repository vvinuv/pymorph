import os
import configparser
import sys
import subprocess
import json
import numpy as np
#from .flagfunc import GetFlag, isset, SetFlag
import traceback
from .errors_warnings import catch_pipeline_issues
from .galfit_config_run import GalfitConfigRunFunc
from .read_fitlog import read_fitlog

def uniform_params(mag, rad, has_bar, number_random, n=4, 
                   delmag=0.5, delr=[0.6, 1.3], deln=[0.5, 1.5],
                   magbar=0.5, rbar=1.2):


    mag = np.random.uniform(mag-delmag, mag+delmag, number_random)
    mag = np.around(mag, 2)


    #Radius
    rad = np.random.uniform(rad * delr[0], rad * delr[1], number_random)
    rad = np.around(rad, 2)

    
    #Sersic index
    n = np.random.uniform(n * deln[0], n * deln[1], number_random)
    n = np.around(n, 2)

   
    if has_bar:
        magbar = np.around(mag-magbar, 2)
        rbar = np.around(rbar * rad , 2)
        nbar = np.where(n > 3, 2.5, n * 0.8)
        nbar = np.around(nbar, 2)
        
        return mag, rad, n, magbar, rbar, nbar

    else:

        return mag, rad, n

class GalfitDetailed:

    def __init__(self, config_file, galfit_conf, target, neighbours):

        config = configparser.ConfigParser()
        config.read(config_file)

        components = config.get('galfit', 'components').split(',')
        self.components = [cm.strip() for cm in components]
        self.GALFIT_PATH = config.get('external', 'GALFIT_PATH')

        self.bbox = config.getfloat("galfit", "bbox")
        self.barbox = config.getfloat("galfit", "barbox")

        self.target = target
        self.neighbours = neighbours
        self.galfit_conf = galfit_conf


    def detailed(self, number_random=50):

        new = '\n'
        xcntr_img = self.target["X_IMAGE"]
        ycntr_img = self.target["Y_IMAGE"]
        axis_ratio = 1 / self.target["ELONGATION"]
        pos_ang = self.target["GALFIT_ANGLE"]
        SexSky = self.target["BACKGROUND"]
        mag_auto = self.target["MAG_AUTO"]
        flux_radius = self.target["FLUX_RADIUS"]
        fstring = self.target["NAME"]
        print('target["X_IMAGE"]', self.target["X_IMAGE"])
        
        objects = {}
        out_info = {}

        axis_ratio = np.around(np.ones(number_random) * axis_ratio, 
                               2).tolist()
        pos_ang = np.around(np.ones(number_random) * pos_ang, 2).tolist()

        if 'bar' in self.components:
            results = uniform_params(mag_auto, 
                                     flux_radius, 
                                     True,
                                     number_random, n=4) 
            magbar = results[3]
            rbar = results[4]
            nbar = results[5]
            objects['bar'] = {'xcntr': xcntr_img, 'ycntr': ycntr_img,
                            'mag': [m for m in magbar],
                            'rad': [r for r in rbar],
                            'n': [n1 for n1 in nbar],
                            'ell': axis_ratio,
                            'pa': pos_ang}
        else:
            results = uniform_params(mag_auto, 
                                         flux_radius, 
                                         False,
                                         number_random, 
                                         n=4)

        mag = results[0]
        rad = results[1]
        n = results[2]

        objects['bulge'] = {'xcntr': xcntr_img,'ycntr': ycntr_img,
                        'mag': [m for m in mag],
                        'rad': [r for r in rad],
                        'n': [n1 for n1 in n],
                        'ell': axis_ratio,
                        'pa': pos_ang}

        objects['disk'] = {'xcntr': xcntr_img, 'ycntr': ycntr_img,
                        'mag': [m for m in mag],
                        'rad': [r for r in rad],
                        'ell': axis_ratio,
                        'pa': pos_ang}


            
        objects['sky'] = SexSky

        #print(objects)
        #print(len(n))
        #input()

        #sys.exit()
        #gcr.GalfitRun()

        f = open(self.galfit_conf, 'r')
        lines = f.readlines()
        f.close()


        galfit_meta = ''

        for line in lines:
            con = line.startswith('A)') or line.startswith('B)') or line.startswith('C)')
            con = con or line.startswith('D)') or line.startswith('E)') or line.startswith('F)') or line.startswith('G)') or line.startswith('H)') or line.startswith('I)') or line.startswith('J)') or line.startswith('K)') or line.startswith('O)') or line.startswith('P)')
            if con:
                galfit_meta += f'{line}'


        for i in range(number_random):
       

            bx = objects['bulge']['xcntr']
            by = objects['bulge']['ycntr']
            bm = mag[i]
            br = rad[i]
            bn = n[i]
            be = axis_ratio[i] 
            bp = pos_ang[i]

            dx = objects['disk']['xcntr']
            dy = objects['disk']['ycntr']
            dm = mag[i]
            dr = rad[i]
            de = axis_ratio[i]
            dp = pos_ang[i]

            sky = objects['sky']

            gconf = galfit_meta
            gconf += f'{new} 0) sersic {new} 1) {bx} {by} 1 1 {new} 3) {bm} 1 {new} 4) {br} 1 {new} 5) {bn} 1 {new}'
            gconf += f' 6) 0.0 0 {new} 7) 0.0 0 {new} 8) 0.0 0 {new} 9) {be} 1 {new}10) {bp} 1 {new} Z) 0 {new}{new}'
            gconf += f'{new}'

            if self.bbox != -99:
                gconf += f'C0) {self.bbox} {new}'
                gconf += f'{new}'

            gconf += f' 0) expdisk {new} 1) {dx} {dy} 1 1 {new} 3) {dm} 1 {new} 4) {dr} 1 {new} 5) 0 0 {new}'
            gconf += f' 6) 0.0 0 {new} 7) 0.0 0 {new} 8) 0.0 0 {new} 9) {de} 1 {new}10) {dp} 1 {new} Z) 0 {new}{new}'


            if 'bar' in self.components:
                bax = objects['bar']['xcntr']
                bay = objects['bar']['ycntr']
                bam = magbar[i]
                bar = rbar[i]
                ban = nbar[i]
                bae = axis_ratio[i]
                bap = pos_ang[i]

                gconf += f' 0) sersic {new} 1) {bax} {bay} 1 1 {new} 3) {bam} 1 {new} 4) {bar} 1 {new} 5) {ban} 1 {new}'
                gconf += f' 6) 0.0 0 {new} 7) 0.0 0 {new} 8) 0.0 0 {new} 9) {bae} 1 {new}10) {bap} 1 {new} Z) 0 {new}{new}'
                
                if self.barbox != 99:
                    gconf += f'C0) {self.barbox} {new}'
                    gconf += f'{new}'

            gconf += f' 0) sky {new} 1) {sky} 1 {new} 2) 0.0 0 {new} 3) 0.0 0 {new} Z) 0'


            f = open('temp.in', 'w')
            f.write(gconf)
            f.close()

            
            subprocess.run([self.GALFIT_PATH, 'temp.in'], check=True)

            #input()
            #if os.path.exists('fit.log'):
            #    os.remove('fit.log')

            if not os.path.exists('fit.log'):
                continue
            #print(2)
            subprocess.call(['mv', 'fit.log', f'fit_{i}_{fstring}.log'])


            basic_info, fit_info, measured_error_bad = read_fitlog(filename = f'fit_{i}_{fstring}.log', yes_bar = 1)

            fb = 10**(-0.4 * (fit_info['bulge']['mag'][0]))
            fd = 10**(-0.4 * (fit_info['disk']['mag'][0]))
            fbar = 0.
            
            if 'bar' in self.components:
                fbar = 10**(-0.4 * (fit_info['bar']['mag'][0]))
                fit_info['BarT'] = round(fbar / (fb+fd+fbar), 2)
                print(f'BarT {fit_info['BarT']}')

            #print('basic_info', basic_info)
            #print('measured_error_bad', measured_error_bad)
            fit_info['BT'] = round(fb / (fb+fd+fbar), 2)
            fit_info['chi2nu'] = basic_info['chi2nu']
            print(f'BT {fit_info['BT']}') 

            out_info[i+1] = fit_info


        out_info["0"] = objects
        #json.dump(objects, open(f"in_{fstring}.json", 'w'))
        json.dump(out_info, open(f"{fstring}_detailed.json", 'w'))
