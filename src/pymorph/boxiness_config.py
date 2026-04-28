import os
import numpy as np
import sys
import json
import subprocess
import matplotlib.pyplot as plt
from astropy.io import fits
import glob
from .read_fitlog import read_fitlog



def detailed(path, input_file, reverse=False):
    new = '\n'
    path = '/Users/vinu/WorkAll/TNG/morphology/decomposition/Aug30_decomposing/default_galfit_conf/'
    #for input_file in glob.glob(path + 'G_bar_*.in'):#[path+'G_bar_109.in']:#, 'G_bar_17.in', 'G_bar_19.in']: #['G_bar_143884.in', 'G_bar_24.in']:#['G_bar_547546.in']:#'G_bar_629409.in']:#, 'G_bar_12.in', 'G_bar_47.in']: 'G_bar_10.in', 'G_bar_12.in', 'G_bar_143884.in', 
        #os.system('/bin/rm -f galfit.*')
        fstring = input_file.split('.')[0].split('G_')[1]        
        print(fstring)
        
        f = open(input_file, 'r')
        lines = f.readlines()
        f.close()

        galfit = ''

        for line in lines:
            con = line.startswith('A)') or line.startswith('B)') or line.startswith('C)')
            con = con or line.startswith('D)') or line.startswith('E)') or line.startswith('F)') or line.startswith('G)') or line.startswith('H)') or line.startswith('I)') or line.startswith('J)') or line.startswith('K)') or line.startswith('O)') or line.startswith('P)')
            if con:
                galfit += f'{line}'

        fit_path = '/Users/vinu/WorkAll/TNG/morphology/decomposition/Aug30_decomposing/fitted_Dec_2025'
        print(f'{fit_path}/fit_{fstring}.log')
        basic_info, fit_info, measured_error_bad = read_fitlog(filename = f'{fit_path}/fit_{fstring}.log', yes_bar = 1)
        
        #print(fit_info, basic_info)     
        bx = fit_info['bulge']['xctr'][0]
        by = fit_info['bulge']['yctr'][0]
        bm = fit_info['bulge']['mag'][0]
        br = fit_info['bulge']['rad'][0]
        bn = fit_info['bulge']['n'][0]
        be = fit_info['bulge']['ell'][0]
        bp = fit_info['bulge']['pa'][0]


        dx = fit_info['disk']['xctr'][0]
        dy = fit_info['disk']['yctr'][0]
        dm = fit_info['disk']['mag'][0]
        dr = fit_info['disk']['rad'][0]
        de = fit_info['disk']['ell'][0]
        dp = fit_info['disk']['pa'][0]

        bax = fit_info['bar']['xctr'][0]
        bay = fit_info['bar']['yctr'][0]
        bam = fit_info['bar']['mag'][0]
        bar = fit_info['bar']['rad'][0]
        ban = fit_info['bar']['n'][0]
        bae = fit_info['bar']['ell'][0]
        bap = fit_info['bar']['pa'][0]

        sky = fit_info['sky']['mag'][0]
        box = 0.3

        galfit += f'{new} 0) sersic {new} 1) {bx} {by} 1 1 {new} 3) {bm} 1 {new} 4) {br} 1 {new} 5) {bn} 1 {new}'
        galfit += f' 6) 0.0 0 {new} 7) 0.0 0 {new} 8) 0.0 0 {new} 9) {be} 1 {new}10) {bp} 1 {new} Z) 0 {new}{new}'

        galfit += f' 0) expdisk {new} 1) {dx} {dy} 1 1 {new} 3) {dm} 1 {new} 4) {dr} 1 {new} 5) 0 0 {new}'
        galfit += f' 6) 0.0 0 {new} 7) 0.0 0 {new} 8) 0.0 0 {new} 9) {de} 1 {new}10) {dp} 1 {new} Z) 0 {new}{new}'

        galfit += f' 0) sersic {new} 1) {bax} {bay} 1 1 {new} 3) {bam} 1 {new} 4) {bar} 1 {new} 5) {ban} 1 {new}'
        galfit += f' 6) 0.0 0 {new} 7) 0.0 0 {new} 8) 0.0 0 {new} 9) {bae} 1 {new}10) {bap} 1 {new} Z) 0 {new}'
        
        box_galfit = galfit + f'C0) {box} {new}{new}'

        galfit += f' 0) sky {new} 1) {sky} 1 {new} 2) 0.0 0 {new} 3) 0.0 0 {new} Z) 0'
        box_galfit += f' 0) sky {new} 1) {sky} 1 {new} 2) 0.0 0 {new} 3) 0.0 0 {new} Z) 0'
        #print(galfit)

        if not os.path.exists('G_bar_{fstring}.in'):
            f = open(f'G_bar_{fstring}.in', 'w')
            #f.write(galfit)
            f.write(box_galfit)
            f.close()

        f = open('temp.in', 'w')
        f.write(box_galfit)
        f.close()

        #sys.exit()
        #input('check constraint')

        
        if 0:
            subprocess.call(['/usr/local/bin/galfit', 'temp.in'])
            #input()
            if not os.path.exists('fit.log'):
                continue
            #print(2)
            subprocess.call(['mv', 'fit.log', f'fit_box_{fstring}.log'])
            #input()
            basic_info, fit_info, measured_error_bad = read_fitlog(filename = f'fit_box_{fstring}.log', yes_bar = 1)

            fb = 10**(-0.4 * (fit_info['bulge']['mag'][0] - mag_zero))
            fd = 10**(-0.4 * (fit_info['disk']['mag'][0] - mag_zero))
            fbar = 10**(-0.4 * (fit_info['bar']['mag'][0] - mag_zero))

            fit_info['BT'] = round(fb/(fb+fd+fbar), 2)
            fit_info['BarT'] = round(fbar/(fb+fd+fbar), 2)
             
            #out_info[i+1] = fit_info

            print(f'BT {fit_info['BT']} BarT {fit_info['BarT']}')
        #input()

    #json.dump(objects, open(f"in_{fstring}.json", 'w'))
    #json.dump(out_info, open(f"out_{fstring}.json", 'w'))
if 1:
    detailed(reverse=True)
sys.exit()

