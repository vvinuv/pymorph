import os
import numpy as np
import sys
import json
import subprocess

def getline(values):
    line = values.pop(0)
    line = line.split()
    return line

def getfit(f):
    values = f.read()
    values = values.replace('\n\n','\n') # remove unnecessary whitespace
    # now replace any unneeded charaters that might affect the I/O process...
    for bad_char in ['(',')','[',']',',','*', '--']:
        values = values.replace(bad_char,' ')

    values = values.split('\n')
    return values


def load_component(data_line, err_line):
    """This function will construct and load dictionaries for a object when passed the already split
    data_line containing fit parameters and err_line containing errors on the fit parameters"""
    is_comp = False # for tracking problems with the err params on galfit 
    # construct object dictionary
    # these are all fitted values, but many may be unused for a particular object type
    obj = {'xctr':[-999.,-999.],
           'yctr':[-999.,-999.],
           'mag':[-999.,-999.],
           'rad':[-999.,-999.],
           'n':[-999.,-999.],
           'ell':[-999.,-999.],
           'pa':[-999.,-999.],
           'boxy':[-999.,-999.]
           }

    # now load data
    obj['xctr'][0]=float(str(data_line[2]))
    obj['yctr'][0]=float(str(data_line[3]))
    obj['mag'][0]=float(data_line[4])

    if data_line[0] == "sky":
        obj['mag'][1]=float(err_line[0])

    else:
        obj['xctr'][1]=float(str(err_line[0]))
        obj['yctr'][1]=float(str(err_line[1]))
        obj['mag'][1]=float(err_line[2])

        if data_line[0] in ['sersic','expdisk']:
            obj['rad'][0]=float(data_line[5])
            obj['rad'][1]=float(err_line[3])

            if data_line[0] == 'sersic':
                obj['n'][0]=float(data_line[6])
                obj['n'][1]=float(err_line[4])
                pos = 7
            else:
                pos = 6


            obj['ell'][0]=float(data_line[pos])
            obj['ell'][1]=float(err_line[pos-2])
            pos +=1

            obj['pa'][0]=float(data_line[pos])
            obj['pa'][1]=float(err_line[pos-2])
            pos +=1

            try:
                obj['boxy'][0]=float(data_line[pos])
                obj['boxy'][1]=float(err_line[pos-2])
            except IndexError: #if Galfit version 3.0 or later, then no boxyness is reported
                obj['boxy'][0]=-999.
                obj['boxy'][1]=-999.

    # now replace any nan or inf parameters
   # for key in obj.keys():
   #     print(key)
    for key in obj.keys():
        if np.isnan(obj[key][1]):
            obj[key][1] = -9999.99
            is_comp = True
        elif np.isinf(obj[key][1]):
            obj[key][1] = -6666.66
            is_comp = True
    return obj, is_comp

#load_component(['expdisk', ':', '87.34', '85.91', '18.88', '7.22', '0.61', '70.49'], ['0.98', '0.59', '0.20', '0.37', '0.04', '6.35'])

def read_fitlog(filename = 'fit.log', yes_bar = 0):
    """ This function will read the fit log and return all the relevant
    information in 2 Dictionaries, 1 with the basic info and one with
    the fit info"""

    measured_error_bad = False

    neighbor = 0
    basic_info = {}
    fit_info = {}
    if os.path.exists(filename):
        f = open(filename,'r')
        values = getfit(f)
        f.close()
        while len(values) > 0:
            line = getline(values)
            #print(line)
            try:
                if(str(line[0]) == 'Input'):
                    basic_info['Input'] = 1
                elif(str(line[0]) == 'Init.'):
                    basic_info['initial_conf'] = str(line[4])
                elif(str(line[0]) == 'Restart'):
                    basic_info['restart_conf'] = str(line[3])
                elif(str(line[0]) == 'Chi^2/nu'):
                    basic_info['chi2nu'] = float(line[2])
                elif(str(line[0]) == 'Chi^2'):
                    basic_info['dof'] = int(line[5])
                    #print('basic_info', basic_info)
                   #sys.exit()
                elif(str(line[0]) == 'Output'):
                    continue

                # for galaxy bulge
                else: # it must be part of the fit...
                    
                    if str(line[0]) == 'sersic':
                        #print(1)
                        if 'bulge' not in fit_info:
                            key = 'bulge'
                        elif ('bar' not in fit_info and yes_bar):
                            key = 'bar'
                        else:
                            neighbor +=1
                            key = 'neighbor' + str(neighbor)
                    elif (str(line[0]) == 'expdisk'):
                        key = 'disk'
                    elif(str(line[0]) == 'psf'):
                        key = 'point'
                    elif(str(line[0]) == 'sky'):
                        key = 'sky'

                    err_line = getline(values)
                    #print(key)
                    #print('err_line', err_line)
                    #print('line', line) 
                    fit_info[key], is_comp= load_component(line, err_line)
                    measured_error_bad += is_comp
                    #print('measured_error_bad, is_comp', measured_error_bad, is_comp)
            except:
                pass

    else:
        print("{} File in {} does not exist!!!!".format(filename, os.getcwd()))

    return basic_info, fit_info, measured_error_bad



def load_sersic(comp, lns):
    obj = {}
    obj['xctr'] = float(lns[0].split()[1])
    obj['yctr'] = float(lns[0].split()[2])
    obj['mag'] = float(lns[1].split()[1])
    obj['rad'] = float(lns[2].split()[1])
    obj['n'] = float(lns[3].split()[1])
    obj['ell'] = float(lns[4].split()[1])
    obj['pa'] = float(lns[5].split()[1])
    return obj

def load_disk(comp, lns):
    obj = {}
    obj['xctr'] = float(lns[0].split()[1])
    obj['yctr'] = float(lns[0].split()[2])
    obj['mag'] = float(lns[1].split()[1])
    obj['rad'] = float(lns[2].split()[1])
    obj['ell'] = float(lns[3].split()[1])
    obj['pa'] = float(lns[4].split()[1])
    return obj

number_random = 1
new = '\n'
mag_zero = 25

for input_file in ['G_bar_629409.in', 'G_bar_12.in', 'G_bar_47.in']:
    fstring = input_file.split('_')[2].split('.')[0]
    f = open(input_file, 'r')
    lines = f.readlines()
    f.close()



    objects = {}
    k = 0
    for i, l in enumerate(lines):
        if '0) sersic' in l and k == 0:
            k += 1
            b = f''.join(lines[i: i+6])
            obj = load_sersic('bulge', lines[i+1: i+7])
            objects['bulge'] = obj
        if '0) sersic' in l and k == 1:
            b = f''.join(lines[i: i+6])
            obj = load_sersic('bar', lines[i+1: i+7])
            objects['bar'] = obj

        if '0) expdisk' in l:
            d = f''.join(lines[i: i+5])
            obj = load_disk('disk', lines[i+1: i+6])
            objects['disk'] = obj

        if '0) sky' in l:
            s = f''.join(lines[i: i+3])
            objects['sky'] = float(lines[i+1].split()[1])
    print(objects)
    #print(b)
    #f = open('b', 'w')
    #f.write(b)
    #f.close()

    magb = np.random.uniform(objects['bulge']['mag']-1, objects['bulge']['mag']+1, number_random)
    magb = np.around(np.concatenate(([objects['bulge']['mag']], magb)), 2)
    magd = np.random.uniform(objects['disk']['mag']-1, objects['disk']['mag']+1, number_random)
    magd = np.around(np.concatenate(([objects['disk']['mag']], magd)), 2)
    magbar = np.random.uniform(1, 0.5, number_random+1) + magb
    magbar = np.around(magbar, 2)
    #magd = np.concatenate(([objects['bar']['mag']], magd))

    rb = np.random.uniform(objects['bulge']['rad'] * 0.4, objects['bulge']['rad'] * 1.5, number_random)
    rb = np.around(np.concatenate(([objects['bulge']['rad']], rb)), 2)
    rd = np.random.uniform(objects['disk']['rad'] * 0.4, objects['disk']['rad'] * 1.5, number_random)
    rd = np.around(np.concatenate(([objects['disk']['rad']], rd)), 2)
    rbar = np.random.uniform(objects['bar']['rad'] * 0.6, objects['bulge']['rad'] * 1.8, number_random)
    rbar = np.around(np.concatenate(([objects['bar']['rad']], rbar)), 2)

    n = np.random.uniform(objects['bulge']['n'] * 0.4, objects['bulge']['n'] * 1.5, number_random)
    n = np.around(np.concatenate(([objects['bar']['n']], n)), 2)
    nbar = np.random.uniform(objects['bar']['n'] * 0.5, objects['bar']['n'] * 1, number_random)
    nbar = np.around(np.concatenate(([objects['bar']['n']], nbar)), 2)

    eb = np.random.uniform(objects['bulge']['ell'] * 0.8, objects['bulge']['ell'] * 1.2, number_random)
    eb[eb > 0.9] = 0.9
    eb[eb < 0.4] = 0.4
    eb = np.around(np.concatenate(([objects['bulge']['ell']], eb)), 2)
    ed = np.random.uniform(objects['disk']['ell'] * 0.8, objects['disk']['ell'] * 1.2, number_random)
    ed[ed > 0.9] = 0.9
    ed[ed < 0.4] = 0.4
    ed = np.around(np.concatenate(([objects['disk']['ell']], ed)), 2)
    ebar = np.random.uniform(objects['bar']['ell'] * 0.6, objects['bar']['ell'] * 1.0, number_random)
    ebar[ebar > 0.5] = 0.5
    ebar[ebar < 0.1] = 0.1
    ebar = np.around(np.concatenate(([objects['bar']['ell']], ebar)), 2)

    objects['bulge']['mag'] = [m for m in magb]
    objects['disk']['mag'] = [m for m in magd]
    objects['bar']['mag'] = [m for m in magbar]
    objects['bulge']['rad'] = [r for r in rb]
    objects['disk']['rad'] = [r for r in rd]
    objects['bar']['rad'] = [r for r in rbar]
    objects['bulge']['n'] = [n1 for n1 in n]
    objects['bar']['n'] = [n1 for n1 in nbar]
    objects['bulge']['ell'] = [e for e in eb]
    objects['disk']['ell'] = [e for e in ed]
    objects['bar']['ell'] = [e for e in ebar]

    print(objects)

    out_info = {}
    for i in range(eb.shape[0]):
        if i == 0:
            subprocess.call(['rm', 'fit.log'])
            subprocess.call(['/home/vinu/github/galfit3.0.7b/galfit', input_file])

        if not os.path.exists('fit.log'):
            continue
        else:
            basic_info, fit_info, measured_error_bad = read_fitlog(filename = 'fit.log', yes_bar = 1)

            fb = 10**(-0.4 * (fit_info['bulge']['mag'][0] - mag_zero))
            fd = 10**(-0.4 * (fit_info['disk']['mag'][0] - mag_zero))
            fbar = 10**(-0.4 * (fit_info['bar']['mag'][0] - mag_zero))

            fit_info['BT'] = round(fb/(fb+fd+fbar), 2)
            fit_info['BarT'] = round(fbar/(fb+fd+fbar), 2)
            
            #print(fit_info)
            out_info[i+1] = fit_info

            f = open(basic_info['initial_conf'], 'r')
            lines = f.readlines()
            f.close()

            new = '\n'
            galfit = ''


            for line in lines:
                #print(line)
                con = line.startswith('A)') or line.startswith('B)') or line.startswith('C)')
                con = con or line.startswith('D)') or line.startswith('E)') or line.startswith('F)') or line.startswith('G)') or line.startswith('H)') or line.startswith('I)') or line.startswith('J)') or line.startswith('K)') or line.startswith('O)') or line.startswith('P)')
                if con:
                    galfit += f'{line}'
            #print(galfit)
            #print(basic_info)
            #print(fit_info['bulge']['xctr'][0])
            #print(measured_error_bad)
            bx = objects['bulge']['xctr']
            by = objects['bulge']['yctr']
            bm = magb[i]
            br = rb[i]
            bn = n[i]
            be = eb[i]
            bp = objects['bulge']['pa']

            dx = objects['disk']['xctr']
            dy = objects['disk']['yctr']
            dm = magd[i]
            dr = rd[i]
            de = ed[i]
            dp = objects['disk']['pa']


            bax = objects['bar']['xctr']
            bay = objects['bar']['yctr']
            bam = magbar[i]
            bar = rbar[i]
            ban = nbar[i]
            bae = ebar[i]
            bap = objects['bar']['pa']

            sm = objects['sky']


            '''
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

            sx = fit_info['sky']['xctr'][0]
            sy = fit_info['sky']['yctr'][0]
            sm = fit_info['sky']['mag'][0]
            '''

            galfit += f'{new} 0) sersic {new} 1) {bx} {by} 1 1 {new} 3) {bm} 1 {new} 4) {br} 1 {new} 5) {bn} 1 {new}'
            galfit += f' 6) 0.0 0 {new} 7) 0.0 0 {new} 8) 0.0 0 {new} 9) {be} 1 {new}10) {bp} 1 {new} Z) 0 {new}{new}'

            galfit += f' 0) expdisk {new} 1) {dx} {dy} 1 1 {new} 3) {dm} 1 {new} 4) {dr} 1 {new} 5) 0 0 {new}'
            galfit += f' 6) 0.0 0 {new} 7) 0.0 0 {new} 8) 0.0 0 {new} 9) {de} 1 {new}10) {dp} 1 {new} Z) 0 {new}{new}'

            galfit += f' 0) sersic {new} 1) {bax} {bay} 1 1 {new} 3) {bam} 1 {new} 4) {bar} 1 {new} 5) {ban} 1 {new}'
            galfit += f' 6) 0.0 0 {new} 7) 0.0 0 {new} 8) 0.0 0 {new} 9) {bae} 1 {new}10) {bap} 1 {new} Z) 0 {new}{new}'

            galfit += f' 0) sky {new} 1) {sm} 1 {new} 2) 0.0 0 {new} 3) 0.0 0 {new} Z) 0'

            #print(galfit)



            f = open('temp.in', 'w')
            f.write(galfit)
            f.close()

            if i > 0:
                subprocess.call(['rm', 'fit.log'])
                subprocess.call(['/home/vinu/github/galfit3.0.7b/galfit', 'temp.in'])

    json.dump(objects, open(f"in_{fstring}.json", 'w'))
    json.dump(out_info, open(f"out_{fstring}.json", 'w'))

