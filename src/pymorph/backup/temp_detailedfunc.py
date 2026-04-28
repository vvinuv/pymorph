import os
import numpy as np
import sys
import json
import subprocess
import matplotlib.pyplot as plt
from .configfunc import GalfitConfigFunc


class Detailed:
    def __init__(self):
        pass
        #self.detailed = detailed(self, xcntr_img, ycntr_img,
        #                         good_object, psffile)
        
    def detailed(self, xcntr_img, ycntr_img, good_object, psffile):
        number_random = 50
        objects = {}
        out_info = {}

        magb, rb, n, eb, magd, rd, ed, magbar, rbar, nbar, ebar = uniform_params(bulge_mag, bulge_rad, bulge_n, bar_n, number_random)

        #magb, rb, n, eb, magd, rd, ed, magbar, rbar, nbar, ebar = normal_params(objects['bulge']['mag'], objects['bulge']['rad'], objects['bulge']['n'], objects['bulge']['ell'], objects['disk']['mag'], objects['disk']['rad'], objects['disk']['ell'], objects['bar']['mag'], objects['bar']['n'], objects['bar']['ell'], number_random)
        objects['bulge']['mag'] = [m for m in magb]
        print(objects['bulge']['mag'])
        objects['bulge']['rad'] = [r for r in rb]
        objects['bulge']['n'] = [n1 for n1 in n]

        #sys.exit()
        objects['disk']['mag'] = [m for m in magd]
        objects['disk']['rad'] = [r for r in rd]

        objects['bar']['mag'] = [m for m in magbar]
        objects['bar']['rad'] = [r for r in rbar]
        objects['bar']['n'] = [n1 for n1 in nbar]

        print(objects)
        #sys.exit()
        for i in range(n.shape[0]):


            CF = GalfitConfigFunc(xcntr_img, ycntr_img, good_object, psffile,
                                  detailed=[magb[i], rb[i], n[i], magd[i], rd[i], 
                                            magbar[i], rbar[i], nbar[i]])        

            Gfile = CF.galfit_conf 
            subprocess.call(['/Users/vinu/github/galfit3.0.7b/galfit', f'{Gfile}'])
            if not os.path.exists('fit.log'):
                continue
            basic_info, fit_info, measured_error_bad = read_fitlog(filename = 'fit.log', yes_bar = 1)

            fb = 10**(-0.4 * (fit_info['bulge']['mag'][0] - mag_zero))
            fd = 10**(-0.4 * (fit_info['disk']['mag'][0] - mag_zero))
            fbar = 10**(-0.4 * (fit_info['bar']['mag'][0] - mag_zero))

            fit_info['BT'] = round(fb / (fb + fd + fbar), 2)
            fit_info['BarT'] = round(fbar / (fb + fd + fbar), 2)
             
            out_info[i+1] = fit_info

            print(f'BT {fit_info['BT']} BarT {fit_info['BarT']}')

        json.dump(objects, open(f"in_{fstring}.json", 'w'))
        json.dump(out_info, open(f"out_{fstring}.json", 'w'))


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
    obj['ell'] = float(lns[7].split()[1])
    obj['pa'] = float(lns[8].split()[1])
    return obj

def load_disk(comp, lns):
    obj = {}
    obj['xctr'] = float(lns[0].split()[1])
    obj['yctr'] = float(lns[0].split()[2])
    obj['mag'] = float(lns[1].split()[1])
    obj['rad'] = float(lns[2].split()[1])
    obj['ell'] = float(lns[7].split()[1])
    obj['pa'] = float(lns[8].split()[1])
    return obj


def uniform_params(mag, rad, n=4, nbar=None, number_random=50,
                   delmag=0.5, delr=[0.6, 1.3], deln=[0.5, 1.5],
                   magbar=0.5, rbar=1.2):


    mag = np.random.uniform(mag-delmag, mag+delmag, number_random)
    mag = np.around(mag, 2)

    magbar = np.around(mag-magbar, 2)

    #Radius
    rad = np.random.uniform(rad * delr[0], rad * delr[1], number_random)

    rbar = np.around(rbar * rad , 2)


    #Sersic index
    n = np.random.uniform(n * deln[0], n * deln[1], number_random)
    n = np.around(n, 2)

    nbar = np.where(n > 3, 2.5, n * 0.8)
    nbar = np.around(nbar, 2)

    return mag, rad, n, nbar


def normal_params(bmag, brad, bn, be, dmag, drad, de,
                   barmag, barn, bare, number_random):

    dist_mod = 36.65
    r_pc_init = 10**(-0.279 * (bmag - dist_mod) - 2.457) 
    magb = np.random.normal(bmag, 0.2, number_random)
    magb = np.around(magb, 2)
    magd = np.random.normal(dmag, 0.2, number_random)
    magd = np.around(magd, 2)
    #Total magnitude of bar
    magbar = np.around(np.random.uniform(0.4, 0.6, number_random) + magb, 2)

    #Radius of bulge
    rb = np.random.normal(brad, brad * 0.3, number_random)
    rb_pc = 10**(-0.279 * (magb - dist_mod) - 2.457) #Quilley & de Lapparent, 2023
    rb = rb_pc / r_pc_init * brad
    rb = np.where(rb < rb * 0.2, rb * 0.2, rb)
    rb = np.around(rb, 2)

    rd = np.random.normal(drad, drad * 0.3, number_random)
    rd_pc = 10**(-0.208 * (magd - dist_mod) - 0.434) #Quilley & de Lapparent, 2023
    rd = rd_pc / r_pc_init * drad
    rd = np.where(rd < rd * 0.2, rd * 0.2, rd)
    rd = np.around(rd, 2)

    #Radius of bar
    rbar = rb * np.random.normal(1.1, 1.6, number_random)
    rbar = np.where(rbar < rbar * 0.2, rbar * 0.2, rbar)
    rbar = np.around(rbar, 2)

    n = np.random.normal(bn, bn * 0.2, number_random)
    n[n < 0.5] = 0.5
    n[n > 8.0] = 8.0
    n = np.around(n, 2)

    nbar = np.random.normal(barn, barn * 0.1, number_random)
    nbar[nbar < 0.5] = 0.5
    nbar[nbar > 8.0] = 8.0
    nbar = np.around(nbar, 2)


    eb = np.random.normal(be, be * 0.2, number_random)
    eb[eb > 0.9] = 0.9
    eb[eb < 0.3] = 0.3
    eb = np.around(eb, 2)

    ed = np.random.normal(de, de * 0.2, number_random)
    ed[ed > 0.9] = 0.9
    ed[ed < 0.3] = 0.3
    ed = np.around(ed, 2)

    ebar = np.random.normal(bare, bare * 0.2, number_random)
    
    return magb, rb, n, eb, magd, rd, eb, magbar, rbar, nbar, ebar

#detailed()
#sys.exit()
if 0:
    magb, rb, n, eb, magd, rd, ed, magbar, rbar, nbar, ebar = uniform_params(16, 5, 4, 0.6, 16, 4, 0.6, 17, 1., 0.2, 100)
    #magb, rb, n, eb, magd, rd, ed, magbar, rbar, nbar, ebar = normal_params(16, 5, 4, 0.6, 16, 4, 0.6, 17, 1., 0.2, 100)
    fig, axes = plt.subplots(ncols=4, nrows=3, figsize=(12, 9))
    axes[0,0].hist(magb)
    axes[0,0].set_xlabel('Bulge Mag')
    axes[0,1].hist(rb)
    axes[0,1].set_xlabel('Bulge rad')
    axes[0,2].hist(n)
    axes[0,2].set_xlabel('Bulge Sersic')
    axes[0,3].hist(eb)
    axes[0,3].set_xlabel('Bulge Elli')
    axes[1,0].hist(magd)
    axes[1,0].set_xlabel('Disk Mag')
    axes[1,1].hist(rd)
    axes[1,1].set_xlabel('Disk rad')
    axes[1,2].hist(ed)
    axes[1,2].set_xlabel('Disk elli')
    axes[2,0].hist(magbar)
    axes[2,0].set_xlabel('Bar Mag')
    axes[2,1].hist(rbar)
    axes[2,1].set_xlabel('Bar rad')
    axes[2,2].hist(nbar)
    axes[2,2].set_xlabel('Bar Sersic')
    axes[2,3].hist(ebar)
    axes[2,3].set_xlabel('Bar Elli')
    plt.show()
