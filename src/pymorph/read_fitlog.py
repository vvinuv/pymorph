import os
import numpy as np
import copy

def getline(values):
    #print('values', values)
    line = values.pop(0)
    #print('line', line)
    line = line.split()
    return line 

def getfit(f):  
    values = f.read()
    #print(values)
    values = values.replace('\n\n','\n') # remove unnecessary whitespace
    # now replace any unneeded charaters that might affect the I/O process...
    for bad_char in ['(',')','[',']',',','*', '--']:
        values = values.replace(bad_char,' ')
           
    values = values.split('\n')
    return values
                    
                    
def load_component(data_line, err_line, boxiness=[-999., -999.]):
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

    #print(boxiness)
    # now load data
    obj['xctr'][0]=float(str(data_line[2]))
    obj['yctr'][0]=float(str(data_line[3]))
    obj['mag'][0]=float(data_line[4])
    #print('obj', obj)
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

            obj['boxy'][0]=boxiness[0]
            obj['boxy'][1]=boxiness[1]
    #print('obj', obj)

    # now replace any nan or inf parameters
    #print('obj', obj)
    #for key in obj.keys():
    #    print(obj[key][0])
    #    print('keys', key)
    for key in obj.keys():
        if np.isnan(obj[key][1]):
            obj[key][1] = -9999.99
            is_comp = True
        elif np.isinf(obj[key][1]):
            obj[key][1] = -6666.66
            is_comp = True
    #print(obj)
    return obj, is_comp


def read_fitlog(filename = 'fit.log', yes_bar = 0):
    """ This function will read the fit log and return all the relevant
    information in 2 Dictionaries, 1 with the basic info and one with
    the fit info"""

    measured_error_bad = False

    boxiness = [-999., -999.]
    neighbor = 0
    basic_info = {}
    fit_info = {}
    if os.path.exists(filename):
        f = open(filename,'r')
        values = getfit(f)
        f.close()
        #print('1 values', values)
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
                            val = copy.deepcopy(values)
                            if val[1].split(':')[0].strip() == 'c0':
                                box = float(val[1].split(":")[1].strip())
                                err = float(val[2].split(':')[-1].strip())
                                boxiness = [box, err]
                                #print(f'Bulge {val[1].split(":")[0].strip()} {box} {err}')

                        elif ('bar' not in fit_info and yes_bar):
                            key = 'bar'
                            val = copy.deepcopy(values)
                            if val[1].split(':')[0].strip() == 'c0':
                                box = float(val[1].split(":")[1].strip())
                                err = float(val[2].split(':')[-1].strip())
                                boxiness = [box, err]
                                #print(f'Bar {val[1].split(":")[0].strip()} {box} {err}')
                        else:
                            neighbor +=1
                            key = 'neighbor' + str(neighbor)
                    elif (str(line[0]) == 'expdisk'):
                        key = 'disk'
                    elif(str(line[0]) == 'psf'):
                        key = 'point'
                    elif(str(line[0]) == 'sky'):
                        key = 'sky'

                    #print(len(values))
                    err_line = getline(values)
                    
                    #print(boxiness)          
                    fit_info[key], is_comp= load_component(line, err_line, boxiness=boxiness)
                    boxiness = [-999., -999.]
                    measured_error_bad += is_comp
                    #print('measured_error_bad, is_comp', measured_error_bad, is_comp)
            except:
                pass

    else:
        print("{} File in {} does not exist!!!!".format(filename, os.getcwd()))

    return basic_info, fit_info, measured_error_bad

if __name__=='__main__':
    path = 'morphology/decomposition/Aug30_decomposing/Boxiness/Bar_Sersic/fitted_Jan8_2026'
    #path = 'morphology/decomposition/Aug30_decomposing/Boxiness/Bulge_Sersic'
    #path = 'morphology/decomposition/Aug30_decomposing/Boxiness/Bulge_Bar_Sersic'
    a, b, c = read_fitlog(filename = f'{path}/fit_box_10.log', yes_bar = 1)
    print(b)
