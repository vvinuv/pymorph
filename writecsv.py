import os
import csv
import datetime
import numpy as np
import fileinput
from .cosmocal import CosmoCal 
import traceback
from .flagfunc import *
from .pymorphutils import RaDegToHMS, DecDegToDMS, output_params
try:
    from .writedbfunc import WriteDB
except:
    print('No mysql database or python mysql.connector module')

class WriteCSV(object):
    """The class which will write html and csv output. This class will also 
       check whether the fit is good or bad using the Chisq and Goodness value
       It will also notify the goodness/badness of fit"""

    def __init__(self):

        self.chi2sq = 1.0
        self.pixelscale = pixelscale

    def writeparams(self, target):


        params = [
                'Name','ra','dec','z','MorphType','mag_auto',
                'magerr_auto','SexHalfRad','num_targets','C',
                'C_err','A','A_err','S','S_err','G','M','magzp',
                'bulge_xctr','bulge_xctr_err','bulge_yctr',
                'bulge_yctr_err','Ie','Ie_err','AvgIe',
                'AvgIe_err','re_pix','re_pix_err','re_kpc',
                're_kpc_err','n','n_err','eb','eb_err','bpa',
                'bpa_err','bboxy','bboxy_err','disk_xctr',
                'disk_xctr_err','disk_yctr','disk_yctr_err',
                'Id','Id_err','rd_pix','rd_pix_err','rd_kpc',
                'rd_kpc_err','ed','ed_err','dpa','dpa_err',
                'dboxy','dboxy_err','BD','BT','BarT','p_xctr',
                'p_xctr_err','p_yctr','p_yctr_err','Ip',
                'Ip_err','Pfwhm','Pfwhm_kpc','bar_xctr',
                'bar_xctr_err','bar_yctr','bar_yctr_err',
                'Ibar','Ibar_err','rbar_pix','rbar_pix_err',
                'rbar_kpc','rbar_kpc_err','n_bar','n_bar_err',
                'ebar','ebar_err','barpa','barpa_err',
                'barboxy','barboxy_err','chi2nu','dof',
                'goodness_upper','goodness_lower','run',
                'SexSky','YetSky','GalSky','GalSky_err',
                'dis_modu','distance','fit','FitFlag','flag',
                'Manual_flag','Comments','Date','Version',
                'Filter','Total_Run','rootname'
                ]

        x = datetime.date.today()
        target['Date'] = f'{x.year}-{x.month}-{x.day}' 
        target['Name'] = f'{target['rootname'}_{target['gal_id']}'
        target['z'] = z
        target['num_targets'] = -999
        target['flag'] = self.flag


        # Reading fit.log
        if target['YES_BAR']:
            basic_info, fit_info, measured_error_bad = read_fitlog(filename='fit.log', 
                                                      yes_bar=1)
        else:
            basic_info, fit_info, measured_error_bad = read_fitlog(filename='fit.log', 
                                                      yes_bar=0)

        #print('basic_info', basic_info)
        if basic_info['dof'] != -999:
            try:
                from scipy import stats
                #sf=1-cdf which gives the probability of large chi2 given dof
                chisqprob_upper = lambda chisq, df: stats.chi2.sf(chisq, df)
                chi2 = basic_info['chi2nu'] * basic_info['dof']
                goodness_upper = chisqprob_upper(chi2, basic_info['dof']) 
                #cdf which gives the probability of sum upto a chi2 given dof
                goodness_lower = stats.chi2.cdf(chi2, basic_info['dof']) 
                goodness_upper = float(f'{goodness_upper:.2E}')
                goodness_lower = float(f'{goodness_lower:.2E}')
            except:
                goodness_upper = -999
                goodness_lower = -999
        else:
            goodness_upper = -999
            goodness_lower = -999

        target['goodness_upper'] = goodness_upper
        target['goodness_lower'] = goodness_lower
        #print('goodness_upper', goodness_upper)
        #print('goodness_lower', goodness_lower)

        #print('fit_info', fit_info)
        #print('measured_error_bad', measured_error_bad, self.flag)
        if measured_error_bad:
            try:
                #set the flag
                target['flag'] = SetFlag(target['flag'], 
                                             GetFlag('ERRORS_FAILED'))
                target['flag'] = SetFlag(target['flag'], 
                                             GetFlag('GALFIT_FAIL'))

                #print(2, target['flag'])
            except:# badflag:
                # the flag is already set
                pass
            
 
        if self.decompose:
            initial_conf = basic_info['initial_conf']
            restart_conf = basic_info['restart_conf']
            #print(restart_conf)
            # move the restart file to a reasonably named output file
            new_outname = initial_conf.replace('in','out')
            try:
                os.rename(restart_conf, new_outname)
            except:
                print("Failed to find restart file!! Galfit may have crashed!!")
            basic_info['restart_conf'] = new_outname

            target['chi2nu'] = basic_info['chi2nu']
            target['dof'] = basic_info['dof']
            # first check all err components and replace if nan or inf

            if 'bulge' in fit_info:
                target['bulge_xctr'] = fit_info['bulge']['xctr'][0]
                target['bulge_yctr'] = fit_info['bulge']['yctr'][0]
                target['Ie'] = fit_info['bulge']['mag'][0]
                target['re_pix'] = fit_info['bulge']['rad'][0]
                target['n'] = fit_info['bulge']['n'][0]
                target['eb'] = fit_info['bulge']['ell'][0]
                target['bpa'] = fit_info['bulge']['pa'][0]
                target['bboxy'] = fit_info['bulge']['boxy'][0]

                target['bulge_xctr_err'] = fit_info['bulge']['xctr'][1]
                target['bulge_yctr_err'] = fit_info['bulge']['yctr'][1]
                target['Ie_err'] = fit_info['bulge']['mag'][1]
                target['re_pix_err'] = fit_info['bulge']['rad'][1]
                target['n_err'] = fit_info['bulge']['n'][1]
                target['eb_err'] = fit_info['bulge']['ell'][1]
                target['bpa_err'] = fit_info['bulge']['pa'][1]
                target['bboxy_err'] = fit_info['bulge']['boxy'][1]

            if 'disk' in fit_info:
                target['disk_xctr'] = fit_info['disk']['xctr'][0]
                target['disk_yctr'] = fit_info['disk']['yctr'][0]
                target['Id'] = fit_info['disk']['mag'][0]
                target['rd_pix'] = fit_info['disk']['rad'][0]
                target['ed'] = fit_info['disk']['ell'][0]
                target['dpa'] = fit_info['disk']['pa'][0]
                target['dboxy'] = fit_info['disk']['boxy'][0]

                target['disk_xctr_err'] = fit_info['disk']['xctr'][1]
                target['disk_yctr_err'] = fit_info['disk']['yctr'][1]
                target['Id_err'] = fit_info['disk']['mag'][1]
                target['rd_pix_err'] = fit_info['disk']['rad'][1]
                target['ed_err'] = fit_info['disk']['ell'][1]
                target['dpa_err'] = fit_info['disk']['pa'][1]
                target['dboxy_err'] = fit_info['disk']['boxy'][1]

            if 'point' in fit_info:
                target['p_xctr'] = fit_info['point']['xctr'][0]
                target['p_yctr'] = fit_info['point']['yctr'][0]
                target['Ip'] = fit_info['point']['mag'][0]

                target['p_xctr_err'] = fit_info['point']['xctr'][1]
                target['p_yctr_err'] = fit_info['point']['yctr'][1]
                target['Ip_err'] = fit_info['point']['mag'][1]

            if 'bar' in fit_info:
                target['bar_xctr'] = fit_info['bar']['xctr'][0]
                target['bar_yctr'] = fit_info['bar']['yctr'][0]
                target['Ibar'] = fit_info['bar']['mag'][0]
                target['rbar_pix'] = fit_info['bar']['rad'][0]
                target['n_bar'] = fit_info['bar']['n'][0]
                target['ebar'] = fit_info['bar']['ell'][0]
                target['barpa'] = fit_info['bar']['pa'][0]
                target['barboxy'] = fit_info['bar']['boxy'][0]

                target['bar_xctr_err'] = fit_info['bar']['xctr'][1]
                target['bar_yctr_err'] = fit_info['bar']['yctr'][1]
                target['Ibar_err'] = fit_info['bar']['mag'][1]
                target['rbar_pix_err'] = fit_info['bar']['rad'][1]
                target['nbar_err'] = fit_info['bar']['n'][1]
                target['ebar_err'] = fit_info['bar']['ell'][1]
                target['barpa_err'] = fit_info['bar']['pa'][1]
                target['barboxy_err'] = fit_info['bar']['boxy'][1]

            if 'sky' in fit_info:
                target['GalSky'] = fit_info['sky']['mag'][0]
                target['GalSky_err'] = fit_info['sky']['mag'][1]

            # Converting fitted params to physical params
            if(z != -999 and z > 0):
                phy_parms = CosmoCal(z, self.H0, self.WM, 
                                    self.WV, self.pixelscale).cal()

                target['dis_modu'] = round(phy_parms[2], 2)

                if 'bulge' in comp:
                    if target['re_pix'] != -999:
                        re_kpc = phy_parms[3] * target['re_pix']
                        target['re_kpc'] = round(re_kpc, 2)
                        re_kpc_err = phy_parms[3] * target['re_pix_err']
                        target['re_kpc_err'] = round(re_kpc_err, 2)

                if 'disk' in comp:
                    if target['rd_pix'] != -999:
                        rd_kpc = phy_parms[3] * target['rd_pix']
                        target['rd_kpc'] = round(rd_kpc, 2)
                        rd_kpc_err = phy_parms[3] * target['rd_pix_err']
                        target['rd_kpc_err'] = round(rd_kpc_err, 2)

                if 'bar' in comp:
                    if target['rbar_pix'] != -999:
                        rbar_kpc = phy_parms[3] * target['rbar_pix']
                        target['rbar_kpc'] = round(rbar_kpc, 2)
                        rbar_kpc_err = phy_parms[3] * target['rbar_pix_err']
                        target['rbar_kpc_err'] = round(rbar_kpc_err, 2)

                if 'point' in comp:
                    target['Pfwhm_kpc'] = 0.5 * phy_parms[3]

            # Finding derived parameters
            if 'bulge' in comp and 'disk' in comp:
                print(target['Ie'] , target['Id'], target['Ibar'])
                fb = 10**(-0.4 * (target['Ie'] - self.mag_zero))
                fd = 10**(-0.4 * (target['Id'] - self.mag_zero))

                if 'point' in comp:
                    fp = 10**(-0.4 * (target['Ip'] - self.mag_zero))
                else:
                    fp = 0.0

                if target['YES_BAR']:
                    print(target['Ibar'])
                    fbar = 10**(-0.4 * (target['Ibar'] - self.mag_zero)) 
                else:
                    fbar = 0.0

                try:

                    target['BD'] = round(fb / fd, 2)
                    target['BT'] = round(fb / (fb + fd + fp + fbar), 2)
                    target['BarT'] = round(fbar / (fb + fd + fp + fbar), 2)

                except:

                    target['BD'] = -999
                    target['BT'] = -999
                    target['BarT'] = -999

            elif 'bulge' in comp:
                target['BT'] = 1.0
                target['BD'] = 0.0
                target['BarT'] = 0.0

            elif 'disk' in comp:
                target['BD'] = 1.0
                target['BT'] = 0.0
                target['BarT'] = 0.0

            if 'bulge' in comp:
                try:
                      pixelscale = self.pixelscale
                except:
                      pixelscale = 1

                try:

                    AvgIe1 = target['Ie']
                    AvgIe2 = 2 * 3.14 * pixelscale * pixelscale
                    AvgIe3 = target['re_pix'] *  target['re_pix']
                    AvgIe4 = np.sqrt(1 -  target['eb']**2.0)
                    AvgIe = AvgIe1 + 2.5 * np.log10(AvgIe2 * AvgIe3 * AvgIe4)
                    AvgIe = round(AvgIe, 2)
                    target['AvgIe'] = AvgIe
                    target['AvgIe_err'] = -999

                except OverflowError:

                    target['AvgIe'] = np.inf
                    target['AvgIe_err'] = np.inf

                    target['flag'] = SetFlag(target['flag'], 
                                             GetFlag('AVGIE_FAILED'))

            # Now test for fitting problems and set flags for analysis
            target['FitFlag'] = 0
            HitLimit = 0

            chi2nu = target['chi2nu']
            if not self.detail:
                if 'bulge' in comp:
                    con1 = abs(target['Ie'] - self.UMag) < 0.2
                    con2 = abs(target['Ie'] - self.LMag) < 0.2
                    if con1 or con2:
                        target['FitFlag'] = SetFlag(target['FitFlag'], 
                                                    Get_FitFlag('IE_AT_LIMIT'))
                    
                    con1 = abs(target['re_pix'] - self.LRe) < 0.1
                    con2 = abs(target['re_pix'] - self.URe) < 1.0
                    if con1 or con2:
                        target['FitFlag'] = SetFlag(target['FitFlag'],
                                                    Get_FitFlag('RE_AT_LIMIT'))

                    con1 = abs(target['n'] - self.LN) < 0.03
                    con2 = abs(target['n'] - self.UN) < 0.5
                    if con1 or con2:
                        target['FitFlag'] = SetFlag(target['FitFlag'],
                                                    Get_FitFlag('N_AT_LIMIT'))

                    con1 = 
                    if abs(target['eb'] - 0.0) < 0.05 or abs(target['eb'] - 0.0) > 0.95:
                        target['FitFlag'] = SetFlag(target['FitFlag'],Get_FitFlag('EB_AT_LIMIT'))
                if 'disk' in comp:
                    if abs(all_params['Id'] - self.UMag) < 0.2 or abs(all_params['Id'] - self.LMag) < 0.2:
                        all_params['FitFlag'] = SetFlag(all_params['FitFlag'],Get_FitFlag('ID_AT_LIMIT'))
                    if abs(all_params['rd_pix'] - self.LRd) < 0.1 or abs(all_params['rd_pix'] - self.URd) < 1.0:
                        all_params['FitFlag'] = SetFlag(all_params['FitFlag'],Get_FitFlag('RD_AT_LIMIT'))
                    if abs(all_params['ed'] - 0.0) < 0.05 or abs(all_params['ed'] - 0.0) > 0.95:
                        all_params['FitFlag'] = SetFlag(all_params['FitFlag'],Get_FitFlag('ED_AT_LIMIT'))

            if not all_params['FitFlag']:
                error_mesg4 = str(error_mesg4) + 'One of the parameters'
                error_mesg5 = str(error_mesg5) + '          hits limit!'

            sigma_percentage = {1:0.682689, 2:0.954499, 3:0.997300, 4:0.999936,
                                5:0.999999}
            sigma_percentage = {1: 0.317311, 2: 0.0455, 3: 0.0027, 4: 6.399e-5, 5: 1e-6}

            print('sigma_percentage', sigma_percentage[self.pymorph_config['goodness_limit']])
            #3sigma gaussian has 0.997 probability. Then probabilty 0.003 shows than the model is not agreeing with the data. This means it is significant  
            if all_params['goodness_upper'] < sigma_percentage[self.pymorph_config['goodness_limit']]: 
                error_mesg2 = str(error_mesg2) + 'Upper Goodness is poor!'
                all_params['FitFlag'] = SetFlag(all_params['FitFlag'], Get_FitFlag('SMALL_UPPER_GOODNESS'))
                error_mesg1 = str(error_mesg1) + 'Chi2nu is large!'

            
            if all_params['goodness_lower'] < sigma_percentage[self.pymorph_config['goodness_limit']]: 
                error_mesg2 = str(error_mesg2) + 'Lower Goodness is poor!'
                all_params['FitFlag'] = SetFlag(all_params['FitFlag'], Get_FitFlag('SMALL_LOWER_GOODNESS'))
                error_mesg1 = str(error_mesg1) + 'Chi2nu is low!'

            if all_params['chi2nu'] > self.pymorph_config['chi2nu_limit']:
                all_params['FitFlag'] = SetFlag(all_params['FitFlag'],Get_FitFlag('LARGE_CHISQ'))

            if abs(all_params['bulge_xctr'] - self.xcntr) > self.center_deviation_limit or \
                   abs(all_params['bulge_yctr'] - self.ycntr) > self.center_deviation_limit or \
                   abs(all_params['disk_xctr'] - self.xcntr) > self.center_deviation_limit or \
                   abs(all_params['disk_yctr'] - self.ycntr) > self.center_deviation_limit:
                if all_params['bulge_xctr'] == -999 or all_params['disk_yctr'] == -999:
                    pass
                else:
                    all_params['FitFlag'] = SetFlag(all_params['FitFlag'], Get_FitFlag('FAKE_CNTR'))
                    error_mesg3 = str(error_mesg3) + 'Fake Center!'

            if all_params['FitFlag'] > 0:
                img_notify = str(self.PYMORPH_PATH) + '/html/goodfit.gif'
                good_fit = 1
            else:
                img_notify = str(self.PYMORPH_PATH) + '/html/badfit.gif'
                good_fit = 0

            # This can be a waste of time if the list is wrong...
            # Finding number of runs in the csv file 
            all_params['run'] = 1
            if os.path.exists(self.final_result_file):
                for line_res in csv.reader(open(self.final_result_file).readlines()[1:]):    
                    if(str(line_res[0]) == self.fstring):
                        all_params['run'] += 1
            #if self.GalSky != 9999:
            #    all_params['GalSky'] = self.GalSky

        #print('params_to_write', params_to_write)
        #print('all_params', all_params)

        #for key in params_to_write.keys():
        #    print(key, all_params[params_to_write[key][0]])

        #Final flag and FitFlag
        self.flag = all_params['flag']
        self.fit_flag = all_params['FitFlag']
        # Writing csv file 
        #print('self.final_result_file', self.final_result_file)
        f_res = open(self.final_result_file, "a")
        writer = csv.writer(f_res)
        writer.writerow([all_params[params_to_write[key][0]] for key in params_to_write.keys()])
        f_res.close()
        #sys.exit()
        # Writing data base
        #print('params_to_write', params_to_write)
        #print('all_params', all_params)
        try:
            WDB = WriteDB()
            WDB._first_db(params_to_write, all_params, self.pymorph_config)
        except Exception as e:
            print('No database can be created!')
            print(e)

        
        



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
                                print(f'Bulge {val[1].split(":")[0].strip()} {box} {err}')

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
            except Exception as e:
                print('Error occured reading fit.log', e)

    else:
        print("{} File in {} does not exist!!!!".format(filename, os.getcwd()))

    return basic_info, fit_info, measured_error_bad
if __name__ == '__main__':
    
    params_to_write = output_params(dbparams=None, decompose=True)

    WF = WriteHtmlCSV('cl1358_1.0', 13.540899999999993, 13.9297, 221.8572889, 8.4680823, 22.9636, 0.0319, 1, 0.656651, 6, 2.546, 25.256, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, ['bulge', 'disk'], True, False, False, 71.0, 0.27, 0.73, 0.045)
    WF.writeparams(params_to_write, 65, 9999, 0.6)
