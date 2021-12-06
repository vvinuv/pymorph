import sys
import datetime
import time
import mysql.connector
import numpy as np

class WriteDB(object):

    def __init__(self):
        pass 
        
    def _first_db(self, param_values, all_params, pymorph_config):
        

        gal_id = all_params['Name']
        try:
            conn = mysql.connector.connect(host=pymorph_config['host'],
                                           port=3306,
                                           user=pymorph_config['user'],
                                           password=pymorph_config['password'])
        except Exception as e:
            print('DB Writing problem')
            print(e)
            sys.exit()


        cursor = conn.cursor()

        pymorph_config['create_db'] = True

        if pymorph_config['create_db']:
            cmd = "CREATE DATABASE IF NOT EXISTS {} DEFAULT CHARACTER SET 'utf8'".format(pymorph_config['database'])
            pymorph_config['create_db'] = False

        cursor.execute(cmd)

        cursor.execute('USE {}'.format(pymorph_config['database']))

        x = datetime.date.today()
        yyyy = x.year
        mm = x.month
        dd = x.day
        curr_date = str(yyyy) + '.' + str(mm) + '.' + str(dd)
        #all_params['Date'] = curr_date
        #all_params['Version'] = 3. #pymorph_config.VERSION
        #all_params['Filter'] = 'A' #pymorph_config.FILTER
        #all_params['Rootname'] = 'R' #pymorph_config.rootname
        #all_params['YetSky'] = 0.0 #c.SkyMin

        print('Writing First DB')


        #param_values = {1: ['Name', 'varchar(500)'], 2: ['ra_gal', 'FLOAT']}
        pymorph_config['create_table'] = True      
        if pymorph_config['create_table']:
            cmd1 = "CREATE TABLE IF NOT EXISTS {} (".format(pymorph_config['table'])
            cmd2 = ', '.join(["{} {}".format(param_values[key][0], \
                    param_values[key][1]) for key in param_values.keys()]) 
            cmd = '{} {} )'.format(cmd1, cmd2)
            #print(cmd)
            cursor.execute(cmd)
            pymorph_config['create_table'] = False


        total_run = 1
        run = 1
        try:
            cmd = "SELECT NAME, ROOTNAME FROM {}".format(pymorph_config['table'])
            cursor.execute(cmd)
            rows = np.array(cursor.fetchall())
            #print(rows)
            if len(rows) != 0:
                run = rows[rows[:, 0] == gal_id].shape[0]
                total_run = rows[(rows[:, 0] == gal_id) & (rows[:, 1] == pymorph_config['rootname'])].shape[0]
                #print(run, total_run)
        except:
            pass


        
        all_params['Total_Run'] = total_run
        all_params['Run'] = run

        cmd = "INSERT INTO {} ({}) VALUES (".format(pymorph_config['table'], ','.join([param_values[key][0] for key in param_values.keys()])) 

        for p in param_values.keys():
            if param_values[p][1] in ('INT', 'BIGINT', 'FLOAT'):
                cmd = cmd + '{}, '.format(all_params[param_values[p][0]])
            else:
                cmd = cmd + '"{}",'.format(all_params[param_values[p][0]])

        cmd = cmd[:-2] + ')'
        #print(cmd)
        cursor.execute(cmd)
        conn.commit()
        cursor.close()
        conn.close()


def WriteDbDetail(Name, param_values_dict, err_dict, SexSky, GalSky, RunNo, flag, FitFlag, Chi2nu, model_type = '', goodness = 9999):
    gal_id = Name
    hst = pymorph_config.host
    dba = pymorph_config.database
    pwd = pymorph_config.pword
    usr = pymorph_config.usr
    if model_type == '':
        tbl = c.table + 'Detailed'
    else:
        tbl = c.table + model_type
    host = c.host
    try:
        Conn = mysql.connect (host = "%s" %hst,
                                user = "%s" %usr,
                                passwd = "%s" %pwd,
                                db = "%s" %dba) 
    except e:
        print("Error %d: %s" % (e.args[0], e.args[1]))
        sys.exit (1)
    cursor = Conn.cursor()
    DictParamWithValue = {}
    DictParamWithType1 = {}
    AllParams = []
    x=datetime.date.today()
    yyyy, mm, dd = x.year, x.month, x.day
    DaTe = str(yyyy) + '.' + str(mm) + '.' + str(dd)
    AllParams.append('Name')
    AllParams.append('Date')
    AllParams.append('rootname')
    DictParamWithType1['Name'] = 'VARCHAR(500)' 
    DictParamWithType1['Date'] = 'VARCHAR(50)' 
    DictParamWithType1['rootname'] = 'VARCHAR(500)'
    DictParamWithValue['Name'] = Name
    DictParamWithValue['Date'] = DaTe
    DictParamWithValue['rootname'] = pymorph_config.rootname
    DictParamWithValue['SexHalfRad'] = c.SexHalfRad
    DictParamWithValue['goodness'] = goodness
    print(c.SexHalfRad)
    
    DictParamWithType2 = {'xb':'FLOAT', 'yb':'FLOAT', \
                          'xd':'FLOAT', 'yd':'FLOAT', \
                          'Ie':'FLOAT','Ie_err':'FLOAT',\
                          're_pix':'FLOAT', 're_err_pix':'FLOAT',\
                          'n':'FLOAT', 'n_err':'FLOAT', \
                          'eb':'FLOAT', \
                          'Id':'FLOAT', 'Id_err':'FLOAT', 'rd_pix':'FLOAT',\
                          'rd_err_pix':'FLOAT', \
                          'ed':'FLOAT', \
                          'BT':'FLOAT', 'chi2nu':'FLOAT', 'run':'INT', \
                          'SexSky':'FLOAT', 'GalSky':'FLOAT', \
                          'flag':'BIGINT','FitFlag':'BIGINT',\
                          'dpa':'FLOAT', 'bpa':'FLOAT',\
                          'ismain':'VARCHAR(10)', 'SexHalfRad':'FLOAT', \
                          'goodness':'FLOAT'}
    ParamToWrite = ['xb', 'yb', 'xd', 'yd','bpa', 'dpa', 'Ie','Ie_err',\
                    're_pix','re_err_pix', 'n', 'n_err', 'eb', \
                    'Id', 'Id_err', 'rd_pix', 'rd_err_pix', \
                    'ed', 'BT', 'chi2nu', 'run', \
                    'SexSky', 'GalSky', 'flag', 'FitFlag', 'ismain',
                    'SexHalfRad', 'goodness']
    ParamType = ['FLOAT', 'FLOAT', 'FLOAT', 'FLOAT', 'FLOAT', 'FLOAT', \
                 'FLOAT', 'FLOAT', 'FLOAT', 'FLOAT', 'FLOAT', \
                 'FLOAT', 'FLOAT', 'FLOAT', 'FLOAT', 'FLOAT', \
                 'FLOAT', 'FLOAT', 'FLOAT', 'FLOAT', \
                 'INT', 'FLOAT', 'FLOAT', 'BIGINT', 'BIGINT',
                 'VARCHAR(10)','FLOAT', 'FLOAT']
    DictParamWithType = {}  #Dictionary with Type
    DictParamWithType.update(DictParamWithType1)
    DictParamWithType.update(DictParamWithType2)
    DictParamWithValue['xb'] = param_values_dict[1][2][0]
    DictParamWithValue['yb'] = param_values_dict[1][2][1]
    DictParamWithValue['Ie'] = param_values_dict[1][3]
    DictParamWithValue['re_pix'] = param_values_dict[1][4]
    DictParamWithValue['n'] = param_values_dict[1][5]
    DictParamWithValue['eb'] = param_values_dict[1][6]
    DictParamWithValue['bpa'] = param_values_dict[1][7]
    DictParamWithValue['ismain'] = param_values_dict[1][11]
    DictParamWithValue['xd'] = param_values_dict[2][2][0]
    DictParamWithValue['yd'] = param_values_dict[2][2][1]
    DictParamWithValue['Id'] = param_values_dict[2][3]
    DictParamWithValue['rd_pix'] = param_values_dict[2][4]
    DictParamWithValue['ed'] = param_values_dict[2][5]
    DictParamWithValue['dpa'] = param_values_dict[2][6]
    DictParamWithValue['Ie_err'] = err_dict[1][2]
    DictParamWithValue['re_err_pix'] = err_dict[1][3]
    DictParamWithValue['n_err'] = err_dict[1][4]
    DictParamWithValue['Id_err'] = err_dict[2][2]
    DictParamWithValue['rd_err_pix'] = err_dict[2][3]
    fb = 10**((pymorph_config.mag_zero - param_values_dict[1][3]) / 2.5)
    fd = 10**((pymorph_config.mag_zero - param_values_dict[2][3]) / 2.5)
    try:
        DictParamWithValue['BT'] = fb / (fb + fd)
    except:
        DictParamWithValue['BT'] = 9999.0
    DictParamWithValue['chi2nu'] = Chi2nu
    DictParamWithValue['run'] = RunNo
    DictParamWithValue['SexSky'] = SexSky
    DictParamWithValue['GalSky'] = GalSky
    DictParamWithValue['flag'] = flag
    DictParamWithValue['FitFlag'] = FitFlag
    for p in ParamToWrite:
        AllParams.append(p)
    cmd = "CREATE TABLE if not exists %s (" % tbl + ','.join(["%s %s" %(p, \
              DictParamWithType[p]) for p in AllParams]) + ")" 
    cursor.execute(cmd)
    cmd = "INSERT INTO %s values (" % tbl 
    for p in AllParams:
        if DictParamWithType[p] in ('INT', 'BIGINT', 'FLOAT'):
            cmd = cmd + str(DictParamWithValue[p]) + ', '
        else:
            cmd = cmd + "'" + str(DictParamWithValue[p]) + "', "
    cmd = str(cmd[:-2]) + ')'
    cursor.execute(cmd)
    cursor.close()
    Conn.close()


if __name__=='__main__':
    param_values = {1: ['Name', 'VARCHAR(500)'], 2: ['ra_gal', 'FLOAT'], 3: ['dec_gal', 'FLOAT'], 4: ['z', 'FLOAT'], 5: ['MorphType', 'INT'], 6: ['mag_auto', 'FLOAT'], 7: ['magerr_auto', 'FLOAT'], 8: ['SexHalfRad', 'FLOAT'], 9: ['num_targets', 'FLOAT'], 10: ['C', 'FLOAT'], 11: ['C_err', 'FLOAT'], 12: ['A', 'FLOAT'], 13: ['A_err', 'FLOAT'], 14: ['S', 'FLOAT'], 15: ['S_err', 'FLOAT'], 16: ['G', 'FLOAT'], 17: ['M', 'FLOAT'], 18: ['magzp', 'FLOAT'], 19: ['bulge_xctr', 'FLOAT'], 20: ['bulge_xctr_err', 'FLOAT'], 21: ['bulge_yctr', 'FLOAT'], 22: ['bulge_yctr_err', 'FLOAT'], 23: ['Ie', 'FLOAT'], 24: ['Ie_err', 'FLOAT'], 25: ['AvgIe', 'FLOAT'], 26: ['AvgIe_err', 'FLOAT'], 27: ['re_pix', 'FLOAT'], 28: ['re_pix_err', 'FLOAT'], 29: ['re_kpc', 'FLOAT'], 30: ['re_kpc_err', 'FLOAT'], 31: ['n', 'FLOAT'], 32: ['n_err', 'FLOAT'], 33: ['eb', 'FLOAT'], 34: ['eb_err', 'FLOAT'], 35: ['bpa', 'FLOAT'], 36: ['bpa_err', 'FLOAT'], 37: ['bboxy', 'FLOAT'], 38: ['bboxy_err', 'FLOAT'], 39: ['disk_xctr', 'FLOAT'], 40: ['disk_xctr_err', 'FLOAT'], 41: ['disk_yctr', 'FLOAT'], 42: ['disk_yctr_err', 'FLOAT'], 43: ['Id', 'FLOAT'], 44: ['Id_err', 'FLOAT'], 45: ['rd_pix', 'FLOAT'], 46: ['rd_pix_err', 'FLOAT'], 47: ['rd_kpc', 'FLOAT'], 48: ['rd_kpc_err', 'FLOAT'], 49: ['ed', 'FLOAT'], 50: ['ed_err', 'FLOAT'], 51: ['dpa', 'FLOAT'], 52: ['dpa_err', 'FLOAT'], 53: ['dboxy', 'FLOAT'], 54: ['dboxy_err', 'FLOAT'], 55: ['BD', 'FLOAT'], 56: ['BT', 'FLOAT'], 57: ['p_xctr', 'FLOAT'], 58: ['p_xctr_err', 'FLOAT'], 59: ['p_yctr', 'FLOAT'], 60: ['p_yctr_err', 'FLOAT'], 61: ['Ip', 'FLOAT'], 62: ['Ip_err', 'FLOAT'], 63: ['Pfwhm', 'FLOAT'], 64: ['Pfwhm_kpc', 'FLOAT'], 65: ['bar_xctr', 'FLOAT'], 66: ['bar_xctr_err', 'FLOAT'], 67: ['bar_yctr', 'FLOAT'], 68: ['bar_yctr_err', 'FLOAT'], 69: ['Ibar', 'FLOAT'], 70: ['Ibar_err', 'FLOAT'], 71: ['rbar_pix', 'FLOAT'], 72: ['rbar_pix_err', 'FLOAT'], 73: ['rbar_kpc', 'FLOAT'], 74: ['rbar_kpc_err', 'FLOAT'], 75: ['n_bar', 'FLOAT'], 76: ['n_bar_err', 'FLOAT'], 77: ['ebar', 'FLOAT'], 78: ['ebar_err', 'FLOAT'], 79: ['barpa', 'FLOAT'], 80: ['barpa_err', 'FLOAT'], 81: ['barboxy', 'FLOAT'], 82: ['barboxy_err', 'FLOAT'], 83: ['chi2nu', 'FLOAT'], 84: ['Goodness', 'FLOAT'], 85: ['Run', 'INT'], 86: ['SexSky', 'FLOAT'], 87: ['YetSky', 'FLOAT'], 88: ['GalSky', 'FLOAT'], 89: ['GalSky_err', 'FLOAT'], 90: ['dis_modu', 'FLOAT'], 91: ['distance', 'FLOAT'], 92: ['fit', 'INT'], 93: ['FitFlag', 'BIGINT'], 94: ['flag', 'BIGINT'], 95: ['Manual_flag', 'BIGINT'], 96: ['Comments', 'VARCHAR(1000)'], 97: ['Date', 'VARCHAR(50)'], 98: ['Version', 'FLOAT'], 99: ['Filter', 'VARCHAR(500)'], 100: ['Total_Run', 'INT'], 101: ['rootname', 'VARCHAR(500)'], 102: ['Cluster', 'VARCHAR(500)'], 103: ['ObsID', 'INT'], 104: ['distance_psf_gal', 'FLOAT']}
    all_params = {'Name': 'cl1358_10.0', 'ra_gal': 221.8584559, 'dec_gal': 8.4736461, 'z': 9999, 'MorphType': -999, 'mag_auto': 20.3232, 'magerr_auto': 0.0076, 'SexHalfRad': 5.191, 'num_targets': 1, 'C': 9999, 'C_err': 9999, 'A': 9999, 'A_err': 9999, 'S': 9999, 'S_err': 9999, 'G': 9999, 'M': 9999, 'magzp': 25.256, 'bulge_xctr': 33.09, 'bulge_xctr_err': 0.09, 'bulge_yctr': 32.07, 'bulge_yctr_err': 0.07, 'Ie': 27.65, 'Ie_err': 0.05, 'AvgIe': 28.399203544535066, 'AvgIe_err': 0.05609498855159012, 're_pix': 12.82, 're_pix_err': 0.9, 're_kpc': -999, 're_kpc_err': -999, 'n': 4.0, 'n_err': 0.0, 'eb': 0.3, 'eb_err': 0.01, 'bpa': 82.67, 'bpa_err': 1.07, 'bboxy': -999.0, 'bboxy_err': -999.0, 'disk_xctr': 31.53, 'disk_xctr_err': 0.08, 'disk_yctr': 32.78, 'disk_yctr_err': 0.05, 'Id': 27.85, 'Id_err': 0.04, 'rd_pix': 2.87, 'rd_pix_err': 0.05, 'rd_kpc': -999, 'rd_kpc_err': -999, 'ed': 0.59, 'ed_err': 0.01, 'dpa': -77.23, 'dpa_err': 2.01, 'dboxy': -999.0, 'dboxy_err': -999.0, 'BD': 1.2022644346174163, 'BT': 0.5459219227804843, 'p_xctr': -999, 'p_xctr_err': -999, 'p_yctr': -999, 'p_yctr_err': -999, 'Ip': -999, 'Ip_err': -999, 'Pfwhm': -999, 'Pfwhm_kpc': -999, 'bar_xctr': -999, 'bar_xctr_err': -999, 'bar_yctr': -999, 'bar_yctr_err': -999, 'Ibar': -999, 'Ibar_err': -999, 'rbar_pix': -999, 'rbar_pix_err': -999, 'rbar_kpc': -999, 'rbar_kpc_err': -999, 'n_bar': -999, 'n_bar_err': -999, 'ebar': -999, 'ebar_err': -999, 'barpa': -999, 'barpa_err': -999, 'barboxy': -999, 'barboxy_err': -999, 'chi2nu': 0.831, 'Goodness': -999, 'Run': 1, 'SexSky': 0.6559478, 'YetSky': -999, 'GalSky': 0.656, 'GalSky_err': 0.0, 'dis_modu': -999, 'distance': -999, 'fit': -999, 'FitFlag': 1, 'flag': 6, 'Manual_flag': -999, 'Comments': -999, 'Date': -999, 'Version': -999, 'Filter': -999, 'Total_Run': -999, 'rootname': -999, 'Cluster': -999, 'ObsID': -999, 'distance_psf_gal': 62.9}
    

    pymorph_config = dict(host='localhost', database='Galaxy', table='cluster', user='vinu', password='kerala', dbparams=['Cluster:HI'], rootname='-999')

    print(pymorph_config['host'])
    WDB = WriteDB()
    WDB._first_db(param_values, all_params, pymorph_config)
