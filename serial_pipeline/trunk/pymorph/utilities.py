import datetime
import config as c
import time
from writehtmlfunc import *

def WriteDb(ParamValues):
    gal_id = ParamValues[0]
    hst = c.host
    dba = c.database
    pwd = c.pword
    usr = c.usr
    tbl = c.table
    host = c.host
    try:
        Conn = mysql.connect (host = "%s" %hst,
                                user = "%s" %usr,
                                passwd = "%s" %pwd,
                                db = "%s" %dba) 
    except mysql.Error, e:
        print "Error %d: %s" % (e.args[0], e.args[1])
        sys.exit (1)
    cursor = Conn.cursor()

    DictParamWithValue = {}
    DictParamWithType1 = {}
    AllParams = []
    for dbparam in c.dbparams:
        DBparam = dbparam.split(':')
        try:
            DictParamWithType1[DBparam[0]] = DBparam[2]
        except:
            DictParamWithType1[DBparam[0]] = 'varchar(500)'
        DictParamWithValue[DBparam[0]] = DBparam[1]
        AllParams.append(DBparam[0])
    Total_Run = 1
    Run = 1
    try:
        cmd = "SELECT Name, rootname from %s" % tbl
        cursor.execute(cmd)
        rows = cursor.fetchall()
        for row in rows:
            if str(row[0]) == gal_id:
                Total_Run += 1
            if str(row[0]) == gal_id and \
               str(row[1]) == c.rootname:
                Run += 1
    except:
        pass
    #    Run = 1
    #    try:
    #        cmd = "SELECT Name, rootname from %s" % tbl
    #        print '1'
    #        cursor.execute(cmd)
    #        print '2'
    #        rows = cursor.fetchall()
    #        print '3'
    #        for row in rows:
    #            print row[0], row[1]
    #            if str(row[0]) == gal_id and \
    #               str(row[1]) == DictParamWithValue['rootname']:
    #                Run += 1
    #    except:
    #        pass
    x=datetime.date.today()
    yyyy, mm, dd = x.year, x.month, x.day
    DaTe = str(yyyy) + '.' + str(mm) + '.' + str(dd)
    AllParams.append('Date')
    AllParams.append('Version')
    AllParams.append('Filter')
    AllParams.append('Total_Run')
    AllParams.append('rootname')
    AllParams.append('YetSky')
    DictParamWithType1['Date'] = 'varchar(50)' 
    DictParamWithType1['Version'] = 'float'
    DictParamWithType1['Filter'] = 'varchar(500)'
    DictParamWithType1['Total_Run'] = 'int'
    DictParamWithType1['rootname'] = 'varchar(500)'
    DictParamWithType1['YetSky'] = 'float'
    DictParamWithValue['Date'] = DaTe
    DictParamWithValue['Version'] = c.VERSION
    DictParamWithValue['Filter'] = c.FILTER
    DictParamWithValue['Total_Run'] = Total_Run
    DictParamWithValue['rootname'] = c.rootname
    DictParamWithValue['YetSky'] = c.SkyMin
    print 'writing db'
    if c.decompose:
        DictParamWithType2 = {'Name':'varchar(500)', 'ra_gal':'float', 
                              'dec_gal':'float',
                              'z':'float', 'mag_auto':'float', 
                              'magerr_auto':'float', 'num_targets':'float',
                              'bulge_xctr':'float',
                              'bulge_xctr_err':'float','bulge_yctr':'float',
                              'bulge_yctr_err':'float', 'Ie':'float','Ie_err':'float',
                              're_pix':'float', 're_err_pix':'float',
                              're_kpc':'float', 're_err_kpc':'float' ,'n':'float', 
                              'n_err':'float', 'AvgIe':'float', 'AvgIe_err':'float',
                              'eb':'float', 'eb_err':'float','bpa':'float', 
                              'bpa_err':'float', 'bboxy':'float', 'bboxy_err':'float',
                              'disk_xctr':'float','disk_xctr_err':'float',
                              'disk_yctr':'float','disk_yctr_err':'float',
                              'Id':'float', 'Id_err':'float', 'rd_pix':'float',
                              'rd_err_pix':'float', 'rd_kpc':'float', 
                              'rd_err_kpc':'float', 'ed':'float', 'ed_err':'float',
                              'dpa':'float', 'dpa_err':'float',
                              'dboxy':'float', 'dboxy_err':'float', 
                              'BD':'float', 'BT':'float', 'Point':'float', 
                              'Point_err':'float', 'Pfwhm':'float', 
                              'Pfwhm_kpc':'float', 'chi2nu':'float', 
                              'Goodness':'float', 'run':'int', 'C':'float', 
                              'C_err':'float', 'A':'float', 'A_err':'float', 
                              'S':'float', 'S_err':'float', 'G':'float', 'M':'float',
                              'SexSky':'float', 'GalSky':'float',
                              'GalSky_err':'float','dis_modu':'float', 
                              'distance':'float', 'fit':'int', 
                              'flag':'bigint', 'SexHalfRad':'float',
                              'bar_xctr':'float','bar_xctr_err':'float',
                              'bar_yctr':'float','bar_yctr_err':'float',
                              'BarMag':'float', 'BarMagErr':'float', 
                              'BarRePix':'float', 'BarRePixErr':'float',
                              'BarReKpc':'float', 'BarReKpcErr':'float',
                              'BarIndex':'float', 'BarIndexErr':'float',
                              'BarEll':'float', 'BarEllErr':'float',
                              'BarBoxy':'float',
                              'Manual_flag':'int', 'MorphType':'int',
                              'Comments':'varchar(1000)'}
        ParamToWrite = ['Name','ra_gal','dec_gal','z', 'mag_auto', 'magerr_auto',
                        'num_targets',
                        'bulge_xctr','bulge_xctr_err','bulge_yctr','bulge_yctr_err',
                        'Ie','Ie_err','re_pix', 're_err_pix','re_kpc', 're_err_kpc' ,
                        'n', 'n_err', 'AvgIe', 'AvgIe_err','eb', 'eb_err',
                        'bpa', 'bpa_err', 'bboxy', 'bboxy_err',
                        'disk_xctr','disk_xctr_err','disk_yctr','disk_yctr_err',
                        'Id', 'Id_err', 'rd_pix','rd_err_pix', 'rd_kpc','rd_err_kpc',
                        'ed', 'ed_err','dpa', 'dpa_err','dboxy', 'dboxy_err', 
                        'BD', 'BT', 'Point', 'Point_err', 'Pfwhm','Pfwhm_kpc',
                        'chi2nu', 'Goodness', 'run', 'C', 'C_err', 'A', 'A_err', 
                        'S', 'S_err', 'G', 'M','SexSky', 'GalSky',
                        'GalSky_err','dis_modu', 
                        'distance', 'fit','flag', 'SexHalfRad',
                        'bar_xctr','bar_xctr_err','bar_yctr','bar_yctr_err',
                        'BarMag', 'BarMagErr', 'BarRePix', 'BarRePixErr',
                        'BarReKpc', 'BarReKpcErr','BarIndex', 'BarIndexErr',
                        'BarEll', 'BarEllErr','BarBoxy',
                        'Manual_flag', 'MorphType', 'Comments']
                        
        ParamType = ['varchar(500)', 'float', 'float', 'float', 'float', 'float',
	             'float',
                     'float', 'float', 'float', 'float', 
                     'float', 'float', 'float', 'float', 'float', 'float',
                     'float', 'float', 'float', 'float', 'float', 'float', 
                     'float', 'float', 'float', 'float', 
                     'float', 'float', 'float', 'float', 
                     'float', 'float', 'float', 'float', 'float', 'float', 
                     'float', 'float', 'float', 'float', 'float', 'float', 
                     'float', 'float', 'float', 'float', 'float', 'float', 
                     'float', 'float', 'int', 'float', 'float', 'float', 'float', 
                     'float', 'float', 'float', 'float', 'float', 'float',
                     'float', 'float'
                     'int', 'bigint', 'int', 'float'
                     'float', 'float', 'float', 'float', 
                     'float', 'float', 'float', 'float', 
                     'float', 'float', 'float', 'float', 
                     'float', 'float', 'float',
                     'int','int', 'varchar(500)']
        
    else:
        DictParamWithType2 = {'Name':'varchar(500)', 'ra':'float', \
                        'dec_':'float',\
                        'z':'float', 'mag_auto':'float', \
			'magerr_auto':'float', 'num_targets':'int',\
			'C':'float', 'C_err':'float', 'A':'float',\
                        'A_err':'float', 'S':'float', 'S_err':'float',\
                        'G':'float', 'M':'float', \
                        'flag':'bigint', 'SexHalfRad':'float',\
                        'Manual_flag':'int', 'MorphType':'int',\
                        'Comments':'varchar(500)'}
        ParamToWrite = ['Name','ra','dec_','z', 'mag_auto', 'magerr_auto',
                        'num_targets',\
	                'C', 'C_err', 'A', 'A_err', 'S', 'S_err', 'G',\
                        'M', 'flag', 'SexHalfRad',\
                        'Manual_flag', 'MorphType', 'Comments']
        ParamType = ['varchar(500)', 'float', 'float', 'float', 'float',\
	             'float', 'float','float','float',\
                     'float', 'float', 'float','float', 'float', 'float',\
                     'bigint', 'float', 'int', 'int', 'varchar(500)']
    DictParamWithType = {}  #Dictionary with Type
    DictParamWithType.update(DictParamWithType1)
    DictParamWithType.update(DictParamWithType2)
    ParamValues.append(9999)
    ParamValues.append(9999)
    ParamValues.append('None')
    ii = 0
    for Param in ParamToWrite:
        DictParamWithValue[Param] = ParamValues[ii]
        ii += 1
    for p in ParamToWrite:
        AllParams.append(p)
    DictParamWithValue['run'] = Run
    if c.FirstCreateDB:
        cmd = "CREATE TABLE if not exists %s (" % tbl + ','.join(["%s %s" %(p, \
              DictParamWithType[p]) for p in AllParams]) + ")" 
        cursor.execute(cmd)
        c.FirstCreateDB = 0
    cmd = "INSERT INTO %s values (" % tbl 
    for p in AllParams:
        if DictParamWithType[p] in ('int', 'bigint', 'float'):
            cmd = cmd + str(DictParamWithValue[p]) + ', '
        else:
            cmd = cmd + "'" + str(DictParamWithValue[p]) + "', "
    cmd = str(cmd[:-2]) + ')'
    cursor.execute(cmd)
    cursor.close()
    Conn.close()


def WriteDbDetail(Name, ParamValuesDict, ErrDict, SexSky, GalSky, RunNo, flag, FitFlag, Chi2nu, model_type = '', goodness = 9999):
    gal_id = Name
    hst = c.host
    dba = c.database
    pwd = c.pword
    usr = c.usr
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
    except mysql.Error, e:
        print "Error %d: %s" % (e.args[0], e.args[1])
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
    DictParamWithType1['Name'] = 'varchar(500)' 
    DictParamWithType1['Date'] = 'varchar(50)' 
    DictParamWithType1['rootname'] = 'varchar(500)'
    DictParamWithValue['Name'] = Name
    DictParamWithValue['Date'] = DaTe
    DictParamWithValue['rootname'] = c.rootname
    DictParamWithValue['SexHalfRad'] = c.SexHalfRad
    DictParamWithValue['goodness'] = goodness
    print c.SexHalfRad
    
    DictParamWithType2 = {'xb':'float', 'yb':'float', \
                          'xd':'float', 'yd':'float', \
                          'Ie':'float','Ie_err':'float',\
                          're_pix':'float', 're_err_pix':'float',\
                          'n':'float', 'n_err':'float', \
                          'eb':'float', \
                          'Id':'float', 'Id_err':'float', 'rd_pix':'float',\
                          'rd_err_pix':'float', \
                          'ed':'float', \
                          'BT':'float', 'chi2nu':'float', 'run':'int', \
                          'SexSky':'float', 'GalSky':'float', \
                          'flag':'bigint','FitFlag':'bigint',\
                          'dpa':'float', 'bpa':'float',\
                          'ismain':'varchar(10)', 'SexHalfRad':'float', \
                          'goodness':'float'}
    ParamToWrite = ['xb', 'yb', 'xd', 'yd','bpa', 'dpa', 'Ie','Ie_err',\
                    're_pix','re_err_pix', 'n', 'n_err', 'eb', \
                    'Id', 'Id_err', 'rd_pix', 'rd_err_pix', \
                    'ed', 'BT', 'chi2nu', 'run', \
                    'SexSky', 'GalSky', 'flag', 'FitFlag', 'ismain',
                    'SexHalfRad', 'goodness']
    ParamType = ['float', 'float', 'float', 'float', 'float', 'float', \
                 'float', 'float', 'float', 'float', 'float', \
                 'float', 'float', 'float', 'float', 'float', \
                 'float', 'float', 'float', 'float', \
                 'int', 'float', 'float', 'bigint', 'bigint',
                 'varchar(10)','float', 'float']
    DictParamWithType = {}  #Dictionary with Type
    DictParamWithType.update(DictParamWithType1)
    DictParamWithType.update(DictParamWithType2)
    DictParamWithValue['xb'] = ParamValuesDict[1][2][0]
    DictParamWithValue['yb'] = ParamValuesDict[1][2][1]
    DictParamWithValue['Ie'] = ParamValuesDict[1][3]
    DictParamWithValue['re_pix'] = ParamValuesDict[1][4]
    DictParamWithValue['n'] = ParamValuesDict[1][5]
    DictParamWithValue['eb'] = ParamValuesDict[1][6]
    DictParamWithValue['bpa'] = ParamValuesDict[1][7]
    DictParamWithValue['ismain'] = ParamValuesDict[1][11]
    DictParamWithValue['xd'] = ParamValuesDict[2][2][0]
    DictParamWithValue['yd'] = ParamValuesDict[2][2][1]
    DictParamWithValue['Id'] = ParamValuesDict[2][3]
    DictParamWithValue['rd_pix'] = ParamValuesDict[2][4]
    DictParamWithValue['ed'] = ParamValuesDict[2][5]
    DictParamWithValue['dpa'] = ParamValuesDict[2][6]
    DictParamWithValue['Ie_err'] = ErrDict[1][2]
    DictParamWithValue['re_err_pix'] = ErrDict[1][3]
    DictParamWithValue['n_err'] = ErrDict[1][4]
    DictParamWithValue['Id_err'] = ErrDict[2][2]
    DictParamWithValue['rd_err_pix'] = ErrDict[2][3]
    fb = 10**((c.mag_zero - ParamValuesDict[1][3]) / 2.5)
    fd = 10**((c.mag_zero - ParamValuesDict[2][3]) / 2.5)
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
        if DictParamWithType[p] in ('int', 'bigint', 'float'):
            cmd = cmd + str(DictParamWithValue[p]) + ', '
        else:
            cmd = cmd + "'" + str(DictParamWithValue[p]) + "', "
    cmd = str(cmd[:-2]) + ')'
    cursor.execute(cmd)
    cursor.close()
    Conn.close()

