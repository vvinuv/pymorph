import datetime
import config as c
import time
import MySQLdb as mysql

from writehtmlfunc import *

def WriteDb(ParamValues, all_params):
    gal_id = all_params['Name']
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

    # add additional keys to params
    for dbparam in c.dbparams:
        DBparam = dbparam.split(':')
        all_params[DBparam[0]] = DBparam[1]
        
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
    all_params['Date'] = DaTe
    all_params['Version'] = c.VERSION
    all_params['Filter'] = c.FILTER
    all_params['Total_Run'] = Total_Run
    all_params['rootname'] = c.rootname
    all_params['YetSky'] = c.SkyMin
    print 'writing db'

    if c.FirstCreateDB:
        cmd = "CREATE TABLE if not exists %s (" % tbl + ','.join(["%s %s" %(ParamValues[key][0], \
              ParamValues[key][1]) for key in ParamValues.keys()]) + ")" 
        cursor.execute(cmd)
        c.FirstCreateDB = 0
    cmd = "INSERT INTO %s (%s) values (" %(tbl, ','.join([ParamValues[key][0] for key in ParamValues.keys()])) 
    for p in ParamValues.keys():
        if ParamValues[p][1] in ('int', 'bigint', 'float'):
            cmd = cmd + str(all_params[ParamValues[p][0]]) + ', '
        else:
            cmd = cmd + "'" + str(all_params[ParamValues[p][0]]) + "', "
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

