import MySQLdb as mysql
import sys
import config as c

def WriteDb(ParamValues):
    dba = c.database
    pwd = 'cluster'
    usr = c.usr
    tbl = c.table
    try:
        Conn = mysql.connect (host = "localhost",
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
    AllParams.append('Version')
    AllParams.append('Filter')
    DictParamWithType1['Version'] = 'float'
    DictParamWithType1['Filter'] = 'varchar(500)'
    DictParamWithValue['Version'] = c.VERSION
    DictParamWithValue['Filter'] = c.FILTER
    if c.decompose:
        DictParamWithType2 = {'Name':'varchar(500)', 'ra':'float', \
                        'dec_':'float',\
                        'z':'float', 'Ie':'float','Ie_err':'float',\
                        're_pixels':'float', 're_err_pixels':'float',\
                        're_kpc':'float', 're_err_kpc':'float' ,'n':'float', \
                       'n_err':'float', 'Avg_Ie':'float', 'Avg_Ie_err':'float',\
                        'eb':'float', 'eb_err':'float', \
                        'Id':'float', 'Id_err':'float', 'rd_pixels':'float',\
                        'rd_err_pixels':'float', 'rd_kpc':'float', \
                        'rd_err_kpc':'float', 'ed':'float', 'ed_err':'float', \
                        'BD':'float', 'BT':'float', 'Point':'float', \
                        'Point_err':'float', 'Pfwhm':'float', \
                        'Pfwhm_kpc':'float', 'chi2nu':'float', \
                        'Goodness':'float', 'run':'int', 'C':'float', \
                        'C_err':'float', 'A':'float', 'A_err':'float', \
                        'S':'float', 'S_err':'float', 'G':'float', 'M':'float',\
                        'distance':'float', 'fit':'int', 'flag':'bigint', \
                        'Comments':'varchar(1000)'}
        ParamToWrite = ['Name','ra','dec_','z', 'Ie','Ie_err','re_pixels',\
                        're_err_pixels', 're_kpc', 're_err_kpc' ,'n', \
                        'n_err', 'Avg_Ie', 'Avg_Ie_err', 'eb', 'eb_err', \
                        'Id', 'Id_err', 'rd_pixels',\
                        'rd_err_pixels', 'rd_kpc', 'rd_err_kpc', \
                        'ed', 'ed_err', 'BD', \
                        'BT', 'Point', 'Point_err', 'Pfwhm', 'Pfwhm_kpc', \
                        'chi2nu', 'Goodness', 'run', 'C', 'C_err', 'A', \
                        'A_err', 'S', 'S_err', 'G', 'M', 'distance', \
                        'fit', 'flag', 'Comments']
        ParamType = ['varchar(500)', 'float', 'float', 'float', 'float',\
                     'float', 'float', 'float', 'float', 'float', 'float',\
                     'float', 'float', 'float', 'float', 'float', 'float',\
                     'float', 'float', 'float', 'float', 'float', 'float',\
                     'float', 'float', 'float', 'float', 'float', 'float',\
                     'float', 'float', 'float', 'int', 'float', 'float',\
                     'float', 'float', 'float', 'float', 'float', 'float',\
                     'float', 'int', 'bigint', 'varchar(500)']
       
    else:
        DictParamWithType2 = {'Name':'varchar(500)', 'ra':'float', \
                        'dec_':'float',\
                        'z':'float', 'C':'float', 'C_err':'float', 'A':'float',\
                        'A_err':'float', 'S':'float', 'S_err':'float',\
                        'G':'float', 'M':'float', 'flag':'bigint', \
                        'Comments':'varchar(500)'}
        ParamToWrite = ['Name','ra','dec_','z', 'C', \
                        'C_err', 'A', 'A_err', 'S', 'S_err', 'G', 'M', \
                        'flag', 'Comments']
        ParamType = ['varchar(500)', 'float', 'float', 'float', 'float',\
                     'float', 'float', 'float','float', 'float', 'float',\
                     'float', 'bigint', 'varchar(500)']
    DictParamWithType = {}  #Dictionary with Type
    DictParamWithType.update(DictParamWithType1)
    DictParamWithType.update(DictParamWithType2)
    ParamValues.append('None')
    ii = 0
    for Param in ParamToWrite:
        DictParamWithValue[Param] = ParamValues[ii]
        ii += 1
    for p in ParamToWrite:
        AllParams.appeind(p)
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
#A = ['EDCSNJ1216490-1200091',184.204,-12.0025277778,0.7863,22.76,0.03,1.55,0.04,0.521182261202,0.0134498648052,1.24,0.12,18.8371116745,0.0316672494062,0.47,0.02,21.23,0.01,8.95,0.12,3.00940725017,0.0403495944157,0.35,0.0,0.244343055269,0.196363096362,9999,9999,9999,9999,1.146,0.664,1,3.63747679805,3.88004211969,0.130618587136,9999,0.365522099184,1.4538,0.546744991269,-2.24971293032,25.1297510444,1,1542]
#WriteDb(A)
#CREATE TABLE IF NOT EXISTS book (name char(40), lastname char(40), petname char (40))



