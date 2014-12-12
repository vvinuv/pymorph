

def InitializeHtmlIndexFile():
    #Initialize index.html
    if exists('index.html'):
        pass
    else:
        indexfile = open('index.html', 'w')
        indexfile.writelines(['<HTML>\n<BODY>\n'])
        indexfile.writelines(['</BODY></HTML>'])
        indexfile.close()

def InitializeResultCSV(ParamToWrite):
    if exists('result.csv'):
        pass
    else:
        f_res = open("result.csv", "ab")
        csvlist = ['%s_%d'%(ParamToWrite[par_key][0], par_key)
                   for par_key in ParamToWrite.keys()]
        print csvlist
        writer = csv.writer(f_res)
        writer.writerow(csvlist)
        f_res.close()

def InitializeFailedList(pnames):
    # writing a input catalogue (restart.cat) for failed objects
    for p in pnames:
        f_failed.writelines(['%s '%p])
    f_failed.writelines(['flag \n'])


def OpenOutcat():
    # Opens files to write output
    f_cat = open(c.out_cata, 'w')
    f_failed = open('restart.cat', 'w')


