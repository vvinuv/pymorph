from os.path import exists


class WriteHtmlFunc:
    """The class which will write html output"""
    def __init__(self, files, distance):
        self.files        = files
        self.distance     = distance
        self.write_params = write_params(files, distance)

def write_params(files, distance):
    outfile = open('result_' + str(files)[6:-4] + 'html','w')
    outfile.writelines(['<HTML> \n <BODY> \n'])
    outfile.writelines(['<TABLE BORDER="0" align="center" cellspacing=1',\
                        ' width=100%> \n'])
    if exists('fit.log'):
        for line in open('fit.log','r'): 
            values = line.split() 
            try: 
                if(str(values[0]) == 'Input'):
                    outfile.writelines(['<TR align="center" bgcolor="#CCFFFF">\
                                        <TH> Image </TH> <TD>', \
                                        str(values[3])[13:-18], \
                                        '</TD> </TR> \n' ])
                    indexfile = open('index.html', 'a+')
                    indexfile.writelines(['<a href="result',\
                                          str(values[3])[5:-18],'.html',\
                                          '"> ', str(values[3])[13:-18],\
                                          ' </a> <br>\n'])
                    indexfile.writelines(['</BODY></HTML>\n'])
                    indexfile.close()
                if(str(values[0]) == 'Init.'):
                    outfile.writelines(['<TR align="center" bgcolor="#99CCFF">\
                                        <TH> Init. par. file </TH> <TD>', \
                                        str(values[4]), '</TD> </TR> \n' ])		
                if(str(values[0]) == 'Restart'):
                    outfile.writelines(['<TR align="center" bgcolor="#CCFFFF">\
                                        <TH> Restart file </TH> <TD>', \
                                        str(values[3]), '</TD> </TR>\n' ])
                if(str(values[0]) == 'Chi^2/nu'):
                    outfile.writelines(['<TR align="center" bgcolor=',\
                                        '"#99CCFF"> <TH> &chi; <sup> 2',\
                                        ' </sup> / &nu; </TH> <TD>', \
                                        str(values[2]),'</TD> </TR>\n' ])
            except:
                pass
    outfile.writelines(['</TABLE> \n'])
    #outfile.writelines(['<CENTER><IMG SRC="plot_', str(files)[6:-4], 'png"></CENTER>'])
    outfile.writelines(['<TABLE BORDER="0" align="center" cellspacing=1',\
                        ' width=100%> \n \
                        <TR align="center"> <TD ROWSPAN=2><IMG SRC="plot_', \
                        str(files)[6:-4], 'png"></TD> <TD  bgcolor="#CCFFFF"> \
                        Separation between <br> psf and image</TD> </TR> \n',\
                        ' <TR align="center" bgcolor="#99CCFF"> <TD> ',\
                        str(round(distance, 3))[:5], \
                        '<br> arc sec</TD> </TR> \n </TABLE> \n'])
    outfile.writelines(['<TABLE BORDER="0" align="center" cellspacing=1',\
                        ' width=100%> \n'])	
    outfile.writelines(['<TR align="center" bgcolor="#CCFFFF"> <TH> Component \
                  <TH> xc <TH> yc <TH> mag <TH> radius <TH> n <TH> ellipticity\
                  <TH> pa <TH> box/disk </TH> </TR> \n'])
    if exists('fit.log'):
        for line in open('fit.log','r'): 
            values = line.split() 
            try: 
                if(str(values[0]) == 'sersic'):
                    outfile.writelines(['<TR align="center" bgcolor=',\
                                        '"#99CCFF"> <TD>', str(values[0]), \
                                        ' </TD> <TD> ', str(values[2])[1:-1],\
                                        ' </TD> <TD> ', str(values[3])[:-1], \
                                        ' </TD> <TD> ', str(values[4]), \
                                        ' </TD> <TD> ', str(values[5]), \
                                        ' </TD> <TD> ', str(values[6]), \
                                        ' </TD> <TD> ', str(values[7]), \
                                        ' </TD> <TD> ', str(values[8]), \
                                        ' </TD> <TD> ', str(values[9]), \
                                        ' </TD>  </TR> \n' ])
                if(str(values[0]) == 'expdisk'):
                    outfile.writelines(['<TR align="center" bgcolor="#99CCFF">\
                                         <TD>', str(values[0]), \
                                        ' </TD> <TD> ', str(values[2])[1:-1],\
                                        ' </TD> <TD> ', str(values[3])[:-1], \
                                        ' </TD> <TD> ', str(values[4]), \
                                        ' </TD> <TD> ', str(values[5]), \
                                        ' </TD> <TD> ', ' ' , \
                                        ' </TD> <TD> ', str(values[6]), \
                                        ' </TD> <TD> ', str(values[7]), \
                                        ' </TD> <TD> ', str(values[8]), \
                                        ' </TD>  </TR> \n'])
                if(str(values[0])[:1] == '('):
                    if(str(a) == 'sersic'):
                        outfile.writelines(['<TR align="center" \
                                            bgcolor="#CCFFFF"> <TD>', ' ' , \
                                            ' </TD> <TD> ',\
                                            str(values[0])[1:-1], \
                                            ' </TD> <TD> ', \
                                            str(values[1])[:-1], \
                                            ' </TD> <TD> ', str(values[2]), \
                                            ' </TD> <TD> ', str(values[3]), \
                                            ' </TD> <TD> ', str(values[4]), \
                                            ' </TD> <TD> ', str(values[5]), \
                                            ' </TD> <TD> ', str(values[6]), \
                                            ' </TD> <TD> ', str(values[7]), \
                                            ' </TD> </TR> \n' ])
                    if(str(a) == 'expdisk'):
                        outfile.writelines(['<TR align="center" \
                                            bgcolor="#CCFFFF"> <TD>', ' ' , \
                                            ' </TD> <TD> ',\
                                            str(values[0])[1:-1], \
                                            ' </TD> <TD> ',\
                                            str(values[1])[:-1], \
                                            ' </TD> <TD> ', str(values[2]), \
                                            ' </TD> <TD> ', str(values[3]), \
                                            ' </TD> <TD> ', ' ' , \
                                            ' </TD> <TD> ', str(values[4]), \
                                            ' </TD> <TD> ', str(values[5]), \
                                            ' </TD> <TD> ', str(values[6]), \
                                            ' </TD> </TR> \n' ])
                a=values[0]				
            except:
                pass
    outfile.write('</TABLE> \n')		
    outfile.write('</BODY></HTML> \n')
    outfile.close()
    outfile1 = open('galfit.html', 'w')
    outfile1.write('<HTML><HEAD> \n')
    outfile1.writelines(['<META HTTP-EQUIV="Refresh" CONTENT="10; \
                    URL=file:///Vstr/vstr/vvinuv/june-july07/',\
                    'pipeline/galfit.html"> \n'])
    outfile1.write('</HEAD> \n')
    outfile = 'result_' + str(files)[6:-4] + 'html'
    try:
        for line in open(outfile, 'r'):
            if(str(line) != '<html>'):
                outfile1.writelines([str(line)])
    except:
        pass
    outfile1.close()

