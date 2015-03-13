SHELL =		/bin/sh
CC =		gcc -g -O3 -m32
CFITSIO =	./cfitsio
CURSES =        /usr/ccs/lib
TERMCAP =	/usr/lib/termcap
SUSELIBS = 	-lcfitsio -lncurses -ltermcap -lm -lnsl -ldl
MACOSXLIBS =    -lcfitsio -lcurses -lm
LIBS =		-lcfitsio -lcurses -ltermcap -lm -lnsl -ldl
SPARC =         -lsocket

.c.o:
		$(CC) -c $<


# ============ description of all targets =============
#       -  <<-- ignore error code

MAIN  = galfit.o readfits.o writefits.o printerr.o \
	read_input.o getline.o clearbuff.o sigma.o spline.o \
	getmu.o nrutil.o search_min.o psfprep.o printpar.o \
	fft.o LevMar.o numrec.o egrid2d.o sigcheck.o \
	sersic.o convolve.o convregion.o nuker.o menu.o \
	clearmodel.o moffat.o gaussian.o exponential.o sky.o subimg.o \
	badpixlist.o assign_pars.o mrqmin.o mrqcof.o \
	work_arrays.o outmodel.o mkmodel.o copy_pars.o objselect.o \
	nuke.o pa.o beta.o constraints.o \
        fprintpar.o outmenu.o keypoll.o chi2calc.o errorcheck.o geterr.o \
	rtbis.o rtsec.o copy_struct.o parabolic.o writeheader.o assign_err.o \
	par_monitor.o read_cons.o strings.o spline_psf.o shift_psf.o \
	fftsg2d.o fftsg.o copymat.o smooth.o psf_fwhm.o get_cgrid.o \
	bcucof.o bcuint.o psf.o king.o

STRUCTS = sigma.o groupfits.o psfprep.o readfits.o menu.o \
	read_input.o search_min.o writefits.o printpar.o LevMar.o \
	pixsample.o convolve.o convregion.o nuker.o clearmodel.o \
	gaussian.o moffat.o exponential.o sky.o subimg.o \
	badpixlist.o assign_pars.o mrqmin.o mrqcof.o sigcheck.o \
	work_arrays.o outmodel.o mkmodel.o copy_pars.o objselect.o \
        keypoll.o outmenu.o chi2calc.o errorcheck.o geterr.o rtbis.o rtsec.o \
	copy_struct.o parabolic.o fprintpar.o writeheader.o assign_err.o \
	read_cons.o strings.o psf_fwhm.o get_cgrid.o bcucof.o bcuint.o psf.o \
	king.o

SERSIC = galfit.o sersic.o

FITSIO = readfits.o writefits.o printerr.o

redhat: $(MAIN) 
	${CC} -c initcurses.c
	${CC} -o galfit $(MAIN) initcurses.o -L$(CFITSIO)/lib $(LIBS) -I$(CFITSIO)/include

suse: $(MAIN)
	${CC} -c initcurses.c
	${CC} -o galfit $(MAIN) initcurses.o -L$(CFITSIO) -L$(TERMCAP) \
						$(SUSELIBS) -I$(CFITSIO)
sun: $(MAIN) 
	${CC} -DSUN=1 -c initcurses.c
	${CC} -o galfit $(MAIN) initcurses.o -L$(CFITSIO) -L$(CURSES) \
						$(LIBS) $(SPARC) -I$(CFITSIO)
jaguar: $(MAIN)
	${CC} -DMACOSX=1 -c initcurses.c
	${CC} -o galfit $(MAIN) initcurses.o -L$(CFITSIO) $(MACOSXLIBS)

panther: $(MAIN)
	${CC} -DMACOSX=1 -DPANTHER=1 -c initcurses.c
	${CC} -o galfit $(MAIN) initcurses.o -L$(CFITSIO) $(MACOSXLIBS)

$(STRUCTS): structs.h debug.h
#	${CC} -L$(CFITSIO) $(LIBS) -I$(CFITSIO) -c $*.c

$(SERSIC): sersic.h structs.h debug.h
#	${CC} -L$(CFITSIO) $(LIBS) -I$(CFITSIO) -c $*.c

$(FITSIO): $(CFITSIO)/include/fitsio.h debug.h
#	${CC} -L$(CFITSIO) $(LIBS) -I$(CFITSIO) -c $*.c

clean:
	rm *.o

distclean: clean
	rm galfit