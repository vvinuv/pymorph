pymorph
=======
This pipeline was written as a part of some projects on galaxy morphology by Vinu Vikraman and Alan Meert. Yogesh Wadadekar, Ajit Kembhavi, G V Vijayagovidan and  Mariangela Bernardi were working on those projects. It REQUIRES SExtractor and GALFIT.  

### BUGS
Current version is buggy and don't use this.

### If you want to really use this

>>> import pymorph
>>> p = pymorph.PyMorph()
>>> p.pymorph()

This assumes that there is a configuration file in your directory. It runs also with default run parameters which are mostly different from the parameters in the configuration file. You can of course change this. Keep on read this


PyMorph is a pipeline which integrates SExtractor and GALFIT for estimating parametric and nonparametric quantities to describe the galaxy morphology. It can be used on a large frames with or without object coordinates of target galaxies or individual target objects. It can use PSF which is nearest to the target galaxy based on SExtractor star-galaxy probability or list of PSFs or individual PSF. It also gives the abilities to select the PSF interactively from a given PSF list or using SExtractor. It can repeat the selection procedure until you are satisfied. It can use SExtractor sky, GALFIT sky or the sky value implemented in PyMorph. It can repeat the fitting if GALFIT fails. It uses different initial parameters for GALFIT during refitting. It calculates nonparameteric quantities such as concentration, asymmetry, clumpness, Gini coefficient and moment parameters. These are implimented in PyMorph. It creates plots and a simple html page. It can write database for the output parameters.  

PyMorph assumes a configuration file in the current directory 

