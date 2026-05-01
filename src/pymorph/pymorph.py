import os
import sys
from string import Template
from pathlib import Path
import datetime
import configparser
import subprocess
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from .maskfunc import MaskGenerator
from .galfit_config_run import GalfitConfigRunFunc
from .casgm import CASGMPipeline
#from writecsv import WriteCSV
from .get_params import GetOutputParams
from .psffunc import PSFPipeline
from .writehtml import generate_galaxy_report
from .plotfunc import PlotFunc
from .errors_warnings import catch_pipeline_issues, PipelineCriticalError
from .errors_warnings import PipelineBase
from .detailed_galfit import GalfitDetailed


class GalaxyPipeline:

    def __init__(self, config_file):

        BASE_DIR = Path(__file__).parent

        config = configparser.ConfigParser()
        config.read(config_file)
        self.imagefile = config.get("imagecata", "imagefile")
        self.whtfile = config.get("imagecata", "whtfile")
        self.sex_cata = config.get("imagecata", "sex_cata")
        self.obj_cata = config.get("imagecata", "obj_cata")
        self.rootname = config.get("imagecata", "rootname")
        self.datadir = config.get("imagecata", "datadir")
        self.imagefile = os.path.join(self.datadir, self.imagefile)
        self.whtfile = os.path.join(self.datadir, self.whtfile)
        self.obj_cata = os.path.join(self.datadir, self.obj_cata)

        self.mask_reg = config.getfloat('mask', 'mask_reg')
        self.thresh_area = config.getfloat('mask', 'thresh_area')
        self.threshold = config.getfloat('mask', 'threshold')

        self.pixelscale = config.getfloat('cosmology', 'pixelscale') 

        size = config.get('size', 'size_list')
        size = [int(s) for s in size.split(',')]
        self.ReSize = size[0]
        self.ImSize = size[1]
        self.FracRad = size[2]
        self.Square = size[3]
        self.FixSize = size[4]

        self.searchrad = config.getfloat('size', 'searchrad')

        self.nearest_neighbour = config.getboolean('size', 'nearest_neighbour')

        self.mag_zero = config.getfloat('mag', 'mag_zero')
        self.photo_filter = config.get('mag', 'photo_filter')

        self.maglim = config.get('findfit', 'maglim').split(',')
        self.maglim = [float(mlim) for mlim in self.maglim]

        components = config.get('galfit', 'components').split(',')
        self.components = [cm.strip() for cm in components]


        self.run_galfit = config.getboolean('modes', 'decompose') 

        self.galcut = config.getboolean("modes", "galcut")
        self.detail = config.getboolean("modes", "detail")
        self.find_and_fit = config.getboolean("modes", "find_and_fit")
    
        self.bbox = config.getfloat("galfit", "bbox") 
        self.barbox = config.getfloat("galfit", "barbox") 
        fitting = config.get("galfit", "fitting")
        fitting = [int(tfit) for tfit in fitting.split(',')]
        self.bulge_cntr = fitting[0]
        self.disk_cntr = fitting[1]
        self.sky_fit = fitting[2]
        self.bar_cntr = fitting[3]

        self.chi2nu_limit = config.getfloat("diagnosis", "chi2nu_limit")
        self.center_deviation_limit = config.getfloat("diagnosis", "center_deviation_limit")

        self.NMag = config.getfloat('params_limit', 'NMag')
        self.NRadius = config.getfloat('params_limit', 'NRadius')

        self.LN = config.getfloat('params_limit', 'LN')
        self.UN = config.getfloat('params_limit', 'UN')

        self.sex_keys = {key.upper(): value
                         for key, value in config['sextractor'].items()
                         }
        self.sex_keys['DATADIR'] = self.datadir
        self.sex_keys['WHTFILE'] = self.whtfile
        self.sex_keys['MAGZERO'] = self.mag_zero
        self.sex_keys['PIXEL_SCALE'] = self.pixelscale
        self.sex_keys['SEx_CATA'] = self.sex_cata
        self.sex_keys['PYMORPH_PATH'] = Path(__file__).parent
        self.sex_keys['PIXEL_SCALE'] = self.pixelscale

        self.sex_catalog = None
        self.obj_catalog = None

        self.neighbour_cols = [
            "X_IMAGE",
            "Y_IMAGE",
            "ALPHA_J2000",
            "DELTA_J2000",
            "FLUX_RADIUS",
            "THETA_IMAGE",
            "ELONGATION",
            "A_IMAGE",
            "MAG_AUTO",
            "ISO0",
            "BACKGROUND"
        ]
        
        self.flags = {}




    # KEEP FLAGS
    def set_flag(self, name, flag):
        if name not in self.flags:
            self.flags[name] = flag   # create list


    def check_bulge(self):
        if 'bulge' in self.components:
            self.set_flag("HAS_BULGE", 3)

    def check_disk(self):
        if 'disk' in self.components:
            self.set_flag("HAS_DISK", 4)

    def check_bar(self):
        if 'bar' in self.components:
            self.set_flag("HAS_BAR", 5)

    def check_bulge_box(self):
        if self.bbox != -99:
            self.set_flag("BULGE_BOXINESS", 6)

    def check_bar_box(self):
        if self.barbox != -99:
            self.set_flag("BAR_BOXINESS", 7)

    def check_sky_fit(self):
        if self.sky_fit:
            self.set_flag("SKY_FIT", 8)

    def check_detail(self):
        if self.detail:
            self.set_flag("DETAILED_FIT", 9)
            
    def check_cutout(self):
        if self.galcut:
            self.set_flag("CUT_OUT_IMAGE", 10)

    def check_find_and_fit(self):
        if self.find_and_fit:
            self.set_flag("FIND_AND_FIT", 11)

    def check_flags(self):
        self.check_bulge()
        self.check_disk()
        self.check_bar()
        self.check_bulge_box()
        self.check_bar_box()
        self.check_sky_fit()
        self.check_detail()
        self.check_cutout()
        self.check_find_and_fit()

    def check_radec(self):
        if self.flag_radec:
            self.set_flag("RA_DEC", 12)


    def check_neighbours(self, number):
        if number > 0:
            self.set_flag("NEIGHBOUR_FIT", 13)

        
    #CHECK THE SIZE OF THE IMAGE. THE FLAG IS WRITTEN IN ERRORS_WARNINGS FILE
    @catch_pipeline_issues()
    def check_image_size(self, size):
        return size


    #CHECK WHETHER FIT.LOG EXISTS. IF NOT WE ASSUME THAT GALFIT DIDN'T RUN
    #THE FLAG IS WRITTEN IN ERRORS_WARNINGS FILE
    @catch_pipeline_issues(file_checker=True)
    def check_file(self, filename):
        return filename


    #PARAMETERS AND ERROR CHECK. IF PARAMETERS IS NOT FOUND IN LEGTH GIVE
    #ERROR. IF THERE IS ERROR IS NONE, INF, NAN, NULL OR NO LIST, IT WILL BE
    #REPLACED BY 9999
    #THE FLAG IS WRITTEN IN ERRORS_WARNINGS FILE
    @catch_pipeline_issues(expected_len=7)
    def check_params_and_errors(self, params, errors):
        return params, errors


    def save(self):
        """Save results + flags per object"""
        data = {
            "id": self.obj_id,
            "flags": self.flags,
        }

        with open(f"result_errors.json", "w") as f:
            json.dump(data, f, indent=2)

    def p(self):
        print(self.flags)

    

    
    # ---------------------------------------------------------
    # Run SExtractor
    # ---------------------------------------------------------
    def run_sextractor(self):

        tpl_file = os.path.join(self.sex_keys["PYMORPH_PATH"], 
                                "SEx/default.sex")

        with open(tpl_file) as f:
            tpl = Template(f.read())
        output = tpl.safe_substitute(self.sex_keys)
        with open(os.path.join(Path.cwd(), 'default.sex'), 'w') as f:
            f.write(output)

        output_cat = f"{self.rootname}_sex.cat"
        sex_out_params = os.path.join(self.sex_keys["PYMORPH_PATH"], 
                                      "SEx/default.param")


        cmd = [
            "sex",
            self.imagefile,
            "-WEIGHT_IMAGE", self.whtfile,
            "-CATALOG_NAME", output_cat,
            "-PARAMETERS_NAME", sex_out_params,
            "-c", "default.sex"
        ]

        subprocess.run(cmd, check=True)

        self.sex_cata = output_cat
        return output_cat



    # ------------------------------------------------
    # Read catalog using header lines (#)
    # ------------------------------------------------
    def read_sex_catalog_header(self, filename):

        colnames = []

        with open(filename) as f:
            for line in f:
                if line.startswith("#"):
                    parts = line.strip().split()
                    if len(parts) >= 3:
                        colnames.append(parts[2])
                else:
                    break

        df = pd.read_csv(
            filename,
            sep=r"\s+",
            comment="#",
            names=colnames
        )

        return df

# ------------------------------------------------
    # Read SExtractor catalog
    # ------------------------------------------------
    def read_sex_catalog(self):
        sex_column_names = ['NUMBER', 'X_IMAGE', 'Y_IMAGE', 'ALPHA_J2000', 
                            'DELTA_J2000', 'FLUX_ISO', 'FLUXERR_ISO', 
                            'MAG_ISO', 'MAGERR_ISO', 'FLUX_RADIUS', 
                            'BACKGROUND', 'THETA_IMAGE', 'ELONGATION', 'ISO0',
                            'A_IMAGE', 'FLAGS', 'CLASS_STAR', 'MAG_AUTO', 
                            'MAGERR_AUTO']

        self.sex_catalog = pd.read_csv(
            self.sex_cata,
            sep=r"\s+",
            comment="#",
            names=sex_column_names
        )


    # ------------------------------------------------
    # Convert HMS/DMS → degrees
    # ------------------------------------------------
    def hms_dms_to_deg(self, ra1, ra2, ra3, dec1, dec2, dec3):

        ra = 15 * (ra1 + ra2/60 + ra3/3600)

        sign = -1 if dec1 < 0 else 1
        dec = sign * (abs(dec1) + dec2/60 + dec3/3600)

        return ra, dec


    # ------------------------------------------------
    # Convert degrees -> HMS/DMS
    # ------------------------------------------------

    def deg_to_hms_dms(self, ra_deg, dec_deg):
        """
        Convert Right Ascension from degrees to HH:MM:SS and
        Convert Declination from degrees to DD:MM:SS
        """
        ra_hours = ra_deg / 15.0

        h = int(ra_hours)
        m = int((ra_hours - h) * 60)
        s = (ra_hours - h - m/60) * 3600

        sign = "+" if dec_deg >= 0 else "-"
        dec_deg = abs(dec_deg)

        d = int(dec_deg)
        m = int((dec_deg - d) * 60)
        s = (dec_deg - d - m/60) * 3600

        ra_hms = f"{h:02d}:{m:02d}:{s:04.1f}"
        dec_dms = f"{sign}{d:02d}:{m:02d}:{s:04.1f}"

        return ra_hms, dec_dms


    # ------------------------------------------------
    # Load object catalog
    # ------------------------------------------------
    def load_obj_catalog(self):

        self.obj_catalog = pd.read_csv(
            self.obj_cata,
            sep=r"\s+",
            comment="#",
        )
        self.prepare_radec()


    # ------------------------------------------------
    # Prepare RA/DEC
    # ------------------------------------------------
    def prepare_radec(self):

        cols = self.obj_catalog.columns

        if "ra" in cols and "dec" in cols:
            return

        if all(c in cols for c in
               ["ra1","ra2","ra3","dec1","dec2","dec3"]):

            ra_list = []
            dec_list = []

            for _, row in self.obj_catalog.iterrows():

                ra, dec = self.hms_dms_to_deg(
                    row["ra1"], row["ra2"], row["ra3"],
                    row["dec1"], row["dec2"], row["dec3"]
                )

                ra_list.append(ra)
                dec_list.append(dec)

            self.obj_catalog["ra"] = ra_list
            self.obj_catalog["dec"] = dec_list


    # ------------------------------------------------
    # Check if object exists in sextractor catalog
    # ------------------------------------------------
    def position_in_sexcat(self, ra, dec, tol_arcsec=1):

        
        sex_coords = SkyCoord(
            self.sex_catalog["ALPHA_J2000"].values * u.deg,
            self.sex_catalog["DELTA_J2000"].values * u.deg
        )

        target = SkyCoord(ra*u.deg, dec*u.deg)

        #print(target)
        #print(sex_coords)
        sep = target.separation(sex_coords)
        #print(sep.shape)
        #print(np.argmin(sep))
        #print(target)

        print('np.min(sep).arcsec', np.min(sep).arcsec)
        return np.min(sep).arcsec < tol_arcsec


    # ------------------------------------------------
    # Find neighbours in SExtractor catalog using RA DEC
    # ------------------------------------------------

    def neighbours_radec(self, ra, dec, radius_arcmic=1):

        radius_arcsec = radius_arcmic * 60

        sex_coords = SkyCoord(
            self.sex_catalog["ALPHA_J2000"].values * u.deg,
            self.sex_catalog["DELTA_J2000"].values * u.deg
        )

        target = SkyCoord(ra * u.deg, dec * u.deg)

        sep = target.separation(sex_coords)

        if sep.arcsec.any() < self.searchrad:
            target_found = True

        # subset dataframe
        neighbours = self.sex_catalog[self.neighbour_cols].copy()
        # corresponding separations

        if len(neighbours) == 0:
            return None, neighbours

        # find closest object (target)
        min_idx = sep.arcsec.argmin()

        target_sex = neighbours.iloc[min_idx]

        # remove target from neighbours
        neighbours = neighbours.drop(neighbours.index[min_idx])

        return target_sex, neighbours.reset_index(drop=True)


    # ------------------------------------------------
    # Find neighbours in SExtractor catalog using X/Y
    # ------------------------------------------------

    def neighbours_xy(self, x, y, dx=50, dy=50):

        # select objects inside box
        mask = (
            (np.abs(self.sex_catalog["X_IMAGE"] - x) < dx) &
            (np.abs(self.sex_catalog["Y_IMAGE"] - y) < dy)
        )

        neighbours = self.sex_catalog[self.neighbour_cols].copy()

        if len(neighbours) == 0:
            return None, neighbours

        # compute distance (Euclidean)
        dx_arr = neighbours["X_IMAGE"].values - x
        dy_arr = neighbours["Y_IMAGE"].values - y

        dist = np.sqrt(dx_arr**2 + dy_arr**2)

        # find closest object
        min_idx = np.argmin(dist)

        target_sex = neighbours.iloc[min_idx]
        
        # remove closest object from neighbours
        neighbours = neighbours.drop(neighbours.index[min_idx])

        return target_sex, neighbours.reset_index(drop=True)

    # ------------------------------------------------
    # Process one target object and its neighbours
    # ------------------------------------------------


    def process_target(self, target, radius_arcmic=0.5):

        self.NAME = f"{self.rootname}_{int(target["GAL_ID"])}"

        ra = target.get("RA", np.nan)
        dec = target.get("DEC", np.nan)

        x = target.get("XIMG", np.nan)
        y = target.get("YIMG", np.nan)


        position = False
        
        self.flag_radec = False

        if not pd.isna(ra) and not pd.isna(dec):

            position = self.position_in_sexcat(ra, dec)

            target_sex, neighbours = self.neighbours_radec(ra, dec,
                                                           radius_arcmic=0.5)
            self.flag_radec = True
            
            
        else:

            target_sex, neighbours = self.neighbours_xy(x, y)

        #print(target_sex)
        #print(position)

        #print(target_sex)
        target = target.to_dict()
        target["GAL_ID"] = int(target["GAL_ID"])
        target.setdefault("Z", -999)
        target.setdefault("SKY", -999)
        if "SKY" in target:
            target["SKY"] = target["SKY"] if target["SKY"] is not None else -999
            #print('target["SKY"]', target["SKY"])
        target['GALFIT_ANGLE'] = target_sex["THETA_IMAGE"] - 90 
        target["POSITION"] = position
        target["ROOTNAME"] = self.rootname
        #NAME = f'{self.rootname}_{target['GAL_ID']}'
        target = {"NAME": self.NAME, **target}
        target["MAG_ZERO"] = self.mag_zero
        target["FILTER"] = self.photo_filter
        
        target['DATE'] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        target.update(target_sex)

        #print(target["ALPHA_J2000"])
        ra_hms, dec_dms = self.deg_to_hms_dms(target["ALPHA_J2000"], 
                                              target["DELTA_J2000"])

        target["RA_HMS"] = ra_hms
        target["DEC_DMS"] = dec_dms

        target.pop("ALPHA_J2000", None)
        target.pop("DELTA_J2000", None)

        #target["RA"] = 9999
        #target["DEC"] = 9999

        if target['SKY'] == -999:
            target.pop('SKY', None)
        else:
            target['BACKGROUND'] = target.pop('SKY', None)


        #print(self.FracRad,  target['FLUX_RADIUS'])

        #sys.exit()
        if self.ReSize:
            size = int(self.FracRad * target['FLUX_RADIUS'])
        else:
            size = self.FixSize 


        #CRITICAL FLAG 
        self.check_image_size(size)

    

        target['IMAGE_SIZE'] = int(size)
        #print('size', target['IMAGE_SIZE'])

        neighbours['GALFIT_ANGLE'] = neighbours["THETA_IMAGE"] - 90
        #neighbours["RA"] = 9999
        #neighbours["DEC"] = 9999

        galaxies = dict()
        galaxies['target'] = target
        galaxies['neighbours'] = neighbours

        return galaxies


    def generate_target_images(self, galaxies):

        img = fits.open(self.imagefile)
        img_data = img[0].data
        wcs = WCS(img[0].header)

        wht = fits.open(self.whtfile)
        wht_data = wht[0].data

        target = galaxies.get("target")
        gal_id = int(target.get("GAL_ID"))
        size = target.get("IMAGE_SIZE") 
        
        #self.check_image_size(size)

        ra = target.get("RA")
        dec = target.get("DEC")

        x = target.get("XIMG") or target.get("X_IMAGE")
        y = target.get("YIMG") or target.get("Y_IMAGE")

        # determine pixel position
        if ra is not None and dec is not None:
            x, y = wcs.world_to_pixel_values(ra, dec)

        position = (x, y)

        #print('position, size', position, size)
        cutout_img = Cutout2D(
            img_data,
            position,
            (size, size),
            wcs=wcs
        )

        cutout_wht = Cutout2D(
            wht_data,
            position,
            (size, size),
            wcs=wcs
        )

        img_name = f"I_{self.rootname}_{gal_id}.fits"
        wht_name = f"W_{self.rootname}_{gal_id}.fits"

        fits.writeto(
            img_name,
            cutout_img.data,
            header=cutout_img.wcs.to_header(),
            overwrite=True
        )

        fits.writeto(
            wht_name,
            cutout_wht.data,
            header=cutout_wht.wcs.to_header(),
            overwrite=True
        )
        
        #sys.exit()

    def transform_to_cutout(self, target, neighbours):
        
        #print(self.FracRad,  target['FLUX_RADIUS'])

        if self.ReSize:
            size = self.FracRad * target['FLUX_RADIUS']
        elif self.ImSize:
            size = fits.open(self.imagefile)[0].data.shape[0]
        else:
            size = self.FixSize

        #print('size', size)
        #sys.exit()
        #CRITICAL FLAG
        #self.check_image_size(size)

        target['IMAGE_SIZE'] = int(size)
        #print('size', target['IMAGE_SIZE'])


        size = target['IMAGE_SIZE']
        half = size // 2

        # Target center
        x0 = target["X_IMAGE"]
        y0 = target["Y_IMAGE"]

        #print('x0', x0)
        # Shift
        dx = x0 - half
        dy = y0 - half
        #print('dx', dx, half)
        
        # --- Transform target ---
        target["X_IMAGE"] = target["X_IMAGE"] - dx
        target["Y_IMAGE"] = target["Y_IMAGE"] - dy
        target["DX"] = dx
        target["DY"] = dy


        #print('x0', x0, target["X_IMAGE"])

        # --- Transform neighbours ---
        neighbours["X_IMAGE"] = neighbours["X_IMAGE"] - dx
        neighbours["Y_IMAGE"] = neighbours["Y_IMAGE"] - dy

        # --- Keep only neighbours inside cutout ---
        mask = (
                (neighbours["X_IMAGE"] >= 0) &
                (neighbours["X_IMAGE"] < size) &
                (neighbours["Y_IMAGE"] >= 0) &
                (neighbours["Y_IMAGE"] < size)
            )

        neighbours = neighbours.loc[mask].reset_index(drop=True)
        
        return target, neighbours


class PyMorph:
    
    def __init__(self, config_file='config.ini'):

        config = configparser.ConfigParser()
        config.read(config_file)
        self.galcut = config.getboolean("modes", "galcut")
        self.detail = config.getboolean("modes", "detail")
        self.number_random = config.getint("modes", "number_random")


        if self.galcut and not self.detail:
            self.run_cutout_model()
        elif not self.galcut and not self.detail:
            self.run()
        elif self.galcut and self.detail:
            self.run_cutout_detailed()
        elif not self.galcut and self.detail:
            self.run_detailed()

    def flags_galfit(self, target, pipeline):

        if target["chi2nu"] > pipeline.chi2nu_limit:
            pipeline.set_flag("LARGE_CHISQ", 16)
        
        if "bulge" in pipeline.components:

            xdiff_sq = (target["X_IMAGE"] - target["bulge_xctr"])**2
            ydiff_sq = (target["Y_IMAGE"] - target["bulge_yctr"])**2
            cntr_diff = np.sqrt(xdiff_sq + ydiff_sq)

            if cntr_diff > pipeline.center_deviation_limit:
                pipeline.set_flag("FAKE_CNTR", 17)

            if abs(target["MAG_AUTO"] - pipeline.NMag - target["bulge_mag"]) < 0.5 or abs(target["MAG_AUTO"] + pipeline.NMag - target["bulge_mag"]) < 0.5: 

                pipeline.set_flag("IB_AT_LIMIT", 18)

            if abs(target["bulge_n"] - pipeline.UN) < 0.5 or abs(target["bulge_n"] - pipeline.LN) < 0.1:

                pipeline.set_flag("N_AT_LIMIT", 19)

            if abs(target["bulge_Re"] -  target["FLUX_RADIUS"] * pipeline.NRadius) < 1 or target["bulge_Re"] < 0.5:

                pipeline.set_flag("RE_AT_LIMIT", 20)

        if "disk" in pipeline.components:

            xdiff_sq = (target["X_IMAGE"] - target["disk_xctr"])**2
            ydiff_sq = (target["Y_IMAGE"] - target["disk_yctr"])**2
            cntr_diff = np.sqrt(xdiff_sq + ydiff_sq)

            if cntr_diff > pipeline.center_deviation_limit:
                pipeline.set_flag("FAKE_CNTR", 17)

            if abs(target["MAG_AUTO"] - pipeline.NMag - target["disk_mag"]) < 0.5 or abs(target["MAG_AUTO"] + pipeline.NMag - target["disk_mag"]) < 0.5: 

                pipeline.set_flag("ID_AT_LIMIT", 21)

            if abs(target["disk_Re"] -  target["FLUX_RADIUS"] * pipeline.NRadius) < 1 or target["disk_Re"] < 0.5:

                pipeline.set_flag("RD_AT_LIMIT", 22)

        if "bar" in  pipeline.components:

            xdiff_sq = (target["X_IMAGE"] - target["bar_xctr"])**2
            ydiff_sq = (target["Y_IMAGE"] - target["bar_yctr"])**2
            cntr_diff = np.sqrt(xdiff_sq + ydiff_sq)

            if cntr_diff > pipeline.center_deviation_limit:
                pipeline.set_flag("FAKE_CNTR", 17)

            if abs(target["MAG_AUTO"] - pipeline.NMag - target["bar_mag"]) < 0.5 or abs(target["MAG_AUTO"] + pipeline.NMag - target["bar_mag"]) < 0.5: 

                pipeline.set_flag("IBAR_AT_LIMIT", 23)

            
            if abs(target["bar_n"] - 2.5) < 0.1 or abs(target["bar_n"] - 0.1) < 0.1:
                pipeline.set_flag("NBAR_AT_LIMIT", 24)
            
            if abs(target["bar_Re"] -  target["FLUX_RADIUS"] * pipeline.NRadius) < 1 or target["bar_Re"] < 0.5:

                pipeline.set_flag("RBAR_AT_LIMIT", 25)


    def critical_errors(self, pipe):
        if "SMALL_IMAGE_SIZE" in pipe.flags:
            pipe.flags = {"SMALL_IMAGE_SIZE": 0}
        if "GALFIT_FAILED" in pipe.flags:
            pipe.flags = {"GALFIT_FAILED": 1}
        if 'INVALID_GALFIT_PARAMS' in pipe.flags:
            pipe.flags = {"INVALID_GALFIT_PARAMS": 2}


    def sub_run(self, pipe, obj):

        pipe.check_flags()

        galaxies = pipe.process_target(obj)

        pipe.check_radec()

        target = galaxies['target']
        neighbours = galaxies['neighbours']

        #NEIGHBOUR FLAG
        pipe.check_neighbours(neighbours.shape[0])

        
        #print(pipe.imagefile, pipe.whtfile)
        #sys.exit()
        #print(galaxies)
        target, neighbours = pipe.transform_to_cutout(galaxies['target'],
                                                   galaxies['neighbours'])
        #print('target', target)
        #print('neighbours', neighbours)
        #print('galaxies', galaxies)

        pipe.generate_target_images(galaxies)

        #print('pipe.flags', pipe.flags)

        #PSF CLASS
        psf_pipe = PSFPipeline("config.ini")
        psf_pipe.process_target(target)

        target.update(psf_pipe.result)

        #print('target', target)

        #MASK CLASS
        mask_gen = MaskGenerator("config.ini", target, neighbours)
        mask_gen.run()

        #GALFIT CONFIG and RUN CLASS
        gcr = GalfitConfigRunFunc("config.ini")
        gcr.write_config(target, neighbours)

        if pipe.run_galfit:
            gcr.GalfitRun()
        #CHECK WHETHER FIT.LOG EXISTS
        fitfile = pipe.check_file("fit.log")

        #CASGM CLASSS
        try:
            casgm_pipe = CASGMPipeline()
            casgm_pipe.compute_CASGM(target)
            casgm_dict = casgm_pipe.result
        except:
            self.set_flag("CASGM_FAILED", 14)
            keys = ['R20', 'R50', 'R80', 'R90', 'C', 'A', 'S', 'C_err',
             'A_err', 'S_err', 'gini', 'gini_err', 'm20', 'm20_err']
            casgm_dict = {}
            for key in keys:
                casgm_dict[key] = 9999

        target.update(casgm_dict)

        #PARSE GALFIT OUTPUT FILE
        g = GetOutputParams("config.ini", pipe)

        g.parse_galfit("fit.log", gcr.components, target['Z'])
        g.flatten(gcr.components)

        target.update(g.result)

        self.flags_galfit(target, pipe)
        summary_flag = sum(2**v for v in pipe.flags.values())
        target["FLAGS"] = summary_flag

        if os.path.exists("fit.log"):
            os.remove("fit.log")

        
        #print(target)
        #SAVE CSV FORMAT
        galaxies = {}
        galaxies["target"] = target
        galaxies["neighbours"] = g.neighbours

        #print(galaxies["target"])

        
        #PLOTTING IMAGES AND SURFACE BRIGHTNESS PROFILE

        #plotter = PlotFunc(f"O_{target["NAME"]}.fits",
        #                      f"EM_{target["NAME"]}.fits")
        #plotter.plot_summary(f"P_{target["NAME"]}.png")


        pd.DataFrame([target]).to_csv(
                                      "results.csv",
                                      mode="a",
                                  header=not os.path.exists("results.csv"),
                                      index=False)

        neighbours = {'Name': target['NAME'], **g.neighbours} 
        pd.DataFrame([neighbours]).to_csv(
                                      "results_neighbours.csv",
                                      mode="a",
                                      header=False,
                                      index=False)


        #generate_galaxy_report(galaxies,
        #                       output_file=f"R_{target["NAME"]}.html",
        #                       image_path=f"P_{target["NAME"]}.png")

        #wcsv = WriteCSV("config.ini")
        #wcsv.writeparams(target)
        #sys.exit()


    def run_cutout_model(self):

        pipe = GalaxyPipeline("config.ini")

        
        pipe.load_obj_catalog()
        obj_catalog = pipe.obj_catalog

        #print(obj_catalog)
        
        for i, obj in obj_catalog.iterrows():

            folder, filename = os.path.split(pipe.imagefile)
            pipe.imagefile =  os.path.join(folder, obj['GIMG'])
            pipe.whtfile =  os.path.join(folder, obj['WIMG'])
            #print(obj['GIMG'], obj['WIMG'])
            #obj['RA'] = 9999
            #obj['DEC'] = 9999
            obj.drop(['GIMG', 'WIMG'], inplace=True)

            pipe.run_sextractor()          # run once

            pipe.read_sex_catalog()
            #galaxies = pipe.process_target(obj)
            #pipe.check_radec()

            try:

                self.sub_run(pipe, obj)

            except PipelineCriticalError as e:
                print("Caught:", e.info["error"], e.info["issue"])
                msg = f"{e.info["error"]}_{e.info["issue"]}"
                print(msg)
                self.critical_errors(pipe)

                continue   # go to next iteration:
            
        #print('pipe.flags', pipe.flags)
 


    def run(self):

        pipe = GalaxyPipeline("config.ini")
        
        pipe.run_sextractor()          # run once

        pipe.read_sex_catalog()

        pipe.load_obj_catalog()

        obj_catalog = pipe.obj_catalog
        for i, obj in obj_catalog.iterrows():
            #print(obj)


            try:

                self.sub_run(pipe, obj)

            except PipelineCriticalError as e:
                print("Caught:", e.info["error"], e.info["issue"])
                msg = f"{e.info["error"]}_{e.info["issue"]}"
                print(msg)
                self.critical_errors(pipe)

                continue   # go to next iteration:

        #print('pipe.flags', pipe.flags)


    def sub_detailed(self, pipe, obj):

        pipe.check_flags()

        galaxies = pipe.process_target(obj)

        pipe.check_radec()

        #print('galaxies 2', galaxies)

        target = galaxies['target']
        neighbours = galaxies['neighbours']

        #print('target 0', target)
        #NEIGHBOUR FLAG
        pipe.check_neighbours(neighbours.shape[0])


        #print(pipe.imagefile, pipe.whtfile)
        #sys.exit()
        #print(galaxies)
        target, neighbours = pipe.transform_to_cutout(galaxies['target'],
                                                   galaxies['neighbours'])       
        #print('target 1', target)
        #print('neighbours', neighbours)
        #print('galaxies', galaxies)

        pipe.generate_target_images(galaxies)

        #print('pipe.flags', pipe.flags)

        #PSF CLASS
        psf_pipe = PSFPipeline("config.ini")
        psf_pipe.process_target(target)

        target.update(psf_pipe.result)

        #print('target', target)

        #MASK CLASS
        mask_gen = MaskGenerator("config.ini", target, neighbours)
        mask_gen.run()

        #GALFIT CONFIG and RUN CLASS
        gcr = GalfitConfigRunFunc("config.ini")
        gcr.write_config(target, neighbours)

        gd = GalfitDetailed('config.ini', gcr.galfit_conf,
                            target, neighbours)
        gd.detailed(self.number_random)

        #except PipelineCriticalError as e:
        #    print("Caught:", e.info["error"], e.info["issue"])
        #    msg = f"{e.info["error"]}_{e.info["issue"]}"
        #    print(msg)
        #    self.critical_errors(pipe)

        #    continue   # go to next iteration:

    def run_detailed(self):


        pipe = GalaxyPipeline("config.ini")


        pipe.run_sextractor()          # run once

        pipe.read_sex_catalog()

        pipe.load_obj_catalog()

        obj_catalog = pipe.obj_catalog
        for i, obj in obj_catalog.iterrows():
            #print(obj)
            
            self.sub_detailed(pipe, obj)


            #except PipelineCriticalError as e:
            #    print("Caught:", e.info["error"], e.info["issue"])
            #    msg = f"{e.info["error"]}_{e.info["issue"]}"
            #    print(msg)
            #    self.critical_errors(pipe)

            #    continue   # go to next iteration:



    def run_cutout_detailed(self):

        pipe = GalaxyPipeline("config.ini")


        pipe.load_obj_catalog()
        obj_catalog = pipe.obj_catalog

        #print(obj_catalog)

        for i, obj in obj_catalog.iterrows():

            folder, filename = os.path.split(pipe.imagefile)
            pipe.imagefile =  os.path.join(folder, obj['GIMG'])
            pipe.whtfile =  os.path.join(folder, obj['WIMG'])
            #print(obj['GIMG'], obj['WIMG'])
            #obj['RA'] = 9999
            #obj['DEC'] = 9999
            obj.drop(['GIMG', 'WIMG'], inplace=True)

            pipe.run_sextractor()          # run once

            pipe.read_sex_catalog()
            #galaxies = pipe.process_target(obj)
            #pipe.check_radec()
            #print('galaxies 1', galaxies)

            #try:

            self.sub_detailed(pipe, obj)

            #except PipelineCriticalError as e:


if __name__ == '__main__':
    pipe = GalaxyPipeline("config.ini")

    pipe.run_sextractor()          # run once

    pipe.read_sex_catalog()

    pipe.load_obj_catalog()

    obj_catalog = pipe.obj_catalog

    for i, obj in obj_catalog.iterrows():
        print(obj)
        pymorph(obj)
