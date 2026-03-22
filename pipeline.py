import os
import sys
import configparser
import subprocess
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from maskfunc import MaskGenerator
from galfit_config_run import GalfitConfigRunFunc
from casgm import CASGMPipeline


class GalaxyPipeline:

    def __init__(self, config_file):

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

        self.mag_zero = config.getfloat('mag', 'mag_zero')
        self.maglim = config.get('mag', 'maglim').split(',')
        self.maglim = [float(mlim) for mlim in self.maglim]


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

    # ---------------------------------------------------------
    # Run SExtractor
    # ---------------------------------------------------------
    def run_sextractor(self, config_file="SEx/default.sex"):

        output_cat = f"{self.rootname}_sex.cat"
        sex_out_params = "SEx/default.param"

        cmd = [
            "sex",
            self.imagefile,
            "-WEIGHT_IMAGE", self.whtfile,
            "-CATALOG_NAME", output_cat,
            "-PARAMETERS_NAME", sex_out_params,
            "-c", config_file
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

        sep = target.separation(sex_coords)

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

        # mask within radius
        mask = sep.arcsec < radius_arcsec

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

        return target_sex.to_dict(), neighbours.reset_index(drop=True)


    # ------------------------------------------------
    # Process one target object and its neighbours
    # ------------------------------------------------


    def process_target(self, target, size, radius_arcmic=0.5):

        ra = target.get("RA", np.nan)
        dec = target.get("DEC", np.nan)

        x = target.get("XIMG", np.nan)
        y = target.get("YIMG", np.nan)

        position = False

        if not pd.isna(ra) and not pd.isna(dec):

            position = self.position_in_sexcat(ra, dec)

            target_sex, neighbours = self.neighbours_radec(ra, dec,
                                                           radius_arcmic=0.5)

        else:

            target_sex, neighbours = self.neighbours_xy(x, y)



        #print(target_sex)
        target = target.to_dict()
        target['IMAGE_SIZE'] = size
        target['GALFIT_ANGLE'] = target_sex["THETA_IMAGE"] - 90 
        target["POSITION"] = position
        target["ROOTNAME"] = self.rootname
        target["NAME"] = f'{self.rootname}_{target['GAL_ID']}'
        target["MAG_ZERO"] = self.mag_zero
        target.update(target_sex)

        neighbours['GALFIT_ANGLE'] = neighbours["THETA_IMAGE"] - 90

        galaxies = dict()
        galaxies['target'] = target
        galaxies['neighbours'] = neighbours

        return galaxies


    def generate_target_images(self, galaxies, size):

        img = fits.open(self.imagefile)
        img_data = img[0].data
        wcs = WCS(img[0].header)

        wht = fits.open(self.whtfile)
        wht_data = wht[0].data

        target = galaxies.get("target")
        gal_id = target.get("GAL_ID")

        ra = target.get("RA")
        dec = target.get("DEC")

        x = target.get("XIMG") or target.get("X_IMAGE")
        y = target.get("YIMG") or target.get("Y_IMAGE")

        # determine pixel position
        if ra is not None and dec is not None:
            x, y = wcs.world_to_pixel_values(ra, dec)

        position = (x, y)

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


    def transform_to_cutout(self, target, neighbours):

        size = target['IMAGE_SIZE']
        half = size // 2

        # Target center
        x0 = target["X_IMAGE"]
        y0 = target["Y_IMAGE"]

        print(x0)
        # Shift
        dx = x0 - half
        dy = y0 - half
        print(dx)
        
        # --- Transform target ---
        target["X_IMAGE"] = target["X_IMAGE"] - dx
        target["Y_IMAGE"] = target["Y_IMAGE"] - dy

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



if __name__ == '__main__':
    size = 200
    pipe = GalaxyPipeline("config.ini")

    pipe.run_sextractor()          # run once

    pipe.read_sex_catalog()

    pipe.load_obj_catalog()

    obj_catalog = pipe.obj_catalog

    print(obj_catalog.iloc[0])
    galaxies = pipe.process_target(obj_catalog.iloc[0], size)

    target = galaxies['target']
    neighbours = galaxies['neighbours']

    #print(galaxies)
    target, neighbours = pipe.transform_to_cutout(galaxies['target'], 
                                                  galaxies['neighbours'])
    print(target, neighbours)
    #print('galaxies', galaxies)
    
    pipe.generate_target_images(galaxies, size)

    mask_gen = MaskGenerator("config.ini", target, neighbours)
    mask_gen.run()

    gcr = GalfitConfigRunFunc("config.ini")
    gcr.write_config(target, neighbours)

    #gcr.GalfitRun()

    casgm_pipe = CASGMPipeline()

    casgm_result = casgm_pipe.compute_CASGM(
                                     target['NAME'],
                                     target['X_IMAGE'],
                                     target['Y_IMAGE'],
                                     target['FLUX_RADIUS'],
                                     target['ELONGATION'],
                                     target['THETA_IMAGE'])
    print(result)
    sys.exit()



