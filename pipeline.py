import os
import configparser
import subprocess
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.nddata import Cutout2D
from astropy.wcs import WCS


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
            "ISO0"
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
        neighbours = self.sex_catalog.loc[mask, self.neighbour_cols].copy()

        # corresponding separations
        sep_subset = sep.arcsec[mask]

        if len(neighbours) == 0:
            return None, neighbours

        # find closest object (object1)
        min_idx = sep_subset.argmin()

        target_sex = neighbours.iloc[min_idx]

        # remove object1 from neighbours
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

        neighbours = self.sex_catalog.loc[mask, self.neighbour_cols].copy()

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


    def process_target(self, target, radius_arcmic=0.5):

        ra = target.get("ra", np.nan)
        dec = target.get("dec", np.nan)

        x = target.get("ximg", np.nan)
        y = target.get("yimg", np.nan)

        position = False

        if not pd.isna(ra) and not pd.isna(dec):

            position = self.position_in_sexcat(ra, dec)

            target_sex, neighbours = self.neighbours_radec(ra, dec, 
                                                           radius_arcmic=0.5)

        else:

            target_sex, neighbours = self.neighbours_xy(x, y)

        #print(target_sex)
        target.loc["position"] = position
        target = target.to_dict()
        target.update(target_sex)
        galaxies = dict()
        galaxies['target'] = target
        galaxies['neighbours'] = neighbours

        return galaxies


    def generate_target_images(self, galaxies, size=120):

        img = fits.open(self.imagefile)
        img_data = img[0].data
        wcs = WCS(img[0].header)

        wht = fits.open(self.whtfile)
        wht_data = wht[0].data

        target = galaxies.get("target")
        gal_id = target.get("gal_id")

        ra = target.get("ra")
        dec = target.get("dec")

        x = target.get("ximg") or target.get("X_IMAGE")
        y = target.get("yimg") or target.get("Y_IMAGE")

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


    def generate_elliptical_mask(df, x0=1000, y0=1000, size=200):

        half = size // 2

        # create empty mask
        mask = np.zeros((size, size), dtype=np.uint8)

        # grid of pixel coordinates
        yy, xx = np.indices((size, size))

        # origin of cutout in global coordinates
        x_origin = x0 - half
        y_origin = y0 - half

        for _, obj in df.iterrows():

            # global coordinates
            xg = obj["X_IMAGE"]
            yg = obj["Y_IMAGE"]

            # convert to local (cutout) coordinates
            xc = xg - x_origin
            yc = yg - y_origin

            # ellipse parameters
            a = obj["A_IMAGE"]
            elong = obj["ELONGATION"]
            b = a / elong

            theta = np.deg2rad(obj["THETA_IMAGE"])

            # shift grid
            x_shift = xx - xc
            y_shift = yy - yc

            # rotate coordinates
            cos_t = np.cos(theta)
            sin_t = np.sin(theta)

            x_rot = x_shift * cos_t + y_shift * sin_t
            y_rot = -x_shift * sin_t + y_shift * cos_t

            # ellipse equation
            ellipse = (x_rot / a) ** 2 + (y_rot / b) ** 2 <= 1

            # set mask pixels
            mask[ellipse] = 1

        return mask
if __name__ == '__main__':
    pipe = GalaxyPipeline("config.ini")

    pipe.run_sextractor()          # run once

    pipe.read_sex_catalog()

    pipe.load_obj_catalog()

    obj_catalog = pipe.obj_catalog
    print(obj_catalog.iloc[0])
    galaxies = pipe.process_target(obj_catalog.iloc[0])

    print('galaxies', galaxies)

    pipe.generate_target_images(galaxies, size=150)
