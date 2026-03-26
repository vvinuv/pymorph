import os
import numpy as np
import configparser
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u


class PSFPipeline:

    def __init__(self, config_file):

        config = configparser.ConfigParser()
        config.read(config_file)

        self.DATADIR = config["imagecata"].get("datadir", None)
        self.obj_cata = config["imagecata"].get("obj_cata", None)
        self.psflist = config["psf"].get("psflist", None)
        


    # --------------------------------------------
    # Convert HMS/DMS → degrees
    # --------------------------------------------
    def hms_dms_to_deg(self, ra1, ra2, ra3, dec1, dec2, dec3):

        ra = 15 * (ra1 + ra2/60 + ra3/3600)

        sign = -1 if dec1 < 0 else 1
        dec = sign * (abs(dec1) + dec2/60 + dec3/3600)

        return ra, dec


    # --------------------------------------------
    # Read PSF list from config
    # --------------------------------------------
    def get_psf_list(self):

        if self.psflist is None:
            return []

        # Case 1: @filelist
        if self.psflist.startswith("@"):

            fname = os.path.join(self.DATADIR, self.psflist[1:])

            with open(fname) as f:
                psf_files = [line.strip() for line in f if line.strip()]

        # Case 2: comma-separated
        else:
            psf_files = [p.strip() for p in self.psflist.split(",")]

        print(psf_files)
        return psf_files


    # --------------------------------------------
    # Extract RA DEC from PSF filename
    # --------------------------------------------
    def parse_psf_filename(self, filename):

        name = filename.split("/")[-1]

        ra1 = float(str(name)[4:6])
        ra2 = float(str(name)[6:8])
        ra3 = float(str(name)[8:10]) + float(str(name)[10]) / 10.0
        dec1 = float(str(name)[11:-10])
        dec2 = float(str(name)[-10:-8])
        dec3 = float(str(name)[-8:-6]) + float(str(name)[-6]) / 10.0

        ra, dec = self.hms_dms_to_deg(
            ra1, ra2, ra3,
            dec1, dec2, dec3
        )

        return ra, dec


    # --------------------------------------------
    # Compute nearest PSF in arcsec
    # --------------------------------------------
    def find_nearest_psf(self, ra, dec):

        
        psf_files = self.get_psf_list()

        if len(psf_files) == 0:
            return None, None

        target = SkyCoord(ra*u.deg, dec*u.deg)

        min_dist = np.inf
        best_psf = None
        best_coord = None
        #print(ra, dec)
        for psf in psf_files:
            try:
                psf = os.path.join(self.DATADIR, psf)
                psf_ra, psf_dec = self.parse_psf_filename(psf)

                psf_coord = SkyCoord(psf_ra*u.deg, psf_dec*u.deg)

                dist = target.separation(psf_coord).arcsec

                if dist < min_dist:
                    min_dist = dist
                    best_psf = psf
                    best_coord = (psf_ra, psf_dec)

            except Exception:
                continue

        return best_psf, min_dist, best_coord


    # --------------------------------------------
    # Update PSF header and save
    # --------------------------------------------
    def update_psf_header(self, psf_file, ra_psf, dec_psf):

        hdu = fits.open(psf_file)

        hdu[0].header["RA_PSF"] = ra_psf
        hdu[0].header["DEC_PSF"] = dec_psf

        hdu.writeto(psf_file, overwrite=True)

        return psf_file


    # --------------------------------------------
    # Main function per target
    # --------------------------------------------
    def process_target(self, target):

        # Case 1: PSF column exists → direct PSF
        psf_file = target.get("PSF", None)

        ra = target.get("RA")
        dec = target.get("DEC")
        if (psf_file and isinstance(psf_file, str)) or ra is None:

            #ra_psf, dec_psf = self.parse_psf_filename(psf_file)

            #output = self.update_psf_header(
            #    psf_file,
            #    ra_psf,
            #    dec_psf,
            #    f"updated_{psf_file}"
            #)

            self.result = {
                           "PSF": psf_file,
                           "distance_psf_arcsec": 0.0
                           }

        else:
            # Case 2: find nearest PSF
            #ra = target.get("RA")
            #dec = target.get("DEC")
            psf_file, dist, (ra_psf, dec_psf) = self.find_nearest_psf(ra, dec)

            #print(psf_file, dist, ra_psf, dec_psf)
            self.update_psf_header(
                psf_file,
                ra_psf,
                dec_psf
            )
            
            self.result = {"PSF": psf_file,
                           "distance_psf_arcsec": int(dist)}

if __name__=='__main__':
    psf_pipe = PSFPipeline("config.ini")

    target = {
        "gal_id": 101,
        "ra": 221.84,
        "dec": 8.64
    }

    result = psf_pipe.process_target(target)

    print(result)
