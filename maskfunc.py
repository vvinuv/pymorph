import numpy as np
from astropy.io import fits
import configparser


class MaskGenerator:

    def __init__(self, config_file, target, neighbours):

        config = configparser.ConfigParser()
        config.read(config_file)

        self.mask_reg = config.getfloat('mask', 'mask_reg')
        self.thresh_area = config.getfloat('mask', 'thresh_area')
        self.threshold = config.getfloat('mask', 'threshold')

        self.target = target
        self.neighbours = neighbours
        self.image_size = target['IMAGE_SIZE']
        self.mask_image = np.zeros(
            (self.image_size, self.image_size), dtype=np.uint8
        )
        self.mask_file = f"M_{target['NAME']}.fits"
        self.elli_mask_file = f"EM_{target['NAME']}.fits"


    # ---------- ELLIPSE ----------
    def create_ellipse_mask(self, x0, y0, a, b, theta):
        """
        Create a single ellipse mask
        """
        ny, nx = self.mask_image.shape
        y, x = np.mgrid[0:ny, 0:nx]

        # shift
        x_shift = x - x0
        y_shift = y - y0

        # rotation
        theta_rad = np.deg2rad(theta)
        cos_t = np.cos(theta_rad)
        sin_t = np.sin(theta_rad)

        x_rot = x_shift * cos_t + y_shift * sin_t
        y_rot = -x_shift * sin_t + y_shift * cos_t

        # ellipse
        return (x_rot**2 / a**2 + y_rot**2 / b**2) <= 1


    def generate(self):
        """
        Generate two masks:
        1. conditional_mask → based on conditions
        2. full_mask → all neighbours
        """

        target = self.target

        # --- initialize two masks ---
        conditional_mask = np.zeros_like(self.mask_image)
        full_mask = np.zeros_like(self.mask_image)

        a_t = target["A_IMAGE"]
        b_t = a_t / target["ELONGATION"]
        area_t = np.pi * a_t * b_t

        for _, neigh in self.neighbours.iterrows():

            a_n = neigh["A_IMAGE"]
            b_n = a_n / neigh["ELONGATION"]
            area_n = np.pi * a_n * b_n

            # distance
            d = np.sqrt(
                (target["X_IMAGE"] - neigh["X_IMAGE"])**2 +
                (target["Y_IMAGE"] - neigh["Y_IMAGE"])**2
            )

            cond1 = area_n < self.thresh_area * area_t
            cond2 = d > self.threshold * (a_t + a_n)

            # --- scale ellipse ---
            a = self.mask_reg * a_n
            b = self.mask_reg * b_n

            ellipse = self.create_ellipse_mask(
                neigh["X_IMAGE"],
                neigh["Y_IMAGE"],
                a,
                b,
                neigh["THETA_IMAGE"]
            )

            # --- FULL MASK (no condition) ---
            full_mask[ellipse] = 1

            # --- CONDITIONAL MASK ---
            if cond1 or cond2:
                conditional_mask[ellipse] = 1

        # store in class
        self.mask_image = conditional_mask
        self.full_mask = full_mask


    # ---------- SAVE ----------
    def save(self):
        """
        Save mask to FITS file
        """
        hdu = fits.PrimaryHDU(self.mask_image)
        hdu.writeto(self.mask_file, overwrite=True)

        hdu = fits.PrimaryHDU(self.full_mask) 
        hdu.writeto(self.elli_mask_file, overwrite=True)

    # ---------- FULL PIPELINE ----------
    def run(self):
        """
        Generate + save mask
        """
        self.generate()
        self.save()
        return self.mask_file
