import numpy as np
from astropy.io import fits
import configparser


class MaskGenerator:

    def __init__(self, target, neighbours):

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


    # ---------- GENERATE MASK ----------
    def generate(self):
        """
        Generate mask image from neighbours
        """
        target = self.target

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

            cond1 = area_n < 0.3 * area_t
            cond2 = d > 2 * (a_t + a_n)

            if cond1 or cond2:
                # scale mask
                a = 3 * a_n
                b = 3 * b_n

                ellipse = self.create_ellipse_mask(
                    neigh["X_IMAGE"],
                    neigh["Y_IMAGE"],
                    a,
                    b,
                    neigh["THETA_IMAGE"]
                )

                self.mask_image[ellipse] = 1

        return self.mask_image


    # ---------- SAVE ----------
    def save(self):
        """
        Save mask to FITS file
        """
        hdu = fits.PrimaryHDU(self.mask_image)
        hdu.writeto(self.mask_file, overwrite=True)


    # ---------- FULL PIPELINE ----------
    def run(self):
        """
        Generate + save mask
        """
        self.generate()
        self.save()
        return self.mask_file
