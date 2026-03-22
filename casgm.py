import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter, rotate, shift


class CASGMPipeline:

    # ------------------------------------------------
    # MAIN FUNCTION (UNCHANGED NAME)
    # ------------------------------------------------
    def compute_CASGM(self, fstring, x0, y0, flux_radius, elongation, theta_image,
                     r_frac=2.0, smooth_sigma=2.0, n_iter=100):

        hdul = fits.open(f'I_{fstring}.fits')
        image = hdul[0].data
        hdul.close()

        hdul = fits.open(f'M_{fstring}.fits')
        mask = hdul[0].data
        hdul.close()

        q = 1.0 / elongation

        # ---------- ELLIPTICAL RADIUS ----------
        def elliptical_radius():
            y, x = np.indices(image.shape)
            x = x - x0
            y = y - y0

            theta = np.deg2rad(theta_image)
            cos_t, sin_t = np.cos(theta), np.sin(theta)

            x_rot = x * cos_t + y * sin_t
            y_rot = -x * sin_t + y * cos_t

            return np.sqrt(x_rot**2 + (y_rot**2) / (q**2))

        r_ell = elliptical_radius()

        mask[r_ell > r_frac * flux_radius] = 1
        image = image * (~mask.astype(bool))

        # ---------- C ----------
        def compute_radii(img):
            r_flat = r_ell.flatten()
            flux = img.flatten()

            valid = np.isfinite(flux) & (flux > 0)
            r_flat = r_flat[valid]
            flux = flux[valid]

            idx = np.argsort(r_flat)
            r_sorted = r_flat[idx]
            flux_sorted = flux[idx]

            cumflux = np.cumsum(flux_sorted)
            frac = cumflux / cumflux[-1]

            def get_r(f):
                i = np.searchsorted(frac, f)
                return r_sorted[i]

            return get_r(0.2), get_r(0.5), get_r(0.8), get_r(0.9)

        R20, R50, R80, R90 = compute_radii(image)
        C = R80 / R20

        # ---------- A ----------
        def rotate_about_center(image, x0, y0):
            ny, nx = image.shape
            cx, cy = nx / 2, ny / 2

            shift_x = cx - x0
            shift_y = cy - y0

            shifted = shift(image, (shift_y, shift_x), order=1)
            rotated = rotate(shifted, 180, reshape=False, order=1)
            unshifted = shift(rotated, (-shift_y, -shift_x), order=1)

            return unshifted

        def asymmetry(image, x0, y0):
            rot = rotate_about_center(image, x0, y0)
            num = np.sum(np.abs(image - rot))
            den = np.sum(np.abs(image))
            return num / den if den != 0 else np.nan

        def asymmetry_minimization(image, x0, y0, R50):
            delta = 0.01 * R50

            points = [
                (x0, y0),
                (x0 - delta, y0 - delta),
                (x0 - delta, y0 + delta),
                (x0 + delta, y0 - delta),
                (x0 + delta, y0 + delta),
                (x0 - delta, y0),
                (x0 + delta, y0),
                (x0, y0 - delta),
                (x0, y0 + delta),
            ]

            A_values = []

            for (x, y) in points:
                try:
                    A_values.append(asymmetry(image, x, y))
                except:
                    A_values.append(np.nan)

            A_values = np.array(A_values)

            A_min = np.nanmin(A_values)
            best_idx = np.nanargmin(A_values)

            return A_min, points[best_idx]

        A, best_center = asymmetry_minimization(image, x0, y0, R50)

        # ---------- S ----------
        smooth = gaussian_filter(image, sigma=smooth_sigma)
        residual = image - smooth

        inner_mask = r_ell < (0.3 * R50)
        residual[inner_mask] = 0

        S = np.sum(np.abs(residual)) / np.sum(np.abs(image))

        # ---------- ERRORS ----------
        noise = np.std(image[mask == 0])

        C_list, A_list, S_list = [], [], []

        for _ in range(n_iter):

            noisy = image + np.random.normal(0, noise, image.shape)

            try:
                r20, r50, r80, r90 = compute_radii(noisy)
                C_list.append(r80 / r20)

                rot = rotate_about_center(noisy, best_center[0], best_center[1])
                A_list.append(np.sum(np.abs(noisy - rot)) / np.sum(np.abs(noisy)))

                sm = gaussian_filter(noisy, sigma=smooth_sigma)
                res = noisy - sm
                res[inner_mask] = 0
                S_list.append(np.sum(np.abs(res)) / np.sum(np.abs(noisy)))

            except:
                continue

        result_CAS = {
            "R20": R20, "R50": R50, "R80": R80, "R90": R90,
            "C": C, "A": A, "S": S,
            "C_err": np.std(C_list),
            "A_err": np.std(A_list),
            "S_err": np.std(S_list),
        }

        result_GM = self.compute_gini_m20_with_error(image)

        result = {}
        result.update(result_CAS)
        result.update(result_GM)

        return result

    # ------------------------------------------------
    def compute_gini(self, image):

        pixels = image.flatten()
        pixels = pixels[pixels > 0]

        if len(pixels) < 2:
            return np.nan

        pixels = np.sort(pixels)
        n = len(pixels)

        index = np.arange(1, n + 1)
        mean_flux = np.mean(pixels)

        return (1 / (mean_flux * n * (n - 1))) * \
               np.sum((2 * index - n - 1) * pixels)

    # ------------------------------------------------
    def compute_m20(self, image):

        ny, nx = image.shape
        y, x = np.mgrid[0:ny, 0:nx]

        flux = image.copy()
        flux[flux < 0] = 0

        total_flux = np.sum(flux)
        if total_flux <= 0:
            return np.nan

        x_c = np.sum(x * flux) / total_flux
        y_c = np.sum(y * flux) / total_flux

        r2 = (x - x_c)**2 + (y - y_c)**2
        M_tot = np.sum(flux * r2)

        pixels = flux.flatten()
        r2_flat = r2.flatten()

        order = np.argsort(pixels)[::-1]

        pixels_sorted = pixels[order]
        r2_sorted = r2_flat[order]

        flux_cumsum = np.cumsum(pixels_sorted)

        threshold = 0.2 * total_flux

        idx = np.where(flux_cumsum <= threshold)[0]

        if len(idx) < len(pixels_sorted):
            idx = np.append(idx, len(idx))

        M20 = np.sum(pixels_sorted[idx] * r2_sorted[idx])

        return np.log10(M20 / M_tot)

    # ------------------------------------------------
    def compute_gini_m20_with_error(self, image, n_iter=100):

        noise_sigma = np.std(image[image != 0])

        gini_vals, m20_vals = [], []

        for _ in range(n_iter):

            noisy = image + np.random.normal(0, noise_sigma, image.shape)
            noisy[noisy < 0] = 0

            g = self.compute_gini(noisy)
            m = self.compute_m20(noisy)

            if not np.isnan(g):
                gini_vals.append(g)

            if not np.isnan(m):
                m20_vals.append(m)

        return {
            "gini": np.mean(gini_vals),
            "gini_err": np.std(gini_vals),
            "m20": np.mean(m20_vals),
            "m20_err": np.std(m20_vals)
        }
