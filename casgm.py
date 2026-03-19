import numpy as np
import numpy.ma as ma
from scipy.ndimage import rotate, shift, gaussian_filter
from astropy.io import fits
from matplotlib.pylab import plt

def rotate_numpy_bilinear(image, x_shift, y_shift, angle_deg):

    
    angle = np.deg2rad(angle_deg)

    ny, nx = image.shape
    cy, cx = ny // 2, nx // 2

    y, x = np.mgrid[0:ny, 0:nx]

    x_shift = x - cx - x_shift
    y_shift = y - cy - y_shift

    cos_t = np.cos(angle)
    sin_t = np.sin(angle)

    x_rot = cos_t * x_shift + sin_t * y_shift + cx
    y_rot = -sin_t * x_shift + cos_t * y_shift + cy

    x0 = np.floor(x_rot).astype(int)
    x1 = x0 + 1
    y0 = np.floor(y_rot).astype(int)
    y1 = y0 + 1

    x0 = np.clip(x0, 0, nx-1)
    x1 = np.clip(x1, 0, nx-1)
    y0 = np.clip(y0, 0, ny-1)
    y1 = np.clip(y1, 0, ny-1)

    Ia = image[y0, x0]
    Ib = image[y1, x0]
    Ic = image[y0, x1]
    Id = image[y1, x1]

    wa = (x1 - x_rot) * (y1 - y_rot)
    wb = (x1 - x_rot) * (y_rot - y0)
    wc = (x_rot - x0) * (y1 - y_rot)
    wd = (x_rot - x0) * (y_rot - y0)

    return wa*Ia + wb*Ib + wc*Ic + wd*Id


def compute_CASGM(fstring, x0, y0, flux_radius, elongation, theta_image,
                r_frac=2.0, smooth_sigma=2.0, n_iter=100):
    """
    Compute Concentration (C), Asymmetry (A), Smoothness (S)

    Parameters
    ----------
    image : 2D array
    x0, y0 : center
    elongation : a/b
    theta_deg : position angle SExtractor THETA_IMAGE
    r_pet_frac : aperture scaling (Petrosian-like)
    smooth_sigma : Gaussian smoothing for S
    n_iter : Monte Carlo for error

    Returns
    -------
    dict with C, A, S and errors
    """

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
    image *= ~mask #ma.array(image, mask=mask)


    #plt.imshow(image)
    #plt.show()
    # ---------- C (Concentration) ----------
    def compute_radii(img):
        r_flat = r_ell.flatten()
        flux = img.flatten()
        #print(r_flat)
        mask = np.isfinite(flux) & (flux > 0)
        r_flat = r_flat[mask]
        flux = flux[mask]
        #print(r_flat, flux)
        idx = np.argsort(r_flat)
        r_sorted = r_flat[idx]
        #print(r_sorted)
        flux_sorted = flux[idx]
        #print(flux_sorted.shape, flux_sorted[np.isfinite(flux_sorted)].shape)
        #print(flux_sorted)
        cumflux = np.cumsum(flux_sorted)
        frac = cumflux / cumflux[-1]
        
        #print(cumflux)
        def get_r(f):
            i = np.searchsorted(frac, f)
            return r_sorted[i]

        return get_r(0.2), get_r(0.5), get_r(0.8), get_r(0.9)

    def compute_radii_myown(img, mask):
        r_flat = r_ell.flatten()
        flux = img.flatten()
        mask = mask.flatten().astype(bool)
        r_flat = r_flat[~mask]
        flux = flux[~mask]
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
    

    # ---------- A (Asymmetry) ----------

    def rotate_about_center(image, x0, y0):
        ny, nx = image.shape
        cx, cy = nx / 2, ny / 2

        shift_x = cx - x0
        shift_y = cy - y0

        shifted = shift(image, (shift_y, shift_x), order=1, mode='nearest')
        rotated = rotate(shifted, 180, reshape=False, order=1)
        unshifted = shift(rotated, (-shift_y, -shift_x), order=1, mode='nearest')

        return unshifted


    def asymmetry(image, x0, y0):
        rot = rotate_about_center(image, x0, y0)
        num = np.sum(np.abs(image - rot))
        den = np.sum(np.abs(image))

        return num / den if den != 0 else np.nan


    def asymmetry_minimization(image, x0, y0, R50):
        """
        Compute asymmetry at 9 positions and return minimum
        """

        delta = 0.01 * R50

        # 9 test points
        points = [
            (x0, y0),  # center

            # corners
            (x0 - delta, y0 - delta),
            (x0 - delta, y0 + delta),
            (x0 + delta, y0 - delta),
            (x0 + delta, y0 + delta),

            # edge midpoints
            (x0 - delta, y0),
            (x0 + delta, y0),
            (x0, y0 - delta),
            (x0, y0 + delta),
        ]

        A_values = []

        for (x, y) in points:
            try:
                A = asymmetry(image, x, y)
                A_values.append(A)
            except:
                A_values.append(np.nan)

        A_values = np.array(A_values)

        # minimum asymmetry
        A_min = np.nanmin(A_values)

        # best center
        best_idx = np.nanargmin(A_values)
        best_center = points[best_idx]

        return {
            "A_min": A_min,
            "A_all": A_values,
            "best_center": best_center,
            "points": points
        }

    result = asymmetry_minimization(
        image,
        x0, y0, R50
    )

    A = result["A_min"]
    best_center = result["best_center"]

    # ---------- S (Smoothness) ----------
    smooth = gaussian_filter(image, sigma=smooth_sigma)
    residual = image - smooth

    # ignore central region (important!)
    inner_mask = r_ell < (0.3 * R50)
    residual[inner_mask] = 0

    S_num = np.sum(np.abs(residual))
    S_den = np.sum(np.abs(image))
    S = S_num / S_den

    # ---------- ERRORS (Monte Carlo) ----------
    noise = np.std(image[mask == 0])
    #print('noise', noise)
    C_list, A_list, S_list = [], [], []

    for _ in range(n_iter):
        noisy = image + np.random.normal(0, noise, image.shape)

        try:
            r20, r50, r80, r90 = compute_radii(noisy, mask)
            C_list.append(r80 / r20)
            #print(C_list)

            rot = rotate_about_center(image, best_center[0], best_center[1])
            A_list.append(np.sum(np.abs(noisy - rot)) / np.sum(np.abs(noisy)))

            sm = gaussian_filter(noisy, sigma=smooth_sigma)
            res = noisy - sm
            res[inner_mask] = 0
            S_list.append(np.sum(np.abs(res)) / np.sum(np.abs(noisy)))

        except:
            continue

    result_CAS = {
        "R20": R20,
        "R50": R50,
        "R80": R80,
        "R90": R90,
        "C": C,
        "A": A,
        "S": S,
        "C_err": np.std(C_list),
        "A_err": np.std(A_list),
        "S_err": np.std(S_list),
    }
    result_GM = compute_gini_m20_with_error(image)
    result_CASGM = {}
    result_CASGM.update(result_CAS)
    result_CASGM.update(result_GM)

    return result_CASGM


def compute_gini(image):

    pixels = image.flatten()
    pixels = pixels[pixels > 0]

    if len(pixels) < 2:
        return np.nan

    pixels = np.sort(pixels)
    n = len(pixels)

    index = np.arange(1, n + 1)
    mean_flux = np.mean(pixels)

    gini = (1 / (mean_flux * n * (n - 1))) * \
           np.sum((2 * index - n - 1) * pixels)

    return gini

def compute_m20(image):

    ny, nx = image.shape
    y, x = np.mgrid[0:ny, 0:nx]

    flux = image.copy()
    flux[flux < 0] = 0

    total_flux = np.sum(flux)
    if total_flux <= 0:
        return np.nan

    # centroid
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

    # important: include boundary pixel
    if len(idx) < len(pixels_sorted):
        idx = np.append(idx, len(idx))

    M20 = np.sum(pixels_sorted[idx] * r2_sorted[idx])

    return np.log10(M20 / M_tot)

def compute_gini_m20_with_error(image, noise_sigma=0.1, n_iter=100):

    noise_sigma = np.std(image[image != 0])

    gini_vals = []
    m20_vals = []

    for _ in range(n_iter):

        noise = np.random.normal(0, noise_sigma, image.shape)
        img_noisy = image + noise

        img_noisy[img_noisy < 0] = 0

        g = compute_gini(img_noisy)
        m = compute_m20(img_noisy)

        if not np.isnan(g):
            gini_vals.append(g)

        if not np.isnan(m):
            m20_vals.append(m)

    gini_vals = np.array(gini_vals)
    m20_vals = np.array(m20_vals)

    return {
        "G": np.mean(gini_vals),
        "G_err": np.std(gini_vals),
        "M20": np.mean(m20_vals),
        "M20_err": np.std(m20_vals)
    }
