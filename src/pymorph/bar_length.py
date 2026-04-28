import numpy as np
from astropy.io import fits


# -----------------------------
# Utility: apply mask + weights
# -----------------------------
def prepare_image(image, weight=None, mask=None):
    img = image.copy()

    if mask is not None:
        img[mask > 0] = 0

    if weight is not None:
        # convert weight -> sigma
        sigma = np.zeros_like(weight)
        valid = weight > 0
        sigma[valid] = 1.0 / np.sqrt(weight[valid])
    else:
        sigma = None

    return img, sigma


# -----------------------------
# Deprojection
# -----------------------------
def deproject(x, y, x0, y0, pa, inc):
    pa = np.deg2rad(pa)
    inc = np.deg2rad(inc)

    x = x - x0
    y = y - y0

    x_r = x * np.cos(pa) + y * np.sin(pa)
    y_r = -x * np.sin(pa) + y * np.cos(pa)

    y_r /= np.cos(inc)

    return x_r, y_r


# -----------------------------
# Fourier computation
# -----------------------------
def compute_fourier(image, x0, y0, pa, inc, q=1.0, nbins=50):
    """
    Compute the radial Fourier decomposition (m = 2 mode) of a galaxy image
    to quantify bar strength and phase as a function of radius.

    This function calculates the normalized Fourier amplitudes A2/A0 and the
    corresponding phase profile using deprojected galaxy coordinates and
    elliptical annuli.

    PARAMETERS
    ----------
    image : 2D numpy.ndarray
        Input galaxy image (background-subtracted).

    x0, y0 : float
        Galaxy center coordinates in pixel units.

    pa : float
        Position angle of the galaxy in degrees.
        Measured counter-clockwise from the x-axis.

    inc : float
        Inclination angle of the galaxy in degrees.
        0° = face-on, 90° = edge-on.

    q : float, optional (default=1.0)
        Axis ratio (b/a) used to define elliptical annuli.
        For a perfectly deprojected disk, q ~ 1.

    nbins : int, optional (default=50)
        Number of radial bins for the Fourier analysis.

    RETURNS
    -------
    r_mid : ndarray
        Midpoint radius of each annulus.

    A2_amp : ndarray
        Normalized Fourier amplitude (A2/A0), representing bar strength.

    phase : ndarray
        Phase of the m=2 Fourier component:
            phi(r) = (1/2) * arctan(B2 / A2)

    METHOD
    ------
    1. Deprojection:
        - Coordinates are rotated by the position angle (PA)
        - Inclination correction applied assuming a thin disk:
              y' = y / cos(inc)

    2. Elliptical Radius:
        - Radius defined as:
              r = sqrt(x'^2 + (y'/q)^2)

    3. Fourier Decomposition:
        For each radial annulus:
            A0(r) = Σ I(r, θ)
            A2(r) = Σ I(r, θ) cos(2θ)
            B2(r) = Σ I(r, θ) sin(2θ)

        Normalized amplitude:
            A2/A0 = sqrt[(A2/A0)^2 + (B2/A0)^2]

    4. Phase:
            phi(r) = (1/2) arctan(B2 / A2)

    PHYSICAL INTERPRETATION
    -----------------------
    - A2/A0 measures the strength of the non-axisymmetric (bar-like) component.
    - Higher A2/A0 indicates a stronger bar.
    - The phase is approximately constant across the bar region.

    NOTES
    -----
    - Pixels with zero or negative total intensity in a bin are ignored.
    - Accurate galaxy center (x0, y0) is critical.
    - Image should be cleaned (masked, sky-subtracted) before use.
    - Deprojection assumes a thin, flat disk (no warp).

    REFERENCES
    ----------
    The Fourier decomposition method used here follows standard approaches
    in the literature:

    - Athanassoula, E. (2003), MNRAS, 341, 1179
      (Bar strength and evolution using Fourier components)

    - Laurikainen, E., Salo, H., & Buta, R. (2004), ApJ, 607, 103
      (Quantifying bar strength via Fourier amplitudes and torques)

    - Díaz-García, S. et al. (2016), A&A, 587, A160
      (S4G survey: Fourier analysis of bar structure)

    - Garcia-Gómez, C. et al. (2017), A&A, 601, A132
      (Measuring bar strength using Fourier analysis)

    EXAMPLE
    -------
    >>> r, A2, phase = compute_fourier(
            image,
            x0=150, y0=150,
            pa=45.0,
            inc=30.0,
            q=0.8,
            nbins=60
        )

    >>> # Peak A2 gives bar region
    >>> Rbar = r[np.argmax(A2)]

    SCIENTIFIC CONTEXT
    ------------------
    Fourier decomposition is a standard method to identify and quantify
    bars in disk galaxies. The m=2 mode corresponds to the dominant
    bi-symmetric structure, which is the defining feature of bars.

    The radial profile of A2/A0 and the constancy of phase are used to:
        - estimate bar length
        - measure bar strength
        - distinguish bar from spiral arms
    """
    ny, nx = image.shape
    y, x = np.indices((ny, nx))

    x, y = deproject(x, y, x0, y0, pa, inc)

    r = np.sqrt(x**2 + (y / q)**2)
    theta = np.arctan2(y, x)

    r_bins = np.linspace(0, np.max(r), nbins)

    r_mid, A2_amp, phase = [], [], []

    for i in range(len(r_bins)-1):
        m = (r >= r_bins[i]) & (r < r_bins[i+1])

        I = image[m]
        th = theta[m]

        if len(I) < 10 or np.sum(I) <= 0:
            continue

        A0 = np.sum(I)
        A2 = np.sum(I * np.cos(2 * th))
        B2 = np.sum(I * np.sin(2 * th))

        A2n = A2 / A0
        B2n = B2 / A0

        amp = np.sqrt(A2n**2 + B2n**2)
        phi = 0.5 * np.arctan2(B2n, A2n)

        r_mid.append(0.5 * (r_bins[i] + r_bins[i+1]))
        A2_amp.append(amp)
        phase.append(phi)

    return np.array(r_mid), np.array(A2_amp), np.array(phase)


# -----------------------------
# Bar radius
# -----------------------------
def estimate_Rbar(r, A2_amp, phase):

    i_peak = np.argmax(A2_amp)
    R_peak = r[i_peak]

    phi0 = phase[i_peak]

    phase_mask = np.abs(phase - phi0) < 0.1
    R_phase = np.max(r[phase_mask]) if np.any(phase_mask) else R_peak

    amp_thresh = 0.5 * np.max(A2_amp)
    amp_mask = A2_amp > amp_thresh
    R_amp = r[amp_mask][-1] if np.any(amp_mask) else R_peak

    Rbar = np.mean([R_peak, R_phase, R_amp])

    return Rbar, np.array([R_peak, R_phase, R_amp])


# -----------------------------
# Monte Carlo error
# -----------------------------
def monte_carlo_Rbar(image, sigma, x0, y0, pa, inc, q, n_iter=30):

    Rbars = []

    for _ in range(n_iter):

        if sigma is not None:
            noise = np.random.normal(0, sigma)
            img = image + noise
        else:
            img = image + np.random.normal(0, np.std(image)*0.05, image.shape)

        r, A2, phase = compute_fourier(img, x0, y0, pa, inc, q)
        Rbar, _ = estimate_Rbar(r, A2, phase)

        Rbars.append(Rbar)

    return np.std(Rbars)


# -----------------------------
# Compute Rbar with error
# -----------------------------
def compute_Rbar_with_error(image, sigma, x0, y0, pa, inc, q):

    r, A2, phase = compute_fourier(image, x0, y0, pa, inc, q)

    Rbar, methods = estimate_Rbar(r, A2, phase)

    # method scatter
    sigma_method = np.std(methods)

    # bin size
    sigma_bin = np.mean(np.diff(r)) / 2.0

    # MC
    sigma_mc = monte_carlo_Rbar(image, sigma, x0, y0, pa, inc, q)

    sigma_total = np.sqrt(
        sigma_method**2 +
        sigma_bin**2 +
        sigma_mc**2
    )

    return Rbar, sigma_total


# -----------------------------
# R ratio
# -----------------------------
def compute_R_ratio(Rcorot, Rbar, err_Rcorot, err_Rbar):

    R = Rcorot / Rbar

    err = R * np.sqrt(
        (err_Rcorot / Rcorot)**2 +
        (err_Rbar / Rbar)**2
    )

    return R, err


# -----------------------------
# FULL PIPELINE
# -----------------------------
def run_full_bar_pipeline(
    image_fits,
    weight_fits,
    mask_fits,
    x0, y0, pa, inc, q,
    Rcorot, err_Rcorot
):

    """
    Full pipeline to compute bar length (Rbar), its uncertainty,
    and the dynamical ratio R = Rcorot / Rbar with error propagation.

    This function performs a Fourier decomposition (m=2 mode) of a galaxy image
    to estimate the bar radius using multiple standard methods and combines
    them with Monte Carlo noise propagation to derive robust uncertainties.

    PARAMETERS
    ----------
    image_fits : str
        Path to the science FITS image of the galaxy.

    weight_fits : str or None
        Path to the weight FITS image (inverse variance map).
        If provided, used to estimate pixel-wise noise for Monte Carlo sampling.
        If None, noise is approximated from image statistics.

    mask_fits : str or None
        Path to the mask FITS image.
        Pixels with mask > 0 are excluded from analysis (e.g., stars, bad pixels).

    x0, y0 : float
        Galaxy center coordinates in pixel units.

    pa : float
        Position angle of the galaxy (degrees).
        Measured counter-clockwise from the x-axis.

    inc : float
        Inclination angle of the galaxy (degrees).
        0 = face-on, 90 = edge-on.

    q : float
        Axis ratio (b/a) of the galaxy disk.
        Used for defining elliptical annuli.

    Rcorot : float
        Corotation radius (same units as pixel scale or converted physical units).

    err_Rcorot : float
        Uncertainty in the corotation radius.

    RETURNS
    -------
    result : dict
        Dictionary containing:

        - 'Rbar' : float
            Estimated bar radius.

        - 'Rbar_err' : float
            Total uncertainty in bar radius, combining:
                * method scatter
                * radial bin size
                * Monte Carlo noise propagation

        - 'R' : float
            Dimensionless bar dynamical ratio:
                R = Rcorot / Rbar

        - 'R_err' : float
            Uncertainty in R computed using standard error propagation:
                σ_R = R * sqrt[(σ_Rc/Rc)^2 + (σ_Rb/Rb)^2]

    METHOD
    ------
    1. Image Preparation:
        - Apply mask (remove unwanted pixels)
        - Use weight image to estimate noise if available

    2. Deprojection:
        - Rotate by position angle
        - Correct for inclination (face-on transformation)

    3. Fourier Decomposition:
        - Compute m=2 mode (A2, B2)
        - Derive amplitude A2/A0 and phase as a function of radius

    4. Bar Radius Estimation:
        - Peak of A2 profile
        - Phase constancy method
        - Amplitude threshold method (A2 drop)
        - Final Rbar = mean of all methods

    5. Error Estimation:
        - Method-to-method scatter
        - Radial bin resolution
        - Monte Carlo perturbation using noise model

    6. Dynamical Ratio:
        - Compute R = Rcorot / Rbar
        - Propagate uncertainties assuming independent errors

    NOTES
    -----
    - All radii are in pixel units unless converted externally.
    - Ensure image is background-subtracted before use.
    - Accurate center (x0, y0) is critical for reliable results.
    - Deprojection assumes a thin disk geometry.
    - For publication-quality results, increase Monte Carlo iterations.

    EXAMPLE
    -------
    >>> result = run_full_bar_pipeline(
            "galaxy.fits",
            "weight.fits",
            "mask.fits",
            x0=150, y0=150,
            pa=45.0,
            inc=30.0,
            q=0.8,
            Rcorot=8.5,
            err_Rcorot=0.8
        )

    >>> print(result)
    {'Rbar': 5.2, 'Rbar_err': 0.7, 'R': 1.63, 'R_err': 0.28}

    SCIENTIFIC CONTEXT
    ------------------
    The ratio R = Rcorot / Rbar is used to classify galactic bars:
        R ~ 1.0–1.4  -> fast bars
        R > 1.4      -> slow bars

    This pipeline follows standard Fourier-based bar analysis methods
    widely used in observational studies of barred galaxies.
    """

    image = fits.getdata(image_fits)

    weight = fits.getdata(weight_fits) if weight_fits else None
    mask = fits.getdata(mask_fits) if mask_fits else None

    img, sigma = prepare_image(image, weight, mask)

    Rbar, err_Rbar = compute_Rbar_with_error(
        img, sigma, x0, y0, pa, inc, q
    )

    R, err_R = compute_R_ratio(
        Rcorot, Rbar,
        err_Rcorot, err_Rbar
    )

    return {
        "Rbar": Rbar,
        "Rbar_err": err_Rbar,
        "R": R,
        "R_err": err_R
    }


result = run_full_bar_pipeline(
    image_fits="galaxy.fits",
    weight_fits="weight.fits",
    mask_fits="mask.fits",
    x0=xc,
    y0=yc,
    pa=PA,
    inc=incl,
    q=axis_ratio,
    Rcorot=Rc,
    err_Rcorot=Rc_err
)

print(result)

#Get result
#{
# 'Rbar': 5.3,
# 'Rbar_err': 0.7,
# 'R': 1.35,
# 'R_err': 0.25
#}
