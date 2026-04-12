# PyMorph v3.1.1(Stable release for Python 3)

This repository is a modernized, Python 3 compatible fork of the original **PyMorph v3.1.1** pipeline which was written in Python 2. Since modern packages are not compatible with Python 2, so the pipeline has been updated to use modern packages compatible with **Python 3**.

## Key Changes Made

1 **Indenation errors:** Fixed indentation inconsistencies of Python 3 and 2.  

2 **String/Byte Handling:** Fixed numerous `TypeError: a bytes-like object is required` errors by updating file I/O operations from binary (`ab`, `rb`) to text mode (`a`, `r`) with proper `newline` handling in `pymorph.py` and `writehtmlfunc.py`.  

3 **Integer Division & Slicing:** Resolved `TypeError: slice indices must be integers` by forcing integer conversion in image masking and convolution routines (`pymconvolve.py` and `ellimaskfunc_easy.py`).  

4 **Legacy Imports:** Replaced deprecated `PyFITS` and `PyWCS` with modern `astropy.io.fits` and `astropy.wcs`.


---

## Prerequisites

* **Python 3.10+**
* **Astropy**
* **NumPy & SciPy**
* **GALFIT 3.0.5** (Must be in your system `$PATH`)
* **SExtractor**

---

## How to Use

1. **Prepare Catalog:** Format your input catalog with these headers:
   `gal_id gimg ra1 ra2 ra3 dec1 dec2 dec3 z mzero star`
    The ra and dec must be broken down in hours, mins and secs from degrees
2. **Inputs:** Import your input image, weight image and psf. If you're not giving an input psf, then you can use the @psflist. The pipeline will automatically use the star that is closest to your object as psf   

3. **Configuration:** Update your local paths and settings in `config.py`.

4. **Execution:**
   ```bash
   python pymorph.py
