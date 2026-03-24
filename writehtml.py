import os
import csv
import datetime
import numpy as np
import fileinput
from cosmocal import CosmoCal 
import traceback


def generate_galaxy_report(d, output_file="report.html", image_path="image.png"):

    target = d["target"]
    neighbours = d.get("neighbours", {})

    name = target.get("NAME", "Galaxy")

    html = f"""
    <html>
    <head>
    <title>Analysis report of {name}</title>
    <style>
        body {{ font-family: Arial; margin: 20px; }}
        table {{ border-collapse: collapse; width: 100%; margin-bottom: 20px; }}
        th, td {{ border: 1px solid black; padding: 6px; text-align: center; }}
        th {{ background-color: #f2f2f2; }}
        .flex {{ display: flex; }}
        .left {{ width: 50%; }}
        .right {{ width: 50%; }}
    </style>
    </head>

    <body>

    <h2>Analysis report of {name}</h2>

    <!-- ================= TOP TABLE ================= -->
    <table>
    <tr>
        <th>RA</th><th>DEC</th>
        <th>Image</th><th>GALFIT Input</th><th>Restart File</th>
    </tr>
    <tr>
        <td>{target.get('RA_HMS')}</td>
        <td>{target.get('DEC_DMS')}</td>
        <td>{target.get('input_image')}</td>
        <td>{target.get('init_file')}</td>
        <td>{target.get('restart_file')}</td>
    </tr>

    <tr>
        <th>Output Image</th><th>Number of DoF</th><th>Chi2nu</th>
        <th>B/T</th><th>Bar/T</th>
    </tr>
    <tr>
        <td>{target.get('output_image')}</td>
        <td>{target.get('ndof')}</td>
        <td>{target.get('chi2nu')}</td>
        <td>{target.get('BT')}</td>
        <td>{target.get('BarT')}</td>
    </tr>

    <tr>
        <th colspan="1">Filter</th>
    </tr>
    <tr>
        <td colspan="1">{target.get('FILTER')}</td>
    </tr>
    </table>

    <!-- ================= IMAGE + SEXTRACTOR ================= -->
    <div class="flex">

        <div class="left">
            <img src="{image_path}" width="90%">
        </div>

        <div class="right">
        <table>
            <tr><th colspan="2">SExtractor Parameters</th></tr>
            <tr><td>X_IMAGE</td><td>{target.get('X_IMAGE')}</td></tr>
            <tr><td>Y_IMAGE</td><td>{target.get('Y_IMAGE')}</td></tr>
            <tr><td>MAG_AUTO</td><td>{target.get('MAG_AUTO')}</td></tr>
            <tr><td>FLUX_RADIUS</td><td>{target.get('FLUX_RADIUS')}</td></tr>
            <tr><td>ELLONGATION</td><td>{target.get('ELONGATION')}</td></tr>
            <tr><td>THETA_IMAGE</td><td>{target.get('THETA_IMAGE')}</td></tr>
            <tr><td>ISO0</td><td>{target.get('ISO0')}</td></tr>
            <tr><td>A_IMAGE</td><td>{target.get('A_IMAGE')}</td></tr>
            <tr><td>BACKGROUND</td><td>{target.get('BACKGROUND')}</td></tr>
        </table>
        </div>

    </div>

    <!-- ================= TARGET COMPONENTS ================= -->
    <table>
    <tr>
        <th>ID</th><th>TYPE</th>
        <th>MAG</th><th>MAG<sub>e</sub></th>
        <th>RE</th><th>RE<sub>e</sub></th>
        <th>N</th><th>N<sub>e</sub></th>
        <th>ε</th><th>ε<sub>e</sub></th>
        <th>PA</th><th>PA<sub>e</sub></th>
        <th>KPC</th><th>Kpc<sub>e</sub></th>
        <th>BOXINESS</th><th>BOXINESS<sub>e</sub></th>
    </tr>
    """

    # ---------- TARGET COMPONENTS ----------
    components = ["bulge", "disk", "bar"]

    for i, comp in enumerate(components, start=1):

        if f"{comp}_mag" not in target:
            continue

        html += f"""
        <tr>
        <td>{i}</td>
        <td>{comp.upper()}</td>
        <td>{target.get(f"{comp}_mag")}</td>
        <td>{target.get(f"{comp}_mag_err")}</td>
        <td>{target.get(f"{comp}_Re")}</td>
        <td>{target.get(f"{comp}_Re_err")}</td>
        <td>{target.get(f"{comp}_n","")}</td>
        <td>{target.get(f"{comp}_n_err","")}</td>
        <td>{target.get(f"{comp}_ell")}</td>
        <td>{target.get(f"{comp}_ell_err")}</td>
        <td>{target.get(f"{comp}_pa")}</td>
        <td>{target.get(f"{comp}_pa_err")}</td>
        <td>{target.get(f"{comp}_kpc","")}</td>
        <td>{target.get(f"{comp}_kpc_err","")}</td>
        <td>{target.get(f"{comp}_boxiness","")}</td>
        <td>{target.get(f"{comp}_boxiness_err","")}</td>
        </tr>
        """

    html += "</table>"

    # ================= SKY =================
    html += f"""
    <table>
    <tr><th colspan="3">Fitted Sky</th></tr>
    <tr><td>X</td><td>Y</td><td>Value</td></tr>
    <tr>
        <td>{target.get('sky_x')}</td>
        <td>{target.get('sky_y')}</td>
        <td>{target.get('sky_val')}</td>
    </tr>
    </table>
    """

    # ================= CASGM =================
    html += f"""
    <table>
    <tr>
        <th>C</th><th>C_ERR</th>
        <th>A</th><th>A_ERR</th>
        <th>S</th><th>S_ERR</th>
        <th>GINI</th><th>GINI_ERR</th>
        <th>M20</th><th>M20_ERR</th>
    </tr>
    <tr>
        <td>{target.get('C')}</td>
        <td>{target.get('C_err')}</td>
        <td>{str(target.get('A'))[:4]}</td>
        <td>{target.get('A_err')}</td>
        <td>{str(target.get('S'))[:4]}</td>
        <td>{target.get('S_err')}</td>
        <td>{target.get('gini')}</td>
        <td>{target.get('gini_err')}</td>
        <td>{target.get('m20')}</td>
        <td>{target.get('m20_err')}</td>
    </tr>
    </table>
    """

    # ================= neighbours =================
    html += """
    <table>
    <tr>
        <th>ID</th><th>TYPE</th>
        <th>X</th><th>Y</th>
        <th>MAG</th><th>RE</th><th>N</th>
        <th>ELL</th><th>PA</th>
    </tr>
    """

    for nid, comp in neighbours.items():
        html += f"""
        <tr>
        <td>{nid}</td>
        <td>{comp.get('type').upper()}</td>
        <td>{comp.get('x')}</td>
        <td>{comp.get('y')}</td>
        <td>{comp.get('mag')}</td>
        <td>{comp.get('Re')}</td>
        <td>{comp.get('n',"")}</td>
        <td>{comp.get('ell')}</td>
        <td>{comp.get('PA')}</td>
        </tr>
        """

    html += "</table>"

    html += "</body></html>"

    # ================= WRITE FILE =================
    with open(output_file, "w") as f:
        f.write(html)

    print(f"Saved HTML report: {output_file}")


if __name__=='__main__':
    d = {'target': {'NAME': 'cl1358_9.0', 'GAL_ID': 9.0, 'RA': 244.9011765, 'DEC': 24.7287144, 'MAG': 20.0, 'Z': 0.1, 'GALFIT_ANGLE': np.float64(-117.26), 'POSITION': np.True_, 'ROOTNAME': 'cl1358', 'MAG_ZERO': 25.256, 'FILTER': 'SDSS_r', 'DATE': '2026-03-24 15:19:55', 'X_IMAGE': np.float64(93.0), 'Y_IMAGE': np.float64(93.0), 'FLUX_RADIUS': np.float64(15.54), 'THETA_IMAGE': np.float64(-27.26), 'ELONGATION': np.float64(1.081), 'A_IMAGE': np.float64(12.332), 'MAG_AUTO': np.float64(18.315), 'ISO0': np.float64(3034.0), 'BACKGROUND': np.float64(0.003969985), 'RA_HMS': '16:43:43.4', 'DEC_DMS': '+24:43:43.4', 'IMAGE_SIZE': 186, 'PSF': '/Users/vinu/github/pymorph_chatgpt/data_large/psf_1447219+0828392.fits', 'distance_psf_arcsec': 98440, 'YES_BAR': True, 'R20': np.float64(5.39), 'R50': np.float64(5.39), 'R80': np.float64(5.39), 'R90': np.float64(5.39), 'C': np.float64(1.0), 'A': np.float32(0.35535), 'S': np.float32(0.17568), 'C_err': np.float64(0.0306), 'A_err': np.float64(0.0055), 'S_err': np.float64(0.0014), 'gini': np.float64(0.4244), 'gini_err': np.float64(0.0019), 'm20': np.float64(-0.8172), 'm20_err': np.float64(0.0092), 'input_image': 'I_cl1358_10.fits', 'init_file': 'G_cl1358_10.in', 'restart_file': 'galfit.24', 'output_image': 'O_cl1358_10.fits', 'sky_x': 35.0, 'sky_y': 35.0, 'sky_err': 0.0, 'sky_val': 0.0101, 'chi2': 5758, 'ndof': 4670, 'chi2nu': 1.233, 'bulge_xctr': 33.75, 'bulge_yctr': 33.6, 'bulge_xctr_err': 0.65, 'bulge_yctr_err': 0.94, 'bulge_mag': 22.42, 'bulge_Re': 5.43, 'bulge_n': 1.0, 'bulge_ell': 0.9, 'bulge_pa': -70.48, 'bulge_mag_err': 1.02, 'bulge_Re_err': 1.34, 'bulge_n_err': 0.25, 'bulge_ell_err': 0.1, 'bulge_pa_err': 47.96, 'bulge_boxiness': -0.54, 'bulge_boxiness_err': 0.49, 'disk_xctr': 35.17, 'disk_yctr': 36.27, 'disk_xctr_err': 0.43, 'disk_yctr_err': 0.39, 'disk_mag': 21.87, 'disk_Re': 4.37, 'disk_ell': 0.9, 'disk_pa': -89.98, 'disk_mag_err': 0.24, 'disk_Re_err': 0.41, 'disk_ell_err': 0.04, 'disk_pa_err': 25.77, 'bar_xctr': 33.82, 'bar_yctr': 34.49, 'bar_xctr_err': 0.13, 'bar_yctr_err': 0.19, 'bar_mag': 22.73, 'bar_Re': 4.55, 'bar_n': 1.55, 'bar_ell': 0.9, 'bar_pa': -65.77, 'bar_mag_err': 1.78, 'bar_Re_err': 4.48, 'bar_n_err': 0.59, 'bar_ell_err': 0.23, 'bar_pa_err': 10.73, 'bar_boxiness': -0.54, 'bar_boxiness_err': 0.49, 'BarT': 0.22, 'BD': 0.6, 'BT': 0.29, 'bulge_Re_kpc': np.float64(0.47), 'bulge_Re_kpc_err': np.float64(0.11), 'disk_Re_kpc': np.float64(0.37), 'disk_Re_kpc_err': np.float64(0.04), 'bar_Re_kpc': np.float64(0.39), 'bar_Re_kpc_err': np.float64(0.38)}, 'neighbours': {4: {'type': 'sersic', 'x': 41.94, 'y': 9.42, 'x_err': 0.07, 'y_err': 0.05, 'mag': 22.99, 'Re': 1.06, 'n': 1.91, 'ell': 0.7, 'PA': -77.41, 'mag_err': 0.94, 'Re_err': 0.93, 'n_err': 2.08, 'ell_err': 0.03, 'PA_err': 5.68, 'c0': -0.54, 'c0_err': 0.49}, 5: {'type': 'sersic', 'x': 13.51, 'y': 10.21, 'x_err': 0.04, 'y_err': 0.06, 'mag': 21.74, 'Re': 1.4, 'n': 3.06, 'ell': 0.59, 'PA': -148.0, 'mag_err': 0.55, 'Re_err': 0.9, 'n_err': 1.96, 'ell_err': 0.02, 'PA_err': 2.36}}}
    generate_galaxy_report(d, "report.html", image_path="galaxy.png")
