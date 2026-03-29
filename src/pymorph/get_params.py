import os
import re
import configparser
from astropy.cosmology import Planck18 as cosmo
import numpy as np
from .errors_warnings import catch_pipeline_issues, PipelineCriticalError

class GetOutputParams:

    def __init__(self, config_file, pipeline):

        self.pipeline = pipeline

        config = configparser.ConfigParser()
        config.read(config_file)


        self.pixelscale = config.getfloat('cosmology', 'pixelscale')
        #self.H0 = config.getfloat('cosmology', 'H0')
        #self.WM = config.getfloat('cosmology', 'WM')
        #self.WV = config.getfloat('cosmology', 'WV')
        #self.redshift = config.getfloat('cosmology', 'redshift')

    # ====================
    # SAFELY RETURN VALUES
    # ====================
    def get_safe(self, arr):
        arr = np.asarray(arr)

        result = []
        for val in arr:
            try:
                if val is None or (isinstance(val, float) and np.isnan(val)):
                    result.append(9999)
                else:
                    result.append(val)
            except Exception:
                result.append(9999)

        return np.array(result)


    def append_and_remove(self, source_file, target_file):
        if os.path.exists(source_file):
            with open(source_file, 'r') as src, open(target_file, 'a') as tgt:
                tgt.write(src.read())

            os.remove(source_file)
            #print(f"{source_file} appended and deleted.")
        else:
            print(f"{source_file} not found.")


    @catch_pipeline_issues(expected_len=7)
    def check_param_length_ser(self, lst):
        return lst

    @catch_pipeline_issues(expected_len=6)
    def check_param_length_exp(self, lst):
        return lst

    @catch_pipeline_issues()
    def check_params_and_errors(self, params, errors):
        return params, errors

    @catch_pipeline_issues()
    def check_file(self, filename):
        return filename

    # =========================================================
    # PARSE GALFIT OUTPUT
    # =========================================================
    def parse_galfit(self, filename, components, z):

        self.components = components
        self.z = z

        self.result = {
            "meta": {},
            "target": {},
            "neighbours": {}
        }

        with open(filename, "r") as f:
            lines = f.readlines()

        #self.append_and_remove("fit.log", "fit2.log")

        i = 0
        comp_id = 0
        current_comp_id = None
        section = "target"

            
        while i < len(lines):
            line = lines[i].strip()

            # -------------------------------
            # HEADER
            # -------------------------------
            if "Input image" in line:
                self.result["meta"]["input_image"] = line.split(":", 1)[1].strip().split("[")[0]

            elif "Init. par. file" in line:
                self.result["meta"]["init_file"] = line.split(":", 1)[1].strip()

            elif "Restart file" in line:
                self.result["meta"]["restart_file"] = line.split(":", 1)[1].strip()

            elif "Output image" in line:
                self.result["meta"]["output_image"] = line.split(":", 1)[1].strip()

            # -------------------------------
            # COMPONENTS
            # -------------------------------
            elif line.startswith(("sersic", "expdisk")):
                comp_id += 1
                current_comp_id = comp_id

                comp_type = line.split()[0]

                nums = re.findall(r"[-+]?\d*\.?\d+(?:e[-+]?\d+)?", line)
                nums = list(map(float, nums))
                 
                # errors
                err_line = lines[i + 1]
                err_nums = re.findall(r"[-+]?\d*\.?\d+(?:e[-+]?\d+)?", err_line)
                err_nums = list(map(float, err_nums))
                

                #print('nums', nums)
                #print('err_nums', err_nums)
                if line.startswith("sersic"):
                    nums = self.check_param_length_ser(nums)
                    #print('sersic1', nums)
                    nums_err_nums = self.check_params_and_errors(nums, err_nums)
                    #print('sersic2', nums_err_nums[0])
                    #print('sersic err', nums_err_nums[1])
                    nums = nums_err_nums[0]
                    err_nums = nums_err_nums[1]

                elif line.startswith("expdisk"):
                    nums = self.check_param_length_exp(nums)
                    nums_err_nums = self.check_params_and_errors(nums, err_nums)
                    nums = nums_err_nums[0]
                    err_nums = nums_err_nums[1]
                    #print('exp2', nums_err_nums[0])
                    #print('exp err', nums_err_nums[1])
                #    err_nums  = self.check_elements_exp(err_nums)
                #    print('exp', nums)
                #    print('exp err', err_nums)



                                #result = self.get_safe(err_nums)
                #print("result", result)

                x, y = nums[0], nums[1]

                comp_dict = {
                    "type": comp_type,
                    "x": x,
                    "y": y,
                    "x_err": err_nums[0],
                    "y_err": err_nums[1],
                }

                if comp_type == "sersic":
                    comp_dict.update({
                        "mag": nums[2],
                        "Re": nums[3],
                        "n": nums[4],
                        "ell": nums[5],
                        "PA": nums[6],

                        "mag_err": err_nums[2],
                        "Re_err": err_nums[3],
                        "n_err": err_nums[4],
                        "ell_err": err_nums[5],
                        "PA_err": err_nums[6],
                    })

                elif comp_type == "expdisk":
                    comp_dict.update({
                        "mag": nums[2],
                        "Re": nums[3],
                        "ell": nums[4],
                        "PA": nums[5],

                        "mag_err": err_nums[2],
                        "Re_err": err_nums[3],
                        "ell_err": err_nums[4],
                        "PA_err": err_nums[5],
                    })

                if section == "target":
                    self.result["target"][comp_id] = comp_dict
                else:
                    self.result["neighbours"][comp_id] = comp_dict

                i += 2
                continue

            # -------------------------------
            # c0
            # -------------------------------
            elif line.startswith("c0"):
                match_val = re.search(r":\s*([-+]?\d*\.?\d+(?:e[-+]?\d+)?)", line)
                c0_val = float(match_val.group(1)) if match_val else None

                err_line = lines[i + 1]
                match_err = re.search(r":\s*([-+]?\d*\.?\d+(?:e[-+]?\d+)?)", err_line)
                c0_err = float(match_err.group(1)) if match_err else None

                if current_comp_id is not None:
                    if section == "target":
                        self.result["target"][current_comp_id]["c0"] = c0_val
                        self.result["target"][current_comp_id]["c0_err"] = c0_err
                    else:
                        self.result["neighbours"][current_comp_id]["c0"] = c0_val
                        self.result["neighbours"][current_comp_id]["c0_err"] = c0_err

                i += 2
                continue

            # -------------------------------
            # SKY
            # -------------------------------
            elif line.startswith("sky"):
                section = "neighbours"
                current_comp_id = None

                brackets = re.findall(r"\[([^\]]+)\]", line)

                pos_vals = list(map(float, brackets[0].split(",")))
                self.result["meta"]["sky_x"], self.result["meta"]["sky_y"] = pos_vals

                self.result["meta"]["sky_err"] = float(brackets[2])

                if len(brackets) > 2:
                    self.result["meta"]["sky_val"] = float(brackets[1])

                i += 2
                continue

            # -------------------------------
            # CHI-SQUARE
            # -------------------------------
            elif "Chi^2 =" in line:
                match = re.search(
                    r"Chi\^2\s*=\s*([\d\.eE+-]+),\s*ndof\s*=\s*(\d+)", line
                )
                if match:
                    self.result["meta"]["chi2"] = int(float(match.group(1))) 
                    self.result["meta"]["ndof"] = int(match.group(2))

            elif "Chi^2/nu" in line:
                match = re.search(r"Chi\^2/nu\s*=\s*([\d\.eE+-]+)", line)
                if match:
                    self.result["meta"]["chi2nu"] = float(match.group(1))

            i += 1

        return self.result


    # =========================================================
    # FLATTEN RESULT
    # =========================================================
    def flatten(self, components):

        self.flat = {}
        self.flat.update(self.result["meta"])

        param_map = {
            "x": "xctr",
            "y": "yctr",
            "Re": "Re",
            "n": "n",
            "ell": "ell",
            "PA": "pa",
            "c0": "boxiness"
        }

        for cid, comp in self.result["target"].items():

            if comp.get("type") == "sersic" and cid == 1:
                prefix = "bulge"
            elif comp.get("type") == "sersic" and cid > 1:
                prefix = "bar"
            elif comp.get("type") == "expdisk":
                prefix = "disk"

            for k, v in comp.items():

                if k == "type":
                    continue

                if k.endswith("_err"):
                    base = k.replace("_err", "")
                    new_k = param_map.get(base, base) + "_err"
                else:
                    new_k = param_map.get(k, k)

                self.flat[f"{prefix}_{new_k}"] = v

        
        self.result["target"] = self.flat
        self.result.pop('meta', None)
        self.get_bt(components)


        self.neighbours = self.result["neighbours"] 
        
        

        for comp in components:
            out_key = f"{comp}_Re_kpc"
            out_key_e = f"{comp}_Re_kpc_err"
            rad = self.result["target"][f"{comp}_Re"]
            rad_e = self.result["target"][f"{comp}_Re_err"]

            if self.z <= 0:
                self.result["target"][out_key] = -9999
                self.result["target"][out_key_e]= -9999
            else:
                self.result["target"][out_key] = self.radius_pix_to_kpc(rad, 
                                                              self.pixelscale,
                                                              self.z)

                self.result["target"][out_key_e] = self.radius_pix_to_kpc(rad_e,
                                                             self.pixelscale,
                                                             self.z)


        self.result = self.result["target"]


    def get_bt(self, components):
        if 'bulge' in components and 'disk' in components:
            fb = 10**(-0.4 * self.result["target"]['bulge_mag'])
            fd = 10**(-0.4 * self.result["target"]['disk_mag'])
            ft = fb + fd
        if 'bar' in components:
            fbar = 10**(-0.4 * self.result["target"]['bar_mag'])
            ft += fbar
            self.result["target"]['BarT'] = round(fbar / ft, 2)
        else:
            self.result["target"]['BarT'] = 9999

        self.result["target"]['BD'] = round(fb / fd, 2)
        self.result["target"]['BT'] = round(fb / ft, 2)

        if 'bulge' in components and "disk" not in components:
            self.result["target"]['BT'] = 1.0
            self.result["target"]['BD'] = 0.0
            self.result["target"]['BarT'] = 0.0
        if 'disk' in components and 'bulge' not in components:
            self.result["target"]['BD'] = 0.0
            self.result["target"]['BT'] = 0.0
            self.result["target"]['BarT'] = 0.0



    def radius_pix_to_kpc(self, r_pix, pixel_scale, z):

        kpc_per_arcsec = cosmo.kpc_proper_per_arcmin(z).value / 60.0
       
        r_kpc = r_pix * pixel_scale * kpc_per_arcsec 
        
        return round(r_kpc, 2)



if __name__=='__main__':
    g = GetOutputParams()

    result = g.parse_galfit("fit.log")
    flat = g.flatten()
