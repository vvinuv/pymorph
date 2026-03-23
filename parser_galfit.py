import re


class GalfitUtils:

    # =========================================================
    # PARSE GALFIT OUTPUT
    # =========================================================
    def parse_galfit_final(self, filename):

        self.result = {
            "meta": {},
            "target": {},
            "neighbors": {}
        }

        with open(filename, "r") as f:
            lines = f.readlines()

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

                x, y = nums[0], nums[1]

                # errors
                err_line = lines[i + 1]
                err_nums = re.findall(r"[-+]?\d*\.?\d+(?:e[-+]?\d+)?", err_line)
                err_nums = list(map(float, err_nums))

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
                    self.result["neighbors"][comp_id] = comp_dict

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
                        self.result["neighbors"][current_comp_id]["c0"] = c0_val
                        self.result["neighbors"][current_comp_id]["c0_err"] = c0_err

                i += 2
                continue

            # -------------------------------
            # SKY
            # -------------------------------
            elif line.startswith("sky"):
                section = "neighbors"
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
                    self.result["meta"]["chi2"] = float(match.group(1))
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
    def flatten(self, name_map=["bulge", "disk", "bar"]):

        flat = {}
        flat.update(self.result["meta"])

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

                flat[f"{prefix}_{new_k}"] = v

        
        self.result["target"] = flat
        self.result.pop('meta', None)
        self.get_bt()

    def get_bt(self):
        comp = ['bulge', 'disk', 'bar']
        if 'bulge' in comp and 'disk' in comp:
            fb = 10**(-0.4 * self.result["target"]['bulge_mag'])
            fd = 10**(-0.4 * self.result["target"]['disk_mag'])
            ft = fb + fd
        if 'bar' in comp:
            fbar = 10**(-0.4 * self.result["target"]['bar_mag'])
            ft += fbar
            self.result["target"]['BarT'] = round(fbar / ft, 2)
        else:
            self.result["target"]['BarT'] = 9999

        self.result["target"]['BD'] = round(fb / fd, 2)
        self.result["target"]['BT'] = round(fb / ft, 2)

        if 'bulge' in comp and "disk" not in comp:
            self.result["target"]['BT'] = 1.0
            self.result["target"]['BD'] = 0.0
            self.result["target"]['BarT'] = 0.0
        if 'disk' in comp and 'bulge' not in comp:
            self.result["target"]['BD'] = 0.0
            self.result["target"]['BT'] = 0.0
            self.result["target"]['BarT'] = 0.0

if __name__=='__main__':
    g = GalfitUtils()

    result = g.parse_galfit_final("fit.log")
    flat = g.flatten()
