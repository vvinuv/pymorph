import os
import re
import numpy as np
import pandas as pd


class GalfitParser:

    def __init__(self, filename=None):
        self.filename = filename
        self.result = {
            "meta": {},
            "target": {},
            "neighbours": {}
        }

    # =========================================================
    # PARSE GALFIT OUTPUT
    # =========================================================
    def parse(self, filename=None):

        if filename:
            self.filename = filename

        if not self.filename:
            raise ValueError("Filename not provided")

        with open(self.filename, "r") as f:
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
                self.result["meta"]["input_image"] = line.split(":", 1)[1].strip()

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

                nums = list(map(float, re.findall(r"[-+]?\d*\.?\d+(?:e[-+]?\d+)?", line)))
                err_nums = list(map(float, re.findall(
                    r"[-+]?\d*\.?\d+(?:e[-+]?\d+)?", lines[i + 1]
                )))

                comp_dict = {
                    "type": comp_type,
                    "x": nums[0],
                    "y": nums[1],
                    "x_err": err_nums[0],
                    "y_err": err_nums[1],
                }

                if comp_type == "sersic":
                    comp_dict.update({
                        "mag": nums[2],
                        "Re": nums[3],
                        "n": nums[4],
                        "ellip": nums[5],
                        "PA": nums[6],
                        "mag_err": err_nums[2],
                        "Re_err": err_nums[3],
                        "n_err": err_nums[4],
                        "ellip_err": err_nums[5],
                        "PA_err": err_nums[6],
                    })

                elif comp_type == "expdisk":
                    comp_dict.update({
                        "mag": nums[2],
                        "Re": nums[3],
                        "ellip": nums[4],
                        "PA": nums[5],
                        "mag_err": err_nums[2],
                        "Re_err": err_nums[3],
                        "ellip_err": err_nums[4],
                        "PA_err": err_nums[5],
                    })

                if section == "target":
                    self.result["target"][comp_id] = comp_dict
                else:
                    self.result["neighbours"][comp_id] = comp_dict

                i += 2
                continue

            # -------------------------------
            # C0 (BOXINESS)
            # -------------------------------
            elif line.strip().startswith("c0"):
                val = re.search(r":\s*([-+]?\d*\.?\d+(?:e[-+]?\d+)?)", line)
                err = re.search(r":\s*([-+]?\d*\.?\d+(?:e[-+]?\d+)?)", lines[i + 1])

                c0_val = float(val.group(1)) if val else None
                c0_err = float(err.group(1)) if err else None

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

                pos = list(map(float, brackets[0].split(",")))
                self.result["meta"]["sky_x"], self.result["meta"]["sky_y"] = pos

                self.result["meta"]["sky_val"] = float(brackets[1])
                self.result["meta"]["sky_err"] = float(brackets[2])

                i += 2
                continue

            # -------------------------------
            # CHI2
            # -------------------------------
            elif "Chi^2 =" in line:
                match = re.search(r"Chi\^2\s*=\s*([\d\.eE+-]+),\s*ndof\s*=\s*(\d+)", line)
                if match:
                    self.result["meta"]["chi2"] = float(match.group(1))
                    self.result["meta"]["ndof"] = int(match.group(2))

            elif "Chi^2/nu" in line:
                match = re.search(r"Chi\^2/nu\s*=\s*([\d\.eE+-]+)", line)
                if match:
                    self.result["meta"]["chi2nu"] = float(match.group(1))

            i += 1
        print(self.result)
        #return self.result

    # =========================================================
    # FLATTEN RESULT
    # =========================================================
    
    def flatten(self):

        flat = {}
        flat.update(self.result["meta"])

        name_map = {1: "bulge", 2: "disk", 3: "bar"}

        param_map = {
            "x": "xctr",
            "y": "yctr",
            "Re": "Re",
            "n": "n",
            "ellip": "ellip",
            "PA": "pa",
            "c0": "boxiness"
        }

        for cid, comp in self.result["target"].items():

            #print(cid, comp)
            prefix = name_map.get(cid, f"comp{cid}")
            #print(prefix)
            if comp.get("type") == "sersic" and cid == 1:
                prefix = "bulge"
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

        print(flat)
        return flat

    # =========================================================
    # WRITE CSV
    # =========================================================
    def to_csv(self, filename="galaxy.csv"):
        flat = self.flatten()
        df = pd.DataFrame([flat])
        df.to_csv(filename, index=False)
        print(f"Saved CSV: {filename}")

if __name__=='__main__':
    parser = GalfitParser('fit.log')
    results = parser.parse()
    #parser.flatten()
    #parser.to_csv('galaxy.csv')
