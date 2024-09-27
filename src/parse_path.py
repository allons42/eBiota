import re
import os
from tqdm import tqdm

# from eBiota_utils import get_metabolite_name_dict

def get_metabolite_name_dict_direct():
    fnames = open("tmp/BFS_result/all_metabolites_names.txt")
    metabolites_names_dict = {}
    lines = fnames.readlines()
    fnames.close()
    
    reg = re.compile('^name=(?P<name>.*) id=(?P<id>.*)')
    
    for line in lines:
        if len(line) < 2:
            break
        line = line.replace("&gt;",">")
        line = line.replace("&amp;","&")
        line = line.replace("&apos;","'")
        regMatch = reg.match(line)
        linebits = regMatch.groupdict()
        name = linebits["name"]
        metID = linebits["id"][2: -2]
        metabolites_names_dict[name] = metID.strip()
        metabolites_names_dict[name.strip()] = metID.strip()
    return metabolites_names_dict


def translate_path():
    d = get_metabolite_name_dict_direct()
    root_in = "tmp/BFS_result/path_from_secretion_to_intake/"
    root_out = "tmp/BFS_result/path_from_secretion_to_intake_simplified/"
    if not os.path.exists(root_out):
        os.makedirs(root_out)
    fs = os.listdir(root_in)
    print(len(fs), "GEMs are included.")
    for f in tqdm(fs, desc="translating paths"):
        if f.endswith(".txt"):
            with open(os.path.join(root_in, f)) as fin:
                lines = fin.readlines()
            with open(os.path.join(root_out, f), "w") as fout:
                fout.write("secretion\tintake\tlevel\n")
                for line in lines:
                    if len(line) < 5:
                        continue
                    line = line.replace("&gt;",">")
                    line = line.replace("&amp;","&")
                    line = line.replace("&apos;","'")
                    x,y = line.split(" linked to ")
                    y,z = y.split(" with level=")
                    fout.write(f'{d[x]}_e\t{d[y]}_e\t{z}')