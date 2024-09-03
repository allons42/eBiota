import cobra
import re
import os
import time
import copy
import pandas as pd
from tqdm import tqdm
import pickle
import numpy as np
from multiprocessing import Pool
import matplotlib.pyplot as plt
from collections import Counter, defaultdict

from eBiota_utils import get_metabolite_name_dict

def get_metabolite_name_dict_direct():
    fnames = open("BFS_result/all_metabolites_names.txt")
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
        metID = linebits["id"][2:-2]
        metabolites_names_dict[name] = metID.strip()
        metabolites_names_dict[name.strip()] = metID.strip()
    return metabolites_names_dict


def trans():
    d = get_metabolite_name_dict()
    root_in = "BFS_result/path_from_secretion_to_intake/"
    root_out = "tmp/"
    fs = os.listdir(root_in)
    print(len(fs), "GEMs are considered")
    for f in fs:
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