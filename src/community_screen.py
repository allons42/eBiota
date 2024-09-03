import cobra
import re
import os
import time
import copy
import pandas as pd
from tqdm import tqdm
import pickle
from multiprocessing import Pool
from collections import Counter, defaultdict
import matplotlib.pyplot as plt

from eBiota_utils import get_good_bacteria, config, get_start_and_linker_metabolites

same_gem_memory = {}
mem_path = "tmp/reaction_id_of_every_bacterium/"

def same_gem(gcf1, gcf2):
    global same_gem_memory, mem_path
    if gcf1 > gcf2:
        gcf1, gcf2 = gcf2, gcf1
    if (gcf1, gcf2) not in same_gem_memory:
        with open(mem_path+gcf1+".pkl", "rb") as f:
            rxn1 = pickle.load(f)
        with open(mem_path+gcf2+".pkl", "rb") as f:
            rxn2 = pickle.load(f)

        if len(rxn1&rxn2)/len(rxn1|rxn2) >= 0.95:
            same_gem_memory[(gcf1, gcf2)] = True
        else:
            same_gem_memory[(gcf1, gcf2)] = False

    return same_gem_memory[(gcf1, gcf2)]


def save_rxn_id(f):
    gem = cobra.io.read_sbml_model(os.path.join(config["path_GEM"], f))
    rxn = set([r.id for r in gem.reactions])
    print("saving reaction id at", f"{mem_path}{f.split('.xml')[0]}.pkl")
    with open(f"{mem_path}{f.split('.xml')[0]}.pkl", "wb") as fout:
        pickle.dump(rxn, fout)


def unique():
    gcf_root = config["path_GEM"]
    fs = os.listdir(gcf_root)
    bac_good = get_good_bacteria()
    df = pd.read_csv("stats/carveme_gem_taxonomy.csv")
    df = df[["assembly", "species"]]
    bac2species = {row["assembly"]:row["species"] for idx,row in df.iterrows()}
    bac_good_pruned = {}

    # save reactions ids for each bacterium
    if not os.path.exists(mem_path):
        os.mkdir(mem_path)
    pool = Pool(config["max_proc"])
    for f in fs:
        pool.apply_async(func=save_rxn_id, args=(f,))
    pool.close()
    pool.join() # wait for all child processes to finish

    # remove similar gems in each path
    for k,v in tqdm(bac_good.items(), desc="remove similar gems in each path"):
        new_list = [v[0]]
        bid, newid = 0, 1
        for gcf,growth,prod in v[1:]: # examine each bacteria
            while bid < newid and prod < new_list[bid][2]*0.95: # find bacteria with similar growth*production
                bid += 1
            for i in range(bid, newid):
                gcf_comp = new_list[i][0]
                if bac2species[gcf_comp] == bac2species[gcf]: # find bacteria with same species
                    if same_gem(gcf_comp, gcf): # find bacteria with similar GEM
                        break
            else:
                new_list.append((gcf, growth, prod))
                newid += 1
        bac_good_pruned[k] = new_list
    
    with open("tmp/bac_good_level30_pruned.pkl", "wb") as f:
        pickle.dump(bac_good_pruned, f)


def get_combination():
    if config["prune"]:
        unique()
        bac_good = get_good_bacteria(pkl="tmp/bac_good_level30_pruned.pkl")
    else:
        bac_good = get_good_bacteria()
    
    start_metabolites, linker_metabolites = get_start_and_linker_metabolites()
    target_metabolites = "h2_e" # default target product
    if config["substrate"] != "default":
        start_metabolites = config["substrate"].split(",")
    if config["intermediate"] != "default":
        linker_metabolites = config["intermediate"].split(",")
    if config["product"] != "default":
        target_metabolites = config["product"].split(",")

    path1 = {} # key: (linker, O2, glc), value: substrate
    for k in bac_good:
        if k[0] in start_metabolites and k[1] in linker_metabolites:
            tmp = (k[1], k[2], k[3])
            if tmp not in path1:
                path1[tmp] = []
            path1[tmp].append(k[0])

    path_all = []
    for cond,substrates in path1.items():
        for k in bac_good:
            if (k[0],k[2],k[3]) == cond:
                for x in substrates:
                    # x,k[0],k[1],k[2],k[3] -> substrate,linker,product,O2,glc
                    if x == k[1] or ('glc__D_e' in [x,k[0]] and k[3]=="without_glucose") or (k[1]=='glc__D_e' and k[3]=='with_glucose'):
                        continue
                    if k[1] in target_metabolites:
                        path_all.append(( x,k[0],k[1],k[2],k[3],len(bac_good[(x,cond[0],cond[1],cond[2])]),len(bac_good[k]) )) 

    path_all_merged = defaultdict(list)
    path_cnt = defaultdict(int)
    for p in path_all:
        # p: substrate,linker,product,O2,glc,num1,num2
        path_all_merged[(p[0],p[2],p[3],p[4])].append(p[1])
        path_cnt[(p[0],p[2],p[3],p[4])] += p[5]*p[6]
    path_all_merged_sorted = sorted(path_all_merged.items(), key=lambda x:path_cnt[x[0]])

    # Determine whether to add glucose based on substrate C source
    met_df = pd.read_excel("stats/exchange_metabolites.xlsx", engine="openpyxl")
    met_carbon = [row["id"] for _,row in met_df.iterrows() if row["biolog C source"]=='C']
    # met_type = {row["id"]:row["Compound Type"] for _,row in met_df.iterrows()}
    path_all_merged_sole_carbon = []
    for key,linker in path_all_merged_sorted:
        # key: substrate,product,O2,glc
        # linker: list of linker metabolites
        sub = key[0]
        # pro = key[1]
        glc = (key[3]=="with_glucose")
        if ((sub in met_carbon)^glc) or (sub=="glc__D_e" and glc):
            #if (met_type[pro] not in ["Amino Acids 20", "Inorganic", "Metal Ions"]) or pro=="h2_e":
            path_all_merged_sole_carbon.append((key,linker))

    with open("tmp/all_path_sole_carbon.pkl", "wb") as f:
        pickle.dump(path_all_merged_sole_carbon, f)

    with open("tmp/all_path_sole_carbon.csv", "w") as f:
        f.write("substrate,product,O2,glucose,intermediates,number_of_pairs\n")
        for key,linker in path_all_merged_sole_carbon:
            f.write(f"{key[0]},{key[1]},{key[2]=='with_O2'},{key[3]=='with_glucose'},\"{linker}\",{path_cnt[key]}\n")
