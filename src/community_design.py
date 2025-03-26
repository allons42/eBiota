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
from eBiota_utils import get_metabolites_list, get_start_and_linker_metabolites, get_good_bacteria, config

def get_metabolite_type():
    df = pd.read_excel("stats/exchange_metabolites.xlsx", engine="openpyxl")
    met_type = {}
    for idx,row in df.iterrows():
        met_type[row["id"]] = row["Compound Type"]
    return met_type


def get_selected_path():
    # path_all: list of ( ( substrate,product,O2,glc),[linkers] )
    with open("tmp/all_path_sole_carbon.pkl", "rb") as f:
        path_all = pickle.load(f)
    
    #met_type = get_metabolite_type()
    #path_all_screened = [p for p in path_all if (met_type[p[0][1]] not in ["Amino Acids 20", "Inorganic", "Metal Ions"] or p[0][1]=="h2_e")]
    return path_all


# return: gem
def read_gem(name, medium=None, o2=True, glc=True):
    try:
        gem = cobra.io.read_sbml_model(os.path.join(config["path_rewrite"], name + config["suffix"]))
    except:
        print(os.path.join(config["path_rewrite"], name + config["suffix"]) + ": not found")
        return None
    
    if medium:
        for rxn in gem.exchanges:
            rxt = rxn.reactants[0]
            if rxt.id not in medium.keys():
                rxn.lower_bound = 0 # Prohibit absorption of substances other than culture medium
            else:
                rxn.lower_bound = -medium[rxt.id]
    
    if "EX_o2_e" in gem.reactions:
        gem.reactions.get_by_id("EX_o2_e").lower_bound = -10 if o2 else 0
    if "EX_glc__D_e" in gem.reactions:
        gem.reactions.get_by_id("EX_glc__D_e").lower_bound = -10 if glc else 0
    return gem

      
def merge_model(names, GEM):
    # find the biomass function of each GEM
    biomass_func=[]
    for x in GEM:
        exp = str(x.objective.expression).split(' ')[0][4:]
        biomass_func.append(exp)
    
    # merge the GEMs
    new_model = GEM[0]# .copy() # a faster version of deepcopy in cobra 
    for i in range(1,len(GEM)):
        existing_rxn = new_model.reactions.query(lambda rxn: rxn.id in GEM[i].reactions)
        for rxn in existing_rxn:
            the_other =  GEM[i].reactions.get_by_id(rxn.id)
            rxn.lower_bound = min(rxn.lower_bound, the_other.lower_bound)
            rxn.upper_bound = max(rxn.upper_bound, the_other.upper_bound)
           
        needcopyrxn = GEM[i].reactions.query(lambda rxn: rxn.id not in new_model.reactions)
        #copyrxn = copy.deepcopy(needcopyrxn)
        new_model.add_reactions(needcopyrxn)
        #copyrxn = copy.deepcopy(GEM[i].reactions) # copy all the reactions
        #new_model.add_reactions([rxn for rxn in copyrxn if not rxn in new_model.reactions])

        interface = new_model.problem
        new_vars = [
            interface.Variable.clone(v)
            for v in GEM[i].variables
            if v.name not in new_model.variables
        ]
        new_model.add_cons_vars(new_vars)
        new_cons = [
            interface.Constraint.clone(c, model=new_model.solver)
            for c in GEM[i].constraints
            if c.name not in new_model.constraints
        ]
        new_model.add_cons_vars(new_cons, sloppy=True)
    
    obj = {}
    for i in range(len(names)):
        biomass=new_model.reactions.get_by_id(biomass_func[i])
        obj[biomass] = 1.0
    new_model.objective = obj
    
    return new_model, biomass_func


def FBA(names, gem, biomass_func, interests):
    try:
        sol = gem.optimize() # faster than pFBA
        for bm in biomass_func:
            gem.reactions.get_by_id(bm).lower_bound = 0.99*sol[bm]
        gem.objective = {gem.reactions.get_by_id(interests[1]): 1}
        sol = cobra.flux_analysis.pfba(gem)
    except:
        return None
    
    res = {}
    for i in range(len(names)):
        res[f"Bac{i+1}"] = names[i]
    for i in range(len(names)):
        res[f"Growth{i+1}"] = sol[biomass_func[i]]
    if "EX_glc__D_e" in gem.reactions:
        interests.append("EX_glc__D_e")
    else:
        for i in range(len(names)):
            res[f"EX_glc__D_e_{i+1}"] = 0
        
    # e.g. interests = ["EX_glc__D_e", "EX_3hpp_e"]
    for i in range(len(interests)):
        met = gem.metabolites.get_by_id(interests[i][3:]) # remove "EX_"
        x_produce = [0 for j in range(len(names))]
        for rxn in met.reactions:
            if rxn.id != interests[i]:
                flux = sol[rxn.id] * rxn.get_coefficient(met) # +:producing, -:consuming
                for j in range(len(names)):
                    if rxn.id.startswith(names[j]):
                        x_produce[j] += flux
        for j in range(len(names)):
            res[interests[i]+"_"+str(j+1)] = x_produce[j]
    
    cross = defaultdict(dict)
    for met in gem.metabolites:
        if met.id.endswith("_e"):
            x_produce = [0 for j in range(len(names))]
            for rxn in met.reactions:
                flux = sol[rxn.id] * rxn.get_coefficient(met) # +:producing, -:consuming
                for j in range(len(names)):
                    if rxn.id.startswith(names[j]):
                        x_produce[j] += flux
            
            if len(names) == 2:
                if x_produce[0] > 1e-8 and x_produce[1] < -1e-8:
                    cross["12"][met.id] = min(x_produce[0], -x_produce[1])
                elif x_produce[1] > 1e-8 and x_produce[0] < -1e-8:
                    cross["21"][met.id] = min(-x_produce[0], x_produce[1])
            elif len(names) == 3:
                if x_produce[0] > 1e-8 and x_produce[1] < -1e-8:
                    cross["12"][met.id] = min(x_produce[0], -x_produce[1])
                if x_produce[0] > 1e-8 and x_produce[2] < -1e-8:
                    cross["13"][met.id] = min(x_produce[0], -x_produce[2])
                if x_produce[1] > 1e-8 and x_produce[0] < -1e-8:
                    cross["21"][met.id] = min(-x_produce[0], x_produce[1])
                if x_produce[1] > 1e-8 and x_produce[2] < -1e-8:
                    cross["23"][met.id] = min(x_produce[1], -x_produce[2])
                if x_produce[2] > 1e-8 and x_produce[0] < -1e-8:
                    cross["31"][met.id] = min(-x_produce[0], x_produce[2])
                if x_produce[2] > 1e-8 and x_produce[1] < -1e-8:
                    cross["32"][met.id] = min(-x_produce[1], x_produce[2])

    for k in ["12", "13", "21", "23", "31", "32"]:
        res[f"cross{k}"] = cross[k]
    return res


def coculture(pair, medium, O2_bool, glucose_bool, rxn_in, rxn_out, intermediate):
    GEMs = []
    for x in pair:
        gem = read_gem(x, medium=medium, o2=O2_bool, glc=glucose_bool)
        if gem is None:
            return None
        GEMs.append(gem)
    
    single = []
    for gem in GEMs:
        single.append(gem.optimize().objective_value)
    names = [gem.id for gem in GEMs]
    gem, biomass_func = merge_model(names, GEMs)
    gem.reactions.get_by_id(rxn_in).lower_bound = -10
    oneres = FBA(names, gem, biomass_func, [rxn_in, rxn_out])
    if oneres:
        oneres["intermediate"] = intermediate
        for i in range(len(pair)):
            oneres[f"Bac{i+1}_single_growth"] = single[i]
    else:
        with open(config["ERR_LOG"], "a") as err:
            err.write(f"Infeasible: {rxn_in[3:]}, {rxn_out[3:]}, {O2_bool}, {glucose_bool}, {pair}\n")
    return oneres


def run_community_design():
    ERR_LOG = config["ERR_LOG"]
    FBA_LOG = config["LOG"]

    start_metabolites, linker_metabolites = get_start_and_linker_metabolites()
    medium = get_metabolites_list()
    print("substrates:", len(start_metabolites))
    print(start_metabolites)
    print("-"*50)
    print("intermediates:", len(linker_metabolites))
    print(linker_metabolites)
    print("-"*50)
    print("medium:", len(medium))
    print(medium)

    if config["prune"]:
        bac_good = get_good_bacteria(pkl=f"tmp/bac_good_{config['target']}_pruned.pkl")
    else:
        bac_good = get_good_bacteria()
    
    path_all = get_selected_path()
    bool_dict = {True:"with_", False:"without_"}
    root_output = config["path_output"]
    if not os.path.exists(root_output):
        os.mkdir(root_output)

    # use coculture to test：
    # e.g. pair=('GCF_000006685.1', 'GCF_000014585.1'), O2_bool=True, glucose_bool=True, 
    #            rxn_in="glc__D_e", rxn_out="h2_e", intermediate=["co2_e", "acald_e"]

    # path_cnt_h2 = [p for p in path_cnt if p[0][1] == "h2_e"]

    # key: substrate, product, O2, glucose
    # value: list of intermediates
    for key,value in tqdm(path_all):
        start_time = time.time()
        substrate, product, O2, glucose = key[0], key[1], key[2], key[3]
        if config["oxygen"] != "default":
            if config["oxygen"] != (O2 == "with_O2"):
                continue
        if config["glucose"] != "default":
            if config["glucose"] != (glucose == 'with_glucose'):
                continue
        
        if f"{substrate}__to__{product}__{O2}__{glucose}.txt" in os.listdir(root_output):
            continue
        
        pairs = defaultdict(list) # dictionary of { (bac1,bac2): [intermediate] }
        for intermediate in value:
            try:
                p1 = bac_good[substrate, intermediate, O2, glucose][:10]
                p2 = bac_good[intermediate, product, O2, glucose][:10]
            except:
                with open(ERR_LOG, "a") as err:
                    err.write(f"This path has no candidate: {key},{intermediate}")
                continue
                
            if config["designated_bacteria"] != "default":
                for bac in config["designated_bacteria"].split(","):
                    for p in bac_good[substrate, intermediate, O2, glucose]:
                        if p[0] == bac:
                            p1.append(p)
                    for p in bac_good[intermediate, product, O2, glucose]:
                        if p[0] == bac:
                            p2.append(p)

            if config["community_size"] == 2:
                for a,_,_ in p1:
                    for b,_,_ in p2:
                        if a == b or (config["designated_bacteria"] != "default" and a not in config["designated_bacteria"].split(",") and b not in config["designated_bacteria"].split(",")):
                            continue
                        pairs[(a,b)].append(intermediate)
            elif config["community_size"] == 3:
                for a,_,_ in p1:
                    for b,_,_ in p2:
                        if a == b or (config["designated_bacteria"] != "default" and a not in config["designated_bacteria"].split(",") and b not in config["designated_bacteria"].split(",")):
                            continue
                        for c,_,_ in p1:
                            if c == a or c == b or (a,c,b) in pairs:
                                continue
                            pairs[(a,c,b)].append(intermediate)
                        for c,_,_ in p2:
                            if c == a or c == b or (a,b,c) in pairs:
                                continue
                            pairs[(a,b,c)].append(intermediate)
        
        rxn_in = "EX_" + substrate
        rxn_out = "EX_" + product
        O2_bool = (O2 == "with_O2")
        glucose_bool = (glucose == 'with_glucose')
        #print(len(pairs))
        
        pool = Pool(config["max_proc"])
        pool_res = []
        for pair,intermediate in pairs.items():
            if pair[0] != pair[1]: # avoid using the same bacteria
                pool_res.append(pool.apply_async(func=coculture, args=(pair,medium,O2_bool,glucose_bool,rxn_in,rxn_out,intermediate)))
        pool.close()
        pool.join() # wait for all child processes to finish
        
        all_res = []
        for x in pool_res:
            oneres = x.get()
            if oneres != None:
                all_res.append(oneres)

        if config["community_size"] == 2:
            if config["target"] =="production":
                res_sorted = sorted(all_res, key=lambda x: x[rxn_out+"_1"]*x["Growth1"]+x[rxn_out+"_2"]*x["Growth2"], reverse=True)
            else:
                res_sorted = sorted(all_res, key=lambda x: x[rxn_in+"_1"]*x["Growth1"]+x[rxn_in+"_2"]*x["Growth2"], reverse=True)
            fout = open(os.path.join(root_output, f"{substrate}__to__{product}__{O2}__{glucose}.txt"), "w")
            fout.write(f"Bac1\tBac2\tGrowth1\tGrowth2\tIntermediate\t{rxn_in}_1\t{rxn_in}_2\tGlucose_absorption_1\tGlucose_absorption_2\t{rxn_out}_1\t{rxn_out}_2\tCross_feeding_forward\tCross_feeding_reverse\tBac1_mono_growth\tBac2_mono_growth\tTotal_{config['target']}\n")
            for tmp_d in res_sorted:
                fout.write(tmp_d["Bac1"] + "\t")
                fout.write(tmp_d["Bac2"] + "\t")
                fout.write(f"{round(tmp_d['Growth1'],8)+0.0:.8f}\t") # Add zero to turn negative zero into positive zero for nicer display
                fout.write(f"{round(tmp_d['Growth2'],8)+0.0:.8f}\t")
                fout.write(",".join(tmp_d["intermediate"]) + "\t")
                fout.write(f"{round(tmp_d[rxn_in+'_1'],8)+0.0:.8f}\t")
                fout.write(f"{round(tmp_d[rxn_in+'_2'],8)+0.0:.8f}\t")
                fout.write(f"{round(tmp_d['EX_glc__D_e_1'],8)+0.0:.8f}\t")
                fout.write(f"{round(tmp_d['EX_glc__D_e_2'],8)+0.0:.8f}\t")
                fout.write(f"{round(tmp_d[rxn_out+'_1'],8)+0.0:.8f}\t")
                fout.write(f"{round(tmp_d[rxn_out+'_2'],8)+0.0:.8f}\t")
                
                fout.write(f"{len(tmp_d['cross12'])}")
                for k,v in tmp_d["cross12"].items():
                    fout.write(f",{k}:{round(v,8)+0.0:.8f}")
                fout.write(f"\t{len(tmp_d['cross21'])}")
                for k,v in tmp_d["cross21"].items():
                    fout.write(f",{k}:{round(v,8)+0.0:.8f}")
                    
                fout.write(f"\t{round(tmp_d['Bac1_single_growth'],8)+0.0:.8f}\t") 
                fout.write(f"{round(tmp_d['Bac2_single_growth'],8)+0.0:.8f}\t")
                if config["target"] =="production":
                    total = tmp_d[rxn_out + '_1'] * tmp_d['Growth1'] + tmp_d[rxn_out + '_2'] * tmp_d['Growth2']
                    fout.write(f"{round(total,8)+0.0:.8f}\t")
                else:
                    total = tmp_d[rxn_in + '_1'] * tmp_d['Growth1'] + tmp_d[rxn_in + '_2'] * tmp_d['Growth2']
                    fout.write(f"{round(total,8)+0.0:.8f}\t")
                fout.write("\n")
            fout.close()
            with open(FBA_LOG, "a") as log:
                log.write(f"{substrate}\t{product}:\t{len(pairs)} pairs are done in {time.time()-start_time} seconds\n")
        
        elif config["community_size"] == 3:
            if config["target"] =="production":
                res_sorted = sorted(all_res, key=lambda x: x[rxn_out+"_1"]*x["Growth1"]+x[rxn_out+"_2"]*x["Growth2"]+x[rxn_out+"_3"]*x["Growth3"], reverse=True)
            else:
                res_sorted = sorted(all_res, key=lambda x: x[rxn_in+"_1"]*x["Growth1"]+x[rxn_in+"_2"]*x["Growth2"]+x[rxn_in+"_3"]*x["Growth3"], reverse=True)
            fout = open(os.path.join(root_output, f"{substrate}__to__{product}__{O2}__{glucose}.txt"), "w")
            fout.write(f"Bac1\tBac2\tBac3\tGrowth1\tGrowth2\tGrowth3\tIntermediate\t{rxn_in}_1\t{rxn_in}_2\t{rxn_in}_3\t")
            fout.write(f"Glucose_absorption_1\tGlucose_absorption_2\tGlucose_absorption_3\t{rxn_out}_1\t{rxn_out}_2\t{rxn_out}_3\tCross_feeding_12\tCross_feeding_13\tCross_feeding_21\tCross_feeding_23\tCross_feeding_31\tCross_feeding_32\tBac1_mono_growth\tBac2_mono_growth\tBac3_mono_growth\tTotal_{config['target']}\n")
            for tmp_d in res_sorted:
                fout.write(tmp_d["Bac1"] + "\t")
                fout.write(tmp_d["Bac2"] + "\t")
                fout.write(tmp_d["Bac3"] + "\t")
                fout.write(f"{round(tmp_d['Growth1'],8)+0.0:.8f}\t")
                fout.write(f"{round(tmp_d['Growth2'],8)+0.0:.8f}\t")
                fout.write(f"{round(tmp_d['Growth3'],8)+0.0:.8f}\t")
                fout.write(",".join(tmp_d["intermediate"]) + "\t")
                fout.write(f"{round(tmp_d[rxn_in+'_1'],8)+0.0:.8f}\t")
                fout.write(f"{round(tmp_d[rxn_in+'_2'],8)+0.0:.8f}\t")
                fout.write(f"{round(tmp_d[rxn_in+'_3'],8)+0.0:.8f}\t")
                fout.write(f"{round(tmp_d['EX_glc__D_e_1'],8)+0.0:.8f}\t")
                fout.write(f"{round(tmp_d['EX_glc__D_e_2'],8)+0.0:.8f}\t")
                fout.write(f"{round(tmp_d['EX_glc__D_e_3'],8)+0.0:.8f}\t")
                fout.write(f"{round(tmp_d[rxn_out+'_1'],8)+0.0:.8f}\t")
                fout.write(f"{round(tmp_d[rxn_out+'_2'],8)+0.0:.8f}\t")
                fout.write(f"{round(tmp_d[rxn_out+'_3'],8)+0.0:.8f}\t")

                for k in ["12", "13", "21", "23", "31", "32"]:
                    fout.write(f"{len(tmp_d['cross'+k])}")
                    for kk,v in tmp_d["cross"+k].items():
                        fout.write(f",{kk}:{round(v,8)+0.0:.8f}")
                    fout.write("\t")

                fout.write(f"{round(tmp_d['Bac1_single_growth'],8)+0.0:.8f}\t") 
                fout.write(f"{round(tmp_d['Bac2_single_growth'],8)+0.0:.8f}\t")
                fout.write(f"{round(tmp_d['Bac3_single_growth'],8)+0.0:.8f}\t")

                if config["target"] =="production":
                    total = tmp_d[rxn_out + '_1'] * tmp_d['Growth1'] + tmp_d[rxn_out + '_2'] * tmp_d['Growth2'] + tmp_d[rxn_out + '_3'] * tmp_d['Growth3']
                    fout.write(f"{round(total,8)+0.0:.8f}\t")
                else:
                    total = tmp_d[rxn_in + '_1'] * tmp_d['Growth1'] + tmp_d[rxn_in + '_2'] * tmp_d['Growth2'] + tmp_d[rxn_in + '_3'] * tmp_d['Growth3']
                    fout.write(f"{round(total,8)+0.0:.8f}\t")
                fout.write("\n")
            fout.close()
            with open(FBA_LOG, "a") as log:
                log.write(f"{substrate}\t{product}:\t{len(pairs)} pairs are done in {time.time()-start_time} seconds\n")
