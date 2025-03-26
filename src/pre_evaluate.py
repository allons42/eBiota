import cobra
import os
import time
from tqdm import tqdm
from multiprocessing import Pool

from eBiota_utils import get_metabolites_list, config


def test_FBA_in_specific_media(gem, st, ed_list, o2=True, glc__D=True, growth_only=False): 
    if "EX_o2_e" in gem.reactions:
        gem.reactions.get_by_id("EX_o2_e").lower_bound = -10 if o2 else 0
            
    if "EX_glc__D_e" in gem.reactions:
        gem.reactions.get_by_id("EX_glc__D_e").lower_bound = -10 if glc__D else 0
    
    if "EX_"+st in gem.reactions:
        gem.reactions.get_by_id("EX_"+st).lower_bound = -10
    elif not growth_only:
        return {}
    
    bool_dict = {True:"with_", False:"without_"}
    #C_number = gem.metabolites.get_by_id("EX_"+st).elements.get("C",0)
    try:
        sol = cobra.flux_analysis.pfba(gem)
    except:
        return {}
    bm = str(gem.objective.expression).split(' ')[0][4:]
    res = {}
    
    if growth_only:
        return sol[bm]
    
    for ed in ed_list:
        bm_max = sol[bm]
        if bm_max < 1e-8:
            label = 'no_growth'
        elif ("EX_" + ed) not in gem.reactions:
            label = 'no_exchange_reaction'
        else:
            with gem:
                if sol["EX_" + st] > -1e-8:
                    label = 'no_absorption'
                elif sol["EX_" + ed] > 1e-8:
                    label = "good"
                else:
                    label = 'normal'
        res_ex_ed = 0 if (label=='no_exchange_reaction' or label=="no_growth") else sol["EX_" + ed]
        res_ex_st = 0 if (label=='no_exchange_reaction' or label=="no_growth") else sol["EX_" + st]
        res[(st, ed, bool_dict[o2]+"O2", bool_dict[glc__D]+"glucose")] = (label, bm_max, res_ex_st, res_ex_ed)
    return res


def pre_eval(name, possible_path):
    path_GEM = config["path_GEM"]
    suffix = config["suffix"]
    if config["target"] == "production":
        output_dir = "tmp/bacteria_label_production/"
    elif config["target"] == "degradation":
        output_dir = "tmp/bacteria_label_degradation/"
    else:
        output_dir = "tmp/bacteria_label_production/"
        
    try:
        gem = cobra.io.read_sbml_model(os.path.join(path_GEM, name + suffix))
    except:
        print(os.path.join(path_GEM, name + suffix) +': not found')
        return
    
    global basic
    for rxn in gem.exchanges:
        rxt = rxn.reactants[0].id
        if rxt not in basic.keys():
            rxn.lower_bound = 0 # Prohibit absorption of substances other than culture medium
        else:
            rxn.lower_bound = -basic[rxt]
    
    # GEM after preprocess
    # print(gem.medium)
    
    res = {}
    try_dict = {}
    for st,ed in possible_path:
        if st==ed or st in basic.keys(): # Not considering substances in the culture medium as starting materials
            continue
        elif st not in try_dict:
            try_dict[st] = [ed]
        else:
            try_dict[st].append(ed)
    
    for st in try_dict:
        with gem:
            if "EX_"+st not in gem.reactions:
                #print("EX_%s not in %s"%(st,gem.id))
                continue
            
            tmp_res = test_FBA_in_specific_media(gem, st, try_dict[st], o2=True, glc__D=True)
            res.update(tmp_res)
            tmp_res = test_FBA_in_specific_media(gem, st, try_dict[st], o2=False, glc__D=True)
            res.update(tmp_res)
            tmp_res = test_FBA_in_specific_media(gem, st, try_dict[st], o2=True, glc__D=False)
            res.update(tmp_res)
            tmp_res = test_FBA_in_specific_media(gem, st, try_dict[st], o2=False, glc__D=False)
            res.update(tmp_res)

    # If there is an exception in the previous section, there will be no file output
    sorted_res = sorted(res, key = lambda k: (k[0], k[1], k[2],k[3]))
    fout = open(output_dir + name + '.txt', 'w')
    fout.write("substrate,product,O2,glucose,growth,absorption,production,label\n")
    for key in sorted_res:
        fout.write("%s,%s,%s,%s,%.4f,%.4f,%.4f,%s\n"%(key[0],key[1],key[2],key[3],res[key][1],res[key][2],res[key][3],res[key][0]))
    fout.close()
    
    
def run_evaluate():
    input_dir = "tmp/BFS_result/path_from_secretion_to_intake_simplified/"
    if config["target"] == "production":
        output_dir = "tmp/bacteria_label_production/"
    elif config["target"] == "degradation":
        output_dir = "tmp/bacteria_label_degradation/"
    else:
        print(f"Warning: target {config['target']} not recognized. Will continue with production.")
        config["target"] = "production"
        output_dir = "tmp/bacteria_label_production/"
        
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    global basic
    basic = get_metabolites_list()
    all_bac = os.listdir(input_dir)
    
    # multi process
    fs = os.listdir(output_dir)
    print(f"{len(all_bac)} bacteria in total, {len(fs)} have been evaluated before.")
    done = set(fs)
    
    start = time.time()
    n_proc = config["max_proc"]
    pool = Pool(n_proc) # Maximum number of process pools
    pbar = tqdm(total=len(all_bac))
    tasks = []

    print(f"Start to evaluate with {n_proc} processors.")

    for gcf in all_bac:
        if gcf not in done:
            gcf_path = []
            with open(os.path.join(input_dir, gcf)) as fin:
                lines = fin.readlines()
                for line in lines:
                    tmp = line.strip().split("\t")
                    gcf_path.append((tmp[1],tmp[0]))
            
            task = pool.apply_async(func=pre_eval, args=(gcf[:-4], gcf_path), callback=lambda _: pbar.update())
            tasks.append(task)
            done.add(gcf)

    for task in tasks:
        task.wait()
    pool.close()
    pool.join()
    pbar.close()

    end = time.time()
    print("evaluation runtime: %.2f seconds"%(end-start))

if __name__ == '__main__':
    run_evaluate()