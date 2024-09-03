import cobra
import re
import os
import time
import copy
from multiprocessing import Pool
import pandas as pd

from eBiota_utils import get_metabolites_list


def test_fva_in_specific_media(gem, st, ed_list, o2=True, glc__D=True, growth_only=False): 
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
            sol2 = sol
            if sol2["EX_" + st] > -1e-8:
                label = 'no_absorption'
            elif sol2["EX_" + ed] > 1e-8:
                label = "good"
            else:
                label = 'normal'
        res_ex_ed = 0 if (label=='no_exchange_reaction' or label=="no_growth") else sol2["EX_" + ed]
        res_ex_st = 0 if (label=='no_exchange_reaction' or label=="no_growth") else sol2["EX_" + st]
        res[(st, ed, bool_dict[o2]+"O2", bool_dict[glc__D]+"glucose")] = (label, bm_max, res_ex_st, res_ex_ed)
    return res


def pre_eval(name, possible_path):
    try:
        gem = cobra.io.read_sbml_model("/data3/jhhou/ebiota/GEM_rewrite/bacteria/Carve_RefSeq/" + name + '.xml.gz')
    except:
        print("/data3/jhhou/ebiota/GEM_rewrite/bacteria/Carve_RefSeq/" + name + '.xml.gz: not found')
        return
    
    global basic, output_dir
    for rxn in gem.exchanges:
        rxt = rxn.reactants[0].id
        if rxt not in basic.keys():
            rxn.lower_bound = 0 # 禁止吸收培养基以外的物质
        else:
            rxn.lower_bound = -basic[rxt]
    
    # 此处的gem为预处理完成的gem
    # print(gem.medium)
    
    res = {}
    try_dict = {}
    for st,ed in possible_path:
        if st==ed or st in basic.keys(): # 不考虑培养基中的物质作为起始物质
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
            
            tmp_res = test_fva_in_specific_media(gem, st, try_dict[st], o2=True, glc__D=True)
            res.update(tmp_res)
            tmp_res = test_fva_in_specific_media(gem, st, try_dict[st], o2=False, glc__D=True)
            res.update(tmp_res)
            tmp_res = test_fva_in_specific_media(gem, st, try_dict[st], o2=True, glc__D=False)
            res.update(tmp_res)
            tmp_res = test_fva_in_specific_media(gem, st, try_dict[st], o2=False, glc__D=False)
            res.update(tmp_res)

    # 若前面部分出现异常，则不会有文件输出
    sorted_res = sorted(res, key = lambda k: (k[0], k[1], k[2],k[3]))
    fout = open(output_dir + name + '.txt','w')
    fout.write("substrate,product,O2,glucose,growth,absorption,production,label\n")
    for key in sorted_res:
        fout.write("%s,%s,%s,%s,%.4f,%.4f,%.4f,%s\n"%(key[0],key[1],key[2],key[3],res[key][1],res[key][2],res[key][3],res[key][0]))
    fout.close()
    
    
def run_evaluate():
    input_dir = "BFS_result/path_from_secretion_to_intake_simplified/"
    output_dir = "tmp/bacteria_label_level30/"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    global basic
    basic = get_metabolites_list()
    all_bac = os.listdir(input_dir)
    
    # multi process
    fs = os.listdir(output_dir)
    print(len(fs), "GEMs already done")
    done = set(fs)
    
    start = time.time()
    n_proc = 10

    while len(done) < len(all_bac):

        cnt = len(done)
        limit = cnt + 1000
        pool = Pool(n_proc) # 进程池数量上限

        for gcf in all_bac:
            if len(done) >= limit:
                print(len(done))
                break
            if gcf not in done:
                gcf_path = []
                with open(os.path.join(input_dir, gcf)) as fin:
                    lines = fin.readlines()
                    for line in lines:
                        tmp = line.strip().split("\t")
                        gcf_path.append((tmp[1],tmp[0]))
                pool.apply_async(func=pre_eval, args=(gcf[:-4], gcf_path))
                done.add(gcf)

        pool.close()
        pool.join()

    end = time.time()
    print("runtime:%.2f seconds"%(end-start))

if __name__ == '__main__':
    run_evaluate()