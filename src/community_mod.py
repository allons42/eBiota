import cobra
import os
import pandas as pd
import pickle
from eBiota_utils import get_metabolites_list, config


def get_metabolite_type():
    df = pd.read_excel("stats/exchange_metabolites.xlsx", engine="openpyxl")
    met_type = {}
    for idx,row in df.iterrows():
        met_type[row["id"]] = row["Compound Type"]
    return met_type

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
                rxn.lower_bound = 0 # 禁止吸收培养基以外的物质
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
        gem.reactions.get_by_id(biomass_func[0]).lower_bound = 0.99*sol[biomass_func[0]]
        gem.reactions.get_by_id(biomass_func[1]).lower_bound = 0.99*sol[biomass_func[1]]
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
        res["EX_glc__D_e_1"] = 0
        res["EX_glc__D_e_2"] = 0
        
    # e.g. interests = ["EX_glc__D_e", "EX_3hpp_e"]
    for i in range(len(interests)):
        met = gem.metabolites.get_by_id(interests[i][3:])
        a_produce, b_produce = 0, 0
        for rxn in met.reactions:
            if rxn.id != interests[i]:
                flux = sol[rxn.id] * rxn.get_coefficient(met) # +:producing, -:consuming
                if rxn.id.startswith(names[0]):
                    a_produce += flux
                elif rxn.id.startswith(names[1]):
                    b_produce += flux
        res[interests[i]+"_1"] = a_produce
        res[interests[i]+"_2"] = b_produce
    
    cross12 = {}
    cross21 = {}
    for met in gem.metabolites:
        if met.id.endswith("_e"):
            a_produce, b_produce = 0, 0
            for rxn in met.reactions:
                flux = sol[rxn.id] * rxn.get_coefficient(met) # +:producing, -:consuming
                if rxn.id.startswith(names[0]):
                    a_produce += flux
                elif rxn.id.startswith(names[1]):
                    b_produce += flux
            if a_produce > 1e-8 and b_produce < -1e-8:
                cross12[met.id] = min(a_produce, -b_produce)
            elif b_produce > 1e-8 and a_produce < -1e-8:
                cross21[met.id] = min(-a_produce, b_produce)
    res["cross12"] = cross12
    res["cross21"] = cross21
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


def gene_mod(input_tsv):
    medium = get_metabolites_list()
    with open('ststs/bigg_reaction_info.pkl', "rb") as f:
        rxn2met = pickle.load(f)
    df = pd.read_csv(input_tsv, sep="\t")
    tmp = input_tsv.split("/")[-1].split(".")[0].split("__to__")
    subs = tmp[0]
    prod = tmp[1].split("__with")[0]
    rxn_in = "EX_" + subs
    rxn_out = "EX_" + prod
    O2_bool = True if "with_O2" in tmp[1] else False
    glc_bool = True if "with_glcose" in tmp[1] else False

    fout = open(".".join(input_tsv.split(".")[:-1]) + "_mod.tsv", "w")
    fout.write("\t".join(df.columns) + "\tBest_knockout\tBest_knockin\n")

    for i,row in df.iterrows():
        bac1, bac2 = row["Bac1"], row["Bac2"]
        if bac1[-2] == "_" and bac2[-2] == "_":
            b1 = bac1[:-2] + "." + bac1[-1]
            b2 = bac2[:-2] + "." + bac2[-1]
        GEM1 = read_gem(b1, medium=medium, o2=O2_bool, glc=glc_bool)
        GEM2 = read_gem(b2, medium=medium, o2=O2_bool, glc=glc_bool)
        gem, biomass_func = merge_model([bac1, bac2], [GEM1, GEM2])
        gem.reactions.get_by_id("EX_" + subs).lower_bound = -10
        bestko, bestflux, = "None", row["Total_production"]
        for rxn in gem.reactions:
            with gem:
                if not rxn.id.startswith("EX_"):
                    rxn.knock_out()
                    oneres = FBA([bac1, bac2], gem, biomass_func, [rxn_in, rxn_out])
                    if oneres:
                        tot = oneres[rxn_in+'_1']*oneres['Growth1']+oneres[rxn_in+'_2']*oneres['Growth2']
                        if tot > bestflux:
                            bestko = rxn.id
                            bestflux = tot
        
        bestki, bestflux, = "None", row["Total_production"]
        add_id = lambda x: x if x.endswith("_e") else (bac1 + "_" + x)
        for ins in rxn2met:
            with gem:
                ins_info = rxn2met[ins]
                ins_mets = {add_id(k): v for k,v in ins_info["mets"].items()}
                R = cobra.Reaction(ins)
                R.name = ins
                R.subsystem = ''
                R.bounds = ins_info["bounds"]
                R.add_metabolites({
                    gem.metabolites.get_by_id(k): v for k,v in ins_mets.items()
                })
                gem.add_reactions([R])
                oneres = FBA([bac1, bac2], gem, biomass_func, [rxn_in, rxn_out])
                if oneres:
                    tot = oneres[rxn_in+'_1']*oneres['Growth1']+oneres[rxn_in+'_2']*oneres['Growth2']
                    if tot > bestflux:
                        bestki = ins
                        bestflux = tot

        for k in df.columns:
            fout.write(str(row[k]) + "\t")
            fout.write(f"{bestko}\t{bestki}\n")