import os
import json
import cobra
import pandas as pd
from tqdm import tqdm
import pickle

def read_config():
    # read config file
    with open('config.json', 'r') as file:
        config = json.load(file)
    return config

config = read_config()

def check_config():
    # check the format of config
    if config["suffix"] != ".xml":
        print("Error: suffix must be '.xml'. For '.xml.gz' files, please run 'gzip -d *.gz' to decompress.")
        exit()

    if not os.path.exists(config["path_GEM"]):
        print(f"Error: path_GEM {config['path_GEM']} not found.")
        exit()
    config["path_rewrite"] = "tmp/rewrite_GEM"

    if os.path.exists(config["path_output"]):
        print(f"Warning: Output directory {config['path_output']} already exists. Please make sure it is correct.")
        # exit()

    if not os.path.exists(config["medium"]):
        print(f"Error: medium {config['medium']} not found.")
        exit()

    if config["target"] not in ["production", "degradation"]:
        print(f"Error: target should be 'production' or 'degradation'.")
        exit()

    if type(config["max_proc"]) != int or config["max_proc"] < 1:
        print("Error: max_proc should be a positive integer.")
        exit()

    if type(config["prune"]) != bool:
        print("Error: prune should be a boolean, true or false.")
        exit()

    if config["community_size"] not in {2, 3}:
        print("Error: currently only support community of size 2 or 3.")
        exit()

    if not os.path.exists("stats/model_weights_2bac.pth"):
        print("Error: Required DeepCooc model files are missing. Please download the necessary files by following the instructions in the README or the eBiota documentation: https://ebiota.readthedocs.io/en/latest/install.html")
        exit()

def update_config(arg):
    if arg.outdir:
        config["path_output"] = arg.outdir

    with open('config.json', 'w') as f:
        json.dump(config, f, indent=4)


# access the medium, return a dict:
#     key: id of exchange metabolites, ends with "_e"
#     value: max uptake rate
def get_metabolites_list(medium="Basic"):
    medium_file = {
        "Basic": "BasicLB_medium",
        "DM38": "DM38_medium",
        "ASF": "ASF_medium"}
    
    df = pd.read_csv(f"stats/{medium_file[medium]}.csv", sep=",")
    basic = {key:values for key, values in zip(df['id'], df['maxFlux'])}
    return basic

# read start and linker metabolites
def get_start_and_linker_metabolites():
    df = pd.read_excel("stats/exchange_metabolites.xlsx", engine="openpyxl")
    substrate = df.dropna(axis=0,subset=["Selected Substrate"])
    substrate = substrate[substrate["Media"] != "Basic"]
    substrate = substrate["id"].tolist()

    linker = df.dropna(axis=0,subset=["Intermediate"])
    linker = linker[linker["Media"] != "Basic"]
    linker = linker["id"].tolist()
    return substrate, linker

    
def get_metabolite_name_dict():
    fnames = open("stats/all_metabolites_info.tsv")
    metabolites_names_dict = {}
    line = fnames.readline()
    while True:
        line = fnames.readline()
        if len(line) < 2:
            break
        line = line.split("\t")
        name = line[1]
        metID = line[0][:-2]
        metabolites_names_dict[metID] = name
    fnames.close()
    return metabolites_names_dict


def rewrite():
    root_old = config["path_GEM"]
    root_new = config["path_rewrite"]
    suffix = config["suffix"]
    if not os.path.exists(root_new):
        os.mkdir(root_new)

    files = os.listdir(root_old)
    files = [x for x in files if x.endswith(suffix)]
    
    for file in tqdm(files, desc="rewriting GEMs for community design"):
        if os.path.exists(os.path.join(root_new, file)):
            continue
        gem = cobra.io.read_sbml_model(os.path.join(root_old, file))
        ID = gem.id
        for x in gem.metabolites:
            if x.id[-1]!='e':
                x.id = ID + '_' + x.id
        for x in gem.reactions:
            if x.id[:3]!='EX_':
                x.id = ID + '_' + x.id
        ## do not edit gene ids, to avoid conflict
        #for x in gem.genes:
        #    x.id = ID + '_' + x.id
        gem.repair()
        cobra.io.write_sbml_model(gem, os.path.join(root_new, file))


def get_good_bacteria(pkl=None):
    # bac_good: bac_good[path] = {list of tuples, e.g.(gcf,growth,growth*production)}
    if pkl:
        with open(pkl, "rb") as f:
            bac_good = pickle.load(f)
        return bac_good
    
    target = config["target"]
    target_row = "production" if target == "production" else "absorption"
    
    if os.path.exists(f"tmp/bac_good_{target}.pkl"):
        with open(f"tmp/bac_good_{target}.pkl", "rb") as f:
            bac_good = pickle.load(f)
        return bac_good

    bac_good = {}
    root_eval = f"tmp/bacteria_label_{target}/"
    fs = os.listdir(root_eval)
    with tqdm(fs) as pbar:
        pbar.set_description("reading bacteria labels")
        for bac in pbar:
            df = pd.read_csv(root_eval + bac)
            good_df = df[df["label"] == "good"]
            for idx,row in good_df.iterrows():
                tmp_path = (row["substrate"], row["product"], row["O2"], row["glucose"])
                if tmp_path not in bac_good:
                    bac_good[tmp_path] = [ (bac[:-4], row["growth"], row["growth"] * row[target_row]) ]
                else:
                    bac_good[tmp_path].append( (bac[:-4], row["growth"], row["growth"] * row[target_row]) )

    with open(f"tmp/bac_good_{target}.pkl", "wb") as f:
        pickle.dump(bac_good, f)
    return bac_good


# return: gem
def read_gem(file, medium=None, o2=True, glc=True):
    try:
        gem = cobra.io.read_sbml_model(file)
    except:
        print(file + ": not found")
        return None
    
    if medium:
        for rxn in gem.exchanges:
            rxt = rxn.reactants[0]
            if rxt.id not in medium.keys():
                rxn.lower_bound = 0 # do not uptake metabolites not in the medium
            else:
                rxn.lower_bound = -medium[rxt.id]
    
    if "EX_o2_e" in gem.reactions:
        gem.reactions.get_by_id("EX_o2_e").lower_bound = -10 if o2 else 0
    if "EX_glc__D_e" in gem.reactions:
        gem.reactions.get_by_id("EX_glc__D_e").lower_bound = -10 if glc else 0
    return gem


# Merge multiple GEMs into one, setting the objective function to the sum of each bacterium's biomass.
# Parameters: GEM_path (root path), names (list of all GEM file names), suffix (either .xml or .xml.gz).
#     If a medium is specified, O2 and glucose must also be defined; if medium is None, all substances are included by default.
# Return value: the merged model and the IDs of each bacterium's biomass function.
def merge_model(GEM_path, names, medium=None, O2=True, glucose=True, base_threshold=0.1):
    GEM = []
    sing = []
    suffix = config["suffix"]
    for name in names:
        try:
            GEM.append(read_gem(os.path.join(GEM_path, name+suffix), medium, o2=O2, glc=glucose))
            sing.append(GEM[-1].optimize().objective_value)
        except:
            print(os.path.join(GEM_path, name+suffix) + ": file not found")
            return None
    
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
        new_model.add_reactions(needcopyrxn)

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
        biomass = new_model.reactions.get_by_id(biomass_func[i])
        biomass.lower_bound = sing[i] * base_threshold
        obj[biomass] = 1.0
    new_model.objective = obj
    
    return new_model, biomass_func


# Draw phase plane
# GEM: cobra format Genome-scale Metabolic Model
# rxn1, rxn2: The reactions corresponding to the two axes of the phase plane, such as growth and production. Format: string.
# ref: The relative value of sampled points on the phase plane, calculated as (y/y_max) / (x/x_max).
# relative_ref: Determines whether ref uses relative or absolute proportions.
# axis_lim: Decides whether to control the axis limits to make (0,0) the origin.
# max_point: Determines whether to additionally mark the maximum value point.
# labels: Custom titles and names for the x and y axes.
# extra_points: A list of lists [x,y,text] used to mark additional points.
# example: draw_phase_plane(Se, "BIOMASS__1", "EX_sucr_e")
def draw_phase_plane(GEM, rxn1, rxn2, ref=[0.001, 0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.5, 2, 3.3, 10, 1000], relative_ref=True, axis_lim=True, max_point=True, labels=None, extra_points=None):
    from cobra import Metabolite
    mu = Metabolite(
        'mu_c',
        formula='X',
        name='mu',
        compartment='c')

    func1 = GEM.reactions.get_by_id(rxn1)
    func2 = GEM.reactions.get_by_id(rxn2)
    x,y = [], []
    
    
    x0, y0 = [], []
    GEM.objective = rxn1
    solution = GEM.optimize()
    x0.append(solution[rxn1])
    y0.append(solution[rxn2])
    
    GEM.objective = rxn2
    solution = GEM.optimize()
    x0.append(solution[rxn1])
    y0.append(solution[rxn2])
    
    x_lim, y_lim = x0[0], y0[1]
    if relative_ref:
        scale = x_lim / y_lim
        ref = [i * scale for i in ref]
    
    for coeff in ref:
        func1.add_metabolites({mu:1})
        func2.add_metabolites({mu:-coeff})
        solution = GEM.optimize()
        x.append(solution[rxn1])
        y.append(solution[rxn2])
        func1.subtract_metabolites({mu:1})
        func2.subtract_metabolites({mu:-coeff})
        
    import matplotlib.pyplot as plt
    if not labels:
        labels = ['Phase Plane', rxn1, rxn2]
    plt.figure(figsize=(5,5))
    plt.title(labels[0],fontsize=20)
    plt.xlabel(labels[1],fontsize=14)
    plt.ylabel(labels[2],fontsize=14)
    if axis_lim:
        plt.xlim((0, 1.1 * x0[0]))
        plt.ylim((0, 1.1 * y[1]))
    if max_point:
        plt.scatter(x0, y0, c="black", s=50)
    if extra_points:
        plt.scatter(extra_points[0], extra_points[1], c="black", s=50)
        for i in range(len(extra_points[0])):
            plt.annotate(extra_points[2][i], xy=(extra_points[0][i]+x_lim/100, extra_points[1][i]+y_lim/100), 
                         xycoords='data', xytext=(+30, +30), textcoords='offset points', 
                         fontsize=13, arrowprops=dict(arrowstyle='->', connectionstyle="arc3,rad=.2"))
    plt.plot(x,y,marker='o', color='blue',
        markersize=5, markeredgecolor='r')
    
    plt.show()