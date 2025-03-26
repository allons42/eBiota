import os
import cobra
import pickle
import copy
import numpy as np
import pandas as pd
from random import shuffle
from tqdm import tqdm
from multiprocessing import Pool
from eBiota_utils import config

if config["USE_CUDA"] and config["CUDA_VISIBLE_DEVICES"] != "default":
    os.environ["CUDA_VISIBLE_DEVICES"] = config["CUDA_VISIBLE_DEVICES"]

from sklearn import metrics
import scipy.stats as stats
import torch
import torch.nn as nn
from itertools import permutations
from torch.utils.data import Dataset, DataLoader
import torch.nn.functional as F


#### parse_data_for_DeepCooc.py ####

# tools
def dump_pkl(obj, pkl_name):
    with open(pkl_name, "wb") as dest:  
        pickle.dump(obj, dest) 

def load_pkl(pkl_name):
    with open(pkl_name, "rb") as src: 
        res = pickle.load(src)
    return res

def list_txt(ls, txt_file):
    ls = set(ls)
    ls = list(ls)
    with open(txt_file, "w") as f:
        for info in ls:
            f.write(info)
            f.write("\n")

def txt_list(txt_file):
    with open(txt_file, 'r') as f:
        spes = f.readlines()
    
    spes = [x.strip('\n') for x in spes]
    return spes

def list_txt_microlist(ls, txt_file):
    with open(txt_file, "w") as f:
        for micro_list in ls:
            f.write(str(micro_list))
            f.write("\n")

def txt_list_microlist(txt_file):
    ls = []
    with open(txt_file, "r") as f:
        for line in f.readlines():
            line = line.strip("\n")
            line = line.strip("]").strip("[").split(",")
            t_ls = [x.strip().strip("'") for x in line]
            ls.append(t_ls)
    return ls


## configs
root_output = config["path_output"]
file_input_txts = [x for x in os.listdir(root_output) if x.endswith(".txt")] ## txt only
micro_num = config["community_size"]
data_taxa_type = config["data_taxa_type"]
feature_index_type = config["feature_index_type"]
temp_dir = os.path.join("tmp", "cooc")
path_GEM = config["path_GEM"]

# root for 2w micro，所有rea名称到反应对象的pairs
stats_dir = './stats/'
rea_table = os.path.join(stats_dir, "reaction_table_carveme_2w.pkl")
rea_table = load_pkl(rea_table)
# 换用transport反应的列表
exchange_reaction_table = f'{stats_dir}/reac_transport_table.tsv'
transportReac_list = pd.read_csv(exchange_reaction_table, sep='\t')
transportReac_list = list(transportReac_list['reaction_nm'])

# carveme table with species_ids: 修改为ncbi上面获取的genus name--1205
carveme_list_speciesID = f"{stats_dir}/carveme_gem_taxonomy_ncbi_1205.csv"
carveme_list_WithspeciesID = pd.read_csv(carveme_list_speciesID)

# 存放每个gcf对应的raaction信息的table
gcfs_table = txt_list(
    f"{stats_dir}/np_rea_table_gcf.txt"
)

## date: 2023.1.3 总的来说分成俩类别：使用met index和使用reaction index的
save_suffix = feature_index_type
if feature_index_type == 'metflux_lsa_256': 
    gem_reaction_index_table = np.load(
        f"{stats_dir}/gcfs_MetList_lsa_256.npy"
    ) # 改成 2w, 256的矩阵了
    single_reac_feature_num = 1
elif feature_index_type == 'metflux_ori_431':
    gem_reaction_index_table = np.load(
        f"{stats_dir}/gcfs_MetList_original_431.npy"
    ) # 改成 2w, 431的矩阵了
    single_reac_feature_num = 1
elif feature_index_type == 'moreReacFeature':
    gem_reaction_index_table = np.load(
        f"{stats_dir}/gems_np.npy"
    ) # 改成 2w, 353, 14的矩阵了
    single_reac_feature_num = 14
else:
    pass # 最开始的处理 0， 1向量

## 主要转换到index的funcs
# 从gcf直接到index
def get_gcf_rea_index(
    gcf="", rea_table=None, extra_gems_root=path_GEM # external gems directory
):
    ## given gcf, return reaction vector or reduced reaction vector using SVD
    try:
        gcf_id = gcfs_table.index(gcf)
        rea_index = gem_reaction_index_table[gcf_id, :]
    except:
        extra_gems_list = os.listdir(extra_gems_root)
        extra_gems_list = [x for x in extra_gems_list if (gcf in x)]
        if len(extra_gems_list) == 0:
            ## miss found, using mean vector to replace
            rea_index = gem_reaction_index_table.mean(axis=0)
            print(f"not found rea index. {gcf}")
        else:
            rea_index = np.zeros((gem_reaction_index_table.shape[1]))
            gem = extra_gems_list[0]
            gem = os.path.join(extra_gems_root, gem)

            g = cobra.io.read_sbml_model(gem)
            for rea in g.reactions:
                ids = rea.id
                rea_ind = rea_table.index(ids)
                rea_index[rea_ind] = 1

    return rea_index

# 2022.12.1: using species_id or genus_nm to covert reactions to index
def genus_species_ave_index(
    species="123",
    genus_species_id_type="genus",
    carveme_list=None,
    rea_table=None,
    species_ave_gems_num=15,
):
    species = str(species)
    carveme_list["species_id"] = carveme_list["species_id"].apply(lambda x: str(x))

    # 根据species或者genus选择
    df_species = carveme_list[carveme_list[genus_species_id_type] == species]
    df_species.reset_index(inplace=True, drop=True)
    species_gems = list(df_species["assembly"])

    # get gems from the species id and trans to gem to get reaction index
    shuffle(species_gems)
    if species_ave_gems_num < len(species_gems):
        species_gems = species_gems[:species_ave_gems_num]

    n_rea =  gem_reaction_index_table.shape[1] # len(rea_table)
    n_gems = len(species_gems)

    # 每一个index不再是一个单纯0 1的 
    index = np.zeros(shape=(n_gems, n_rea, single_reac_feature_num), dtype=np.float32)
    for i, gcf in enumerate(species_gems):
        rea_index = get_gcf_rea_index(gcf)
        if len(rea_index.shape) < 2:
            rea_index = rea_index[: , np.newaxis]
        index[i, :] = rea_index

    ## average and flatten
    index = index.mean(axis=0)
    index = index[np.newaxis, :]
    if np.isnan(index).sum() > 0:
        default_ind = np.zeros(shape=(n_rea, single_reac_feature_num), dtype=np.float32)  # gem_reaction_index_table.mean(axis=0)
        default_ind = default_ind[np.newaxis, :]
        index = default_ind # 返回一个default的参数，这个genus没有找到
    return index


def str_com(community="['mycoplasma', 'finegoldia']"):
    # 11.14:
    str = eval(community)
    return str


def com_index(
    community=["mycoplasma", "finegoldia"],
    genus_species_id_type="genus",
    carveme_list=None,
    rea_table=None,
    multi_process_ind=0,
):
    str = community
    coms = []
    for genus_species_id in str:
        if genus_species_id_type == 'gcf':
            # 直接用gcf作为输入情况
            ind = get_gcf_rea_index(
            gcf=genus_species_id,
            rea_table=rea_table,
            )
            ind = ind[np.newaxis, :]
        else:
            # 输入时species或者genus名称
            ind = genus_species_ave_index(
                species=genus_species_id,
                genus_species_id_type=genus_species_id_type,
                carveme_list=carveme_list,
                rea_table=rea_table,
            )
        coms.append(ind)

    return np.concatenate(coms, axis=0), multi_process_ind

# 主要的class: 训练时候为了扩增数据增加的功能都去掉了
class EMPdataset:
    def __init__(
        self,
        label: int,
        csv_file=None,
        cut_size=2,  ## 这一项参数输入micro_num
        dataset_from=0,
        part_set=-1,  ## -1代表不开启part模式
        data_taxa_type="gcf", ## "genus"
    ):
        self.label = label

        # dataset, element of data is [gem1, gem2]
        self.micro_num = cut_size
        self.dataset_from = dataset_from
        self.data_taxa_type = data_taxa_type
        # 读入数据，后续根据数据格式再修改
        self.data_df = pd.read_csv(csv_file, sep='\t', index_col=False)
        # self.data_df.reset_index(inplace=True, drop=True)
        
        self.data = []
        for i in range(self.data_df.shape[0]):
            t_bac_com = []
            for micro_n in range(self.micro_num):
                # 为了处理n菌的数据
                t_bac = self.data_df.loc[i, f'Bac{micro_n+1}']
                t_bac = t_bac[: -2] + '.' + t_bac[-1]
                t_bac_com.append(t_bac)
            self.data.append(t_bac_com)

        # 对于非常大的数据量，开启并行
        single_part_samples = len(self.data) // 15

        if part_set >= 0:
            self.data = self.data[
                part_set
                * single_part_samples : min(
                    (part_set + 1) * single_part_samples, len(self.data)
                )
            ]  # 加上限制

        self.micro_list = self.data.copy()

        print(
            f"This is a type {str(label)} dataset. There are {len(self.data)} micro samples.\nA peek of data is {self.data[0]}."
        )

    def gems(
        self,
    ):
        gems_ls = []
        for micro_com in self.data:
            gems_ls += micro_com

        gems_ls = set(gems_ls)
        gems_ls = list(gems_ls)

        self.gems_ls = gems_ls
        return gems_ls

    def data_tokenize(self, reaction_table=rea_table[0], pool_num=32):
        if pool_num > 1:
            pool = Pool(processes=pool_num)

        # if self.dataset_from == 0:
            # deprecated
            # self.data = [gem2index(x, reaction_table) for x in tqdm(self.data)]

        # 目前的最优方法--从已经做好的table里面直接读取rea的index
        # 每一个sample会搞出来多个样本点
        data_res = [-1] * len(self.data)
        for i, x in enumerate(tqdm(self.data)):
            
            # 原本的方法：按照genus水平：carveme_list_WithspeciesID 里面也包含了genus这一列
            task = pool.apply_async(
                com_index,
                (
                    x,
                    self.data_taxa_type,
                    carveme_list_WithspeciesID,
                    transportReac_list,  #rea_table[0],
                    i,
                ),
            )

            data_np, data_ind = task.get(timeout=10000)
            data_res[data_ind] = data_np
        self.data = data_res

## 一些主要调用的函数
def main_dataset(np_pkl="", save_pkl="size2_dataset.pkl", emps=(), micro_num=3, temp_dir = './temp'):
    # 加入功能：保存
    if len(np_pkl) > 1:
        # 这个时候是已经存储完成，只是读取
        try:
            total_data, labels = load_pkl(np_pkl)
        except:
            dir, basename = os.path.split(np_pkl)
            total_data, labels = np.load(os.path.join(dir, "x_" + basename)), np.load(
                os.path.join(dir, "y_" + basename)
            )
        return total_data, labels

    # 第一次构建数据集，并且保存为pkl
    total_data = []
    labels = []

    micro_genus = []

    for emp in emps:
        t_data = emp.data  # already array
        t_label = [emp.label] * len(t_data)
        t_micro_genus = emp.micro_list

        total_data.append(t_data)
        labels += t_label
        micro_genus += t_micro_genus

    total_data = np.concatenate(
        total_data, axis=0
    )  # np.array(total_data, dtype=np.int)# 为啥是int？
    labels = np.array(labels, dtype=np.float32)

    print(f"Data shape is {total_data.shape}.")
    
    # 保存数据，这里新建一个temp文件夹作为中间结果的保存地方
    pkl_micro_now_root = temp_dir
    os.makedirs(pkl_micro_now_root, exist_ok=True)
    np.save(os.path.join(pkl_micro_now_root, "x_" + save_pkl), total_data)
    np.save(os.path.join(pkl_micro_now_root, "y_" + save_pkl), labels)

    # list_txt_microlist(
    #     micro_genus, os.path.join(pkl_micro_now_root, save_pkl.replace(".npy", ".txt"))
    # ) 不再保存这个信息
    return total_data, labels

## 主要接口函数
def build_data_for_DLmodel(file_input, part, micro_num, data_taxa_type, temp_dir='./temp'):
    data_ = EMPdataset(
            label=0,
            csv_file=file_input,
            cut_size=micro_num,
            dataset_from=1,
            part_set=part,
            data_taxa_type=data_taxa_type,
        )
    
    data_.data_tokenize()
    main_dataset(
            save_pkl=
            os.path.basename(file_input).split('.txt')[0]
            + f"_{data_taxa_type}_dataset_{save_suffix}_micro_{micro_num}.npy",
            emps=(data_, ),
            micro_num=micro_num,
            temp_dir=temp_dir
        )
    # 20240527: 改成相应的文件名打头


#### dataset.py ####
# 数据归一化: 已经弃用
gcfs_all_gem_metList_table = np.load(
    f"{stats_dir}/gcfs_MetList_original_431.npy"
)
mu = np.mean(gcfs_all_gem_metList_table, axis=0)
sigma = np.std(gcfs_all_gem_metList_table, axis=0)

class DatasetFolder(Dataset):
    def __init__(
        self,
        micro_pairs_path="",
        micro_pairs_co_occur_path="",
        data_process="",
        process_transformer=None,
        val_split=0.05,
        exchange_reaction_table=None,
        logger=None,
        micro_pairs_path_metInd=None,
    ):
        # 注意此时输入的path是已经分好的用于train的数据路径
        self.val_split = val_split
        self.logger = logger
        self.total_data, self.labels = np.load(micro_pairs_path), np.load(
            micro_pairs_co_occur_path
        )
        if micro_pairs_path_metInd is not None:
            # 第二路加载met相关的特征
            self.total_data_metInd = np.load(micro_pairs_path_metInd)

        else:
            self.total_data_metInd = None

        try:
            self.total_data_metInd = np.squeeze(self.total_data_metInd, axis=3)
            # 如果只有一维特征，在这里会拉直
            self.total_data = np.squeeze(self.total_data, axis=3)
        except:
            pass

        if len(self.total_data.shape) == 4:
            # 简单只取前面4维度的特征：针对reac index的简单处理
            # date：1.4 不再拉直，对每一维度特征单独处理
            self.total_data = self.total_data[:, :, :, :4]
            
        if self.total_data_metInd is not None:
            # met的index需要标准化
            self.total_data_metInd = (self.total_data_metInd - mu) / sigma
            self.met_dim = self.total_data_metInd.shape[2]
            self.met_num = 1
        else:
            self.met_dim = 0
            self.met_num = 0

        # 提取模型架构相关的信息
        self.input_dim = self.total_data.shape[2]
        self.micro_num = self.total_data.shape[1]
        self.reaction_num = self.total_data.shape[3]

        # 仅仅保留exchange反应，这一步不需要了，因为输出的data已经只有ex反应
    
        # To tensor
        self.total_data, self.labels = torch.Tensor(self.total_data),  torch.Tensor(
            self.labels
        )
        if self.total_data_metInd is not None:
            self.total_data_metInd = torch.Tensor(self.total_data_metInd)
        
        print(
            f"Dataset loaded, shape of X is {self.total_data.size()} and y is {self.labels.shape}."
        )
        if logger is not None:
            logger.info(
                f"Dataset loaded, shape of X is {self.total_data.size()} and y is {self.labels.shape}."
            )

    def __getitem__(self, index):
        # 随机sample micro的顺序
        idx_micro = list(range(self.micro_num))

        x_merged = []
        x = self.total_data[index]
        y = self.labels[index]
        for idx_order in permutations(idx_micro, self.micro_num):
            idx_order = np.array(list(idx_order))
            x_t = x[
                idx_order,
            ]
            x_merged.append(x_t[np.newaxis, :, :])

        x_merged = np.concatenate(x_merged)

        if self.total_data_metInd is not None:
            # 额外输出met信息
            x_merged_metInd = []
            x = self.total_data_metInd[index]
            for idx_order in permutations(idx_micro, self.micro_num):
                idx_order = np.array(list(idx_order))
                x_t = x[
                    idx_order,
                ]
                x_merged_metInd.append(x_t[np.newaxis, :, :])

            x_merged_metInd = np.concatenate(x_merged_metInd)
        else:
            x_merged_metInd = x_merged

        return x_merged, x_merged_metInd, y

    def __len__(self):
        return self.total_data.shape[0]

    def split_validation(
        self,
    ):
        n_samples = self.total_data.shape[0]
        n_val = int(self.val_split * n_samples)
        # idx shuffle
        idx_full = np.arange(n_samples)
        np.random.shuffle(idx_full)

        val_idx = idx_full[:n_val]
        train_idx = idx_full[n_val:]
        # 复制val data
        val_dataset = copy.deepcopy(self)
        val_dataset.total_data = val_dataset.total_data[val_idx, :, :]
        if self.total_data_metInd is not None:
            val_dataset.total_data_metInd = val_dataset.total_data_metInd[val_idx, :, :]
        val_dataset.labels = val_dataset.labels[val_idx]

        self.total_data = self.total_data[train_idx, :, :]
        if self.total_data_metInd is not None:
            self.total_data_metInd = self.total_data_metInd[train_idx, :, :]
        self.labels = self.labels[train_idx]

        self.logger.info(val_dataset.total_data.shape, self.total_data.shape)
        return val_dataset


# 获取dataloader 主函数
def get_dataloader(
    micro_pairs_path="",
    micro_pairs_co_occur_path="",
    mode="test",
    batch_size=10,
    data_process="",
    process_transformer=None,
    exchange_reaction_table=None,
    logger=None,
    micro_pairs_path_metInd=None,
):
    if mode == "test":
        dataset = DatasetFolder(
            micro_pairs_path=micro_pairs_path,
            micro_pairs_co_occur_path=micro_pairs_co_occur_path,
            data_process=data_process,
            process_transformer=process_transformer,
            exchange_reaction_table=exchange_reaction_table,
            micro_pairs_path_metInd=micro_pairs_path_metInd,
        )
        train_dataloader = DataLoader(
            dataset,
            batch_size=batch_size,
            shuffle=False,
            num_workers=8,
            pin_memory=True,
        )
        return (
            train_dataloader, 
            dataset.input_dim,
            dataset.reaction_num,
            dataset.met_dim,
            dataset.met_num,
            )
    else:
        dataset = DatasetFolder(
            micro_pairs_path=micro_pairs_path.replace(".npy", "_training.npy"),
            micro_pairs_co_occur_path=micro_pairs_co_occur_path.replace(
                ".npy", "_training.npy"
            ),
            data_process=data_process,
            process_transformer=process_transformer,
            exchange_reaction_table=exchange_reaction_table,
            logger=logger,
            micro_pairs_path_metInd=micro_pairs_path_metInd.replace(
                ".npy", "_training.npy"
            )if micro_pairs_path_metInd is not None else None,
        )

        # val_dataset = dataset.split_validation()
        val_dataset = DatasetFolder(
            micro_pairs_path=micro_pairs_path.replace(".npy", "_validation.npy"),
            micro_pairs_co_occur_path=micro_pairs_co_occur_path.replace(
                ".npy", "_validation.npy"
            ),
            data_process=data_process,
            process_transformer=process_transformer,
            exchange_reaction_table=exchange_reaction_table,
            logger=logger,
            micro_pairs_path_metInd=micro_pairs_path_metInd.replace(
                ".npy", "_validation.npy"
            ) if micro_pairs_path_metInd is not None else None,
        )

        train_dataloader = DataLoader(
            dataset,
            batch_size=batch_size,
            shuffle=True,
            num_workers=16,
            pin_memory=True,
        )
        val_dataloader = DataLoader(
            val_dataset,
            batch_size=batch_size,
            shuffle=False,
            num_workers=1,
            pin_memory=True,
        )
        return (
            train_dataloader,
            val_dataloader,
            dataset.process_transformer,
            dataset.input_dim,
            dataset.reaction_num,
            dataset.met_dim,
            dataset.met_num,
        )


#### metric_loss ####
def accuracy(output, target):
    with torch.no_grad():
        # pred = torch.argmax(output, dim=1)
        pred = (output > 0.5).detach().type(torch.int32)

        assert pred.shape[0] == len(target)
        correct = 0
        correct += torch.sum(pred == target).item()
    return correct / len(target)

def auc_score(output, target):
    output, target = output.cpu().detach().numpy(), target.cpu().detach().numpy()
    fpr, tpr, _ = metrics.roc_curve(target, output)
    auc_res = metrics.auc(fpr, tpr)

    return auc_res

def focal_loss(output, target):
    alpha = 1.0
    gamma = 3.0

    ce_loss = F.cross_entropy(
        output, target, reduction="none"
    )  # important to add reduction='none' to keep per-batch-item loss
    pt = torch.exp(-ce_loss)  # ce loss 实际是log(softmax)所以直接exp就可以了
    Focal_loss = (alpha * (1 - pt) ** gamma * ce_loss).mean()  # mean over the batch

    return Focal_loss


def cross_entropy(output, target):
    return F.cross_entropy(
        output,
        target,
    )  # ignore_index=PAD_IDX)


#### model.py ####
class GEMreactionNetFc(nn.Module):
    def __init__(self, input_dim=4291, fc_out_dim=128, dropout_rate=0.5) -> None:
        super().__init__()
        self.dropout = nn.Dropout(p=dropout_rate)
        self.linears = nn.Sequential(
            nn.Linear(input_dim, 4096),  # 4096
            nn.Dropout(p=dropout_rate),
            nn.LeakyReLU(negative_slope=0.1),
            nn.Linear(4096, 2048),
            nn.Dropout(p=dropout_rate),
            nn.LeakyReLU(negative_slope=0.1),
            nn.Linear(2048, 1024),
            nn.Dropout(p=dropout_rate),
            nn.LeakyReLU(negative_slope=0.1),
            nn.Linear(1024, fc_out_dim * 2),
            nn.Dropout(p=dropout_rate),
            nn.LeakyReLU(negative_slope=0.1),
            nn.Linear(fc_out_dim * 2, fc_out_dim),
        )

    def forward(self, x):
        x = torch.squeeze(x)
        out = self.linears(x)
        return out

class MicroCoexistenceNetCore(nn.Module):
    # 修改合并3channel的方式
    def __init__(
        self,
        micro_num=3,
        input_dim=4291,
        single_out_dim=128,
        model_mode="FC",
        dropout_rate=0.5,
    ):
        # 这个函数会把x=（-1, micro_num, feature_dim）的输入映射到y=(-1, single_out_dim)的输出
        super().__init__()

        self.micro_num = micro_num
        self.single_out_dim = single_out_dim

        if model_mode.startswith("CNN"):
            pass
        elif model_mode.startswith("FC"):
            self.single_path_net = GEMreactionNetFc(
                input_dim=input_dim,
                fc_out_dim=single_out_dim,
                dropout_rate=dropout_rate,
            )

        self.dropout = nn.Dropout(p=dropout_rate)

        # merge three channels of microbiome
        self.merge_cnn = nn.Conv1d(
            in_channels=self.micro_num, out_channels=1, kernel_size=1, stride=1
        )

    def forward(self, x):
        x_zeros = torch.zeros(
            size=(x.size(0), self.micro_num, self.single_out_dim), device=x.device
        )
        for micro_now in range(self.micro_num):
            x_input = x[:, micro_now, :]
            x_input = torch.unsqueeze(x_input, dim=1)
            x_single_path_out = self.single_path_net(x_input)

            x_zeros[:, micro_now, :] = x_single_path_out

        x_zeros = self.dropout(x_zeros)
        x_merged = self.merge_cnn(x_zeros)
        x_merged = torch.squeeze(x_merged, dim=1)
        return x_merged

class MicroCoexistenceNetv3(nn.Module):
    # 会一次性输入met和reac两个通道
    def __init__(
        self,
        micro_num=3,
        reaction_dim=4291,
        reaction_num=4,
        met_dim=431,
        met_num=0,
        single_out_dim=128,
        model_mode="FC",
        dropout_rate=0.5,
    ):
        super().__init__()

        self.micro_num = micro_num
        self.single_out_dim = single_out_dim
        self.met_num = met_num
        self.reaction_num = reaction_num
        # reaction相关的net
        if reaction_num > 0:
            self.reac_nets = []
            for i in range(reaction_num):
                reac_net_now = MicroCoexistenceNetCore(
                    micro_num=micro_num,
                    input_dim=reaction_dim,
                    single_out_dim=single_out_dim,
                    model_mode=model_mode,
                    dropout_rate=dropout_rate,
                )
                self.reac_nets.append(reac_net_now)
            self.reac_nets = nn.ModuleList(self.reac_nets)

        if met_num > 0:
            self.met_net = MicroCoexistenceNetCore(
                micro_num=micro_num,
                input_dim=met_dim,
                single_out_dim=single_out_dim,
                model_mode=model_mode,
                dropout_rate=dropout_rate,
            )

        self.fc_to_class_label = nn.Sequential(
            nn.Linear(single_out_dim * (reaction_num + met_num), 1024),  # 4096
            nn.Dropout(p=dropout_rate),
            nn.LeakyReLU(negative_slope=0.1),
            nn.Linear(1024, 128),
            nn.Dropout(p=dropout_rate),
            nn.LeakyReLU(negative_slope=0.1),
            nn.Linear(128, 1),
        )
        self.dropout = nn.Dropout(p=dropout_rate)

        # 输出层
        self.out_layer = nn.Sigmoid()

        # merge three channels of microbiome
        self.merge_cnn = nn.Conv1d(
            in_channels=self.micro_num, out_channels=1, kernel_size=1, stride=1
        )

        # initialize model
        for m in self.modules():
            if isinstance(m, (nn.Conv1d, nn.Linear)):
                nn.init.xavier_normal_(m.weight)
                # nn.init.kaiming_normal_(m.weight, mode='fan_out', nonlinearity='relu')
            elif isinstance(m, (nn.BatchNorm1d, nn.GroupNorm)):
                nn.init.constant_(m.weight, 1)
                nn.init.constant_(m.bias, 0)

    def forward(self, x_reac, x_met):
        x_merged = []
        # 先暂时只加入reaction的特征
        for i in range(self.reaction_num):
            x_out_now = self.reac_nets[i](x_reac[:, :, :, i])
            x_merged.append(x_out_now)
        # 加入met的特征融合
        if self.met_num > 0:
            x_out_met = self.met_net(x_met)
            x_merged.append(x_out_met)
        
        x_merged = torch.cat(x_merged, 1)

        x_out = self.fc_to_class_label(x_merged)
        x_out = self.out_layer(x_out)
        x_out = x_out.squeeze(dim=-1)
        return x_out

class MicroCoexistenceNet_hy(nn.Module):
    # 修改合并3channel的方式
    def __init__(self, micro_num=3, main_model=None) -> None:
        super().__init__()
        self.single_permutaion_model = main_model
        self.micro_num = micro_num

    def forward(self, x, x_met):
        # x shape: (n, 6, 3, 400)
        permutation_num = x.size(1)
        x_out_temp = []
        for per_ord in range(permutation_num):
            x_per_ord = x[:, per_ord, :, :]
            x_met_per_ord = x_met[
                :,
                per_ord,
                :,
            ]

            x_out_per_ord = self.single_permutaion_model(x_per_ord, x_met_per_ord)
            x_out_per_ord = torch.unsqueeze(x_out_per_ord, dim=1)
            x_out_temp.append(x_out_per_ord)
        x_out_temp = torch.cat(x_out_temp, dim=1)

        # 分别计算均值和MSE为
        x_coexist_res = x_out_temp.mean(dim=1)
        x_permutation_sim_res = x_out_temp - torch.unsqueeze(x_coexist_res, dim=1)
        x_permutation_sim_res = (x_permutation_sim_res * x_permutation_sim_res).sum(
            dim=1
        )

        return x_coexist_res, x_permutation_sim_res

# 10.5: BCE loss
# 获取conv之后的data 维度的函数
def get_out_data_dim_2(data_dim, micro_num, model):
    test_x = torch.rand(10, micro_num, data_dim)
    test_out = model(test_x)
    return test_out.size(2) * test_out.size(1)


#### other_filter_label.py ####
# induce six types of interaction 
p_threshold = 0.1
interaction_type = ['competition (−/−)', 'mutualism (+/+)', 'neutralism (0/0)', 'parasitism (±)', 'commensalism (+/0)', 'amensalism (-/0)']
def check_single_pair(single_od_list, pair_od_list, threshold = 0.1, hyposis_test = False):
    if hyposis_test:
        stat,p_value = stats.mannwhitneyu(single_od_list,pair_od_list,alternative='less')
        if p_value < p_threshold:
            return 1
        
        stat,p_value = stats.mannwhitneyu(single_od_list,pair_od_list,alternative='greater')
        if p_value < p_threshold:
            return -1
    
        return 0
            
            
    single_od = np.mean(single_od_list)
    pair_od = np.mean(pair_od_list)
    ## 报除零错误
    if single_od:
        diff = pair_od / single_od - 1
    else:
        # divide zero
        diff = pair_od - single_od
        
    if diff > threshold:
        return 1
    elif diff < - threshold:
        return -1
    else:
        return 0

def get_interaction_type(od1 = ([0.05, 0.05], [0.1, 0.1]), od2 = ([0.05, 0.05], [0.1, 0.1]), threshold = 0.1, hyposis_test = False):
    a, b = check_single_pair(* od1, threshold=threshold, hyposis_test = hyposis_test), check_single_pair(* od2, threshold=threshold, hyposis_test = hyposis_test)
    
    if a < 0 and b < 0:
        return interaction_type[0], 0
            
    if a > 0 and b > 0:
        return interaction_type[1], 1
    
    if a == 0 and b == 0:
        return interaction_type[2], 2
    
    if int(a * b) == -1:
        return interaction_type[3], 3
        
    if (a > 0 and b == 0) or (a == 0 and b > 0):
        return interaction_type[4], 4
    
    if (a < 0 and b == 0) or (a == 0 and b < 0):
        return interaction_type[5], 5
    
    return 'unknown_type', -1

def handle_single_inputfile(file_in = '', file_out = ''):
    ## 注意：这里的interaction之类的函数只能针对2个bac的情况使用
    base_name = os.path.basename(file_in)
    base_name = base_name.split('.')[0]
    df_main = pd.read_csv(file_in, sep='\t', index_col=False)
    
    # 后处理6种interaction
    def get_6_interaction_type(x):
        spe1_list, spe2_list = (x['Bac1_mono_growth'], x['Growth1']), (x['Bac2_mono_growth'], x['Growth2'])
        interact_name, _ = get_interaction_type(spe1_list, spe2_list, threshold=0.01, hyposis_test = False)
        return interact_name
    
    if config["community_size"] == 2:
        df_main['Interaction_type'] = df_main[['Growth1', 'Growth2', 'Bac1_mono_growth', 'Bac2_mono_growth']].apply(lambda x: get_6_interaction_type(x), axis=1)
    
    # 后处理eveness
    def eveness_index(x):
        p1 = x['Growth1'] + 1e-6 
        p2 = x['Growth2'] + 1e-6
        H = - (p1 * np.log(p1) + p2 * np.log(p2))
        eveness = H / np.log(2)
        return eveness
    # df_main['eveness_index'] = df_main[['Growth1', 'Growth2', 'Bac1_mono_growth', 'Bac2_mono_growth']].apply(lambda x: eveness_index(x), axis=1)
    
    df_main.to_csv(file_in if len(file_out) < 1 else file_out, sep='\t', index=False)
    # shutil.rmtree(txt_dir)


# 工具应用时的主函数
def DeepCooc():
    for file_input in file_input_txts:
        # 每个txt文件都跑一下预处理
        build_data_for_DLmodel(
            file_input = f'{root_output}/{file_input}', 
            part = -1, # 这个参数不用了
            micro_num = int(micro_num), 
            data_taxa_type = data_taxa_type,
            temp_dir = temp_dir
            )
        
        # 运行deepCooc进行计算
        file_nm = file_input.split('.')[0]
        test_x = f'{temp_dir}/x_{file_nm}_{data_taxa_type}_dataset_moreReacFeature_micro_{micro_num}.npy'
        test_y = f'{temp_dir}/y_{file_nm}_{data_taxa_type}_dataset_moreReacFeature_micro_{micro_num}.npy'
        dl_model_path = f"./stats/model_weights_{micro_num}bac.pth"
        # training strategy
        batch_size = 128
        device = torch.device("cuda" if config["USE_CUDA"] and torch.cuda.is_available() else "cpu")
        zero_class_ratio = 0.5

        # model structure
        exchange_reaction_table = None
        arch = "" 
        model_mode = 'FC_mergedCNN'
        single_out_dim = 256
        dropout_rate = 0.5

        model_save_root = root_output
        os.makedirs(model_save_root, exist_ok=True) 
        (
            test_data_loader,
            input_dim,
            reaction_num,
            met_dim,
            met_num,
        ) = get_dataloader(
            test_x,
            test_y,
            "test",
            batch_size=batch_size,
            data_process=arch,
            exchange_reaction_table=exchange_reaction_table,
            logger=None,
            micro_pairs_path_metInd=None,
        )

        # build model architecture, then print to console
        model_main = MicroCoexistenceNetv3(
            micro_num=micro_num,
            single_out_dim=single_out_dim,
            reaction_dim=input_dim,
            reaction_num=reaction_num,
            met_dim=met_dim,
            met_num=met_num,
            model_mode=model_mode,
            dropout_rate=dropout_rate,
        )
        model = MicroCoexistenceNet_hy(micro_num=micro_num, main_model=model_main)

        # model、loss、metric
        model = model.to(device)
        model.load_state_dict(torch.load(dl_model_path))

        # 运行模型
        size = len(test_data_loader)
        y_pred = []
        for batch, (X, X_metInd, y) in enumerate(test_data_loader):
            print(f'{batch}/{size}.....')
            
            # Compute prediction and loss
            X, X_metInd, y = X.to(device), X_metInd.to(device), y.to(device).long()
            y_pred_t, _ = model(X, X_metInd)
            y_pred_np = y_pred_t.cpu().detach().numpy()
            
            y_pred.append(y_pred_np)
        y_pred = np.concatenate(y_pred, axis=0)
        
        # 保存结果
        input_df = pd.read_csv(f'{root_output}/{file_input}', sep='\t', index_col=False)
        input_df['DeepCooc_Co_occurrence'] = y_pred
        input_df['DeepCooc_Co_occurrence'] = input_df['DeepCooc_Co_occurrence'].apply(lambda x: 1 if x > 0.5 else 0)
        out_file = os.path.join(model_save_root, os.path.basename(f'{root_output}/{file_input}').replace('.txt', '.tsv'))
        input_df.to_csv(out_file, sep='\t', index=False)
        
        # 删除中间文件
        os.remove(test_x)
        os.remove(test_y)
        os.system(f'rm -r {root_output}/{file_input}') # rm the original file
        
        # 获取其他指标
        handle_single_inputfile(out_file, out_file)
        

if __name__ == '__main__':
    DeepCooc()