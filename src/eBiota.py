import os
import sys
import json
import pandas as pd
import argparse
from tqdm import tqdm
import subprocess

from parse_path import trans
from eBiota_utils import rewrite
from pre_evaluate import run_evaluate
from community_screen import get_combination
from community_design import run_community_design

__author__ = "Jiaheng Hou, Haoyu Zhang, Yulin Liao"
__copyright__ = "Copyright (c) 2024 Zhulab"
__license__ = "The MIT License (MIT)"
__version__ = "1.0"


def parse_arg():
    args = argparse.ArgumentParser(description="Ebiota parameters")
    args.add_argument(
        "--Function", "-F",
        default=None,
        type=str, # filter_prod_and_substrate
        help="Task mode choice.",
    )
    
    args.add_argument(
        "--input_csv",
        default=None,
        type=str,
        help="Path to csv file with all parameters ebiota may need."
    )
    args.add_argument(
        "--output",
        default=None,
        type=str,
        help="Path of the output root."
    )
    args = args.parse_args()
    return args

def main_production_target():
    # 用户给定底物+产物，设计菌群
    pass

def main(args):
    function_mode = args.Function
    with open('config.json', 'r') as file:
        config = json.load(file)
    if args.Function=="design":
        cmd = "perl find_path.pl " +  config["path_GEM"]
        suc = subprocess.call(cmd, shell=True)
        if suc == 0:
            print("CoreBFS completed!")
        else:
            print("An error occurred in CoreBFS. Please check the parameters!")
            return
        
        trans()
        rewrite()
        run_evaluate()
        get_combination()
        run_community_design()


if __name__ == "__main__":
    args = parse_arg()
    main(args)
    