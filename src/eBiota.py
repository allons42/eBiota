import argparse
import subprocess
import os
import shutil

from parse_path import translate_path
from eBiota_utils import rewrite, config, update_config
from pre_evaluate import run_evaluate
from community_screen import get_combination
from community_design import run_community_design
# from community_check import DeepCooc

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
        "--outdir",
        default=None,
        type=str,
        help="Path of the output root."
    )
    args = args.parse_args()
    return args


def main(args):
    function_mode = args.Function
    if function_mode == "design":
        print("Start to analyze the input GEMs...")
        cmd = "perl find_path.pl " +  config["path_GEM"]
        suc = subprocess.call(cmd, shell=True)
        if suc == 0:
            print("CoreBFS completed!")
        else:
            print("An error occurred in CoreBFS. Please check the parameters!")
            return
        
        translate_path()
        print("All path detected! Start to evaluate the GEMs...")

        rewrite()
        run_evaluate()
        print("Evaluation completed! Start to build communities...")

        get_combination()
        run_community_design()
        #DeepCooc()
        print("All done!")
    
    if args.Function == "test":
        get_combination()
        run_community_design()

if __name__ == "__main__":
    args = parse_arg()
    update_config(args)
    """
    if os.path.exists("tmp"):
        shutil.rmtree("tmp")
    os.mkdir("tmp")
    """
    main(args)
    