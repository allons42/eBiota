import argparse
import subprocess
import os
import re

from parse_path import run_CoreBFS, translate_path
from eBiota_utils import rewrite, config, check_config, update_config
from pre_evaluate import run_evaluate
from community_screen import get_combination
from community_design import run_community_design
from community_mod import gene_mod

__author__ = "Jiaheng Hou, Haoyu Zhang, Yulin Liao"
__copyright__ = "Copyright (c) 2024 Zhulab"
__license__ = "The MIT License (MIT)"
__version__ = "1.1"


def parse_arg():
    args = argparse.ArgumentParser(description="Ebiota parameters")
    args.add_argument(
        "--Function", "-F",
        default=None,
        type=str, # filter_prod_and_substrate
        help="Task mode choice.",
    )
    
    args.add_argument(
        "--outdir", "-o",
        default=None,
        type=str,
        help="Path of the output directory."
    )

    args.add_argument(
        "--input_tsv", "-i",
        default=None,
        type=str,
        help="Path of the input TSV file, used for gene modification."
    )
    args = args.parse_args()
    return args


def main(args):
    if not os.path.exists("tmp"):
        os.makedirs("tmp")
    function_mode = args.Function

    if function_mode == "doall":
        run_CoreBFS()
        translate_path()
        print("Start to evaluate the GEMs...")
        rewrite()
        run_evaluate()
        print("Start to build communities...")
        get_combination()
        run_community_design()

        from community_check import DeepCooc
        DeepCooc()
        print("All done!")

    elif function_mode == "evaluate":
        run_CoreBFS()
        translate_path()
        print("Start to evaluate the GEMs...")
        rewrite()
        run_evaluate()

    elif function_mode == "design":
        get_combination()
        run_community_design()
        print("Design is completed!")
    
    elif function_mode == "cooc":
        from community_check import DeepCooc
        DeepCooc()
        print("DeepCooc co-occurrence prediction is done!")

    elif args.Function == "genemod":
        print("Start to analyze possible gene modifications. This may take a while...")
        gene_mod(args.input_tsv)
        print("Gene modification is completed!")
        
    elif args.Function == "test":
        passflag = True
        try:
            # use subprocess.run to execute Perl
            result = subprocess.run(['perl', '-v'], stdout=subprocess.PIPE, universal_newlines=True)
            version_match = re.search(r'(v\d+\.\d+\.\d+)', result.stdout).group(1)
            print("Perl version:", version_match)
        except Exception as e:
            print(f"Error: Failed to get Perl version: {e}")
            passflag = False
        try:
            import torch
            print("Pytorch version:", torch.__version__)
            if config["USE_CUDA"]:
                if torch.cuda.is_available():
                    print("CUDA version:", torch.version.cuda)
                    if config["CUDA_VISIBLE_DEVICES"] != "default":
                        device_id = int(config["CUDA_VISIBLE_DEVICES"])
                    else:
                        device_id = 0
                    print("GPU:", torch.cuda.get_device_name(device_id))
            else:
                print("Disable GPU acceleration and use CPU for computation.")
        except:
            print("Error: Pytorch is not correctly installed.")
            passflag = False
        if passflag:
            print("All tests passed.")

    else:
        print("Invalid function mode. Please check the parameters.")


if __name__ == "__main__":
    print("Welcome to use eBiota!")
    print("Usage: python eBiota.py -F [Function mode] -o [Output directory] -i [Input tsv file for modification]")
    print("Function mode: doall, evaluate, design, cooc, genemod, test\n")
    print("Examples:")
    print("python eBiota.py -F doall -o ./results")
    print("python eBiota.py -F evaluate -o ./results")
    print("python eBiota.py -F design -o ./results")
    print("python eBiota.py -F cooc -o ./results")
    print("python eBiota.py -F genemod -i ./results/glc__D_e__to__h2_e__with_O2__with_glucose.tsv")
    print("python eBiota.py -F test\n")

    args = parse_arg()
    update_config(args)
    # if os.path.exists("tmp"):
    #     shutil.rmtree("tmp")
    # os.mkdir("tmp")
    check_config()
    main(args)
    