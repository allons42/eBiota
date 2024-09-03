# Installation

## **Conda**

To ensure project-specific dependencies and avoid modifying system environment, we strongly recommend using conda to create virtual environments. All dependencies can be installed from conda channels without super user rights using the following steps:

**1. Install Mini-/Anaconda**

Follow the instructions provided by conda to Install [Anaconda/Miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

**2. Create conda environment for eBiota and adding package sources**

```bash
# Download latest version of eBiota from github
git clone https://github.com/allons42/e-Biota.git
cd e-Biota/src/

# Create and activate a conda environment "ebiota_env"
conda env create -n ebiota_env --file ebiota_env.yml
conda activate ebiota_env

# install Carveme for GEM rebuild, according to https://carveme.readthedocs.io/
pip install carveme
conda install -c bioconda diamond
```

**3. Test the installation**

```bash
python eBiota.py --test
```

## **Pytorch**

The DeepCooc module, used for co-occurrence analysis, utilizes deep learning and requires the installation of the PyTorch package. Please refer to the [PyTorch](https://pytorch.org/get-started/locally/) official website to choose the appropriate CUDA version or CPU version for installation. We are using the 1.10.1+cu111 version, and other versions should be feasible. If there are any problems, feel free to raise an issue on our [GitHub](https://github.com/allons42/e-Biota).

```bash
# the version we used
pip install torch==1.10.1+cu111 torchvision==0.11.2+cu111 torchaudio==0.10.1 -f https://download.pytorch.org/whl/cu111/torch_stable.html
```

## **Optimizer for Linear Programming**

We recommend using *Gurobi* as LP-solver, as it is usually faster than the default optimizer *glpk*. The *Gurobi* solver is free for academic use ([see here](https://www.gurobi.com/features/academic-named-user-license/)). Please follow the instructions to install *Gurobi*.

## Perl

Our algorithm requires Perl, which is usually already installed in computer system.

```bash
# check the version and availability of Perl
perl -v
```

If Perl is not availble, follow the instructions on the [official website](https://www.perl.org/get.html) to install.

## Collect e-Biota database

The complete dataset can be downloaded from Zenodo：link, including the following results:

1. eBiota-GEM: 21,514 Genome-Scale Metabolic Models (GEMs) constructed using CarveMe based on RefSeq complete genomes.
2. baterial evaluation: The evaluation of the ability to uptake substrates and secret productions for all 21514 GEMs.
3. community design: The result calculated from eBiota-GEM includes various combinations for two-bacterial consortia, covering strain IDs, substrates, products, yields, dual-bacterial growth, single-bacterial growth, co-occurrence predictions, and interactions. 

## **Troubleshooting**
