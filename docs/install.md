# Installation

## **Conda**

To ensure project-specific dependencies and avoid modifying system environment, we strongly recommend using conda to create virtual environments. All dependencies can be installed from conda channels without super user rights using the following steps:

**1. Install Mini-/Anaconda**

Follow the instructions provided by conda to Install [Anaconda/Miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

**2. Create conda environment for eBiota and adding package sources**

```bash
# Download latest version of eBiota
git clone https://github.com/allons42/eBiota.git
cd eBiota/src/

# Create and activate a conda environment "ebiota_env"
conda create -n ebiota_env python=3.9
conda activate ebiota_env
pip install -r requirements.txt

# [Optional] install CarveMe for GEM rebuild, according to https://carveme.readthedocs.io/
pip install carveme
conda install -c bioconda diamond
```

**3. Install DeepCooc**

The DeepCooc module, used for co-occurrence analysis, utilizes deep learning and requires the installation of the PyTorch package. Please refer to the [PyTorch](https://pytorch.org/get-started/locally/) official website to choose the appropriate CUDA version or CPU version for installation. We used the 1.12.0+cu113 version, and other versions should be feasible.

```bash
# Install PyTorch with conda or wheel on Linux
# Conda
conda install pytorch==1.12.0 torchvision==0.13.0 torchaudio==0.12.0 cudatoolkit=11.3 -c pytorch

# Wheel
pip install torch==1.12.0+cu113 torchvision==0.13.0+cu113 torchaudio==0.12.0 --extra-index-url https://download.pytorch.org/whl/cu113
```

After installing PyTorch, you need to download necessary files from Zenodo ([doi: 10.5281/zenodo.13895108](https://doi.org/10.5281/zenodo.13895108)):

```bash
wget -O DeepCooc_files.tar.gz https://zenodo.org/records/13895108/files/DeepCooc_files.tar.gz?download=1
tar -xzvf DeepCooc_files.tar.gz -C stats/
rm DeepCooc_files.tar.gz
```

**4. Test the installation**

```bash
python eBiota.py --test
```

## Perl

Our algorithm requires Perl, which is usually already installed in computer system.

```bash
# check the version and availability of Perl
perl -v
```

If Perl is not installed, follow the instructions on the [official website](https://www.perl.org/get.html) to install.

## **[Optional] Install Gurobi Optimizer**

We recommend using *Gurobi* as Linear Programming solver, as it is usually faster than the default optimizer *glpk*. The *Gurobi* solver is free for academic use ([see here](https://www.gurobi.com/features/academic-named-user-license/)). Please follow the instructions to install *Gurobi*.

## [Optional] Collect eBiota database

The complete dataset can be downloaded from Zenodo ([doi: 10.5281/zenodo.13895108](https://doi.org/10.5281/zenodo.13895108)), including the following results:

1. **GEM.tar.gz**: The eBiota-GEM dataset, containing 21,514 Genome-Scale Metabolic Models (GEMs) constructed using CarveMe based on RefSeq complete genomes.
2. **Baterial_evaluation.tar.gz**: The evaluation of the ability to uptake substrates and secret productions for all 21,514 GEMs.
3. **Community_results.tar.gz**: The results calculated from eBiota-GEM includes various combinations for two-bacterial consortia, covering strain IDs, substrates, products, yields, dual-bacterial growth, single-bacterial growth, co-occurrence predictions, and interactions.
4. **DeepCooc_files.tar.gz**: The parameters of DeepCooc, required by eBiota platform.
