# eBiota: a platform for *ab initio* designing microbial communities with desired functions


[![Documentation Status](https://readthedocs.org/projects/e-biota/badge/?version=latest)](https://e-biota.readthedocs.io/en/latest/index.html)
DOI

We developed eBiota, a platform for the *ab initio* design of artificial microbial communities with desired functions from a large bacterial seed pool, including maximal production or degradation efficiency of specified compounds. To achieve this, three novel algorithms are embedded in eBiota: CoreBFS, ProdFBA, and DeepCooc. CoreBFS is a graph-based, breadth-first search (BFS) algorithm for rapidly determining whether a bacterium has metabolic pathways from substrates to intermediates and/or from intermediates to products. Noting that *ab initio* designing artificial microbial communities from such a large seed pool is computationally expensive, CoreBFS allows efficient searches at acceptable computational cost through pre-stored molecules and metabolic pathways. ProdFBA is a multi-step FBA algorithm that infers metabolic flux distributions for microbial communities with potential optimum production (or degradation) efficiency while ensuring biomass. Compared with classical FBA, ProdFBA provides more artificial microbial communities with the potential to generate target products. DeepCooc is a deep learning algorithm for predicting the co-occurrence of microbial communities, which provides insight into microbial co-existence.

![workflow](img/Fig1.png)

## Installation

We recommend installing eBiota in an virtual environment with conda.

```bash
# Download latest version of eBiota
git clone https://github.com/allons42/e-Biota.git
cd e-Biota/src/

# Create and activate a conda environment "ebiota_env"
conda env create -n ebiota_env --file ebiota_env.yml
conda activate ebiota_env

# install Carveme for GEM rebuild, according to https://carveme.readthedocs.io/
pip install carveme
conda install -c bioconda diamond

# check if the installation is successful
python eBiota.py --test
```

## Quickstart

eBiota supports various functions for microbial community design. For a quickstart, we provide an example to design communities that utilize glucose and produce hydrogen.

```bash
python eBiota.py --Function design --substrate glc__D_e --product h2_e
```

There are plenty of configurations to custom your communities, detailed in `config.json`. For more usage and tutorials, see the [documentation](https://e-biota.readthedocs.io/en/latest/index.html).

## Main Results

The results mentioned in our paper is quite large, and can be downloaded from zenodo: link.

The following results are included:

1. eBiota-GEM: 21,514 Genome-Scale Metabolic Models (GEMs) constructed using CarveMe based on RefSeq complete genomes.
2. baterial evaluation: The evaluation of the ability to uptake substrates and secret productions for all 21514 GEMs.
3. community design result: The result calculated from eBiota-GEM includes various combinations for two-bacterial consortia, covering strain IDs, substrates, products, yields, dual-bacterial growth, single-bacterial growth, co-occurrence predictions, and interactions. 
4. cultivation information: The mapping of eBiota-GEM to BacDive database, including hierarchical category,  culture medium, temperature, and oxygen tolerance. The data is accessed at Oct 15 2022.
5. supplementary tables in papers.

## Citation

eBiota: a platform for *ab initio* designing microbial communities with desired functions
