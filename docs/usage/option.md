Usage Options
================

eBiota supports various functions for microbial community design. For a quickstart, we provide an example to design communities that utilize glucose and produce hydrogen. The whole process would take about 10 minutes.

```bash
python eBiota.py --Function doall --outdir result
```

There are plenty of configurations to customize your communities, detailed in `config.json`. Here list some important configurations:

| Configuration  | meaning                                                      |
| -------------- | ------------------------------------------------------------ |
| suffix         | suffix of GEM files, usually ".xml" or ".xml.gz"             |
| path_GEM       | the path to GEM files                                        |
| path_output    | the path to store final results                              |
| medium         | the path to medium file                                      |
| max_proc       | maximum process number in parallel                           |
| prune          | whether to prune similar GEMs in the results                 |
| target         | "production" or "degradation"                                |
| **substrate**  | the chosen substrate, set to "default" to enumerate all possible substrates (only in production mode) |
| intermediate   | the chosen intermediate, set to "default" to enumerate all possible intermediates |
| **product**    | the chosen product, set to "default" to enumerate all possible products (only in degradation mode) |
| oxygen         | whether the medium contains oxygen, set to "default" to consider both conditions |
| glucose        | whether the medium contains glucose, set to "default" to consider both conditions |
| community_size | the size of microbial communities                            |


Building GEM model from genomes
------------------

We use the open-source tool [CarveMe](https://carveme.readthedocs.io/) to build Genome-scale Metabolic Models (GEMs), which can be reconstructed with existing faa file or RefSeq GCF id. Alternatively, you can directly use existing GEMs from databases such as BiGG.

```bash
# use Carveme to build GEM (.xml or .xml.gz) with faa file
carve genome.faa -o model.xml.gz
# use Carveme to build GEM with RefSeq GCF id
carve --refseq GCF_000005845.2 -o ecoli_k12_mg1655.xml
```

Bacteria evaluation
-------------------

```bash
python eBiota.py --Function evaluate
```
This function will evaluate the GEMs in the `GEM` directory and store the results in the `tmp` directory. The evaluation results include the growth rate and production rate of each bacterium in the medium determined by `config.json`.


Community design
----------------

```bash
python eBiota.py --Function design --outdir result
```
This is the main function to design microbial communities. The program will output the designed communities in the designated directory. All reletive parameters are determined by `config.json`, including the medium, the number of bacteria in the community, the optimization goal, etc.
Notice that the bacteria used for community design should be processed by “Function: evaluate” before. 

Co-occurrence prediction with DeepCooc
--------------------------------------
   
   ```bash 
   python eBiota.py --Function cooc --outdir result
   ```
This function will predict the co-occurrence relationship between bacteria in communities. This function requires the output of "Function: design" as input.