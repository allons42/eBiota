Usage Options
================

eBiota supports various functions for microbial community design. For a quickstart, we provide an example to design communities that utilize glucose and produce hydrogen. The whole process would take about 10 minutes.

```bash
python eBiota.py --Function doall --outdir result
```

All community customization options are available in `config.json`. Below is a detailed explanation of each parameter: 

| Configuration        | meaning                                                      | default                  |
| -------------------- | ------------------------------------------------------------ | ------------------------ |
| path_GEM             | the path of the GEM seed pool                                | GEM                      |
| path_output          | the path to store final results                              | result                   |
| medium               | the medium file for microbial communities                    | stats/BasicLB_medium.csv |
| max_proc             | maximum process number in parallel                           | 10                       |
| USE_CUDA             | whether to enable GPU acceleration (via CUDA) for DeepCooc computations | true                     |
| CUDA_VISIBLE_DEVICES | Specifies the GPU device ID(s) to use. This setting has no effect if `USE_CUDA=false` | 0                        |
| prune                | whether to prune similar GEMs in the results                 | true                     |
| target               | "production" or "degradation"                                | production               |
| substrate            | the chosen substrate, set to "default" to enumerate all possible substrates (only in production mode) | glc__D_e                 |
| intermediate         | the chosen intermediate, set to "default" to enumerate all possible intermediates | default                  |
| product              | the chosen product, set to "default" to enumerate all possible products (only in degradation mode) | h2_e                     |
| designated_bacteria  | selected bacteria in the community                           | default                  |
| oxygen               | whether the medium contains oxygen, set to "true", "false", or "default" to consider both conditions | default                  |
| glucose              | whether the medium contains glucose, set to "true", "false", or "default" to consider both conditions | default                  |
| community_size       | the size of microbial communities, supporting 2 or 3         | 2                        |


Building GEMs from genomes
-------------------

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

It is also feasible to specify or replace specific microbiota, by simply setting the `designated_bacteria` parameter in `config.json` to the desired bacterial IDs.

Notice that the bacteria used for community design should be processed by “Function: evaluate” before. 

Co-occurrence prediction with DeepCooc
--------------------------------------

   ```bash 
   python eBiota.py --Function cooc --outdir result
   ```
This function will predict the co-occurrence relationship between bacteria in communities. This function requires the output of "Function: design" as input.

Gene modification
--------------------------------------

   ```bash 
python eBiota.py --Function genemod --input_tsv glc__D_e__to__h2_e__with_O2__with_glucose.tsv
   ```

This function will enumerate possible gene knockouts or gene insertions for better production. The input TSV file is anticipated to be the output of previous functions, including the bacteria IDs and the culture conditions. This function may be significantly slower than the main workflow, so we recommend calculating only the top communities.

