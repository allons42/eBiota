Growth medium
=====

## Basic LB medium

We started our medium design with Luria-Bertani (LB) medium, one of the most basic media. To allow most bacteria to grow, we figured out 215 essential metabolites that are necessary for the growth of some bacteria, mainly including ions, amino acids, nucleotides, vitamins and various sugars. Afterwards, we removed those labeled as carbon sources according to the Biolog phenotype array to reduce selective utilization. Finally, we added 94 essential substrate to the medium, with limited uptake rate at 0.1 mmol/(gDW·h). Additionally, glucose was added when the substrate was not a carbon source. The configuration of this medium can be found at github.

## Custom medium

eBiota accepts custom medium in CSV format, as demonstrated in the Usage section. The input file for the medium must include two essential columns: **id** and **maxFlux**.

### id

The **id** must be consistent with the metabolite id in the GEM. It typically ends with "_e", representing an exchange metabolite. Notice that there are two commonly used naming system: [BiGG ID](http://bigg.ucsd.edu/universal/metabolites) and [ModelSEED ID](https://modelseed.org/biochem/compounds). The default setting of eBiota is BiGG ID, which is also used in Carveme GEM and BiGG GEM. To use ModelSEED ID, we provide a “SEED Compound” column in Basic LB medium, which can also be searched in the aforementioned website.

### maxFlux

**maxFlux** indicates the maximum uptake rate of the metabolite for the community, with a default maximum value of 1000 mmol/(gDW·h). Typically the uptake rate of main carbon source is set at 10, and that of neccesary amino acids are set at 0.1 or 1.
