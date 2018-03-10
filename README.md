# tRNAMultiEnv

Data analysis for fitness landscape across multiple environments. Please refer to our publication in Nature Ecology and Evolution for details.

## Overview
A fitness landscape (FL) describes the genotype-fitness relationship in a given environment.  To explain and predict evolution, it is imperative to measure the FL in multiple environments because the natural environment changes frequently.  Using a high-throughput method that combines precise gene replacement with next-generation sequencing, we determine the in vivo FL of a yeast tRNA gene comprising over 23,000 genotypes in four environments.  (More details will be updated upon publication of the manuscript)

## Contact
If you have any questions or comments, you could reach me at lichuan [at] umich [dot] edu.

## FILE MANIFEST
### Reads processing
Perl scripts for barcode counts from fastq files. Only perfectly matched fully overlapping paired-end reads were used in estimating genotype frequencies. 
### Calculating fitness
R scripts for calculating fitness. The change in relative genotype frequency during the competition was used to estimate the fitness of each genotype relative to the wild-type.  The fitness of a genotype is calculated by averaging the fitness among biological replicates.  To ensure relatively accurate fitness estimation, we focused on 23,284 genotypes with read counts ≥ 100 before the competition unless otherwise mentioned.  
### Other analysis
Analysis for epistasis, GxE and GxGxE. Epistasis is defined by ε=f_AB-f_A f_B, where fA and fB are the fitness values of two single mutants and fAB is the fitness of the corresponding double mutant.  We similarly tested G×G×E by examining if epistasis is significantly different between two environments, using a t-test at a nominal P-value of 0.05.  We also tested G×G×E using an alternative definition of epistasis.

Other codes in the paper are available upon requests.
