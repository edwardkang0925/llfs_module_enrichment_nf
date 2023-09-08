#!/bin/bash

# $1 is the trait name (ex)ABI which should also be the name of gene score file, for example  ABI.csv)
# $2 is the total number of tests (sum of number of genes across categories)
# $3 is the name of the pipeline (cma, twas, gwas, staar)
#SBATCH -J MEA_nf
#SBATCH --mem=4G
#SBATCH -o logs/MEA_%J.out
#SBATCH -D /scratch/mblab/edwardkang/llfs_module_enrichment_nf/

# set up singularity
eval $(spack load --sh singularityce@3.8.0)
export SINGULARITY_CACHEDIR="/scratch/mblab/edwardkang/singularity/cache"

# set up nextflow environment
eval $(spack load --sh nextflow@22.10.4)

# run nextflow
nextflow run mea_slurm.nf --trait $1 --numTests $2 --pipeline $3 --pvalFileName "/app/data/pvals/$1/$1.csv" -c conf/mea.config