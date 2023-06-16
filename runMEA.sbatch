#!/bin/bash

#SBATCH -J MEA_nf
#SBATCH --mem=4G
#SBATCH -o ./MEA_%J.out
#SBATCH -D /scratch/mblab/edwardkang/llfs_module_enrichment_nf/

# set up singularity
eval $(spack load --sh singularityce@3.8.0)
export SINGULARITY_CACHEDIR="/scratch/mblab/edwardkang/singularity/cache"

# set up nextflow environment
eval $(spack load --sh nextflow@22.10.4)

# run nextflow
nextflow run mea.nf -c conf/mea.config