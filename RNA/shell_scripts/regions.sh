#!/bin/bash
#SBATCH --account=def-kmj477
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --mem=10000M
#SBATCH --time=01:00:00

module load nixpkgs/16.09 intel/2018.3
module load bamtools/2.4.1
module load gcc/5.4.0 freebayes/1.2.0

bamtools coverage -in /home/biomatt/scratch/splittail_cigar/merged_cigar_long.bam \
| /home/biomatt/walleye_rna/coverage_to_regions.py /home/biomatt/projects/def-kmj477/biomatt/splittail/lace_transcriptome/splittail_lace.fasta.fai 10000 >/home/biomatt/scratch/splittail_ref.fa.10000.regions