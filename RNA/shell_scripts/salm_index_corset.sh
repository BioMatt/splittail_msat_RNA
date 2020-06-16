#!/bin/bash
#SBATCH --time=06:00:00
#SBATCH --account=def-kmj477
#SBATCH --mem=24000M
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL


module load nixpkgs/16.09 gcc/7.3.0 openmpi/3.1.2
module load salmon/0.11.3

salmon index -t /home/biomatt/projects/def-kmj477/biomatt/splittail/transcriptome/bwa_ref.fasta -i salm_index_corset quasi -k 31
