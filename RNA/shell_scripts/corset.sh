#!/bin/bash
#SBATCH --account=def-kmj477
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --time=12:00:00
#SBATCH --constraint=broadwell
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=0

module load nixpkgs/16.09 intel/2018.3
module load corset/1.07


corset -D 99999999999 -g 1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,5,5,5,5,5,5,6,6,6,6,6,6 \
-n C.0.4B.1,C.0.4B.2,C.0.4B.4,C.0.4B.5,C.HD3.2,C.HD3.3,C.HD3.4,C.HD3.5,C.HD3.6,C.HD3.8,\
C.HD7.2,C.HD7.3,C.HD7.4,C.HD7.7,C.HD7.8,C.HD7.9,S.0.4B.1,S.0.4B.2,S.0.4B.3,S.0.4B.5,S.HD3.1,\
S.HD3.2,S.HD3.3,S.HD3.4,S.HD3.5,S.HD3.7,S.HD7.2,S.HD7.3,S.HD7.4,S.HD7.5,S.HD7.8,S.HD7.9 \
-i salmon_eq_classes /home/biomatt/projects/def-kmj477/biomatt/splittail/salm_quants_corset/salm_output_*/aux_info/eq_classes.txt \
-p /home/biomatt/scratch/splittail/splittail_broadwell