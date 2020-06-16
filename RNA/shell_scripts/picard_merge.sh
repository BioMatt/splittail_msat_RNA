#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --account=def-kmj477
#SBATCH --mem=80000M
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL

module load nixpkgs/16.09  intel/2018.3
module load samtools/1.9
module load picard/2.18.9

java -jar $EBROOTPICARD/picard.jar MergeSamFiles \
I=/home/biomatt/scratch/splittail_cigar/C-0-4B-1_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/C-0-4B-2_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/C-0-4B-4_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/C-0-4B-5_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/C-HD3-2_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/C-HD3-3_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/C-HD3-4_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/C-HD3-5_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/C-HD3-6_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/C-HD3-8_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/C-HD7-2_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/C-HD7-3_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/C-HD7-4_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/C-HD7-7_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/C-HD7-8_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/C-HD7-9_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/S-0-4B-1_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/S-0-4B-2_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/S-0-4B-3_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/S-0-4B-5_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/S-HD3-1_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/S-HD3-2_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/S-HD3-3_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/S-HD3-4_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/S-HD3-5_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/S-HD3-7_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/S-HD7-2_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/S-HD7-3_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/S-HD7-4_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/S-HD7-5_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/S-HD7-8_cigar.bam \
I=/home/biomatt/scratch/splittail_cigar/S-HD7-9_cigar.bam \
O=/home/biomatt/scratch/splittail_cigar/merged_cigar_long.bam