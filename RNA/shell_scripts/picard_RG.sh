#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --account=def-kmj477
#SBATCH --mem=80000M
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --array=1-32

module load nixpkgs/16.09  intel/2018.3
module load samtools/1.9
module load picard/2.18.9

indir="/home/biomatt/scratch/splittail_cleaned_bam/"

list=splittail_samples.txt
string="sed -n "$SLURM_ARRAY_TASK_ID"p ${list}"
splittail_id=$($string)

echo "$splittail_id"

outdir="/home/biomatt/scratch/splittail_rg/"

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=${indir}${splittail_id}_cleaned.bam O=${outdir}${splittail_id}_rg_sorted.bam \
SO=coordinate \
RGID=${splittail_id} \
RGLB=${splittail_id} \
RGPL=ILLUMINA \
RGPU=${splittail_id} \
RGSM=${splittail_id}
