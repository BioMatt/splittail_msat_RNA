#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --account=def-kmj477
#SBATCH --mem=80000M
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --array=1-32

module load nixpkgs/16.09 intel/2018.3
module load gatk/4.0.8.1

echo "yes"

indir="/home/biomatt/scratch/splittail_markdup/"

list=splittail_samples.txt
string="sed -n "$SLURM_ARRAY_TASK_ID"p ${list}"
splittail_id=$($string)

echo "$splittail_id"

outdir="/home/biomatt/scratch/splittail_cigar/"

java -jar /cvmfs/soft.computecanada.ca/easybuild/software/2017/Core/gatk/4.0.8.1/gatk-package-4.0.8.1-local.jar \
SplitNCigarReads --reference /home/biomatt/projects/def-kmj477/biomatt/splittail/lace_transcriptome/splittail_lace.fasta \
--input ${indir}${splittail_id}_pic_markdup.bam --output ${outdir}${splittail_id}_cigar.bam
