#!/bin/bash
#SBATCH --time=06:00:00
#SBATCH --account=def-kmj477
#SBATCH --mem=120000M
#SBATCH --array=1-32
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL

module load nixpkgs/16.09  intel/2018.3 
module load star/2.7.0a

starindex="/home/biomatt/scratch/splittail_star_index"

indir="/home/biomatt/projects/def-kmj477/biomatt/splittail/clean_reads/"

workingdir="/home/biomatt/scratch/splittail_star_2passmode/"

list=splittail_read_list.txt
string="sed -n "$SLURM_ARRAY_TASK_ID"p ${list}" 
str=$($string) 

var=$(echo $str | awk -F"\t" '{print $1, $2}') 
set -- $var 
c1=$1 
c2=$2

echo "$c1" 
echo "$c2"

splittail=$(echo ${c1} | cut -d _ -f 1)

echo "$splittail"

mkdir ${workingdir}${splittail}

outdir=${workingdir}${splittail}

STAR --runMode alignReads --genomeDir $starindex --readFilesCommand zcat \
--twopassMode Basic \
--readFilesIn ${indir}${c1} ${indir}${c2} \
--outFileNamePrefix ${outdir}/${splittail} \
--outSAMtype BAM Unsorted \
--runThreadN 48
