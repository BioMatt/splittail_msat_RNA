#!/bin/bash
#SBATCH --account=def-kmj477
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --mem=20000M
#SBATCH --array=1-32
#SBATCH --time=08:00:00


module load nixpkgs/16.09 gcc/5.4.0 openmpi/2.1.1
module load salmon/0.12.0

list=splittail_read_list.txt
string="sed -n "$SLURM_ARRAY_TASK_ID"p ${list}"
str=$($string)

var=$(echo $str | awk -F"\t" '{print $1, $2}')
set -- $var
c1=$1
c2=$2

echo "$c1"
echo "$c2"

sample=$(echo ${c1} | cut -d"." -f1)

salmon quant --dumpEq --validateMappings --rangeFactorizationBins 4 --seqBias --gcBias \
-i /home/biomatt/projects/def-kmj477/biomatt/splittail/salm_index_corset -l IU \
-1 /home/biomatt/projects/def-kmj477/biomatt/splittail/clean_reads/${c1} \
-2 /home/biomatt/projects/def-kmj477/biomatt/splittail/clean_reads/${c2} \
-o /home/biomatt/projects/def-kmj477/biomatt/splittail/salm_quants_corset/salm_output_${sample}