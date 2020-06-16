#!/bin/bash
#SBATCH --time=06:00:00
#SBATCH --account=def-kmj477
#SBATCH --mem=0
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --constraint=skylake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48

module load nixpkgs/16.09 intel/2018.3
module load star/2.7.0a

starindex="/home/biomatt/scratch/splittail_star_index2"

cd /home/biomatt/scratch/

STAR --runMode genomeGenerate --runThreadN 2 --genomeDir $starindex \
--genomeFastaFiles /home/biomatt/projects/def-kmj477/biomatt/splittail/lace_transcriptome/splittail_lace.fasta \
--sjdbGTFfile /home/biomatt/projects/def-kmj477/biomatt/splittail/lace_transcriptome/splittail_lace_trans.gff \
--outTmpDir /home/biomatt/scratch/star_2passtemp \
--sjdbFileChrStartEnd /home/biomatt/scratch/splittail_star_1pass/C-0-4B-1/C-0-4B-1SJ.out.tab /home/biomatt/scratch/splittail_star_1pass/C-0-4B-2/C-0-4B-2SJ.out.tab \
/home/biomatt/scratch/splittail_star_1pass/C-0-4B-4/C-0-4B-4SJ.out.tab /home/biomatt/scratch/splittail_star_1pass/C-0-4B-5/C-0-4B-5SJ.out.tab \
/home/biomatt/scratch/splittail_star_1pass/C-HD3-2/C-HD3-2SJ.out.tab /home/biomatt/scratch/splittail_star_1pass/C-HD3-3/C-HD3-3SJ.out.tab /home/biomatt/scratch/splittail_star_1pass/C-HD3-4/C-HD3-4SJ.out.tab \
/home/biomatt/scratch/splittail_star_1pass/C-HD3-5/C-HD3-5SJ.out.tab /home/biomatt/scratch/splittail_star_1pass/C-HD3-6/C-HD3-6SJ.out.tab /home/biomatt/scratch/splittail_star_1pass/C-HD3-8/C-HD3-8SJ.out.tab \
/home/biomatt/scratch/splittail_star_1pass/C-HD7-2/C-HD7-2SJ.out.tab /home/biomatt/scratch/splittail_star_1pass/C-HD7-3/C-HD7-3SJ.out.tab /home/biomatt/scratch/splittail_star_1pass/C-HD7-4/C-HD7-4SJ.out.tab \
/home/biomatt/scratch/splittail_star_1pass/C-HD7-7/C-HD7-7SJ.out.tab /home/biomatt/scratch/splittail_star_1pass/C-HD7-8/C-HD7-8SJ.out.tab /home/biomatt/scratch/splittail_star_1pass/C-HD7-9/C-HD7-9SJ.out.tab \
/home/biomatt/scratch/splittail_star_1pass/S-0-4B-1/S-0-4B-1SJ.out.tab /home/biomatt/scratch/splittail_star_1pass/S-0-4B-2/S-0-4B-2SJ.out.tab /home/biomatt/scratch/splittail_star_1pass/S-0-4B-3/S-0-4B-3SJ.out.tab \
/home/biomatt/scratch/splittail_star_1pass/S-0-4B-5/S-0-4B-5SJ.out.tab /home/biomatt/scratch/splittail_star_1pass/S-HD3-1/S-HD3-1SJ.out.tab /home/biomatt/scratch/splittail_star_1pass/S-HD3-2/S-HD3-2SJ.out.tab \
/home/biomatt/scratch/splittail_star_1pass/S-HD3-3/S-HD3-3SJ.out.tab /home/biomatt/scratch/splittail_star_1pass/S-HD3-4/S-HD3-4SJ.out.tab /home/biomatt/scratch/splittail_star_1pass/S-HD3-5/S-HD3-5SJ.out.tab \
/home/biomatt/scratch/splittail_star_1pass/S-HD3-7/S-HD3-7SJ.out.tab /home/biomatt/scratch/splittail_star_1pass/S-HD7-2/S-HD7-2SJ.out.tab /home/biomatt/scratch/splittail_star_1pass/S-HD7-3/S-HD7-3SJ.out.tab \
/home/biomatt/scratch/splittail_star_1pass/S-HD7-4/S-HD7-4SJ.out.tab /home/biomatt/scratch/splittail_star_1pass/S-HD7-5/S-HD7-5SJ.out.tab /home/biomatt/scratch/splittail_star_1pass/S-HD7-8/S-HD7-8SJ.out.tab \
/home/biomatt/scratch/splittail_star_1pass/S-HD7-9/S-HD7-9SJ.out.tab \
--sjdbOverhang 99 \
--genomeSAindexNbases 1 \
--limitGenomeGenerateRAM 180046865450
