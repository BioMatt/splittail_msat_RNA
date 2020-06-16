#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --account=def-kmj477
#SBATCH --mem=64000M
#SBATCH --mail-user=thorstem@myumanitoba.ca
#SBATCH --mail-type=ALL

module load nixpkgs/16.09 intel/2018.3
module load star/2.7.0a

starindex="/home/biomatt/scratch/splittail_star_index"

STAR --runMode genomeGenerate --runThreadN 4 --genomeDir $starindex \
--genomeFastaFiles /home/biomatt/projects/def-kmj477/biomatt/splittail/lace_transcriptome/splittail_lace.fasta \
--sjdbGTFfile /home/biomatt/projects/def-kmj477/biomatt/splittail/lace_transcriptome/splittail_lace_trans.gff \
--sjdbOverhang 99 \
--limitGenomeGenerateRAM 64000000000 \
--limitIObufferSize 16000000000 \
--genomeChrBinNbits 12