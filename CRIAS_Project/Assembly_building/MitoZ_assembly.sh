#!/bin/bash
#SBATCH --job-name=paleo

mkdir -p /scratch/mahautgr/INTERNSHIP_Elis/ISLAMIC_Samples/6_Assembly_MitoZ

DirMitoZ=/scratch/mahautgr/INTERNSHIP_Elis/ISLAMIC_Samples/6_Assembly_MitoZ/
DirDdup=/home/genouest/cnrs_umr6553/mahautgr/INTERNSHIP_Elis/ISLAMIC_Samples/5_Dedup_Fastp/

    source /local/env/envconda.sh
    conda activate /home/genouest/cnrs_umr6553/mahautgr/my_env/MitoZ

    echo -e "\n_________________________ Sample P27761_1057 _________________________"
#### OPTIONS
# --fq1 <file>   fastq 1 file
# --outprefix <STR>     output prefix
#   --fastq_read_length <INT> read length of fastq reads, used by mitoAssemble. [150]
# --genetic_code <INT>  which genetic code table to use? 'auto' means determined by '--clade' option.
# --clade {Chordata,Arthropoda,Echinodermata,Annelida-segmented-worms,Bryozoa,Mollusca,Nematoda,Nemertea-ribbon-worms,Porifera-sponges}
#           which clade does your species belong to? [Arthropoda]
# --memory <INT>        memory size limit for spades/megahit, no enough memory will make the two programs halt
#                       or exit [50]
# --kmes_megahit The default value is --kmers_megahit 21 29 39 59 79 99 119 141, 
#               which can cause wrong mitogenomes sometimes, therefore, I would recommend trying larger kmers first, and if they do not work, you can try smaller kmers, say --kmers_megahit 39 59 79 99 119 141
# Finally, set a required taxonomic name (e.g. order, class, family, genus) for filtering out non-target-taxon mitochondrial sequences: --requiring_taxa Chordata. 
### try adding those if doesn't work 21 29
    mitoz assemble --outprefix "P27761_1057" --thread_number 16 --clade Chordata --genetic_code 'auto' --fq1 $DirDdup"P27761_1057_rmduplicates.fastq.gz" --fastq_read_length 25 --assembler megahit --kmers_megahit 21 29 39 59 79 99 119 141 --memory 80 --requiring_taxa Chordata
    
    echo -e "\n_________________________ Sample P27761_1058 _________________________"
    
    mitoz assemble --outprefix "P27761_1058" --thread_number 16 --clade Chordata --genetic_code 'auto' --fq1 $DirDdup"P27761_1058_rmduplicates.fastq.gz" --fastq_read_length 25 --assembler megahit --kmers_megahit 21 29 39 59 79 99 119 141 --memory 80 --requiring_taxa Chordata