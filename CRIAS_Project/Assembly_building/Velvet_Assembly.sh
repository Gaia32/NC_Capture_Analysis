#!/bin/bash
#SBATCH --job-name=paleo

####################################################################
# This is a mitogenome analysis for 2 unidentified sequences       #
# This part the e novo assembly method                      	   #
####################################################################

####################################################################
# Velvet requires an index file to be built before the assembly takes place. 
# We must choose a k- mer value for building the index. Longer k- mers result in a more stringent assembly, 
# at the expense of coverage. There is no definitive value of k for any given project. 
# However, there are several absolute rules:
#       k must be less than the read length
#       it should be an odd number.
# Here the k-mer size is 25 because we discarded reads shorter than this

# Here is the documentation of Velvet: https://angus.readthedocs.io/en/2016/week3/LN_assembly.html

date

	echo -e "\n_________________________ Sample P27761_1057 _________________________"
    source /local/env/envconda.sh
    conda activate /home/genouest/cnrs_umr6553/mahautgr/my_env/VelvetOptimiser

    VelvetOptimiser.pl --f '-shortPaired -fastq.gz P27761_1057_rmduplicates.fastq.gz' \
    -s 25 --p "P27761_1057" --t 8

    echo -e "\n_________________________ Sample P27761_1058 _________________________"
    VelvetOptimiser.pl --f '-shortPaired -fastq.gz P27761_1058_rmduplicates.fastq.gz' \
    -s 25 --p "P27761_1058" --t 8

date