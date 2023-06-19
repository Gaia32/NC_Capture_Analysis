#!/bin/bash
#SBATCH --job-name=paleo

date

##########################################################################
# This is a mitogenome analysis for 2 unidentified sequences             #
# This part is the consensus building for the inference in the tree      #
# We build the consensus from teh dog as well as the wolf to see if      #
# there is any impact using the different reerence genomes               #
##########################################################################

DirProjet=(/scratch/mahautgr/INTERNSHIP_Elis/ISLAMIC_Samples/7_Consensus)
cd $DirProjet
mkdir -p MAP_DAMAGE ANGSD REPORT HAPLOTYPE

DirRaw=(/scratch/mahautgr/INTERNSHIP_Elis/ISLAMIC_Samples/6_Multimapping/SAMTOOLS/)
DirDdup=(/home/genouest/cnrs_umr6553/mahautgr/INTERNSHIP_Elis/ISLAMIC_Samples/5_Dedup_Fastp/)
DirMamD=$DirProjet/MAP_DAMAGE/
DirAng=$DirProjet/ANGSD/
DirHapl=$DirProjet/HAPLOTYPE/
DirREPORT=$DirProjet/REPORT/

REF_DOG=/scratch/mahautgr/REF_GENOMES/MITOCHONDRION/DOG_Ref_mitochondrion/NC_002008.4.fasta
REF_GREY_WOLF=/scratch/mahautgr/REF_GENOMES/MITOCHONDRION/GREY_WOLF_Ref_mitochondrion/NC_008092.1.fasta

for file in ${DirRaw}*_GREY_WOLF_aln_samse_mappedQ25.sorted.bam; do # <------- Change the name here DOG or GREY WOLF
	name=${file##*/}
	sample_tag=${name%_G*} # <----- Change the last letter here

    echo "___________ 7bis - Rescale ___________"
    source /local/env/envconda.sh
    conda activate /groups/Paleogenomics/ENV/MapDamage-2.2.1 # <-------- Change the names here too (-i ; $REF ; -o)

    mapDamage -i $DirRaw${sample_tag}_GREY_WOLF_aln_samse_mappedQ25.sorted.bam -r $REF_GREY_WOLF --rescale > $DirMamD${sample_tag}_GREY_WOLF_mappedQ25_try.sorted.rescaled.bam

    echo "_______________________ 8 - AngSD _________________________" # <------- Change the name here (-i)
    
    source /local/env/envconda.sh
    conda activate /home/genouest/cnrs_umr6553/mahautgr/my_env/Angsd
    angsd -doFasta 2 -out $DirAng${sample_tag}_GREY_WOLF_consensus_true -explode 1 -doCounts 1 -minQ 20 -minMapQ 30 -i $DirMamD${sample_tag}_GREY_WOLF_mappedQ25.sorted.rescaled.bam -setMinDepth 3 -nThreads 20 -doDepth 1
    # angsd -doFasta 2 -out $DirAng${sample_tag}_GREY_WOLF_consensus -explode 1 -doCounts 1 -minQ 20 -minMapQ 30 -i $DirMamD${sample_tag}_GREY_WOLF_aln_samse_mappedQ25.sorted.rescaled.bam -setMinDepth 3 -nThreads 20 -doDepth 1
    # angsd -doFasta 2 -out $DirAng${sample_tag}_DOG_consensus -explode 1 -doCounts 1 -minQ 20 -minMapQ 30 -i $DirMamD${sample_tag}_DOG_mappedQ25.sorted.rescaled.bam -setMinDepth 3 -nThreads 20 -doDepth 1
    
done