#!/bin/bash
#SBATCH --job-name=paleo

date

##########################################################################
# This is a mitogenome analysis for 2 unidentified sequences             #
# This part is the consensus building for the inference in the tree 	 #
##########################################################################

DirProjet=(/scratch/mahautgr/INTERNSHIP_Elis/ISLAMIC_Samples/8_Alignment)
cd $DirProjet
mkdir -p 

	echo " _______________________ 10 - Mafft _________________________"
    source /local/env/envmafft-7.475.sh
#### 1057
	mafft --add P27761_1057_consensus.fa --reorder ISLAMIC8SAMPLES_TREE_BASE_ALIGNED_V2_DLOOP_WORK.fst > Alignementginsi_1057.fst
	ginsi Alignementginsi_1057.fst > New_alignement_ginsi_1057.fst

#### 1058
	mafft --add P27761_1058_consensus.fa --reorder New_alignement_ginsi_1057.fst > Alignementginsi_1057_1058.fst
	ginsi Alignementginsi_1057_1058.fst > New_alignement_ginsi_1057_1058.fst

####################################################
# ETAPE 10bis - Cleaning de l'alignement
## Outil - Hmmcleaner
	echo " _________________________ 10bis - Hmmcleaner / trimal / prequal_______________________"

    source /local/env/envconda.sh
    conda activate /groups/Paleogenomics/ENV/trimal 
    trimal -in New_alignement_ginsi_1057_1058.fst -out New_alignement_ginsi_curred_1057_1058.fst -gt 0.9 # <-- on prend celui là

# Options choisies: -resoverlap 0.60 -seqoverlap 60
# * - gt 0.9 --> si il y a des gap dans 90% des positions, on supprime la colonne
# * -resoverlap 0.60 --> trimAI calcul un score à chaque position par rapport à l'élément le plus commun à hauteur de 60%
# * -seqoverlap 60 --> Si la séquence a 60% de ses positions qui ont un score en dessous du resoverlap, la séquence est supprimée.
# --> il faut que 60% de ses positiosn passent le score pour que la seq soit gardée


####################################################
# ETAPE 11 - Génération des arbres
## Outil - Iqtree
	# echo " _________________________ 11 - Iqtree _______________________"

    # source /local/env/envconda.sh
    # conda activate /groups/Paleogenomics/ENV/iqtree

    # -m GTR2+G4+F -alrt 1000

# Options choisies:
# * -s --> pour spécifier le nom du fichier d'alignement en input
# * -m MFP --> model finder as a substitution model
# * -B 1000 --> nbr of boostrap replicates for UFboot
# * -b 100 --> nbr of nonparametric boostrap 
# * -m TIM2+I+G --> run UFBoot ultrafastbootstrap
# * -T 8 specifies the nbr og GPU 

