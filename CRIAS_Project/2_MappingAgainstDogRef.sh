#!/bin/bash
#SBATCH --job-name=paleo

############################################################################################
# This is a mitogenome analysis for 2 unidentified sequences                               #
# This part is the mapping against the dog reference genome                                #
# We mapped against the whole reference to recover more coverage than just the mitogenome  #
# This script was run in the Genouest cluster usinga slurm job                             #
############################################################################################

DirProjet=(/scratch/mahautgr/INTERNSHIP_Elis/ISLAMIC_Samples/6_Mapping_DOG) # <------------- + faire attention à ne pas mettre de / à la fin
cd $DirProjet
mkdir -p BWA SAMTOOLS MAP_DAMAGE FASTA ANGSD REPORT

DirRaw=(/home/genouest/cnrs_umr6553/mahautgr/INTERNSHIP_Elis/ISLAMIC_Samples/RAW_Data/Pig_data/)
DirDdup=(/home/genouest/cnrs_umr6553/mahautgr/INTERNSHIP_Elis/ISLAMIC_Samples/5_Dedup_Fastp/)
DirMap=$DirProjet/BWA/
DirSam=$DirProjet/SAMTOOLS/
DirMamD=$DirProjet/MAP_DAMAGE/
DirFast=$DirProjet/FASTA/
DirHTS=$DirProjet/ANGSD/
DirREPORT=$DirProjet/REPORT/

REF_DOG=/groups/dog/data/canFam3/sequence/bwa_index/Canis_familiaris.CanFam3.1.72.dna_sm.toplevel.fa

for file in ${DirRaw}*_R1.fastq.gz; do 
	name=${file##*/}
	sample_tag=${name%_R*}

	echo -e "\n_________________________ Sample "$sample_tag" _________________________"

	# echo "___________ 6 - BWA ___________" # ETAPE 6 - Alignement sur le génome de référence
	# gunzip -c $DirDdup${sample_tag}_rmduplicates.fastq.gz > $DirDdup${sample_tag}_rmduplicates.fastq
	# sed -n '1~4s/^@/>/p;2~4p' $DirDdup${sample_tag}_rmduplicates.fastq > $DirDdup${sample_tag}_rmduplicates.fasta # on met au format fasta
	# source /local/env/envbwa-0.7.17.sh
	# bwa aln -t 8 -l 1024 -o 2 -n 0.01 $REF_DOG $DirDdup${sample_tag}_rmduplicates.fasta > $DirMap${sample_tag}_DOG_aln.sai
	# bwa samse $REF_DOG $DirMap${sample_tag}_DOG_aln.sai $DirDdup${sample_tag}_rmduplicates.fastq.gz > $DirMap${sample_tag}_DOG_aln_samse.sam

	echo "___________ 7 - SAMTOOLS ___________" # ETAPE 7 - Extraction des reads mitochondriaux
    source /local/env/envsamtools-1.15.sh

    # # on filtre les unmapped
	samtools view -bT $DirDdup${sample_tag}_rmduplicates.fasta -F 4 $DirMap${sample_tag}_DOG_aln_samse.sam > $DirSam${sample_tag}_DOG_aln_samse_mapped.sam 
    # on filtre les reads d'une qualité < 25
    samtools view -h -q 25 $DirSam${sample_tag}_DOG_aln_samse_mapped.sam > $DirSam${sample_tag}_DOG_aln_samse_mappedQ25.sam 
    # sam to bam mapped Q25
    samtools view -S -b $DirSam${sample_tag}_DOG_aln_samse_mappedQ25.sam > $DirSam${sample_tag}_DOG_aln_samse_mappedQ25.bam 
    #sort par read_name
    samtools sort -n --threads 6 $DirSam${sample_tag}_DOG_aln_samse_mappedQ25.bam -o $DirSam${sample_tag}_DOG_sort.ba
    #fixmate
    samtools fixmate -m --threads 6 $DirSam${sample_tag}_DOG_sort.bam $DirSam${sample_tag}_DOG_sort_fixmate.bam
    #sort per coordinates
    samtools sort --threads 6 -o $DirSam${sample_tag}_DOG_sort_fixmate.sorted.bam $DirSam${sample_tag}_DOG_sort_fixmate.bam
    ##flagstat du markdup
    samtools flagstat $DirSam${sample_tag}_DOG_sort_fixmate.sorted.bam > $DirSam${sample_tag}_DOG_sort_fixmate_flagstat.sorted.txt ### attention à l'extention là
    ##index le fichier markdup.bam
    samtools index -@ 8 $DirSam${sample_tag}_DOG_sort_fixmate.sorted.bam > $DirSam${sample_tag}_DOG_sort_fixmate.sorted.bam.bai

    ##extraire les MT
    samtools view $DirSam${sample_tag}_DOG_sort_fixmate.sorted.bam MT -o $DirSam${sample_tag}_DOG_sort_fixmate_MT.sorted.bam
    samtools view $DirSam${sample_tag}_DOG_sort_fixmate_MT.sorted.bam > $DirSam${sample_tag}_DOG_sort_fixmate_MT_view.sorted.bam
    
    samtools flagstat $DirSam${sample_tag}_DOG_sort_fixmate_MT.sorted.bam > $DirSam${sample_tag}_DOG_sort_fixmate_MT_flagstat.sorted.txt

    samtools coverage $DirSam${sample_tag}_DOG_sort_fixmate_MT.sorted.bam > $DirOut${sample_tag}_DOG_samcoverage.txt
    samtools depth $DirSam${sample_tag}_DOG_sort_fixmate_MT.sorted.bam > $DirOut${sample_tag}_DOG_samdepth.txt

    echo "___________ 7bis - Rescale ___________"
    source /local/env/envconda.sh
    conda activate /groups/Paleogenomics/ENV/MapDamage-2.2.1

    mapDamage -i $DirSam${sample_tag}_DOG_sort_fixmate_MT.sorted.bam -r $REF_DOG --rescale 

    echo "_______________________ 8 - AngSD _________________________" #en in on prend le fichier .bam des mapped à q25
    #en out on obtient une séquence consensus du génome mitochondrial pour chaque échantillon
    
    source /local/env/envconda.sh
    conda activate /home/genouest/cnrs_umr6553/mahautgr/my_env/Angsd

    angsd -doFasta 2 -r nomeChrMito: -out prefix.out -explode 1 -doCounts 1 -minQ 20 -minMapQ 30 -i bamfile -setMinDepth 10 -nThreads 20 -doDepth 1
    angsd -doFasta 2  -out prefix.out -explode 1 -doCounts 1 -minQ 20 -minMapQ 30 -i rescaledbamfile -setMinDepth 3 -nThreads 20 -doDepth 1

done

