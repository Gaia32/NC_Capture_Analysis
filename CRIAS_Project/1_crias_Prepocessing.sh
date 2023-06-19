#!/bin/bash
#SBATCH --job-name=paleo
# Thx Deborah Diquelou for some parts of her script.

date

####################################################################
# This is a mitoREF_DOG analysis for 2 unidentified sequences       #
# After consensus building, we run a blast to idenitify the specie #
####################################################################

DirProjet=(/home/genouest/cnrs_umr6553/mahautgr/INTERNSHIP_Elis/ISLAMIC_Samples) # <------------- + faire attention à ne pas mettre de / à la fin
cd $DirProjet

mkdir -p OUT 1_Cutadapt 2_Rmadapt_Fastp 3_Merging_Fastp 4_Filt_Fastp 5_Dedup_Fastp 6_Mapping 7_Samtools 8_Fasta 9_Htsbox

DirRaw=(/home/genouest/cnrs_umr6553/mahautgr/INTERNSHIP_Elis/ISLAMIC_Samples/RAW_Data/Pig_data/) # <----------- Ici le path des données brutes
DirCut=$DirProjet/1_Cutadapt/
DirRmad=$DirProjet/2_Rmadapt_Fastp/
DirMerg=$DirProjet/3_Merging_Fastp/
DirFilt=$DirProjet/4_Filt_Fastp/
DirDdup=$DirProjet/5_Dedup_Fastp/
DirMap=$DirProjet/6_Mapping/
DirSam=$DirProjet/7_Samtools/
DirMamD=$DirProjet/7bis_Map_Damage/
DirFast=$DirProjet/8_Fasta/
DirHTS=$DirProjet/9_Htsbox/
DirOut=$DirProjet/OUT/

for file in ${DirRaw}*_R1.fastq.gz; do 
	name=${file##*/}
	sample_tag=${name%_R*}

	echo "\n$sample_tag"
	echo "_________________________ Sample "$sample_tag" _________________________"

	echo "___________ 1 - RM ADAPT - CUTADAPT ___________" # ETAPE 1 - Elimination des adaptateurs sur R1 et R2
	source /local/env/envcutadapt-1.15.sh
	cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -O 3 -e 0.5 -m 1 -o $DirCut${sample_tag}_R1_cutadapt.fastq.gz --untrimmed-output=$DirCut${sample_tag}_R1_untrim.fastq.gz $DirRaw${sample_tag}_R1.fastq.gz > $DirOut${sample_tag}_R1_cutadapt_report.txt 2>&1
	cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -O 3 -e 0.5 -m 1 -o $DirCut${sample_tag}_R2_cutadapt.fastq.gz --untrimmed-output=$DirCut${sample_tag}_R2_untrim.fastq.gz $DirRaw${sample_tag}_R2.fastq.gz > $DirOut${sample_tag}_R2_cutadapt_report.txt 2>&1

	echo "___________ 2 - RM ADAPT - fastp ___________" ## ETAPE 2 - $fastp_env
	source /local/env/envconda.sh
    fastp_env=/home/genouest/cnrs_umr6553/mahautgr/my_env/Fastp/fastp_0_23_1

	$fastp_env -i $DirCut${sample_tag}_R1_cutadapt.fastq.gz -o $DirRmad${sample_tag}_R1_rmadapt.fastq --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --disable_quality_filtering --disable_length_filtering --disable_trim_poly_g --dont_eval_duplication -h $DirRmad${sample_tag}_R1_rmadapt.html>$DirOut${sample_tag}_R1_rmadapt_report.txt 2>&1
	$fastp_env -i $DirCut${sample_tag}_R2_cutadapt.fastq.gz -o $DirRmad${sample_tag}_R2_rmadapt.fastq --adapter_sequence AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --disable_quality_filtering --disable_length_filtering --disable_trim_poly_g --dont_eval_duplication -h $DirRmad${sample_tag}_R2_rmadapt.html>$DirOut${sample_tag}_R2_rmadapt_report.txt 2>&1

    gzip $DirRmad${sample_tag}_R1_rmadapt.fastq
    gzip $DirRmad${sample_tag}_R2_rmadapt.fastq

	echo "___________ 3 - PAIRING - MERGING - BBTOOLS / fastp ___________" # ETAPE 3 - Merging et filtration
	source /local/env/envconda.sh
	conda activate /home/genouest/cnrs_umr6553/ddiquelou/bbtools
    
	echo -e "REPAIR"
    repair.sh in1=$DirRmad${sample_tag}_R1_rmadapt.fastq.gz in2=$DirRmad${sample_tag}_R2_rmadapt.fastq.gz --overwrite=true out=$DirMerg${sample_tag}_R1_repair.fastq.gz out2=$DirMerg${sample_tag}_R2_repair.fastq.gz outs=$DirMerg${sample_tag}_repair_singleton.fastq.gz>$DirOut${sample_tag}_repair_report.txt 2>&1

    echo -e "\nMERGE"
    $fastp_env --in1 $DirMerg${sample_tag}_R1_repair.fastq.gz --in2 $DirMerg${sample_tag}_R2_repair.fastq.gz -A --disable_adapter_trimming --disable_length_filtering --disable_quality_filtering --dont_eval_duplication --overlap_len_require 25 -m --merged_out $DirMerg${sample_tag}_merged.fastq.gz --out1 $DirMerg${sample_tag}_R1_unmerged.fastq.gz --out2 $DirMerg${sample_tag}_R2_unmerged.fastq.gz --unpaired1 $DirMerg${sample_tag}_R1_singleton.fastq.gz --unpaired2 $DirMerg${sample_tag}_R2_singleton.fastq.gz --failed_out $DirMerg${sample_tag}_merged_failed.fastq.gz>$DirOut${sample_tag}_merging_report.txt 2>&1

# on concatène tout, y compris les singletons issus du pairing
	cat $DirMerg${sample_tag}_merged.fastq.gz $DirMerg${sample_tag}_R1_unmerged.fastq.gz $DirMerg${sample_tag}_R2_unmerged.fastq.gz $DirMerg${sample_tag}_R1_singleton.fastq.gz $DirMerg${sample_tag}_R2_singleton.fastq.gz $DirMerg${sample_tag}_repair_singleton.fastq.gz > $DirMerg${sample_tag}_merging_all.fastq.gz

	echo "___________ 4 - FILTRATION - fastp ___________" # ETAPE 4 - Filtration + longueur
    $fastp_env -i $DirMerg${sample_tag}_merging_all.fastq.gz -o $DirFilt${sample_tag}_filt.fastq.gz -A --disable_adapter_trimming --dont_eval_duplication -l 25 --average_qual 20 --trim_poly_g --trim_poly_x --low_complexity_filter --cut_right --cut_right_mean_quality 20 >$DirOut${sample_tag}_filt_report.txt 2>&1

	echo "___________ 5 - DEDUP - fastp ___________"#  ETAPE 5 - Suppression des duplicats de PCR
	$fastp_env -i $DirFilt${sample_tag}_filt.fastq.gz --dedup -o $DirDdup${sample_tag}_rmduplicates.fastq.gz -A --disable_adapter_trimming --disable_quality_filtering --disable_length_filtering -h $DirDdup${sample_tag}_rmduplicate.html>$DirOut${sample_tag}_rmduplicates_report.txt 2>&1

	echo "___________ 6 - BWA ___________" # ETAPE 6 - Alignement sur le génome de référence
	sed -n '1~4s/^@/>/p;2~4p' $DirDdup${sample_tag}_rmduplicates.fastq.gz > $DirDdup${sample_tag}_rmduplicates.fasta # on met au format fasta

###### irst we try with the dog ref genome
	REF_DOG=/groups/dog/data/canFam3/sequence/bwa_index/Canis_familiaris.CanFam3.1.72.dna_sm.toplevel.fa

	source /local/env/envbwa-0.7.17.sh

	bwa aln -t 8 -l 1024 -o 2 -n 0.01 $REF_DOG $DirDdup${sample_tag}_rmduplicates.fasta > $DirMap${sample_tag}_aln.sai
	bwa samse $REF_DOG $DirMap${sample_tag}_aln.sai $DirDdup${sample_tag}_rmduplicates.fastq.gz > $DirMap${sample_tag}_aln_samse.sam

	echo "___________ 7 - SAMTOOLS ___________" # ETAPE 7 - Extraction des reads mitochondriaux
    source /local/env/envsamtools-1.15.sh

    # on filtre les unmapped
	samtools view -bT $DirDdup${sample_tag}_rmduplicates.fasta -F 4 $DirMap${sample_tag}_aln_samse.sam > $DirSam${sample_tag}_aln_samse_mapped.sam 
    
    # on filtre les reads d'une qualité < 25
    samtools view -h -q 25 $DirSam${sample_tag}_aln_samse_mapped.sam > $DirSam${sample_tag}_aln_samse_mappedQ25.sam 
    
    # sam to bam mapped Q25
    samtools view -S -b $DirSam${sample_tag}_aln_samse_mappedQ25.sam > $DirSam${sample_tag}_aln_samse_mappedQ25.bam 
    
    #sort par read_name
    samtools sort -n --threads 6 $DirSam${sample_tag}_aln_samse_mappedQ25.bam -o $DirSam${sample_tag}_sort.bam

    #fixmate
    samtools fixmate -m --threads 6 $DirSam${sample_tag}_sort.bam $DirSam${sample_tag}_sort_fixmate.bam

    #sort per coordinates
    samtools sort --threads 6 -o $DirSam${sample_tag}_sort_fixmate.sorted.bam $DirSam${sample_tag}_sort_fixmate.bam

    ##flagstat du markdup
    samtools flagstat $DirSam${sample_tag}_sort_fixmate.sorted.bam > $DirSam${sample_tag}_sort_fixmate_flagstat.sorted.bam ### attention à l'extention là

    ##index le fichier markdup.bam
    samtools index -@ 8 $DirSam${sample_tag}_sort_fixmate.sorted.bam > $DirSam${sample_tag}_sort_fixmate.sorted.bam.bai

    ##extraire les MT
    samtools view $DirSam${sample_tag}_sort_fixmate.sorted.bam MT -o $DirSam${sample_tag}_sort_fixmate_MT.sorted.bam
    samtools view $DirSam${sample_tag}_sort_fixmate_MT.sorted.bam > $DirSam${sample_tag}_sort_fixmate_MT_view.sorted.bam
    
    samtools flagstat $DirSam${sample_tag}_sort_fixmate_MT.sorted.bam > $DirSam${sample_tag}_sort_fixmate_MT_flagstat.sorted.txt

    samtools coverage $DirSam${sample_tag}_sort_fixmate_MT.sorted.bam > $DirOut${sample_tag}_samcoverage.txt
    samtools depth $DirSam${sample_tag}_sort_fixmate_MT.sorted.bam > $DirOut${sample_tag}_samdepth.txt

    echo "___________ 8 - Rescale ___________"
    source /local/env/envconda.sh
    conda activate /groups/Paleogenomics/ENV/MapDamage-2.2.1

    mapDamage -i $DirSam${sample_tag}_sort_fixmate_MT.sorted.bam -r $REF_DOG --rescale

    echo "_______________________ 9 - AngSD _________________________" #en in on prend le fichier .bam des mapped à q25
    #en out on obtient une séquence consensus du génome mitochondrial pour chaque échantillon
    
    source /local/env/envconda.sh
    conda activate /home/genouest/cnrs_umr6553/mahautgr/my_env/Angsd

#Extrair consensus mito de um bam com genoma completo

    # angsd -doFasta 2 -r nomeChrMito: -out prefix.out -explode 1 -doCounts 1 -minQ 20 -minMapQ 30 -i bamfile -setMinDepth 10 -nThreads 20 -doDepth 1
    # angsd -doFasta 2  -out prefix.out -explode 1 -doCounts 1 -minQ 20 -minMapQ 30 -i rescaledbamfile -setMinDepth 3 -nThreads 20 -doDepth 1

done

date