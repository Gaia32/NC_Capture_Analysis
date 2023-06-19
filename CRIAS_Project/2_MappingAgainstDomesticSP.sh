#!/bin/bash
#SBATCH --job-name=paleo

date

#############################################################################################
#	 	This is a mitogenome analysis for 2 unidentified sequences          				#
# 		This part is the mapping against 12 differents species (only mitogenome): 		    #
# dog, human, cow, sheep, iberian wolf, grey wolf, goat, pig, deer, chicken, horse, donkey  #
# 		Then we compare the coverage between the different species 							#
#############################################################################################

DirProjet=(/home/genouest/cnrs_umr6553/mahautgr/INTERNSHIP_Elis/ISLAMIC_Samples) # path to your project
cd $DirProjet

DirRaw=$DirProjet/RAW_Data/Pig_data/
DirMulti=$DirProjet/6_Multimapping/
DirDdup=$DirProjet/5_Dedup_Fastp/
DirOut=$DirProjet/OUT/
DirREF=/scratch/mahautgr/REF_GENOMES/MITOCHONDRION/ #path to the refreenc mitogenomes 
ANIMALS=("CHICKEN" "COW" "DEER" "DOG" "GOAT" "GREY_WOLF" "HUMAN" "IBERIAN_WOLF" "PIG" "SHEEP" "HORSE" "DONKEY")

for file in ${DirRaw}*_R1.fastq.gz; do 
	name=${file##*/}
	sample_tag=${name%_R*}
	echo "_________________________ Sample "$sample_tag" _________________________"

	for animal in ${ANIMALS[@]}; do
		REF_GEN=$DirREF${animal}"_Ref_mitochondrion/"*.fasta
		echo $DirDdup${sample_tag}_rmduplicates.fasta
		echo $animal

		echo "___________ 6 - Multimapping ___________"
		source /local/env/envbwa-0.7.17.sh
		bwa aln -t 8 -l 1024 -o 2 -n 0.01 $REF_GEN $DirDdup${sample_tag}_rmduplicates.fasta > $DirMulti"BWA/"${sample_tag}"_"${animal}_aln.sai
		bwa samse $REF_GEN $DirMulti"BWA/"${sample_tag}"_"${animal}_aln.sai $DirDdup${sample_tag}_rmduplicates.fastq.gz > $DirMulti"BWA/"${sample_tag}"_"${animal}_aln_samse.sam

		echo "___________ 7 - SAMTOOLS ___________"
		source /local/env/envsamtools-1.15.sh

		# on filtre les unmapped
		samtools view -bT $DirDdup${sample_tag}_rmduplicates.fasta -F 4 $DirMulti"BWA/"${sample_tag}"_"${animal}_aln_samse.sam > $DirMulti"SAMTOOLS/"${sample_tag}"_"${animal}_aln_samse_mapped.sam 
		# on filtre les reads d'une qualité < 25
		samtools view -h -q 25 $DirMulti"SAMTOOLS/"${sample_tag}"_"${animal}_aln_samse_mapped.sam  > $DirMulti"SAMTOOLS/"${sample_tag}"_"${animal}_aln_samse_mappedQ25.sam 
		# sam to bam mapped Q25
		samtools view -S -b $DirMulti"SAMTOOLS/"${sample_tag}"_"${animal}_aln_samse_mappedQ25.sam  > $DirMulti"SAMTOOLS/"${sample_tag}"_"${animal}_aln_samse_mappedQ25.bam 
		# sort par read_name
		samtools sort --threads 6 $DirMulti"SAMTOOLS/"${sample_tag}"_"${animal}_aln_samse_mappedQ25.bam  -o $DirMulti"SAMTOOLS/"${sample_tag}"_"${animal}_aln_samse_mappedQ25.sorted.bam 
		##flagstat du markdup
    	samtools flagstat $DirMulti"SAMTOOLS/"${sample_tag}"_"${animal}_aln_samse_mappedQ25.sorted.bam > $DirMulti"REPORT/"${sample_tag}"_"${animal}_flagstat.sorted.txt
		# coverage 
		samtools coverage $DirMulti"SAMTOOLS/"${sample_tag}"_"${animal}_aln_samse_mappedQ25.sorted.bam > $DirMulti"REPORT/"${sample_tag}"_"${animal}_mappedQ25_coverage.txt


	done
done

date

# to count the number of reads with 3x or 4x
# cat OUT/P27761_1057_DOG_samdepth.txt | awk 'FRN=1 {print $3}' | grep -e "3" -e "4" | wc -l


