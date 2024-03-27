#!/bin/bash
#SBATCH --job-name=raw2vcf
#SBATCH --output=raw2vcf.out
#SBATCH --mem=80G
#SBATCH --cpus-per-task=10

date

############################################################################
#       This is a script for nuclear genome analysis of 4000 SNIPs         #
#            This is all the steps from raw data to variants calling       #
# this script allows to process the rawdata fromt target-enrichment assay   #
# to the variant calling format file using BWA aln as the aligner          #
# all softwares must be installed prior to the processing                  #
############################################################################

DirProjet=(path/to/the/project/directory)
cd $DirProjet

mkdir -p OUT 1_Cutadapt 2_Rmadapt_Fastp 3_Merging_Fastp 4_Filt_Fastp 5_Dedup_Fastp 6_Mapping_BWA 7_Samtools 8_VCF

DirRaw=(path/to/the/raw/data) 
DirCut=$DirProjet/1_Cutadapt/
DirRmad=$DirProjet/2_Rmadapt_Fastp/
DirMerg=$DirProjet/3_Merging_Fastp/
DirFilt=$DirProjet/4_Filt_Fastp/
DirDdup=$DirProjet/5_Dedup_Fastp/
DirMap=$DirProjet/6_Mapping_BWA/
DirSam=$DirProjet/7_Samtools/
DirVcf=$DirProjet/8_VCF/
DirOut=$DirProjet/OUT/

GENOME=(path/to/the/reference/genome)
fastp_env=(path/to/the/fastp_env/fastp_0_23_1) #this version was used

for file in ${DirRaw}*_R1.fastq.gz; do #change here if the file doesn't end with _R1.fastq.gz
	name=${file##*/}
	sample_tag=${name%_R*}
	echo -e "\n_________________________ Sample "$sample_tag" _________________________"

	echo "___________ 1 - RM ADAPT - CUTADAPT ___________" 
	source /local/env/envcutadapt-1.15.sh
	cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -O 3 -e 0.5 -m 1 -o $DirCut${sample_tag}_R1_cutadapt.fastq.gz --untrimmed-output=$DirCut${sample_tag}_R1_untrim.fastq.gz $DirRaw${sample_tag}_R1.fastq.gz > $DirOut${sample_tag}_R1_cutadapt_report.txt 2>&1
	cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -O 3 -e 0.5 -m 1 -o $DirCut${sample_tag}_R2_cutadapt.fastq.gz --untrimmed-output=$DirCut${sample_tag}_R2_untrim.fastq.gz $DirRaw${sample_tag}_R2.fastq.gz > $DirOut${sample_tag}_R2_cutadapt_report.txt 2>&1

	echo "___________ 2 - RM ADAPT - fastp ___________" 
    

	$fastp_env -i $DirCut${sample_tag}_R1_cutadapt.fastq.gz -o $DirRmad${sample_tag}_R1_rmadapt.fastq.gz --adapter_fasta "/groups/Paleogenomics/DOG/Adapters_&_SNIP/"adapters.fasta --disable_quality_filtering --disable_length_filtering --disable_trim_poly_g --dont_eval_duplication -h $DirRmad${sample_tag}_R1_rmadapt.html>$DirOut${sample_tag}_R1_rmadapt_report.txt 2>&1
	$fastp_env -i $DirCut${sample_tag}_R2_cutadapt.fastq.gz -o $DirRmad${sample_tag}_R2_rmadapt.fastq.gz --adapter_fasta "/groups/Paleogenomics/DOG/Adapters_&_SNIP/"adapters.fasta --disable_quality_filtering --disable_length_filtering --disable_trim_poly_g --dont_eval_duplication -h $DirRmad${sample_tag}_R2_rmadapt.html>$DirOut${sample_tag}_R2_rmadapt_report.txt 2>&1

	echo "___________ 3 - PAIRING - MERGING - BBTOOLS / fastp ___________" 
	source /local/env/envconda.sh
	conda activate /path/to/bbtools_env # version 37.62
    
	echo -e "REPAIR"
     	repair.sh in1=$DirRmad${sample_tag}_R1_rmadapt.fastq.gz in2=$DirRmad${sample_tag}_R2_rmadapt.fastq.gz --overwrite=true out=$DirMerg${sample_tag}_R1_repair.fastq.gz out2=$DirMerg${sample_tag}_R2_repair.fastq.gz outs=$DirMerg${sample_tag}_repair_singleton.fastq.gz>$DirOut${sample_tag}_repair_report.txt 2>&1

    	echo -e "\nMERGE"
    	$fastp_env --in1 $DirMerg${sample_tag}_R1_repair.fastq.gz --in2 $DirMerg${sample_tag}_R2_repair.fastq.gz -A --disable_adapter_trimming --disable_length_filtering --disable_quality_filtering --dont_eval_duplication --overlap_len_require 6 -m --merged_out $DirMerg${sample_tag}_merged.fastq.gz --out1 $DirMerg${sample_tag}_R1_unmerged.fastq.gz --out2 $DirMerg${sample_tag}_R2_unmerged.fastq.gz --unpaired1 $DirMerg${sample_tag}_R1_singleton.fastq.gz --unpaired2 $DirMerg${sample_tag}_R2_singleton.fastq.gz --failed_out $DirMerg${sample_tag}_merged_failed.fastq.gz>$DirOut${sample_tag}_merging_report.txt 2>&1

    	# on concatène tout, SANS les singletons issus du pairing
	cat $DirMerg${sample_tag}_merged.fastq.gz $DirMerg${sample_tag}_R1_unmerged.fastq.gz $DirMerg${sample_tag}_R2_unmerged.fastq.gz > $DirMerg${sample_tag}_merged_umerged.fastq.gz

	echo "___________ 4 - FILTRATION - fastp ___________" 
    	$fastp_env -i $DirMerg${sample_tag}_merged_umerged.fastq.gz -o $DirFilt${sample_tag}_filt.fastq.gz -A --disable_adapter_trimming --dont_eval_duplication -l 25 --average_qual 20 --trim_poly_g --trim_poly_x --low_complexity_filter --cut_right --cut_right_mean_quality 20 >$DirOut${sample_tag}_filt_report.txt 2>&1

	echo "___________ 5 - DEDUP - fastp ___________"#  
	$fastp_env -i $DirFilt${sample_tag}_filt.fastq.gz --dedup -o $DirDdup${sample_tag}_rmduplicates.fastq.gz -A --disable_adapter_trimming --disable_quality_filtering --disable_length_filtering -h $DirDdup${sample_tag}_rmduplicate.html>$DirOut${sample_tag}_rmduplicates_report.txt 2>&1

    	echo -e "\n___________ 6 - BWA ___________"

    	gunzip -c $DirDdup${sample_tag}_rmduplicates.fastq.gz > $DirDdup${sample_tag}_rmduplicates.fq
	sed -n '1~4s/^@/>/p;2~4p' $DirDdup${sample_tag}_rmduplicates.fq > $DirDdup${sample_tag}_rmduplicates.fasta
	
    	source /local/env/envbwa-0.7.17.sh
	bwa aln -t 8 -l 1024 -o 2 -n 0.01 $GENOME $DirDdup${sample_tag}_rmduplicates.fasta > $DirAln${sample_tag}_aln.sai
	bwa samse $GENOME $DirAln${sample_tag}_aln.sai $DirDdup${sample_tag}_rmduplicates.fq > $DirAln${sample_tag}_aln_samse.sam
    	rm $DirDdup${sample_tag}_rmduplicates.fq

    	echo -e "\n___________ 7 - SAMTOOLS ___________"
    	source /local/env/envsamtools-1.15.sh
    	samtools view -b -h -q 25 -F 4 $DirAln${sample_tag}_aln_samse.sam > $DirAln${sample_tag}_aln_samse_mappedQ25.bam 
    	samtools sort -o $DirAln${sample_tag}_aln_samse_mappedQ25.sorted.bam $DirAln${sample_tag}_aln_samse_mappedQ25.bam
    	samtools index $DirAln${sample_tag}_aln_samse_mappedQ25.sorted.bam

   	echo -e "\n___________ 8 - VARIANT CALLING ___________"
	source /local/env/envbcftools-1.9.sh
	time bcftools mpileup -T "/groups/Paleogenomics/DOG/Adapters_&_SNIP/"4003.txt -C0 -B -A -f $GENOME --threads 8 $DirAln${sample_tag}_aln_samse_mappedQ25.sorted.bam | bcftools call -m -o $DirVcf${sample_tag}_mappedQ25.sorted.vcf
    	#keep only 3X depth VCFs
	awk '/DP=0/||/DP=1;/||/DP=2;/ {next;} {print}' $DirVcf${sample_tag}_mappedQ25.sorted.vcf > $DirVcf${sample_tag}_mappedQ25.sorted.3x.vcf

done

date
