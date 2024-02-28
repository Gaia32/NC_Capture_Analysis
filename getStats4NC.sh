#!/bin/bash
#SBATCH --job-name=Report
#SBATCH --output=Report.out
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4

####################################################################
# This script allows to create a table (columns separated by ",")  #
# This table contains the number of reads filtered at each step    #
# of the pipeline, from raw fastq to the number of SNIPs extracted #
# The _suffixes at each step are completely subjectives, therefor  #
# make sure the suffixes and names of folders are well written.    #
# Also this was used to compare two aligners, novoalign and bwa aln#
# you can remove the part on novoalign if it was not used. 		   #
####################################################################

DirProjet=(pathToYourProject)
DirRaw=(pathToYourRawFastq/)
DirCut=$DirProjet/1_Cutadapt/
DirRmad=$DirProjet/2_Rmadapt_Fastp/
DirMerg=$DirProjet/3_Merging_Fastp/
DirFilt=$DirProjet/4_Filt_Fastp/
DirDdup=$DirProjet/5_Dedup_Fastp/
DirMapBWA=$DirProjet/6_BWA_aln/
DirMapNOV=$DirProjet/6_Novo_Aln/
DirVcfBWA=$DirProjet/7_BWA_VCF/
DirVcfNOV=$DirProjet/7_Novo_VCF/
DirOut=$DirProjet/OUT/

for file in ${DirRaw}*.fastq.gz; do 
	name=${file##*/}
	sample_tag=${name%_Mt*}

	echo "Sample_Id,"$sample_tag
	echo -n "Total reads R1," ; zgrep -c "@" $DirRaw${sample_tag}*_R1_001.fastq.gz
	echo -n "Total reads R2," ; zgrep -c "@" $DirRaw${sample_tag}*_R2_001.fastq.gz
	echo -n "R1 trimmed," ; zgrep -c "@" $DirCut${sample_tag}_R1_cutadapt.fastq.gz
	echo -n "R1 untrimmed," ; zgrep -c "@" $DirCut${sample_tag}_R1_untrim.fastq.gz
	echo -n "R2 trimmed," ; zgrep -c "@" $DirCut${sample_tag}_R2_cutadapt.fastq.gz
	echo -n "R2 untrimmed," ; zgrep -c "@" $DirCut${sample_tag}_R2_untrim.fastq.gz

    echo -n "Nbr reads out Rmadapt R1," ; zgrep -c "@" $DirRmad${sample_tag}_R1_rmadapt.fastq.gz
    echo -n "Nbr reads out Rmadapt R2," ; zgrep -c "@" $DirRmad${sample_tag}_R2_rmadapt.fastq.gz
    
	echo -n "R1 repair," ; zgrep -c "@" $DirMerg${sample_tag}_R1_repair.fastq.gz
    echo -n "R2 repair,"; zgrep -c "@" $DirMerg${sample_tag}_R2_repair.fastq.gz
    
    echo -n "Singleton repair," ; zgrep -c "@" $DirMerg${sample_tag}_repair_singleton.fastq.gz
	echo -n "merged," ; zgrep -c "@" $DirMerg${sample_tag}_merged.fastq.gz
	echo -n "unmerged R1," ; zgrep -c "@" $DirMerg${sample_tag}_R1_unmerged.fastq.gz
	echo -n "unmerged R2," ; zgrep -c "@" $DirMerg${sample_tag}_R2_unmerged.fastq.gz
	echo -n "singleton R1," ; zgrep -c "@"  $DirMerg${sample_tag}_R1_singleton.fastq.gz
	echo -n "singleton R2," ; zgrep -c "@" $DirMerg${sample_tag}_R2_singleton.fastq.gz
	echo -n "failed merging," ; zgrep -c "@" $DirMerg${sample_tag}_merged_failed.fastq.gz
    echo -n "merged and unmerged," ; zgrep -c "@" $DirMerg${sample_tag}_merged_umerged.fastq.gz
	cat $DirOut${sample_tag}_filt_report.txt | sed -n -e "
	/total reads:/p
	/Read pairs merged:/p
	/reads failed due to low quality:/p
	/reads failed due to too many N:/p
	/reads failed due to low complexity:/p
	/reads failed due to too short:/p
	/reads with polyX in 3' end:/p
    /reads corrected by overlap analysis:/p"| sed -r 's/(.*) /\1,/;'
	echo -n "Reads all," ; zgrep -c "@" $DirFilt${sample_tag}_filt.fastq.gz
	echo -n "Nbr of reads left after deduplication," ; zgrep -c "@" $DirDdup${sample_tag}_rmduplicates.fastq.gz
    cat $DirOut${sample_tag}_rmduplicates_report.txt | sed -n -e "
	/total reads:/p
    /Duplication rate /p" | sed -r 's/(.*) /\1,/;'
    
	echo -n "mapped BWA q25," ; awk 'FNR==1 {print $1}' $DirMapBWA${sample_tag}_aln_samse_mappedQ25.sorted.txt
  echo -n "mapped sam NOVO," ;awk 'FNR==7 {print $1}' $DirMapNOV${sample_tag}_novo_mapped.sorted.txt
  
	echo -n "Nb of variants called with BWA," ; grep -c "DP=" $DirVcfBWA${sample_tag}_mappedQ25.sorted.vcf
	echo -n "Nb of variants called with 3x depth with BWA," ; grep -c "DP=" $DirVcfBWA${sample_tag}_mappedQ25.sorted.3x.vcf
  echo -n "Nb of variants called with NOVO," ; grep -c "DP=" $DirVcfNOV${sample_tag}_novo_mapped.sorted.vcf
	echo -n "Nb of variants called with 3x depth with NOVO," ; grep -c "DP=" $DirVcfNOV${sample_tag}_novo_mapped.sorted.3x.vcf

done

