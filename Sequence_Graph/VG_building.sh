#!/bin/bash
#SBATCH --job-name=VG_build
#SBATCH --output=VG_construct.out
#SBATCH --mem=150G
#SBATCH --cpus-per-task=16

###########################################
# Use of the chr 1 (301 varianst for the 4000 snips) for testing the sequence graph
###########################################

GENOME=/groups/dog/data/canFam3/sequence/bwa_index/Canis_familiaris.CanFam3.1.72.dna_sm.toplevel.fa
SNIPs=/scratch/mahautgr/REF_GENOMES/SNIPs/722g.990.SNP.INDEL.chrAll.vcf.gz

source /local/env/envconda.sh
conda activate /home/genouest/cnrs_umr6553/mahautgr/my_env/Vg

  vg construct -C -R chr1 -r $GENOME -v $SNIPs >1.vg
  vg index -x 1.xg -g 1.gcsa 1.vg
  
########################################### 
# Paired end reads interleaved in a single FASTQ
# mapping paired reads
# then sortng and indexing reads and calling variants using bcftools
###########################################
  vg map -x 1.xg -g 1.gcsa -f ${sample_tag}_rmduplicates.fq -i > ${sample_tag}_mapped.gam
  
  #sort per coordinates + index
  samtools sort -o 1822_vg_mapped.sorted.bam 1822_vg_mapped.bam
  samtools index 1822_vg_mapped.sorted.bam

	##flagstat des mapped Q25
  samtools flagstat 1822_vg_mapped.sorted.bam > 1822_vg_flagstat_mapped.txt
	samtools coverage 1822_vg_mapped.sorted.bam > 1822_vg_samcoverage_mappedQ25.txt
	samtools depth 1822_vg_mapped.sorted.bam > 1822_vg_samdepth_mapped.txt

	echo -e "\n___________ 8 - VARIANT CALLING ___________"
	source /local/env/envbcftools-1.9.sh
	time bcftools mpileup -T "/groups/Paleogenomics/DOG/Adapters_&_SNIP/"4003.txt -C0 -B -A -f "/groups/dog/data/canFam3/sequence/bwa_index/Canis_familiaris.CanFam3.1.72.dna_sm.toplevel.fa" --threads 16 1822_vg_mapped.sorted.bam | bcftools call -m -o 1822_vg_mapped.sorted.vcf

  
###########################################
# This script allows to convert the gam file obtained by the seuquence graph 
# into a bam to be able to call the variants using bcftools
###########################################
#surject the alignments back into the reference space of sequence "x", yielding a BAM file

  vg surject -x 1.xg -b ${sample_tag}_mapped.gam > ${sample_tag}_mapped.bam
  
###########################################
# calling variants
# giving sme ideas but I am not sure of which one is the best
###########################################
  vg augment -a pileup -Z ${sample_tag}.trans -S ${sample_tag}.support 1.vg ${sample_tag}_mapped.gam > ${sample_tag}.aug.vg
  vg augment 1.vg ${sample_tag}_mapped.gam -i -S > ${sample_tag}_mapped_aug.vg
  vg call ${sample_tag}_mapped_aug.vg > ${sample_tag}.vcf
  
  vg augment 1.vg ${sample_tag}_mapped.gam -A ${sample_tag}_aug.gam > ${sample_tag}.aug.vg
  vg call -z ${sample_tag}.trans -s ${sample_tag}.support ${sample_tag}.aug.vg > ${sample_tag}.vcf
  
  vg pack -x 1.xg -g ${sample_tag}_mapped.gam -Q 5 -o ${sample_tag}_mapped.pack
  vg call 1.xg -k ${sample_tag}_mapped.pack > ${sample_tag}_mapped.vcf
  
  # Compute the read support from the gam
  # -Q 5: ignore mapping and base qualitiy < 5
  # -s 5: ignore first and last 5bp from each read
  vg pack -x 1.xg -g ${sample_tag}_mapped.gam -o ${sample_tag}_mapped_nofilter.pack

  # Generate a VCF from the support.  
  vg call 1.xg -k ${sample_tag}_mapped_nofilter.pack > ${sample_tag}_graph_calls_nofilter.vcf
