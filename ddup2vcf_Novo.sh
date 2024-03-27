#!/bin/bash
#SBATCH --job-name=NovProcess
#SBATCH --output=NovoProcessing.out
#SBATCH --mem=100G
#SBATCH --cpus-per-task=16

#########################################
# this script allows to map reads against a modified reference after deduplication step to variant calling
# it also allows to create the modified reference using the variants from the plassais 2019 722g 
# NovoAlign needs a licence to be used
#########################################

DirProjet=(path/to/the/project/directory)
cd $DirProject
mkdir -p OUT 6_Novo_Aln 7_Novo_VCF

DirDdup=$DirProjet/5_Dedup_Fastp/
DirOut=$DirProjet/OUT/
DirAln=$DirProjet/6_Novo_Aln/
DirVcf=$DirProjet/7_Novo_VCF/

# genome de ref
GENOME=(path/to/the/reference/genome/directory)
DirSNIPS=(path/to/the/VCF/with/variants/directory)
novocraft=(path/to/the/novoalign)
novoutil=${novocraft}novoutil
novoindex=${novocraft}novoindex
novoalign=${novocraft}novoalign
novosort=${novocraft}novosort

#step 1) modify reference genome with all 4000 snips
  #setp 2) index ref genome
  #step 3) map against IUPAC genome for the 4000
  #step 4) sort the mapped reads
  #step 4) call variants using either bcftools 

#_______________________________________________
#_______________________________________________
  
  $novoutil iupac ${DirSNIPS}722g.990.SNP.INDEL.chrAll.vcf.gz ${DirGENOME}Canis_familiaris.CanFam3.1.72.dna_sm.toplevel.fa > /scratch/mahautgr/REF_GENOMES/CanFam3.1.72.IUPAC_TopCommonSNIPs.fa
  $novoindex /scratch/mahautgr/REF_GENOMES/CanFam3.1.72.IUPAC_TopCommonSNIPs.nix /scratch/mahautgr/REF_GENOMES/CanFam3.1.72.IUPAC_TopCommonSNIPs.fa

#_______________________________________________
#_______________________________________________

for file in ${DirDdup}*_rmduplicates.fastq.gz; do
	name=${file##*/}
	sample_tag=${name%_rmdup*}
	echo -e "\n_________________________ Sample "$sample_tag" _________________________"
 
    echo -e "\n___________ 6 - NovoAlign ___________"
    $novoalign -d /scratch/mahautgr/REF_GENOMES/CanFam3.1.72.IUPAC_TopCommonSNIPs.nix -f $DirDdup${sample_tag}_rmduplicates.fastq.gz -k --addM5 off -o BAM > $DirAln${sample_tag}_novo_mapped.bam
    #OPTIONS --addM5 off -o BAM
    #-k to allow quality recalibration
    #if only one file in -f ==> SE mode

    $novosort -m 12G -i -o $DirAln${sample_tag}_novo_mapped.sorted.bam $DirAln${sample_tag}_novo_mapped.bam
    #Minimum Memory Usage (-m)
    # sort by coordinates

    echo -e "\n___________ 8 - VARIANT CALLING ___________"
	source /local/env/envbcftools-1.9.sh
	time bcftools mpileup -T "/groups/Paleogenomics/DOG/Adapters_&_SNIP/"4003.txt -C0 -B -A -f $GENOME --threads 8 $DirAln${sample_tag}_novo_mapped.sorted.bam | bcftools call -m -o $DirVcf${sample_tag}_novo_mapped.sorted.bcftools.vcf
  
    #keep only 3X depth VCFs
	awk '/DP=0/||/DP=1;/||/DP=2;/ {next;} {print}' $DirVcf${sample_tag}_novo_mapped.sorted.bcftools.vcf > $DirVcf${sample_tag}_novo_mapped.sorted.bcftools.3x.vcf

done
