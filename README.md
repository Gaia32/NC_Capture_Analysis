# NC_Capture_Analysis

____*This folder contains scripts which where used to process the nuclear enrichment captures.*____


Description of the scripts & folders:

########################## **CRIAS_Project Folder**

The CRIAS projet was a project with Ana Elisabete Pires and archezoologists. 
The project consisted of analysing two sequences from unknown origin which where supposedly dogs or pigs.
Inside this folder, you can find scripts which where used to find out whihc of those two they where.

########################## **Sequence_Graph Folder**

This folder contains a script that I have used to create a pangenome for the dog using the canfam3 reference genome and the VCF with variants from the 722g project. *Plassais, J. et al, Whole genome sequencing of canids reveals genomic regions under selection and variants influencing morphology, Nature Communications, Volume 10, Article number: 1489 (2019) [Nature ].*

____*Then we have a 3 scripts used in the analysis of the NC enrichment assays analysis.*____

__________________________ **getStats4NC.sh**

This script allows top retrieve a table with stat data (eg. nb of reads, nb opf SNPs) 
at each step of the pipeline.

__________________________ **raw2vcf.sh**

This script contains all steps to process the raw data untill the VCF creation for the variant analysis.

__________________________ **ddup2vcf_novo.sh**

This cript uses novalign aligner to align the sequences against a modified reference genome and then creates the vcfs.
The Novoaligner uses IUPAC letters inside the reference genome to acount for polymorphism.

