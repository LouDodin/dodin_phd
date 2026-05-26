### DATA ###
Lucien Weymiens' experiment
Coevolution T0 host DNA sequences
Obtained with Short-read sequencing technology (double direction)

4 files (from two different sequencing)
Mc_CTAGGTTG-GCTCGAAT-A8DA1FO2H_L002_R1.fastq.gz
Mc_CTAGGTTG-GCTCGAAT-A8DA1FO2H_L002_R2.fastq.gz
Mc_CTAGGTTG-GCTCGAAT-A9NF1TB3R_L002_R1.fastq.gz
Mc_CTAGGTTG-GCTCGAAT-A9NF1TB3R_L002_R2.fastq.gz


#10/12/2025

Put files on Genotoul

Concatenate files
---> Mc_cat_R1.fastq.gz and Mc_cat_R2.fastq.gz

Find reference genome of Micromonas on ncbi -> genome -> Micromonas Commoda
---> GCA_000090985.2_ASM9098v2_genomic.fna



#10/12/2025
@Map on reference genome

1) bwa : Index reference genome

2) bwa : Map reads
map_reads.sh
=====
#!/bin/bash
#SBATCH -p workq
#SBATCH -t 01-00:00:00
#SBATCH --mem=8G

module load bioinfo/bwa-mem2/2.2.1

bwa-mem2 mem -t 8 -R '@RG\tID:Mc_cat\tSM:Mc_cat\tLB:lib1\tPL:ILLUMINA\tPU:unit1' GCA_000090985.2_ASM9098v2_genomic.fna Mc_cat_R1.fastq.gz Mc_cat_R2.fastq.gz > Mc_cat_mem2.sam
=====
---> Mc_cat_mem2.sam

3) samtools : Transform .sam into .bam
sam_to_bam.sh
=====
#!/bin/bash
#SBATCH -p workq
#SBATCH -t 01-00:00:00
#SBATCH --mem=8G

module load bioinfo/samtools/1.21

samtools view Mc_cat_mem2.sam -b -o Mc_cat_mem2.bam
=====
---> Mc_cat_mem2.bam

Check if .bam is good:
% samtools view -H Mc_cat_mem2.bam | grep '^@RG'
---> @RG	ID:Mc_cat	SM:Mc_cat	LB:lib1	PL:ILLUMINA	PU:unit1
It's good !

4) samtools : Sort : order reads by position on genome
% samtools sort Mc_cat_mem2.bam -o Mc_cat_mem2_sorted.bam
---> Mc_cat_mem2_sorted.bam

5) Index genome
% samtools index Mc_cat_mem2_sorted.bam

6) 
% samtools depth Mc_cat_mem2_sorted.bam -a -o Mc_cat_mem2_sorted_deptha

% awk '$3>0' <file> > <output_file>



#17/12/2025
@Obtain coverage per window

1) Download multifasta_lengths.pl (file written by Sheree Yau)

2) Compute it
% perl multifasta_lengths.pl <.fasta>
---> .bed file   | sequencing | length |

3) 
% bedtools coverage -a <file.bam> -b <file_windows1000> > <output_file>
---> Coverage per window


#17/12/2025
@Plot coverage per window
Done in Python : host_coevol_mapping_t0/plot_mapping_t0.py

#09/01/2026
@Plot coverage per window
Done using IGV

I downloaded IGV
Reference genome = GCA_000090985.2_ASM9098v2_genomic.fna
Track = Mc_cat_mem2_sorted.bam and Mc_cat_mem2_sorted.bam.bai

#09/01/2026
@Plot coverage per window
I don't manage to do it in R


#16/01/2026
@SNP calling from .bam
@Using GATK

HaplotypeCaller
https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller

1) Generate the missing reference indexes
% samtools faidx GCA_000090985.2_ASM9098v2_genomic.fna
---> GCA_000090985.2_ASM9098v2_genomic.fna.fai
% gatk CreateSequenceDictionary -R GCA_000090985.2_ASM9098v2_genomic.fna
---> GCA_000090985.2_ASM9098v2_genomic.dict

2) Search SNPs
variant.sh
=====
#!/bin/bash
#SBATCH -p workq
#SBATCH -t 01-00:00:00
#SBATCH --mem=8G

module load devel/python/Python-3.11.1
module load devel/java/17.0.6
module load statistics/R/4.3.1
module load bioinfo/GATK/4.6.2.0

gatk --java-options "-Xmx4g" HaplotypeCaller -R GCA_000090985.2_ASM9098v2_genomic.fna -I Mc_cat_mem2_sorted.bam -O Mc_cat_mem2_sorted_variants.vcf
=====

#21/01/2026
3) Analysis
% bcftools stats Mc_cat_mem2_sorted_variants.vcf > Mc_cat_mem2_sorted_variants_stats.txt
---> Mc_cat_mem2_sorted_variants_stats.txt

Send it to my Macbook
% get Mc_cat_mem2_sorted_variants_stats.txt

Try to plot it but don't manage to do it :
% brew install bcftools
% python3 -m pip install matplotlib
% plot-vcfstats -p Mc_cat_mem2_sorted_variants_plots Mc_cat_mem2_sorted_variants_stats.txt
---> Folder Mc_cat_mem2_sorted_variants_plots/
	index.html
	summary.txt
	.png
	.dat

#23/01/2026
Obtain number of SNPs per chromosome

% bcftools view -v snps Mc_cat_mem2_sorted_variants.vcf > | bcftools query -f '%CHROM\n' > | sort | uniq -c
  58472 CP001323.1
  57311 CP001324.1
  50778 CP001325.1
  50074 CP001326.1
  46320 CP001327.1
  45651 CP001328.1
  43436 CP001329.1
  39552 CP001330.1
  35607 CP001331.1
  28804 CP001332.1
  29054 CP001333.1
  22524 CP001334.1
   2771 CP001335.1
  35597 CP001574.1
  39275 CP001575.1
  40304 CP001576.1
  38044 CP001577.1
   1344 FJ858267.1
    584 FJ859351.1

I reordonnate manually :
Ch : nb_SNP lg_ch pourc_SNP
Ch 1 : 35597 2053059 0.017339
Ch 2 : 58472 1914325 0.030544
Ch 3 : 57311 1759951 0.032564
Ch 4 : 50778 1584431 0.032048
Ch 5 : 50074 1518631 0.032973
Ch 6 : 46320 1431126 0.032366
Ch 7 : 45651 1394110 0.032746
Ch 8 : 39275 1276107 0.030777
Ch 9 : 43436 1260462 0.034460
Ch 10 : 40304 1160640 0.034726
Ch 11 : 39552 1145873 0.034517
Ch 12 : 38044 1084119 0.035092
Ch 13 : 35607 1011177 0.035213
Ch 14 : 28804 832468 0.034601
Ch 15 : 29054 739136 0.039308
Ch 16 : 22524 608929 0.036990
Ch 17 : 2771 214782	0.012901
Mitochondrie : 584 47425 0.012314
Chloroplaste : 1344 72585 0.018516