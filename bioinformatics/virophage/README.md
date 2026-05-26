# 04/03/26
cat *.fasta > all_virophages.fasta
bwa-mem2 mem all_virophages.fasta Mc_cat_R1.fastq.gz Mc_cat_R2.fastq.gz | samtools view -bS - | samtools sort -o all_virophages_mapped.sorted.bam
samtools view -b -F 4 all_virophages_mapped.sorted.bam | samtools sort -o all_virophages_mapped_only.sorted.bam
samtools index all_virophages_mapped_only.sorted.bam
samtools depth all_virophages_mapped_only.sorted.bam > all_virophages_depth.txt

awk '{cov[$1]+=$3; len[$1]++} END {for (c in cov) print c, cov[c]/len[c]}' all_virophages_depth.txt

Résultats :
virophage couverture_moyenne_par_contig
MH920636.1 2.06977
NC_015230.1 545.6
MH919296.1 2






# 05/03/26
## Compresser et indexer le vcf
bgzip Mc_cat_mem2_sorted_variants.vcf
tabix -p vcf Mc_cat_mem2_sorted_variants.vcf.gz

## Générer le concensus
gatk 
bcftools consensus -f GCA_000090985.2_ASM9098v2_genomic.fna Mc_cat_mem2_sorted_variants.vcf.gz
> consensus.fasta

## Générer le fichier qui contient la database de protéines qu'on cherche (ici virophage capsid / major capsid protein virus)
Recherche des protéines sur UniProt
makeblastdb -in capsid_proteins.fasta -dbtype prot -out capsid_db

## tblastx
blastx -query consensus.fasta -db capsid_db -out results.txt -outfmt 6







# 09/03/26
seqtk seq -a Mc_cat_R1.fastq.gz  Mc_cat_R2.fastq.gz > Mc_cat.fa

makeblastdb -in virophage.faa -dbtype prot -out virophage_db

blastx -query Mc_cat.fa -db virophage_db -out blast_results.txt -evalue 1e-5 -outfmt 6 -num_threads 4

makeblastdb -in Mc_cat_top.fa -dbtype nucl -out database

tblastn -query test.faa -db database -task tblastn-fast -evalue 1e-5 -num_threads 8 -outfmt 6 -out mcp_vs_database -max_target_seqs 200000





makeblastdb -in Mc_cat.fa -dbtype nucl -out database_complete

tblastn -query proteins.faa -db database_complete -task tblastn-fast -num_threads 8 -outfmt 6 -out proteins_vs_database_complete -max_target_seqs 2000
#-evalue 1e-5

sort -k3,3nr proteins_vs_database_complete > proteins_vs_database_complete_sorted


seqkit grep -p "" Mc_cat.fa