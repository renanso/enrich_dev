#!/bin/bash
#SBATCH --job-name=seq_xtr	# Job name
#SBATCH --nodes=1	        # N of nodes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem="8G"			# Memory per node; by default using M as unit
#SBATCH --time=1:00:00     # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --export=ALL        #export user’s explicit environment variables to compute node
#SBATCH --output=%x_%j.out		# Standard output log
#SBATCH --error=%x_%j.err		# Standard error log
#SBATCH --partition=normal

# Load module
ml cluster/bedtools/2.28.0

#directories and inputs
#reference genome
REF=/cluster/home/rsouza/enrichment/Silphium_integrifolium_var_Bad_Astra.mainGenome2.fasta
#Input VCF
VCF=/cluster/home/rsouza/enrichment/snpeff_files/vcfs/out_all_S99_called_srt_fltr_2M_ann.vcf.gz
#Annotations
ANN=/cluster/home/rsouza/enrichment/snpeff_files/vcfs/sintegrifolium.ann2.tsv


# 1- extract the information from the annotated vcf

#high impact
cat ${VCF} | grep 'HIGH' | sed 's/|/\t/g' |awk -v OFS='\t' '{ print $11, $1, $2, $2-250, $2+250, $3, "REF="$4,"ALT="$5, $9, $10}' > high_impact.txt

#moderate impact
cat ${VCF} | grep 'MODERATE' | sed 's/|/\t/g' |awk -v OFS='\t' '{ print $11, $1, $2, $2-250, $2+250, $3, "REF="$4,"ALT="$5, $9, $10}' > moderate_impact.txt

#low impact
cat ${VCF} | grep 'LOW' | sed 's/|/\t/g' |awk -v OFS='\t' '{ print $11, $1, $2, $2-250, $2+250, $3, "REF="$4,"ALT="$5, $9, $10}' > low_impact.txt

#upstream

#downstream


#2 get the annotations 
# gene lists
# in case there is need to fix some the genes names
sed -r 's/(\B[0-9]{2}g)([0-9]{1}\b)/\10\2/g' low_impact.txt > test2.txt
sed -r 's/(\B[0-9]{2}g)([0-9]{2}\b)/\10\2/g' test2.txt > test3.txt
sed -r 's/(\B[0-9]{2}g)([0-9]{3}\b)/\10\2/g' test3.txt > test4.txt
sed -r 's/(\B[0-9]{2}g)([0-9]{4}\b)/\10\2/g' test4.txt > high_impact2.txt

#sort list
cat high_impact2.txt | sort -nk 1 > high_impact3.txt

#simplify the gene file
cat ${ANN} | sed 's/\.t[0-20]//g'  | awk -v OFS='\t' '{ print $1, $6, $7, $8, $9, $10}' > sintegrifolium.ann3.tsv

#simplify the gene list

cat high_impact3.txt | cut -f1 > high_impact3_list.txt

# Get the high impact genes from annotation [The -f option instructs grep to read the patterns to look for from file2.txt]

grep -f high_impact3_list.txt sintegrifolium.ann3.tsv > high_impact_SNP2M_gene_ann_04-29.txt

#sort
cat high_impact_SNP2M_gene_ann_04-29.txt | sort -nk 1 > high_impact_SNP2M_gene_ann_04-29b.txt

#3- join snp data with annotation file

join high_impact3.txt high_impact_SNP2M_gene_ann_04-29b.txt > test_high_impact_genes.txt


# remove missing annotations and duplicate genes

cat test_high_impact_genes.txt | grep -v ' \- ' >  test_high_impact_genes2.txt

awk '!a[$1]++' test_high_impact_genes2.txt > test_high_impact_genes3.txt

#4- prepare the final bed file (maybe use pipe instead of underscore)

cat test_high_impact_genes3.txt  | awk -v OFS='\t' '{ print $2, $4, $5, $6"|"$7"|"$8"|"$10"|"$1"|"$11"-"$12"-"$13"-"$14"-"$15}' > SNP2M_gene_info_low.bed

#5- extract the sequences 250 up and down with bedtools (this step can be also done on R)

bedtools getfasta -fi ${REF} -bed SNP2M_gene_info_low.bed -name > SNP2M_250bp_low.fasta

#6- delete intermediate files

rm test_high_impact_genes*
rm high_impact*
rm moderate_impact*
rm low_impact*
rm test*
rm sintegrifolium.ann3.tsv


## next step is the sliding window with 120bp with 1 bp sliding window
## using fasta_windows_v1.1.sh from https://github.com/kdillmcfarland/sliding_windows
## bash fasta_windows_v1.1.sh SNP2M_250bp_high.fasta 120 1