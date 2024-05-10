#!/bin/bash
#SBATCH --job-name=seq_xtr	# Job name
#SBATCH --nodes=1	        # N of nodes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem="32G"			# Memory per node; by default using M as unit
#SBATCH --time=24:00:00     # Time limit hrs:min:sec or days-hours:minutes:seconds
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
VCF=/cluster/home/rsouza/enrichment/snpeff_files/vcfs/out_all_S99_called_srt_fltr_2M_ann_v2.vcf.gz
#Annotations
ANN=/cluster/home/rsouza/enrichment/snpeff_files/vcfs/sintegrifolium.ann.tsv


# 1- extract the information from the annotated vcf

#selecting targets

cat ${VCF} | grep -E 'HIGH|MODERATE|LOW|MODIFIER' | sed 's/|/\t/g' | awk -v OFS='\t' '{ print $14, $1, $2, $2-250, $2+250, $3, "REF="$4,"ALT="$5, $9, $10}' | sed 's/CHR_START-//g' |  sed 's/-CHR_END//g' | sed 's/\.t[1-9]//g' > targets.txt
#ignoring intergenic candidates
#cat ${VCF} | grep -v 'intergenic' | grep -E 'HIGH|MODERATE|LOW|MODIFIER' | sed 's/|/\t/g' | awk -v OFS='\t' '{ print $14, $1, $2, $2-250, $2+250, $3, "REF="$4,"ALT="$5, $9, $10}' | sed 's/CHR_START-//g' |  sed 's/-CHR_END//g' | sed 's/\.t[1-9]//g' > targets.txt
#sort list
cat targets.txt | sort -nk 1 > targets_srt.txt

#simplify the gene file
cat ${ANN} | sed 's/\.t[1-9]//g'  | awk -v OFS='\t' '{ print $1, $6, $7, $8, $9, $10}' > sintegrifolium.ann2.tsv

#simplify the gene list

cat targets_srt.txt | cut -f1 > gene_list.txt

# filter the annotation file to get the meaninful annotations. Removing 'consensus disorder'

cat sintegrifolium.ann2.tsv |  grep -v 'consensus' > sintegrifolium.ann3.tsv

# Get the genes from annotation [The -f option instructs grep to read the patterns to look for from file2.txt]

grep -f gene_list.txt sintegrifolium.ann3.tsv > gene_ann.txt

#sort
cat gene_ann.txt | sort -nk 1 > gene_ann_srt.txt

#3- join snp data with annotation file

join targets_srt.txt gene_ann_srt.txt > gene_targets.txt

# remove missing annotations and duplicate genes

cat gene_targets.txt | grep -v ' \- ' >  gene_targets2.txt

awk '!a[$1]++' gene_targets2.txt > gene_targets3.txt

#4- prepare the final bed file (maybe use pipe instead of underscore)

cat gene_targets3.txt  | awk -v OFS='\t' '{ print $2, $4, $5, $6"|"$7"|"$8"|"$10"|"$1"|"$11"-"$12"-"$13"-"$14"-"$15}' > final_sites.bed

#5- extract the sequences 250 up and down with bedtools (this step can be also done on R)

bedtools getfasta -fi ${REF} -bed final_sites.bed -name > final_seq_250bp.fasta

#6- delete intermediate files

rm gene*
rm targets*
rm sintegrifolium.ann2.tsv
rm sintegrifolium.ann3.tsv

#7- run sliding window with 120bp and 10 bp window
## using fasta_windows_v1.1.sh from https://github.com/kdillmcfarland/sliding_windows

bash fasta_windows_v1.1.sh  final_seq_250bp.fasta 120 10

#8-run blast on the result
# blast needs to be in the PATH
blastn -db sintegrifolium -query windows_final_seq_250bp.fasta -out blast_results.txt -outfmt 6