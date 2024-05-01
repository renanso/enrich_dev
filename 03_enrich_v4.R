library("Biostrings")
DNAStringSet()

dna1<- readDNAStringSet("./windows_test.fasta", format="fasta",
                        nrec=-1L, skip=0L, seek.first.rec=FALSE,
                        use.names=TRUE, with.qualities=FALSE)
dna1[1]

#################
#####BLAST#######

#BiocManager::install(version = "devel")
#BiocManager::install("rBLAST")
library('rBLAST')

# check if makeblastdb is correctly installed. If not need to install BLAST in the system
Sys.which("makeblastdb")

# 1. make database
# Will build databases for each chromosome because of the genome size
#makeblastdb(file = "/Users/renansouza/Documents/sequences/Chr01.fa", db_name = "/Users/renansouza/Documents/sequences/db/small", dbtype = "nucl")
  
# 3. open database
db <- blast("/Users/renansouza/Documents/sequences/db/small")
db
  
## 4. perform search 
blast_results<-predict(db, dna1, BLAST_args = "-perc_identity 99")

## creating a matching pattern to filter probes with specific match
library(stringr)

##chromosome test

original_chr<- str_match(blast_results$qseqid, "_S\\s*(.*?)\\_")
blast_results$original_chr<- as.numeric(original_chr[,2])
blast_results$target_chr<-as.numeric(substr(blast_results$sseqid, 4, 5))
blast_results$chr_match<- blast_results$original_chr==blast_results$target_chr

## filter matching results
library(dplyr)

blast_results2<-filter(blast_results, chr_match == "TRUE")

##position test
##extract original SNP position
position<-str_match(blast_results2$qseqid, "S[0-1]{2}_\\s*(.*?)\\|")
position2<-as.numeric(position[,2])

#logical test
blast_results2$test<-position2 - blast_results2$sstart

blast_results2$pos_match<- blast_results2$test < 1000 ##assuming a gross error margin

##filter position
blast_results3<-filter(blast_results2, pos_match == "TRUE")

##penalize probes that match multiple locations

##################
### Thermodinamics
##################

#prepare data set adding sequence to blast results filtered
seq_name = names(dna1)
sequence = paste(dna1)
df <- data.frame(seq_name, sequence)
colnames(df)<-c('qseqid','sequence')

probes <- (merge(df, blast_results3, by = 'qseqid'))

#devtools::install_github("jensenlab/primer3")

library("primer3")

##Tm and hairpin
#iterate on the column (apply function to each row in a column)

probes$tm<-as.numeric(sapply(probes[,2],calculate_tm))

probe_hairpin<-lapply(probes[,2],calculate_hairpin)
probe_hairpin[1]

#calculate_hairpin(df[1,2])

probe_homodimer<-lapply(df[,2],calculate_homodimer)

# GC content (make this a function)
library(stringr)

gc_fun <- function(x){ 
  num_g <- str_count(x, "G")
  num_c <- str_count(x, "C")
  return ((num_g + num_c) / str_length(x) * 100 )} 

probes$gc<-as.numeric(sapply(probes[,2],gc_fun))

##add filter parameters
par(mfrow=c(1,2))
hist(probes$gc, main = "GC content")
hist(probes$tm, main = "Tm")

##Apply filters

probes2<- probes %>% filter(gc > 35, gc < 45) %>% filter(tm > 65, gc < 75)

par(mfrow=c(1,2))
hist(probes2$gc, main = "GC content")
hist(probes2$tm, main = "Tm")

##final list with probes overlapping by 60 bp (select based on the window position + SNP locus)

## plot location of probes in the genome. Density

## filter to 1 target each 7M.
