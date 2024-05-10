rm(list=ls())
library(data.table)
#blast_results1<-read.table("results_high_v2.txt")
blast_results1 <- fread("blast_results.txt", showProgress = FALSE)

names<-c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")

colnames(blast_results1)<- names

## Applying filters

library(dplyr)
## Filter 1: % of identity and fragment length

probes1<- blast_results1 %>% filter(pident > 90) %>% filter(length > 110)

## Filter 2: Filter out probes that binds to other parts of the genome with more than 90% match

# Return names which have only a single row of data
probes2<-probes1 %>% 
  group_by(qseqid) %>% 
  filter(n()==1)

#Filter 3: Matching pattern to filter probes with specific match to target region

library(stringr)
##chromosome match test

original_chr<- str_match(probes2$qseqid, "_S\\s*(.*?)\\_")
probes2$original_chr<- as.numeric(original_chr[,2])
probes2$target_chr<-as.numeric(substr(probes2$sseqid, 4, 5))
probes2$chr_match<- probes2$original_chr==probes2$target_chr

## filter matching results

probes3<-filter(probes2, chr_match == "TRUE")

##position test
##extract original SNP position
position<-str_match(probes3$qseqid, "S[0-7]{2}_\\s*(.*?)\\|")
position2<-as.numeric(position[,2])

#logical test
probes3$position<- position2
probes3$test<-position2 - probes2$sstart

probes3$pos_match<- probes3$test < 500 ##assuming the maximum distance base on the design

##filter position
probes4<-filter(probes3, pos_match == "TRUE")


#drop the intial dataframe to free up memory

rm(blast_results1, probes1, probes2, probes3)
gc()

##################
### Thermodinamics
##################

#prepare data set adding sequence to blast results filtered

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Biostrings")


library("Biostrings")

dna1<- readDNAStringSet("./windows_final_seq_250bp.fasta", format="fasta",
                        nrec=-1L, skip=0L, seek.first.rec=FALSE,
                        use.names=TRUE, with.qualities=FALSE)
dna1[1]

seq_name = names(dna1)
sequence = paste(dna1)
df <- data.frame(seq_name, sequence)
colnames(df)<-c('qseqid','sequence')

probes5 <- (merge(df, probes4, by = 'qseqid'))

#devtools::install_github("jensenlab/primer3")

library("primer3")

##Tm
#iterate on the column (apply function to each row in a column)

probes5$tm<-round(as.numeric(sapply(probes5[,2],calculate_tm)),1)

# GC content 

gc_fun <- function(x){ 
  num_g <- str_count(x, "G")
  num_c <- str_count(x, "C")
  return ((num_g + num_c) / str_length(x) * 100 )} 

probes5$gc<-round(as.numeric(sapply(probes5[,2],gc_fun)),1)

##checking distributions
par(mfrow=c(1,2))
hist(probes5$gc, main = "GC content")
hist(probes5$tm, main = "Tm")

##Apply filters
str(probes5)
probes6<- probes5 %>% filter(gc > 35, gc < 45) %>% filter(tm > 65, tm < 75)

par(mfrow=c(1,2))
hist(probes6$gc, main = "GC content")
hist(probes6$tm, main = "Tm")

#####################
#secondary structure

##Note that the maximum length of ``seq`` is 60 bp. This is a cap suggested 
#by the Primer3 team as the longest reasonable sequence length for which a 
#two-state NN model produces reliable results.

# We can still analyze the probes by creating fragments up to 60bp.

#frag1 = first 60bp
#frag2 = last 60bp
#frag3 = first 30bp + last 30 bp
#frag4 = middle 60bp

probes6$frag1<-substr(probes6$sequence, 1, 60)
probes6$frag2<-substr(probes6$sequence, 61, 120)
probes6$frag3<-paste0(substr(probes6$sequence, 1, 30), substr(probes6$sequence, 91, 120))
probes6$frag4<-substr(probes6$sequence, 31, 90)


##hairpin Tm calculation
probes6_hairpin_frag1<-lapply(probes6$frag1,calculate_hairpin)
probes6_hairpin_frag2<-lapply(probes6$frag2,calculate_hairpin)
probes6_hairpin_frag3<-lapply(probes6$frag3,calculate_hairpin)
probes6_hairpin_frag4<-lapply(probes6$frag4,calculate_hairpin)

probes6$hairpin_temp_frag1<-sapply(probes6_hairpin_frag1,"[[",2)
probes6$hairpin_temp_frag2<-sapply(probes6_hairpin_frag2,"[[",2)
probes6$hairpin_temp_frag3<-sapply(probes6_hairpin_frag3,"[[",2)
probes6$hairpin_temp_frag4<-sapply(probes6_hairpin_frag4,"[[",2)

#probes6_hairpin[1]
#structure_found<-sapply(probes6_hairpin,"[[",1)
#ds<-sapply(probes6_hairpin,"[[",3)
#dh<-sapply(probes6_hairpin,"[[",4)
#dg<-sapply(probes6_hairpin,"[[",5)
#align_end_1<-sapply(probes6_hairpin,"[[",6)
#align_end_2<-sapply(probes6_hairpin,"[[",7)

probe6_homodimer_frag1<-lapply(probes6$frag1,calculate_homodimer)
probe6_homodimer_frag2<-lapply(probes6$frag2,calculate_homodimer)
probe6_homodimer_frag3<-lapply(probes6$frag3,calculate_homodimer)
probe6_homodimer_frag4<-lapply(probes6$frag4,calculate_homodimer)

probes6$homodimer_temp_frag1<-sapply(probes6_homodimer_frag1,"[[",2)
probes6$homodimer_temp_frag2<-sapply(probes6_homodimer_frag2,"[[",2)
probes6$homodimer_temp_frag3<-sapply(probes6_homodimer_frag3,"[[",2)
probes6$homodimer_temp_frag4<-sapply(probes6_homodimer_frag4,"[[",2)


## filter for hairpin Tm

probes7<- probes6 %>% filter(hairpin_temp_frag1 < 50) %>% 
  filter(hairpin_temp_frag2 < 50) %>%
  filter(hairpin_temp_frag3 < 50) %>%
  filter(hairpin_temp_frag4 < 50)


## filter for homodimer Tm

probes8<- probes7 %>% filter(homodimer_temp_frag1 < 50) %>% 
  filter(homodimer_temp_frag2 < 50) %>%
  filter(homodimer_temp_frag3 < 50) %>%
  filter(homodimer_temp_frag4 < 50)

# Thermodynamics calculations. The most important is the Structure Tm. If the Tm is lower than the reaction temperature,
# there will be no issues with the probe. 

# Other parameters: structure found, ds - change in entropy, dh - change in entalpy, dg - Gibbs free energy


##Filter to keep only one probe per site

probes8$probe_sites<-data.frame(str_match(probes8[,1], "S\\s*(.*?)\\|"))[,2]
probes8$probe_type<-str_sub(data.frame(str_match(probes8[,1], "\\|[A-Z]+\\|"))[,1],2,-2)

probes9<- probes8 %>% distinct(position, .keep_all = TRUE)


## Filter for one probe per gene

probes9$gene<-str_match(probes9[,1], "[a-z]{6}[0-9]{2}[a-z][0-9]{5}")[,1]

probes10<- probes9 %>% distinct(gene, .keep_all = TRUE)

write.table(probes10, "probes_v1.txt")

###############################
######Visualize probe positions
##plot Positions
str(probes10)

## making the plot

# Libraries
library(ggplot2)


##bins every 5M
bins<-seq(from=1, to=0.8e10, by=0.5e8)

##adding the bins to dataframe
probes10$bins<-findInterval(probes10$position, bins)

# extracting the info needed for plots
df2<-probes10[,c(14,17,23,36)]
df2$bins2<-paste0(df2$original_chr,"-",df2$bins)
df3<-data.frame(table(df2$bins2))
colnames(df3)<-c("bins2","freq")
df4 <- merge(df2,df3,by="bins2")
##keep unique bins
df5<- df4 %>% distinct(bins2, .keep_all = TRUE)

ggplot(probes10,aes(x=position, y= original_chr, colour = probe_type)) +
geom_point() +
scale_y_continuous(breaks = seq(min(probes10$original_chr), max(probes10$original_chr), by = 1)) +
theme_bw() +
theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
annotate("segment", x = 1, xend = c(1.84e9,1.39e9,1.1e9,1.1e9,0.9e9,0.75e9,0.76e9), y = 1:7, yend = 1:7,
         colour = "black") +
    theme(
    legend.position = c(0.92, 0.92),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )


ggplot(df5,aes(x=position, y= original_chr)) + 
  geom_count(aes(size=freq)) + 
  scale_y_continuous(breaks = seq(min(probes10$original_chr), max(probes10$original_chr), by = 1)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  annotate("segment", x = 1, xend = c(1.84e9,1.39e9,1.1e9,1.1e9,0.9e9,0.75e9,0.76e9), y = 1:7, yend = 1:7,
           colour = "black") +
  annotate("text", x = 1.7e9, y = 4.5, label = "Bin= 5M") +
 theme(
  legend.position = c(0.92, 0.92),
  legend.justification = c("right", "top"),
  legend.box.just = "right",
  legend.margin = margin(6, 6, 6, 6)
)

###########################
##Formating the final files
###########################

##ordering by probe type and selecting the high and moderate impact
#probes10_srt<-probes10[order(probes10$probe_type),,drop=FALSE] %>% filter(probe_type != "LOW")

#only reorder
probes10_srt<-probes10[order(probes10$probe_type),,drop=FALSE]

##Applying a more stringent filter for GC and Tm

probes10_srt_fltr<- probes10_srt %>% filter(gc > 40, gc < 45) %>% filter(tm > 69, tm < 72.5)

#selecting 1000 probes

probes10_1k<- probes10_srt_fltr[1:1000,]

##1000 probes summary

probes10_1k_sum<- probes10_1k[,c(1,2,5,3,10,11,22,23,20,21)]

write.table(probes10_1k_sum, "final_probes_1k_sum.txt")

##visualizing the probes again

#tm and gc
par(mfrow=c(1,2))
hist(probes10_1k$gc, main = "GC content", breaks = 2)
hist(probes10_1k$tm, main = "Tm", breaks = 2)

# extracting the info needed for plots
df2<-probes10_1k[,c(14,17,23,36)]
df2$bins2<-paste0(df2$original_chr,"-",df2$bins)
df3<-data.frame(table(df2$bins2))
colnames(df3)<-c("bins2","freq")
df4 <- merge(df2,df3,by="bins2")
##keep unique bins
df5<- df4 %>% distinct(bins2, .keep_all = TRUE)

ggplot(df5,aes(x=position, y= original_chr)) + 
  geom_count(aes(size=freq)) + 
  scale_y_continuous(breaks = seq(min(probes10$original_chr), max(probes10$original_chr), by = 1)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  annotate("segment", x = 1, xend = c(1.84e9,1.39e9,1.1e9,1.1e9,0.9e9,0.75e9,0.76e9), y = 1:7, yend = 1:7,
           colour = "black") +
  annotate("text", x = 1.7e9, y = 4.5, label = "Bin= 5M") +
  theme(
    legend.position = c(0.92, 0.92),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )

final_probes<- cbind(paste0(probes10_1k$qseqid,"|Tm=", probes10_1k$tm,"|GC=", probes10_1k$gc), probes10_1k$sequence)

##changing table to fasta format
#add ">" to headers
final_probes[,1] <- paste0(">",final_probes[,1])

#bind rows of headers ans seqs
probes_fasta <- c(rbind(final_probes[,1], final_probes[,2]))
probes_fasta[1:10]

write.table(probes_fasta, "final_probes_1k.txt", row.names=FALSE,sep="\t", quote = FALSE, col.names = FALSE)

