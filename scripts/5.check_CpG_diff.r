library(SynMut)
library(Biostrings)
library(tidyverse)
library(parallel)

# read data
load("scripts/helper/df_orf.rdata")
source("scripts/helper/make_regioned_data.r")

files_seq_mut_2a <- list.files("data/Mutant sequences in fasta format_updated/insertion of mutant_2a_os_intra_nonnsp sequences into WT SARS-2 backbo", "\\.fa$", full.names=T)
files_seq_mut_2b <- list.files("data/Mutant sequences in fasta format_updated/insertion of mutant_2b_4 sequences into WT SARS-2 backbone", "\\.fa$", full.names=T)
files_seq_mut_2c <- list.files("data/Mutant sequences in fasta format_updated/insertion of mutant_2c sequences into WT SARS-2 backbone", "\\.fa$", full.names=T)

name_seq_mut_2a <- gsub("2a.+", "a", basename(files_seq_mut_2a))
name_seq_mut_2b <- gsub("2b.+", "b", basename(files_seq_mut_2b))
name_seq_mut_2c <- gsub("2c.+", "c", basename(files_seq_mut_2c))

files_seq_mut_all <- c(files_seq_mut_2a, files_seq_mut_2b, files_seq_mut_2c)
name_seq_mut_all <- c(name_seq_mut_2a, name_seq_mut_2b, name_seq_mut_2c)

seqs_mut <- lapply(files_seq_mut_all, function(x) {
	readDNAStringSet(x)
})
seqs_mut <- do.call(c, seqs_mut)
names(seqs_mut) <- name_seq_mut_all

### 1. quantify the codon usage difference
list_genes_orf1ab <- df_orf$sequence[grepl("^nsp", df_orf$sequence)]
list_genes_spike <- "S"

seq_mut_orf_orf1ab <- lapply(seqs_mut, function(x) {
	get_orf_seq(list_genes_orf1ab, seq=x)
})
seq_mut_orf_orf1ab <- DNAStringSet(seq_mut_orf_orf1ab)
seq_mut_orf_spike <- subseq(seqs_mut, df_orf$start[df_orf$sequence=="S"], df_orf$stop[df_orf$sequence=="S"])

cb1_s_muatant_spike_seq <- DNAStringSet(
"atgtttgtttttcttgttttattgccactagtctctagtcagtgtgttaatcttacaaccagaactcaattaccccctgcatacactaattctttcacacgtggtgtttattaccctgacaaagttttcagatcctcagttttacattcaactcaggacttgttcttacctttcttttccaatgttacttggttccatgctatacatgtctctgggaccaatggtactaagaggtttgataaccctgtcctaccatttaatgatggtgtttattttgcttccactgagaagtctaacataataagaggctggatttttggtactactttaGATTCGAAAACACAATCACTTTTAATTGTGAATAATGCTACAAACGTCGTCATTAAGGTTTGCGAATTTCAATTTTGTAACGATCCTTTTTTGGGTGTTTATTACCATAAGAACAATAAGTCATGGATGGAATCTGAATTCAGAGTATATAGTTCTGCAAATAACTGTACCTTTGAATATGTATCTCAACCTTTCCTAATGGATCTGGAGGGAAAACAAGGCAACTTTAAGAATTTGCGTGAATTTGTCTTTAAGAACATTGATGGCTATTTTAAGATATATAGTAAACACACTCCAATTAACCTTGTTCGTGATCTGCCACAAGGATTTTCTGCTTTAGAGCCACTCGTTGATCTACCTATTGGTATTAATATTACAAGGTTCCAGACTCTTTTAGCCCTTCATCGTTCATATCTAACACCTGGTGATTCATCTTCTGGCTGGACCGCTGGTGCTGCTGCTTATTATGTAGGTTATCTCCAACCTAGAACTTTTCTACTTAAATATAATGAAAATGGGACCATAACTGATGCAGTAGACTGTGCATTAGATCCACTTTCTGAAACAAAATGCACTCTCAAATCTTTTACAGTGGAAAAAGGTATTTACCAAACTTCAAACTTCCGTGTACAACCAACTGAGTCTATTGTTAGGTTTCCTAATATTACAAACTTATGTCCTTTTGGTGAAGTTTTTAATGCTACAAGATTTGCATCTGTATATGCCTGGAATAGAAAAAGAATCTCAAACTGTGTTGCCGATTATTCAGTGCTTTATAACAGTGCTTCTTTCAGTACTTTTAAGTGTTACGGAGTTTCCCCTACTAAATTAAATGACCTTTGTTTTACAAATGTTTACGCTGATAGCTTTGTTATTAGAGGCGACGAGGTTAGACAGATAGCTCCTGGACAAACTGGCAAAATTGCAGACTATAATTACAAATTACCTGATGATTTTACTGGATGCGTCATCGCTTGGAACTCCAATAATTTAGATTCTAAGGTGGGTGGAAACTATAATTACTTGTACAGATTATTTCGTAAATCTAACCTTAAACCCTTTGAAAGGGACATTTCAACTGAAATTTATCAAGCCGGTTCTACACCATGTAACGGTGTAGAAGGTTTTAATTGCTACTTCCCTCTTCAATCATATGGCTTTCAACCCACTAATGGTGTTGGTTATCAACCTTACCGTGTTGTTGTTTTGAGTTTTGAATTACTACATGCTCCAGCAACTGTTTGTGGACCTAAGAAGTCTACAAACCTCGTTAAGAACAAGTGTGTTAACTTTAATTTCAATGGTCTTACAGGTACTGGAGTATTGACAGAAAGTAATAAGAAATTTTTACCTTTCCAGCAATTCGGAAGAGACATTGCAGACACAACAGATGCAGTTAGGGACCCACAAACTCTTGAAATACTTGATATAACCCCCTGTTCTTTTGGTGGTGTCAGTGTCATTACTCCTGGTACTAATACTTCTAATCAAGTCGCTGTTTTGTATCAAGACGTGAATTGCACAGAAGTTCCTGTCGCAATCCATGCAGATCAGTTAACTCCTACGTGGCGCGTTTACTCCACAGGTTCCAATGTCTTTCAAACTAGAGCTGGTTGCTTAATCGGTGCAGAACATGTCAATAATTCATATGAATGTGATATACCAATTGGTGCTGGTATTTGCGCGAGTTATCAAACACAAACAAATTCACCTCGGCGGGCACGGAGTGTTGCAAGTCAGTCTATAATCGCATATACTATGTCTTTAGGTGCTGAAAATTCTGTAGCTTATTCTAACAATTCTATAGCTATTCCTACTAATTTCACAATTTCAGTTACTACCGAGATTCTTCCAGTTTCAATGACAAAAACCTCCGTAGATTGTACTATGTATATATGCGGCGACTCAACCGAGTGTAGTAACCTCCTTTTGCAGTATGGATCCTTTTGTACACAGCTGAACCGTGCTCTTACAGGTATTGCAGTTGAACAAGATAAAAACACACAGGAAGTTTTTGCGCAGGTAAAACAAATATATAAAACACCTCCAATTAAAGACTTTGGTGGTTTTAATTTTTCACAGATATTGCCTGATCCATCCAAACCTAGTAAGAGGAGTTTCATTGAAGATCTATTGTTTAACAAGGTAACTTTGGCAGACGCCGGTTTTATAAAACAGTACGGCGATTGTTTGGGTGATATTGCTGCAAGGGATCTTATTTGTGCCCAGAAATTTAACGGTTTGACTGTGTTACCTCCATTACTTACCGACGAGATGATTGCTCAATACACATCTGCTCTCCTTGCTGGAACTATTACCTCGGGTTGGACTTTTGGTGCTGGGGCTGCACTTCAAATACCATTTGCTATGCAGATGGCTTACAGGTTTAATGGTATAGGCGTTACACAAAATGTACTTTATGAGAATCAAAAACTTATTGCAAATCAATTTAATTCAGCAATTGGCAAGATTCAAGATAGCCTTTCATCTACGGCTAGCGCACTTGGTAAGTTACAAGATGTTGTTAATCAAAACGCTCAAGCTTTAAATACTTTAGTCAAGCAACTCAGCTCCAATTTTGGTGCAATTTCTTCTGTTCTTAACGACATTCTCAGCAGACTCGATAAAGTTGAGGCTGAGGTCCAAATTGATAGATTGATCACTGGTAGACTTCAATCACTTCAAACATATGTTACACAACAATTGATCAGGGCAGCTGAAATAAGAGCTAGTGCAAACTTAGCCGCTACTAAAATGTCAGAATGTGTTTTAGGGCAATCCAAAAGAGTGGACTTTTGCGGAAAAGGTTATCATCTTATGTCTTTTCCACAATCAGCACCACATGGTGTTGTGTTCCTACACGTTACATACGTACCAGCCCAAGAGAAAAATTTTACTACAGCACCTGCTATCTGCCATGATGGTAAGGCTCACTTTCCAAGGGAAGGTGTTTTTGTCTCTAACGGCACACATTGGTTTGTCACTCAACGTAACTTTTATGAACCCCAGATTATCACTACTGATAATACTTTTGTTTCTGGCAATTGTGATGTGGTCATTGGCATTGTGAATAATACAGTTTATGATCCTTTGCAACCAGAACTTGATTCCTTTAAAGAAGAATTAGACAAATATTTCAAAAACCATACATCTCCAGACGTCGACTTAGGTGATATTAGTGGTATAAATGCTTCCGTCGTTAACATTCAGAAAGAAATCGACCGTCTTAATGAGGTGGCTAAAAATTTGAATGAATCTCTCATTGATTTACAGGAGCTTGGCAAATATGAACAATATATTAAGTGGCCTTGGTATATCTGGCTTGGATTTATTGCTGGTTTAATCGCAATCGTTATGGTTACAATTATGCTTTGTTGTATGACTAGTTGTTGTTCTTGTTTGAAAGGTTGCTGCTCTTGTGGTTCATGTTGTAAATTTGATGAAGATGATTCTGAGCCAGTCTTAAAGGGCGTCAAACTCCATTATACTTAA")

meta_data <- read_csv("data/meta_omsn_tidy_2020_02_20_pan.csv")
names(meta_data)[c(6, 9, 14, 20, 21)] <-
    c("Gene", "virus_species_ori", "Host_ori", "Host", "Virus Species")
cds_omsn <- readDNAStringSet("data/cleaned_cds_osmn_2020_02_20_pan.fasta")
meta_data$width <- width(cds_omsn)
###delete guangdong pangolin
cds_omsn <- cds_omsn[!grepl("pangolin/Guangdong/P2S/2019", meta_data$cds_name)]
meta_data <- meta_data[!grepl("pangolin/Guangdong/P2S/2019", meta_data$cds_name),]

wt_s_muatant_spike_seq <- cds_omsn[which(meta_data$id == "BetaCoV/Wuhan-Hu-1/2019|EPI_ISL_402125" & meta_data$Gene=="spike")]

get_du(wt_s_muatant_spike_seq)
get_du(cb1_s_muatant_spike_seq)

# lapply(1:10, function(x){
#   tmp <- codon_random(wt_s_muatant_spike_seq, keep = TRUE)
#   get_cu(tmp)==get_cu(wt_s_muatant_spike_seq)
# }) # These are all true

CG_shuffled <- mclapply(1:1000000, function(x){
  get_du(codon_random(wt_s_muatant_spike_seq, keep = TRUE))[7] # get CpG
}, mc.cores=16)

## plot the distribution of the CpG content for the shuffled sequences, and label the CpG content of the CB1 mutant, using ggplot
## Pre caculate the histogram data and use geom_col to plot, to avoid the calculation of histogram during plotting which takes a long time

df_CpG <- tibble(CpG=unlist(CG_shuffled))
median(df_CpG$CpG)
mean(df_CpG$CpG)

hist_data <- hist(df_CpG$CpG, breaks=seq(0, max(df_CpG$CpG)+1, by=1), plot=FALSE)
df_hist <- tibble(CpG=hist_data$mids, Frequency=hist_data$counts)

p <- ggplot(df_hist, aes(x=CpG, y=Frequency)) +
  geom_col(fill="lightblue", color="black") +
  geom_vline(xintercept = get_du(cb1_s_muatant_spike_seq)[7], color="red", linetype="dashed", size=1) +
  geom_label(aes(x=get_du(cb1_s_muatant_spike_seq)[7], 
                 y=max(Frequency), 
                 label=paste0("CB1 mutant CpG: ", round(get_du(cb1_s_muatant_spike_seq)[7], 4))),
             color="red", 
             vjust=1.5, 
             hjust=1.1,
             size=5) +
  ## also show the median and mean
  geom_vline(xintercept = median(df_CpG$CpG), color="blue", linetype="dotted", size=1) +
  geom_vline(xintercept = mean(df_CpG$CpG), color="green", linetype="dotted", size=1) +
  geom_label(aes(x=median(df_CpG$CpG), 
                 y=max(Frequency)/2, 
                 label=paste0("Median: ", round(median(df_CpG$CpG), 4))),
             color="blue", 
             vjust=-0.5, 
             hjust=-0.1,
             size=5) +
  geom_label(aes(x=mean(df_CpG$CpG), 
                 y=max(Frequency)/3, 
                 label=paste0("Mean: ", round(mean(df_CpG$CpG), 4))),
             color="green", 
             vjust=-0.5, 
             hjust=-0.1,
             size=5) +
  labs(title="Distribution of CpG content in synonymous codon reshuffled sequences",
       x="CpG content",
       y="Frequency") +
  theme_minimal()
ggsave("results/CpG_content_distribution_cb1_mutant_spike.pdf", plot=p, width=8, height=6)

