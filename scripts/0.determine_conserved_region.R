# Aim: identified highly conserved sequences that are shared by all known huaman and animal sarbecoviruses
# Input: aligned amino acid sequences of different genes
# Output: hihgly conserved regions of different genes in table format

library(tidyverse)

# 1. Generate lists of Sarbecoviruses
df_seqs_ve <- readxl::read_excel("../data/data_seqs_ve.xlsx")
df_seqs_nat <- read_csv("../data/Lytras-etal_nCoV_origins-suppTableS1.csv")

colnames(df_seqs_nat)[1:5] <- names(df_seqs_ve)[c(1,3,4,2,5)]
df_plus <- df_seqs_ve[!df_seqs_ve$`Accession number` %in% df_seqs_nat$`Accession number`,]
df_plus <- df_plus %>% filter(`Accession number` != "NC_045512")

df_seqs <- bind_rows(df_seqs_nat, df_plus)
df_seqs <- df_seqs %>% arrange(`Virus name`)
writeLines(df_seqs$`Accession number`, "../data/fasta/ac_id.txt")
df_seqs$`Host species` <- gsub(" ", "_", df_seqs$`Host species`)
df_seqs <- df_seqs[,c(1,2,4)]
df_seqs <- left_join(df_seqs, df_seqs_ve[,c(2, 4, 5)])

df_seqs <- left_join(df_seqs, df_seqs_nat[,c(3, 4:6)], "Accession number")

write_csv(df_seqs, "../data/data_seqs_all.csv") # manual fill to xlsx

# 2. Alignment of the sequences
library(Biostrings)

seq_gisaid <- readDNAStringSet("../data/fasta/gisaid_hcov-19_2021_09_16_08.fasta")
seq_genbank <- readDNAStringSet("../data/fasta/sequence_nt.fasta")

