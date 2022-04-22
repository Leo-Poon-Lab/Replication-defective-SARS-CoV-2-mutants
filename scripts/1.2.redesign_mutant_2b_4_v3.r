# 2022-04-19
# the partial segment of mutant 2b_4 cannot be synthesized (positions 17335 to 21892), therefore we reshuffled some codons in this segment ( 17341 to 21555 in ORF1b), to test if the new segment can be sythesized.

library(Biostrings)
library(tidyverse)
library(SynMut)

seq_mut_2b_4_full <- readDNAStringSet("../results/Final design for Syntethesis/Mutant_2b_4_seq.fasta")
region <- c(17662, 18930) # please refer to "../data/Mutant_2b_4_toxic.fasta"
seq_mut_2b_4 <- subseq(seq_mut_2b_4_full, region[1], region[2])

df_mutable <- read_csv("../results/df_mutatable_region.csv")
df_mutable_2b_4 <- df_mutable[region[1]:region[2],]
vec_mutatble <-  sapply(seq_len(nrow(df_mutable_2b_4)/3), function(n) {
	all(df_mutable_2b_4$modifiable[(3*n-2):(3*n)])
})	

rd_mut_2b_4 <- input_seq(seq_mut_2b_4, data.frame(vec_mutatble))
set.seed(2022)
rd_mut_2b_4_new <- codon_random(rd_mut_2b_4, keep = TRUE)
seq_mut_2b_4_new <- get_dna(rd_mut_2b_4_new)
seq_mut_2b_4_full_new <- seq_mut_2b_4_full
subseq(seq_mut_2b_4_full_new, region[1], region[2]) <- seq_mut_2b_4_new
names(seq_mut_2b_4_full_new) <- gsub("Mutant_2b", "Mutant_2b_v3", names(seq_mut_2b_4_full_new))
writeXStringSet(seq_mut_2b_4_full_new, "../results/Final design for Syntethesis/Mutant_2b_4_seq_v3_20220422.fasta")

source("./helper/compare_seqs_to_table.r")
df_cmp <- compare_seqs(seq_int=seq_mut_2b_4_full_new, seq_ref=seq_mut_2b_4_full)
write_csv(df_cmp, "../results/Final design for Syntethesis/Mutant_2b_4_info_v3_20220422.csv")
