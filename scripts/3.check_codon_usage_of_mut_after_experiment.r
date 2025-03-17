library(SynMut)
library(Biostrings)
library(tidyverse)

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

seq_ref_orf_orf1ab <- get_orf_seq(list_genes_orf1ab, seq=seq_ref)
seq_ref_orf_spike <- subseq(seq_ref, df_orf$start[df_orf$sequence=="S"], df_orf$stop[df_orf$sequence=="S"])

df_compare <- tibble(sample=names(seqs_mut), 
	codon_dist_orf1ab=codon_dist(seq_mut_orf_orf1ab, seq_ref_orf_orf1ab),
	codon_dist_spike=codon_dist(seq_mut_orf_spike, seq_ref_orf_spike),
	dinu_dist_orf1ab=dinu_dist(seq_mut_orf_orf1ab, seq_ref_orf_orf1ab),
	dinu_dist_spike=dinu_dist(seq_mut_orf_spike, seq_ref_orf_spike)
)
writexl::write_xlsx(df_compare, "results/df_mutant_compare_real_mutants.xlsx")

### 2. visualize the codon usage difference
source("scripts/helper/CA_analysis.r")
CA_analysis(seq_mut=seq_mut_orf_orf1ab[grepl("^[ABCDE]_[abc]$", names(seq_mut_orf_orf1ab))], gene="orf1ab", output_suffix="real_mutants", mutant_pattern="^[ABCDEF]_[abc]$", resecured_list = c("F_a", "F_b", "F_c", "D_c"), color_by_design = TRUE)
CA_analysis(seq_mut=seq_mut_orf_spike[grepl("^[EF]_[abc]$", names(seq_mut_orf_spike))], gene="spike", output_suffix="real_mutants", mutant_pattern="^[ABCDEF]_[abc]$", resecured_list = c("F_a", "F_b", "F_c", "D_c"), color_by_design = TRUE)
