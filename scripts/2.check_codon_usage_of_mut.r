library(SynMut)
library(Biostrings)
library(tidyverse)

# read data
load("./helper/df_orf.rdata")
source("./helper/make_regioned_data.r")

# files_seq_mut <- list.files("../results/mutants/", "fasta", full.names=T)
files_seq_mut <- list.files("../results/Final design for Syntethesis/", "fasta", full.names=T)
name_seq_mut <- gsub(".fasta$", "", list.files("../results/Final design for Syntethesis/", "fasta"))
name_seq_mut <- gsub("_seq_v2", "_v2_seq", name_seq_mut)
name_seq_mut <- gsub("_seq_v3", "_v3_seq", name_seq_mut)
name_seq_mut <- gsub("_seq.*", "", name_seq_mut)
name_seq_mut <- gsub("_os.+", "", name_seq_mut)

seqs_mut <- lapply(files_seq_mut, function(x) {
	readDNAStringSet(x)
})
seqs_mut <- do.call(c, seqs_mut)
names(seqs_mut) <- stringr::str_to_sentence(name_seq_mut)

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
writexl::write_xlsx(df_compare, "../results/df_mutant_compare.xlsx")

### 2. visualize the codon usage difference
source("./helper/CA_analysis.r")
CA_analysis(seq_mut=seq_mut_orf_orf1ab, gene="orf1ab")
CA_analysis(seq_mut=seq_mut_orf_spike, gene="spike")
