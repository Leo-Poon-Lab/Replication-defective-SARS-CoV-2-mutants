# Aim: generating sequences of different mutants
# Input: reference seq and mutable region
# Output: sequences of mutants
## Three groups of mutants

library(SynMut)
library(Biostrings)
library(tidyverse)
dir.create("../results/mutants")

seq_ref <- readDNAStringSet("../data/reference.fasta")
df_orf <- read_csv("../data/ORF_SCoV2.csv")
df_mutable <- read_csv("../results/df_mutatable_region.csv")

# helper function
QC_gene <- function(gene){
	gene_corrected <- sapply(gene, function(x){
		tmp <- df_orf$sequence[tolower(df_orf$sequence) ==tolower(x)]
		stopifnot("Please make sure gene id are correct"=length(tmp)>0)
		return(tmp)
	})
	return(gene_corrected)
}

Check_TRS <- function(x) {} ##TODO

# Mutant group 1: we will mutate the first two lysine codons of selected ORFs (not on conserved regions) to an amber codon (UAG)
## Mutant 1A: amber codons into the E and M genes.

make_amber_mutant <- function(ref=seq_ref, gene=c("e", "m"), target_aa="K", pos_mutable_nt=df_mutable, mut_codon="tag", num_of_mut_each_gene=2){
	stopifnot(!is.na(gene) & is.character(gene))
	gene_corrected <- QC_gene(gene)
	seq_ori <- seq_ref
	df_mut_log_all <- lapply(gene_corrected, function(gene_t){
		pos_start <- df_orf$start[df_orf$sequence==gene_t]
		pos_stop <- df_orf$stop[df_orf$sequence==gene_t]
		seq_mod <- subseq(seq_ori, pos_start, pos_stop)

		check <- strsplit(as.character(translate(seq_mod)), "")[[1]]==target_aa
		pos_target_aa <- which(check)
		if(length(pos_target_aa)<num_of_mut_each_gene){
			warning(paste0("There are only ", length(pos_target_aa), " codon(s) coding for ", target_aa, " ", "on gene ", gene_t, "!"))
		}

		pos_mutable_nt_t <- pos_mutable_nt[pos_start:pos_stop,]
		df_mut_log <- lapply(pos_target_aa, function(pos_mut_aa) {
			pos_mut_nt_i <- pos_mut_aa*3-2
			pos_mut_nt_j <- pos_mut_aa*3
			check_mutable <- all(pos_mutable_nt_t$modifiable[pos_mut_nt_i:pos_mut_nt_j])
			codon_ori <- as.character(subseq(seq_mod, pos_mut_nt_i, pos_mut_nt_j))
			return(tibble(Gene = gene_t, codon_ori=codon_ori, codon_mut=mut_codon, mutable=check_mutable, choose_for_mutation=FALSE, pos_mut_nt_i=pos_mut_nt_i, pos_mut_nt_j=pos_mut_nt_j))
		})
		df_mut_log <- bind_rows(df_mut_log)
		df_mut_log$pos_full_genome_i <- df_mut_log$pos_mut_nt_i + pos_start -1
		df_mut_log$pos_full_genome_j <- df_mut_log$pos_mut_nt_j + pos_start -1
		num_mutatble <- sum(df_mut_log$mutable)
		if(num_mutatble<num_of_mut_each_gene){
			warning(paste0("There are only ", num_mutatble, " codon(s) coding for ", target_aa, " ", "on gene ", gene_t, " which is modifiable!"))
		}
		if(num_mutatble>0){
			df_mut_log$choose_for_mutation[df_mut_log$mutable][seq_len(min(num_mutatble, num_of_mut_each_gene))] <- TRUE
		}		
		return(df_mut_log)
	})
	df_mut_log_bind <- bind_rows(df_mut_log_all)

	sapply(df_mut_log_all, function(df_t) {
		if(nrow(df_t)==0){return(NA)}
		df_t <- df_t %>% filter(choose_for_mutation)
		sapply(seq_len(nrow(df_t)), function(x) {
			subseq(seq_ori, df_t$pos_full_genome_i[x], df_t$pos_full_genome_j[x]) <<- DNAString(toupper(mut_codon))
		})		
	})
	return(list(df_mut_log_bind, seq_ori))	
}

mutant_1a <- make_amber_mutant(ref=seq_ref, gene=c("e", "m"), target_aa="K", pos_mutable_nt=df_mutable, mut_codon="tag", num_of_mut_each_gene=2)
mutant_1a_log <- mutant_1a[[1]]
mutant_1a_seq <- mutant_1a[[2]]
write_csv(mutant_1a_log, "../results/mutants/mutant_1a_info.csv")
writeXStringSet(mutant_1a_seq, "../results/mutants/mutant_1a_seq.fasta")

## mutant 1B: Amber codons into ORFs: 3a, 6, 7a, 7b and 10
df_orf$sequence
mutant_1b <- make_amber_mutant(ref=seq_ref, gene=c(paste0("ORF", c("3a", "6", "7a", "7b", "10"))), target_aa="K", pos_mutable_nt=df_mutable, mut_codon="tag", num_of_mut_each_gene=2)
mutant_1b_log <- mutant_1b[[1]]
mutant_1b_seq <- mutant_1b[[2]]
write_csv(mutant_1b_log, "../results/mutants/mutant_1b_info.csv")
writeXStringSet(mutant_1b_seq, "../results/mutants/mutant_1b_seq.fasta")
