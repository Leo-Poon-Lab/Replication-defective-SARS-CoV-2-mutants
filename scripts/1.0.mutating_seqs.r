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
df_orf$length <- df_orf$stop - df_orf$start + 1
save(seq_ref, df_orf, file="./helper/df_orf.rdata")
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
mutant_1b <- make_amber_mutant(ref=seq_ref, gene=c(paste0("ORF", c("3a", "6", "7a", "7b", "10"))), target_aa="K", pos_mutable_nt=df_mutable, mut_codon="tag", num_of_mut_each_gene=2)
mutant_1b_log <- mutant_1b[[1]]
mutant_1b_seq <- mutant_1b[[2]]
names(mutant_1b_seq) <- "Mutant_1a"
write_csv(mutant_1b_log, "../results/mutants/mutant_1b_info.csv")
writeXStringSet(mutant_1b_seq, "../results/mutants/mutant_1b_seq.fasta")

# Mutant group 2: modify codon usages.
## Mutant 2A: reshuffle codons while keeping the original codon usage
### create a full-genome orf sequence
source("./helper/make_full_genome.r")
source("./helper/compare_seqs_to_table.r")
source("./helper/make_regioned_data.r")

#### mutant 2a_all_inter : reshuffle codons in all genes, allow inter-gene reshuffle
rd_orf_full <- make_rd(list_genes=df_orf$sequence)
set.seed(2022)
rd_2a_all_inter <- codon_random(rd_orf_full, keep = T)
stopifnot(translate(get_dna(rd_2a_all_inter)) == translate(seq_orf_full)) # All are silent mutations
seq_2a_all_inter <- make_full_genome(seq_orf=get_dna(rd_2a_all_inter), list_genes=df_orf$sequence)
names(seq_2a_all_inter) <- "Mutant_2A_all_inter"
writeXStringSet(seq_2a_all_inter, "../results/mutants/mutant_2a_all_inter_seq.fasta")
df_2a_all_inter <- compare_seqs(seq_int=seq_2a_all_inter, seq_ref=seq_ref)
write_csv(df_2a_all_inter, "../results/mutants/mutant_2a_all_inter_info.csv")

#### mutant 2a_all_intra : reshuffle codons in all genes, *avoid* inter-gene reshuffle
orf_2a_all_intra <- sapply(df_orf$sequence, function(gene_t) {
	rd_orf_gene_t <- make_rd(gene_t)
	set.seed(2022)
	rd_2a_all_intra_gene_t <- codon_random(rd_orf_gene_t, keep = T)
	get_dna(rd_2a_all_intra_gene_t)
})
orf_2a_all_intra <- do.call(xscat, orf_2a_all_intra)
seq_2a_all_intra <- make_full_genome(seq_orf=orf_2a_all_intra, list_genes=df_orf$sequence)
names(seq_2a_all_intra) <- "Mutant_2A_all_intra"
writeXStringSet(seq_2a_all_intra, "../results/mutants/mutant_2a_all_intra_seq.fasta")
df_2a_all_intra <- compare_seqs(seq_int=seq_2a_all_intra, seq_ref=seq_ref)
write_csv(df_2a_all_intra, "../results/mutants/mutant_2a_all_intra_info.csv")

#### mutant 2a_os_inter : reshuffle codons in ORF1ab and Spike, allow inter-gene reshuffle
list_genes_os <- df_orf$sequence[grepl("^S$", df_orf$sequence) | grepl("^nsp", df_orf$sequence)]
rd_orf_os <- make_rd(list_genes=list_genes_os)
set.seed(2022)
rd_2a_os_inter <- codon_random(rd_orf_os, keep = T)
seq_2a_os_inter <- make_full_genome(seq_orf=get_dna(rd_2a_os_inter), list_genes=list_genes_os)
names(seq_2a_os_inter) <- "Mutant_2a_os_inter"
writeXStringSet(seq_2a_os_inter, "../results/mutants/mutant_2a_os_inter_seq.fasta")
df_2a_os_inter <- compare_seqs(seq_int=seq_2a_os_inter, seq_ref=seq_ref)
write_csv(df_2a_os_inter, "../results/mutants/mutant_2a_os_inter_info.csv")

#### mutant 2a_os_intra : reshuffle codons in ORF1ab and Spike, avoid inter-gene reshuffle
orf_2a_os_intra <- sapply(list_genes_os, function(gene_t) {
	rd_orf_gene_t <- make_rd(gene_t)
	set.seed(2022)
	rd_2a_os_intra_gene_t <- codon_random(rd_orf_gene_t, keep = T)
	get_dna(rd_2a_os_intra_gene_t)
})
orf_2a_os_intra <- do.call(xscat, orf_2a_os_intra)
seq_2a_os_intra <- make_full_genome(seq_orf=orf_2a_os_intra, list_genes=list_genes_os)
names(seq_2a_os_intra) <- "Mutant_2a_os_intra"
writeXStringSet(seq_2a_os_intra, "../results/mutants/mutant_2a_os_intra_seq.fasta")
df_2a_os_intra <- compare_seqs(seq_int=seq_2a_os_intra, seq_ref=seq_ref)
write_csv(df_2a_os_intra, "../results/mutants/mutant_2a_os_intra_info.csv")

#### mutant 2a_os_intra_nonnsp : reshuffle codons in ORF1ab (full ORF1ab, no nsp) and Spike, allow inter-gene reshuffle
list_genes_orf1ab <- df_orf$sequence[grepl("^nsp", df_orf$sequence)]
list_genes_spike <- "S"
rd_orf_o <- make_rd(list_genes=list_genes_orf1ab)
rd_orf_s <- make_rd(list_genes=list_genes_spike)
set.seed(2022)
rd_2a_o_intra <- codon_random(rd_orf_o, keep = T)
rd_2a_s_intra <- codon_random(rd_orf_s, keep = T)
seq_orf_combined <- xscat(get_dna(rd_2a_o_intra), get_dna(rd_2a_s_intra))
seq_2a_os_intra_nonnsp <- make_full_genome(seq_orf=seq_orf_combined, list_genes=c(list_genes_orf1ab, list_genes_spike))
names(seq_2a_os_intra_nonnsp) <- "Mutant_2a_os_intra_nonsnp"
writeXStringSet(seq_2a_os_intra_nonnsp, "../results/mutants/mutant_2a_os_intra_nonnsp.fasta")
df_2a_os_intra_nonnsp <- compare_seqs(seq_int=seq_2a_os_intra_nonnsp, seq_ref=seq_ref)
write_csv(df_2a_os_intra_nonnsp, "../results/mutants/mutant_2a_os_intra_nonnsp_info.csv")

## Mutant 2B: 
### This part was done based on our previous work: https://academic.oup.com/ve/article/6/1/veaa032/5837024?login=true
### Using the previous alignment file, We tried to mutate/mimic the codon usage pattern to "Human_Caen1", "Human_KFMC-7", "Bat_KJ473821", "Rodent_SJHM", "Pangolin_P1E", "Swine_PHEV"
### We try to mutate two ORFs: ORF1ab and Spike
seq_betacov_orf1ab <- readDNAStringSet("../data/betacov_orf1ab.fasta")
seq_betacov_orf1ab <- seq_betacov_orf1ab[-1]
seq_betacov_spike <- readDNAStringSet("../data/betacov_spike.fasta")
seq_betacov_spike <- seq_betacov_spike[-1]

list_genes_orf1ab <- df_orf$sequence[grepl("^nsp", df_orf$sequence)]
list_genes_spike <- "S"

#### mutant 2b
sapply(seq_along(seq_betacov_orf1ab), function(i) {
	seq_bc_orf1ab_i <- seq_betacov_orf1ab[i]
	seq_bc_spike_i <- seq_betacov_spike[i]
	name_mut_i <- paste0("Mutant_2b_", names(seq_bc_orf1ab_i))
	name_mut_sim_i <- paste0("Mutant_2b_", i)
	
	rd_2b_orf1ab_i <- input_seq(get_orf_seq(list_genes_orf1ab), get_mutable(list_genes_orf1ab))
	set.seed(2022)
	rd_2b_orf1ab_mut_i <- codon_mimic(rd_2b_orf1ab_i, alt=seq_bc_orf1ab_i)
	
	rd_2b_spike_i <- input_seq(get_orf_seq(list_genes_spike), get_mutable(list_genes_spike))
	set.seed(2022)
	rd_2b_spike_mut_i <- codon_mimic(rd_2b_spike_i, alt=seq_bc_spike_i)
	
	seq_orf_combined <- xscat(get_dna(rd_2b_orf1ab_mut_i), get_dna(rd_2b_spike_mut_i))
	seq_2b_i <- make_full_genome(seq_orf=seq_orf_combined, list_genes=c(list_genes_orf1ab, list_genes_spike))

	names(seq_2b_i) <- name_mut_i
	writeXStringSet(seq_2b_i, paste0("../results/mutants/", name_mut_sim_i, "_seq.fasta"))
	df_2b_i <- compare_seqs(seq_int=seq_2b_i, seq_ref=seq_ref)
	write_csv(df_2b_i, paste0("../results/mutants/", name_mut_sim_i, "_info.csv"))
})

## Mutant 2C:
### We aim to generate a mutant 2C which uses most unfavorable codons.
#### In the initial version, we try to substitute all codons into their most under-represented synonymous counterparts, in ORF1ab and Spike respectively.
list_genes_os <- df_orf$sequence[grepl("^S$", df_orf$sequence) | grepl("^nsp", df_orf$sequence)]
list_genes_orf1ab <- df_orf$sequence[grepl("^nsp", df_orf$sequence)]
list_genes_spike <- "S"

find_underrepresented_codons <- function(seq_orf) {
	cu_orf <- get_cu(seq_orf)
	tmp <- sapply(names(AMINO_ACID_CODE), function(aa_i) {
		codon_i <- names(GENETIC_CODE[GENETIC_CODE==aa_i])
		idx_i <- colnames(cu_orf) %in% codon_i
		cu_i <- cu_orf[,idx_i]
		if(length(cu_i)==1){return(NULL)}else{return(names(cu_i)[which.min(cu_i)])}
	})
	return(unlist(tmp))
}

seq_orf_orf1ab <- get_orf_seq(list_genes_orf1ab)
codon_ur_orf1ab <- find_underrepresented_codons(seq_orf_orf1ab)
rd_2c_orf1ab_mut_i <- input_seq(seq_orf_orf1ab, get_mutable(list_genes_orf1ab))
set.seed(2022)
sapply(codon_ur_orf1ab, function(codon_i) {
	print(codon_i)
	rd_2c_orf1ab_mut_i <<- codon_to(rd_2c_orf1ab_mut_i, max.codon=codon_i)
	return("")
})

seq_orf_spike <- get_orf_seq(list_genes_spike)
codon_ur_spike <- find_underrepresented_codons(seq_orf_spike)
rd_2c_spike_mut_i <- input_seq(seq_orf_spike, get_mutable(list_genes_spike))
set.seed(2022)
sapply(codon_ur_spike, function(codon_i) {
	print(codon_i)
	rd_2c_spike_mut_i <<- codon_to(rd_2c_spike_mut_i, max.codon=codon_i)
	return("")
})

seq_orf_combined <- xscat(get_dna(rd_2c_orf1ab_mut_i), get_dna(rd_2c_spike_mut_i))
seq_2c <- make_full_genome(seq_orf=seq_orf_combined, list_genes=c(list_genes_orf1ab, list_genes_spike))

names(seq_2c) <- "Mutant_2c"
writeXStringSet(seq_2c, paste0("../results/mutants/mutant_2c_seq.fasta"))
df_seq_2c <- compare_seqs(seq_int=seq_2c, seq_ref=seq_ref)
write_csv(df_seq_2c, paste0("../results/mutants/mutant_2c_info.csv"))
