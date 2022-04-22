require(Biostrings)
require(tidyverse)

get_codon <- function(seq, pos_nt, pos_codon) {
	stopifnot(length(pos_nt)==length(pos_codon))
	vapply(seq_along(pos_codon), function(i) {
		pos_codon_i <- pos_codon[i]
		check_multiple <- grepl("|", pos_codon_i, fixed = T)
		if(check_multiple){
			pos_codon_i <- strsplit(pos_codon_i, "|", fixed=T)[[1]]
			pos_codon_i <- as.numeric(pos_codon_i)
		} else{
			pos_codon_i <- as.numeric(pos_codon_i)
		}
		codon_i <- vapply(pos_codon_i, function(pos_codon_ii) {
			pos_adj <- 1-pos_codon_ii
			pos_start <- pos_nt[i]+pos_adj
			as.character(subseq(seq, pos_start, pos_start+2))
		}, character(1), USE.NAMES = F)
		if(check_multiple){return(paste0(codon_i, collapse = "|"))}else{return(codon_i)}
	}, character(1), USE.NAMES = F)
}

translate_codon <- function(x) {
	tmp <- strsplit(x, "|", fixed=T)[[1]]
	aa_x <- sapply(tmp, function(y) {
		as.character(translate(DNAString(y)))
	})
	paste0(aa_x, collapse = "|")
}

compare_seqs <- function(seq_int, seq_ref) {
	tmp <- compareStrings(seq_int, seq_ref)
	mut_sites_t <- which(strsplit(tmp, "")[[1]]=="?")
	df_mut_t <- df_mutable %>% filter(position_nt %in% mut_sites_t)
	stopifnot(all(df_mut_t$modifiable))
	df_mut_t <- df_mut_t %>% select(position_nt:orf)
	df_mut_t$nt_ref <- strsplit(as.character(seq_ref), "")[[1]][mut_sites_t]
	df_mut_t$nt_mut <- strsplit(as.character(seq_int), "")[[1]][mut_sites_t]
	df_mut_t$codon_ref <- get_codon(seq_ref, df_mut_t$position_nt, df_mut_t$position_codon)
	df_mut_t$codon_mut <- get_codon(seq_int, df_mut_t$position_nt, df_mut_t$position_codon)
	df_mut_t$aa_ref <- vapply(df_mut_t$codon_ref, translate_codon, character(1))
	df_mut_t$aa_mut <- vapply(df_mut_t$codon_mut, translate_codon, character(1))
	return(df_mut_t)
}
