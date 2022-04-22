require(Biostrings)
require(tidyverse)

source('./helper/substitute_gene.r')

make_full_genome <- function(seq_orf, list_genes) {
	seq_out <- seq_ref
	seq_orf_t <- seq_orf
	tmp <- lapply(list_genes, function(x) {
		pos_start_t <- df_orf$start[df_orf$sequence==x]
		pos_stop_t <- df_orf$stop[df_orf$sequence==x]
		length_orf <- pos_stop_t-pos_start_t+1
		print(x)
		print(c(subseq(seq_out, pos_start_t, pos_stop_t), subseq(seq_orf_t, 1, length_orf)))
		seq_out <<- substitute_gene(gene=x, seq_tmplate=seq_out, seq_sub=subseq(seq_orf_t, 1, length_orf))
		seq_orf_t <<- subseq(seq_orf_t, length_orf+1, width(seq_orf_t))
	})
	return(seq_out)
}

