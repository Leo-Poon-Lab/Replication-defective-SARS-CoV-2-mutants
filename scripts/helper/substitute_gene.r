load(file="./helper/df_orf.rdata")
substitute_gene <- function(gene=NA, seq_tmplate=NA, seq_sub=NA) {
	stopifnot(all(length(gene)==1, length(seq)==1))
	stopifnot(gene %in% df_orf$sequence)
	id_t <- df_orf$sequence==gene
	stopifnot(width(seq_sub)==df_orf$length[id_t])
	stopifnot(width(seq_tmplate)==29903)

	pos_start_t <- df_orf$start[df_orf$sequence==gene]
	pos_stop_t <- df_orf$stop[df_orf$sequence==gene]
	subseq(seq_tmplate, pos_start_t, pos_stop_t) <- seq_sub
	return(seq_tmplate)
}
