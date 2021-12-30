require(SynMut)
load(file="./helper/df_orf.rdata")


get_orf_seq <- function(list_genes) {
	seq_orf_full <- lapply(list_genes, function(x) {
		pos_start_t <- df_orf$start[df_orf$sequence==x]
		pos_stop_t <- df_orf$stop[df_orf$sequence==x]
		subseq(seq_ref, pos_start_t, pos_stop_t)
	})
	seq_orf_full <- do.call(xscat, seq_orf_full)
	return(seq_orf_full)
}

get_mutable <- function(list_genes) {
	check_mutable_orf_full <- lapply(list_genes, function(x) {
		pos_start_t <- df_orf$start[df_orf$sequence==x]
		pos_stop_t <- df_orf$stop[df_orf$sequence==x]
		df_tmp <- df_mutable[pos_start_t:pos_stop_t,]
		sapply(seq_len(nrow(df_tmp)/3), function(n) {
			all(df_tmp$modifiable[(3*n-2):(3*n)])
		})	
	})
	check_mutable_orf_full <- do.call(c, check_mutable_orf_full)
	return(data.frame(check_mutable_orf_full))
}

make_rd <- function(list_genes) {
	seq_orf_full <- get_orf_seq(list_genes)
	check_mutable_orf_full <- get_mutable(list_genes)

	rd_orf_full <- input_seq(seq_orf_full, check_mutable_orf_full)
	return(rd_orf_full)
}

