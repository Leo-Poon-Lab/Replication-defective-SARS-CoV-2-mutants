# Aim: identified highly conserved sequences that are shared by all known huaman and animal sarbecoviruses
# Input: aligned amino acid sequences of different genes
# Output: highly conserved regions of different genes in table format

# 1. Alignment of the sequences
library(Biostrings)
library(tidyverse)
library(SynMut)
library(readxl)
df_seqs <- readxl::read_excel("../data/data_seqs_all.xlsx")

## merge sequences
seq_gisaid <- readDNAStringSet("../data/fasta/gisaid_hcov-19_2021_09_16_08.fasta")
seq_genbank <- readDNAStringSet("../data/fasta/sequence_nt-2.fasta")

names(seq_genbank) <- sapply(names(seq_genbank), function(x){
	tmp <- strsplit(x, " ")[[1]][1]
	strsplit(tmp, ".", fixed=T)[[1]][1]
}, USE.NAMES = F)
names(seq_gisaid) <- sapply(names(seq_gisaid), function(x){
	strsplit(x, "|", fixed=T)[[1]][2]
})

seqs_all <- c(seq_genbank, seq_gisaid)
writeXStringSet(seqs_all[names(seqs_all)=="MN908947"], "../data/reference.fasta")
writeXStringSet(seqs_all[names(seqs_all)!="MN908947"], "../data/seqs_all.fasta")

## align with MAFFT
system("mafft --thread -1 --localpair --maxiterate 1000 --keeplength --add ../data/seqs_all.fasta ../data/reference.fasta > ../data/seqs_all_aln.fasta")

# 2. finding conserved regions
## reference: The coding capacity of SARS-CoV-2, https://doi.org/10.1038/s41586-020-2739-1
## In addition to the canonical viral open reading frames (ORFs), an addtionnal 23 unannotated viral ORFs were indentfied.

seq_aln <- readDNAStringSet("../data/seqs_all_aln.fasta")
df_orf <- read_csv("../data/ORF_SCoV2.csv")
df_orf <- df_orf %>% arrange(start)
df_orf$start<=c(0, df_orf$stop[-length(df_orf$stop)]) # only nsp12 had overlapping reading frame
(df_orf$stop-df_orf$start+1) %% 3 # CDS OK

df_orf_new <- read_tsv("../data/ORF_SCoV2_new.bed", col_names = F)
names(df_orf_new)[1:4] <- c("Reference", "start", "stop", "sequence")
(df_orf_new$stop-df_orf_new$start+1) %% 3 # the start position is zero-base in bed file
df_orf_new$start <- df_orf_new$start+1 # change to one-base

df_orf_all <- bind_rows(df_orf, df_orf_new %>% select(sequence, start, stop))

df_masking <- tibble(position_nt=1:29903, coding=FALSE, position_aa=NA, position_codon=NA, orf=NA, conserved_nt=FALSE, conserved_codon=FALSE, deletion_nt=FALSE, deletion_codon=FALSE, enzyme_binding=FALSE, multiple_orf=FALSE, frameshift_orfs=NA, voc_mut_codon=FALSE, frameshift_orfs_found_in_VOC=NA)
dir.create("../data/alignment_gene")

df_con_sum <- lapply(seq_along(df_orf_all$sequence), function(i){
	gene_i <- df_orf_all$sequence[i]
	seq_aln_i <- subseq(seq_aln, df_orf_all$start[i], df_orf_all$stop[i])
	writeXStringSet(seq_aln_i, paste0("../data/alignment_gene/seq_aln_", gene_i, ".fasta"))
	con_mat_i <- consensusMatrix(seq_aln_i)
	df_conserved_nt <- tibble(pos=seq_len(width(seq_aln_i)[1]))
	pos_conserved_nt <- which(apply(con_mat_i, 2, function(x){sum(x[1:4]>0)==1}))
	df_conserved_nt$conserved_nt <- df_conserved_nt$pos %in% pos_conserved_nt
	pos_indel <- which(apply(con_mat_i, 2, function(x){x[rownames(con_mat_i)=="-"]>0}))
	df_conserved_nt$deletion_nt <- df_conserved_nt$pos %in% pos_indel
	write_csv(df_conserved_nt, paste0("../data/alignment_gene/pos_conserved_nt_", gene_i, ".csv"))
	index_i <- df_orf_all$start[i]:df_orf_all$stop[i]
	df_masking$coding[index_i] <<- TRUE
	df_masking$position_aa[index_i] <<- paste(df_masking$position_aa[index_i], seq_along(index_i), sep="|")
	df_masking$position_aa[index_i] <<- gsub("NA|", "", df_masking$position_aa[index_i], fixed = T)
	df_masking$position_codon[index_i] <<- paste(df_masking$position_codon[index_i], rep(1:3, length(index_i)/3), sep="|")
	df_masking$position_codon[index_i] <<- gsub("NA|", "", df_masking$position_codon[index_i], fixed = T)	
	df_masking$orf[index_i] <<- paste(df_masking$orf[index_i], gene_i, sep="|")
	df_masking$orf[index_i] <<- gsub("NA|", "", df_masking$orf[index_i], fixed = T)
	# T | F
	# T & F
	df_masking$conserved_nt[index_i] <<- df_conserved_nt$conserved_nt | df_masking$conserved_nt[index_i]
	df_masking$deletion_nt[index_i] <<- df_conserved_nt$deletion_nt | df_masking$deletion_nt[index_i]

	df_conserved_codon <- tibble(pos=seq_len(width(seq_aln_i)[1]/3))
	df_conserved_codon$conserved_codon <- sapply(df_conserved_codon$pos, function(x){
		all(df_conserved_nt$conserved_nt[(3*x-2):(3*x)])
	})
	df_conserved_codon$deletion_codon <- sapply(df_conserved_codon$pos, function(x){
		any(df_conserved_nt$deletion_nt[(3*x-2):(3*x)])
	})
	df_masking$conserved_codon[index_i] <<- rep(df_conserved_codon$conserved_codon, each=3) | df_masking$conserved_codon[index_i]
	df_masking$deletion_codon[index_i] <<- rep(df_conserved_codon$deletion_codon, each=3) | df_masking$deletion_codon[index_i]
		
	df_summary <- tibble(gene=gene_i)
	df_summary$percentage_conserved <- round(sum(df_conserved_codon$conserved_codon)/nrow(df_conserved_codon)*100,2)
	df_summary$percentage_deltion <- round(sum(df_conserved_codon$deletion_codon)/nrow(df_conserved_codon)*100,2)
	print(gene_i)
	return(df_summary)
})
df_con_sum <- bind_rows(df_con_sum)
write_csv(df_con_sum, "../results/summary_conserved_gene.csv")

# 3. label enzyme binding sites
df_enzyme <- read_xlsx("../data/SARS-CoV-2 sequences that should not be edited.xlsx")
df_masking$enzyme_binding <- sapply(df_masking$position_nt, function(x){
	# print(x)
	if(any(x<=df_enzyme$`end (WUHAN)` & x>=df_enzyme$`strart (WUHAN)`, na.rm=T)){
		return(TRUE)
	} else {
		return(FALSE)
	}
})

# 4. locating overlapping ORFs and frameshift ORFs
## we note that the frameshifting site at nsp12 is conserved over multiple viruses.
df_masking$multiple_orf <- grepl("|", df_masking$orf, fixed = T)

df_masking$frameshift_orfs <- sapply(df_masking$position_codon, function(x){
	if(is.na(x)){return(NA)}
	tmp <- strsplit(x, "|", fixed=T)[[1]]
	if(length(unique(tmp))==1){
		return(FALSE)
	}else {
		return(TRUE)
	}
})

# 5. check if the mutations on VOCs exist within frameshift ORFs
df_snps <- read_csv("https://raw.githubusercontent.com/Leo-Poon-Lab/VOC-SNPs-checking/main/results/df_defining_snps_nt.csv")
df_snps_voc <- df_snps %>% filter(grepl("(Alpha|Beta|Gamma|Delta)", variant)) 
df_snps_voc <- df_snps_voc %>% filter(!grepl("K417N", variant))

df_masking$voc_mut_codon <- sapply(df_masking$position_nt, function(x) {
	any(x>=df_snps_voc$pos_nt_start & x<=df_snps_voc$pos_nt_stop)
})
df_masking$frameshift_orfs_found_in_VOC <- df_masking$voc_mut_codon & df_masking$frameshift_orfs
frameshift_orf_exclude <- unique(df_masking$orf[df_masking$frameshift_orfs_found_in_VOC]) 
frameshift_orf_exclude <- frameshift_orf_exclude[!is.na(frameshift_orf_exclude)]
frameshift_orf_exclude # We identified 5 overlapping orf regions which seem to be modifiable, as mutations at these regions were found in different VOCs.

# 6. check whether in the  transcription regulatory sequence (TRS) region
## https://www.frontiersin.org/files/Articles/641445/fgene-12-641445-HTML-r1/image_m/fgene-12-641445-t001.jpg
df_TRS_region <- tibble(start=c(70, 21556, 25385, 26237, 26473, 27041, 27388, 27888, 28260))
df_TRS_region$stop <- df_TRS_region$start+5
df_masking$in_trs_region <- sapply(df_masking$position_nt, function(x) {
	any((x<=df_TRS_region$stop) & (x>=df_TRS_region$start))
})


# label Modifiable regions
# the modifiable regions should be: 1. is coding region; 2.not conserved in codon; 3.not enzyme binding sites; 4. not in regions with overlapped orfs, except mutations within such overlapped orf regions which were previously found in VOCs; 5. not in TRS regions.

df_masking$modifiable <- df_masking$coding & !df_masking$conserved_codon & !df_masking$enzyme_binding & (!df_masking$frameshift_orfs | (df_masking$orf %in% frameshift_orf_exclude)) & !df_masking$in_trs_region

write_csv(df_masking, "../results/df_mutatable_region.csv")
df_masking <- read_csv("../results/df_mutatable_region.csv")
