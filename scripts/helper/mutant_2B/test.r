codon_mimic(seq_orf1ab, fasta_orf1ab[6])

cu.target = oligonucleotideFrequency(fasta_orf1ab[6], width = 3, step = 3)[1, ]
object = seq_orf1ab
codon.mimic(cu.target,
            dna.seq = object@dnaseq,
            region = object@region)

convert_to_seq <- function(dna.seq){
    seq <- lapply(as.character(dna.seq), function(x) {
        splitseq(s2c(x))
    })
    return(seq)
}

codon.count <- function(x) {
    base::split(x, GENETIC_CODE[order(names(GENETIC_CODE))])
}

freq <- function(x) {
    lapply(codon.count(x), function(x) {
        x / sum(x)
    })
}


mut.assign.back <- function(seq.region, mut.need){
    seq.mut <- lapply(seq_along(seq.region), function(i) {
        seq.tmp <- seq.region[[i]]
        seq.tmp.aa <- seqinr::translate(s2c(c2s(seq.tmp)))
        aa.names <- names(table(seq.tmp.aa))
        mut.need.tmp <- mut.need[[i]]
        for (j in seq_along(aa.names)) {
            for.sample <- as.character(which(seq.tmp.aa == aa.names[j]))
            pos.tmp <- as.numeric(sample(for.sample))
            mut.cd.tmp <- round(mut.need.tmp[[aa.names[j]]])
            if (!any(is.na(mut.cd.tmp)) & any(mut.cd.tmp>1)) {
                suppressWarnings(seq.tmp[pos.tmp] <-
                        rep(names(mut.cd.tmp), mut.cd.tmp))
            }
        }
        return(seq.tmp)
    })
}


codon.mimic.region <- function(cu.target, dna.seq, region){
    seq <- convert_to_seq(dna.seq)
    seq.region <- mapply(function(x, y) {
        return(x[y])
    }, seq, region, SIMPLIFY = FALSE)
    seq.fixed <- mapply(function(x, y) {
        return(x[!y])
    }, seq, region, SIMPLIFY = FALSE)

    cu.fixed <- get_cu(DNAStringSet(unlist(lapply(seq.fixed, c2s))))
    cu.ori <- get_cu(dna.seq)
    freq.target <- freq(cu.target)

    mut.need <- lapply(seq_along(seq), function(i) {
        count.ori <- codon.count(cu.ori[i,])
        count.fixed <- codon.count(cu.fixed[i,])
        mut.usage <- mapply(function(x, y) {
            sum(x) * y
        }, count.ori, freq.target, SIMPLIFY = FALSE)
        mut.need <- mapply(function(x, y) {
            x - y
        }, mut.usage, count.fixed, SIMPLIFY = FALSE)
        mut.need <- lapply(mut.need, function(x) {
            #fix the negative needs
            x.negative <- x[x < 0]
            x.positive <- x[x > 0]
            if (length(x.negative) > 0) {
                return((1 + sum(abs(
                    x.negative
                )) / sum(x.positive)) * x.positive)
            } else{
                return(x)
            }
        })
    })
    seq.mut <- mut.assign.back(seq.region, mut.need)
    names(seq.mut) <- names(seq.region)
    return(seq.mut)
}
