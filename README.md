# Designing replication-defective SARS-CoV-2 mutant vaccines with biased codon usages

*Authors*

---

## Abstract
This repository hosts the scripts and data used for designing genomic sequences of replication-defective SARS-CoV-2 mutants.

## Details
`./scripts/0.determine_mutable_region.R`: This R code is used for determining mutable codons in different ORFs in SARS-CoV-2 genome. The modifiable regions should be: 1. is coding region; 2.not conserved in codon; 3.not enzyme binding sites; 4. not in regions with overlapped orfs, except mutations within such overlapped orf regions were previously found in VOCs.

`./scripts/1.mutating_seqs.r`: This code includes functions for generating different mutants in this study.

## Environment
```
R version 4.3.1 (2023-06-16)
Platform: x86_64-apple-darwin20 (64-bit)
Running under: macOS 15.3.2

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Asia/Hong_Kong
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
[1] compiler_4.3.1 cli_3.6.2      jsonlite_1.8.8 rlang_1.1.2
```

---
