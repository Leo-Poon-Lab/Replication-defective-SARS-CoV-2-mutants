# Designing replication-defective SARS-CoV-2 mutant vaccines with biased codon usages

*Authors*

---

## Abstract
This repository hosts the scripts and data used for designing genomic sequences of replication-defective SARS-CoV-2 mutants.

## Details
`./scripts/0.determine_mutable_region.R`: This R code is used for determining mutable codons in different ORFs in SARS-CoV-2 genome. The modifiable regions should be: 1. is coding region; 2.not conserved in codon; 3.not enzyme binding sites; 4. not in regions with overlapped orfs, except mutations within such overlapped orf regions were previously found in VOCs.
`./scripts/1.mutating_seqs.r`: This code includes functions for generating different mutants in this study.

## Environment


---
