# PAMfinder
PAMfinder identifies novel “NGG” protospacer adjacent motifs (PAMs) from tumor-normal (T-N) subtracted variant call files (VCFs) generated from next generation sequencing. Somatic NGG PAMs can serve as tumor-specific CRISPR-Cas9 targets.  

## Background
CRISPR-Cas9 requires a specific protospacer adjacent motif (PAM) as a binding signal for Cas9 activity, as the PAM sequence evolutionally allows the bacterial adaptive immune system to differentiate between self (contains the gRNA sequence but lacks the PAM) and non-self DNA (contains both the gRNA sequence and the PAM). Analogously, CRISPR-Cas9 can distinguish between cancer cells (containing PAMs from somatic mutations) and normal cells (lacking PAMs) for selective cancer killing. The most commonly used CRISPR-Cas9 system, derived from Streptococcus pyogenes, recognizes a PAM sequence of “NGG” (any nucleotide followed by two guanine (G) bases).  

## Description
When supplied with a tumor-normal (T-N) subtracted VCF (ideally one generated using mutect2), PAMfinder performs the following:
1)	Filter out variants that don’t pass the quality check
    - Only accepts variants that "PASS" under “Filter”
2)	Filter for variants that produce only one single base substitution (SBS)
    - Indels are excluded
    - SBSs with more than 1 alternative allele are excluded
3)	Filter for variants that generate a novel C or G
4)	Filter for variants based on read depth for both tumor and normal
    - Currently set at >= 18 for whole genome sequencing data
5)	Filter for variants based on variant allele frequency (VAF)
    - Currently set at 30% for cell lines
6)	Sequence reconstruction using alternative allele
7)	Filter for variants that produce novel NGG PAM
8)	Design sgRNA

## Reference
a) Selina Shiqing K Teh, Kirsten Bowland, Eitan Halper-Stromberg, Akhil Kotwal, Alexis Bennett, Alyza Skaist, Jacqueline Tang, Fidel Cai, Antonella Macoretta, Hong Liang, Hirohiko Kamiyama, Sarah Wheelan, Ming-Tseh Lin, Ralph H Hruban, Chien-Fu Hung, Michael Goldstein, Robert B Scharpf, Nicholas J Roberts, James R Eshleman, CRISPR-Cas9 for selective targeting of somatic mutations in pancreatic cancers, NAR Cancer, Volume 6, Issue 2, June 2024, zcae028, https://doi.org/10.1093/narcan/zcae028
b) 
