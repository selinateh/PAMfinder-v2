# Load necessary libraries
library(vcfR)
library(VariantAnnotation)
library(dplyr)
library(tidyr)
library(BSgenome.Hsapiens.UCSC.hg38)

# Specify the VCF file path
vcf_file_path <- "mutect2_test.vcf"  ## Change this to your actual VCF file path!!!

# Step 1: Use vcfR to read the VCF file
vcf_data <- read.vcfR(vcf_file_path)
vcf_header <- vcf_data@meta

# Step 2: Read the VCF file using VariantAnnotation
vcf_variant_data <- readVcf(vcf_file_path)

# Step 3: Extract fixed fields
chromosomes <- seqnames(vcf_variant_data)
positions <- start(vcf_variant_data)
refs <- ref(vcf_variant_data)
alts <- alt(vcf_variant_data)
filters <- rowData(vcf_variant_data)$FILTER

# Create a data frame from the fixed information
fixed_df <- data.frame(
  CHROM = chromosomes,
  POS = positions,
  REF = refs,
  ALT = sapply(alts, function(x) paste(x, collapse = ",")),  ## Keep ALT as is
  FILTER = filters,
  stringsAsFactors = FALSE
)

# Rename columns to desired format
colnames(fixed_df) <- c("Chr", "Pos", "Ref", "Alt", "Filter")

# Step 4: Extract read depth and VAF
ad_data <- extract.gt(vcf_data, element = "AD")
ad_df <- as.data.frame(ad_data)
af_data <- extract.gt(vcf_data, element = "AF")
VAF <- as.data.frame(af_data[, "Panc1328_T"])  ## Change column name based on VCF corresponding to tumor sample!!!
colnames(VAF) <- c("VAF")

# Calculate read depth from Tumor AD
ad_extracted <- ad_df %>%
  mutate(Tumor_Sample = as.character(ad_data[, "Panc1328_T"])) %>%  ## Change column name based on VCF corresponding to tumor sample!!!
  separate(Tumor_Sample, into = c("RefC", "AltC"), sep = ",") %>%  ## Separate AD into RefC and AltC
  mutate(T_DP = as.numeric(RefC) + as.numeric(AltC))  ## Calculate total depth T_DP

# Calculate read depth from Normal AD
ad_extracted <- ad_extracted %>%
  mutate(Normal_Sample = as.character(ad_data[, "Panc1328_N"])) %>%  ## Change column name based on VCF corresponding to normal sample!!!
  separate(Normal_Sample, into = c("RefC_N", "AltC_N"), sep = ",") %>%  ## Separate AD into RefC and AltC for normal
  mutate(N_DP = as.numeric(RefC_N) + as.numeric(AltC_N))  ## Calculate total depth N_DP

dp_df <- ad_extracted %>%
  select(T_DP, N_DP)

# Combine read depth and VAF data with fixed_df
combined_df <- cbind(fixed_df, dp_df, VAF)
combined_df$VAF <- as.numeric(combined_df$VAF)

# Step 5: Initial PASS filter
combined_df_pass <- combined_df %>%
  filter(Filter == "PASS")

# Step 6: Exclude REF or ALT > 1bp and those with more than 1 alternative allele
valid_indices <- which(nchar(combined_df_pass$Ref) == 1 & 
                         sapply(strsplit(combined_df_pass$Alt, ","), function(x) all(nchar(x) == 1) && length(x) == 1))

# Filter the fixed fields based on valid indices
combined_df_filtered <- combined_df_pass[valid_indices, ]

# Step 7: Filter for specific ALT alleles (G or C)
combined_df_GC <- combined_df_filtered %>%
  filter(Alt %in% c("G", "C"))

# Step 8: Filter for read depth >=18 for both tumor and normal
combined_df_DP <- combined_df_GC %>%
  filter(N_DP >= 18 & T_DP >= 18)

# Step 9: Filter VAF >= 0.3
combined_df_VAF <- combined_df_DP %>%
  filter(VAF >= 0.300)

# Step 10: Reconstruct sequences
pm=combined_df_VAF
pm$Chr <- paste0("chr", pm$Chr)
pm=pm %>% rowwise() %>% mutate(seq=as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, Chr, Pos, Pos+22)))
pm2=pm %>% rowwise() %>% mutate(fwdseqreconstruct= gsub(paste("^",Ref,sep=""),Alt,seq))
pm3=pm2 %>% rowwise() %>% mutate(behind_seq=as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, Chr, Pos-22,Pos-1)))
pm4=pm3 %>% rowwise() %>% mutate(reconstructbothdirections=paste(behind_seq,fwdseqreconstruct,sep=''))
pm5=pm4 %>% rowwise() %>% mutate(reconstructbothdirections_rc=as.character(reverseComplement(DNAString(as.character(reconstructbothdirections)))))
pm6=pm5

# Add a strand column based on the ALT values
pm6 <- pm6 %>%
  mutate(strand = case_when(
    Alt == "G" ~ "+",
    Alt == "C" ~ "-",
    TRUE ~ NA_character_  ## Default case if ALT is not G or C
  ))

# Full sequence reconstruction

pm6$reconstructed <- 0

for (i in 1:nrow(pm6)) {
  pm6$reconstructed[i] <- if(grepl('-', pm6$strand[i])) substr(pm6$reconstructbothdirections_rc[i], 1, 24) else if(grepl('+', pm6$strand[i])) substr(pm6$reconstructbothdirections[i], 1, 24)
}

# Step 11: Design sgRNA sequence

PAMcheck <- function(x, y) {
  # Initialize an empty vector to store results
  result <- character(length(x))
  
  for (i in seq_along(x)) {
    if (substr(x[i], nchar(x[i]) - 1, nchar(x[i])) == y) {
      result[i] <- "PAM"
    } else if (substr(x[i], nchar(x[i]) - 2, nchar(x[i]) - 1) == y) {
      result[i] <- "prePAM"
    } else {
      result[i] <- "noPAM"
    }
  }
  
  return(result)
}

pm7=pm6
pm7$PAMcheck <- PAMcheck(pm7$reconstructed, 'GG')
pm8=pm7

# Initialize the sgRNA column
pm8$sgRNA <- NA  ## Make sure the sgRNA column is initialized

# Loop through each row of pm8
for (j in 1:nrow(pm8)) {
  if (pm8$PAMcheck[j] == "PAM") {
    pm8$sgRNA[j] <- substr(pm8$reconstructed[j], 2, 21) 
  } else if (pm8$PAMcheck[j] == "prePAM") {
    pm8$sgRNA[j] <- substr(pm8$reconstructed[j], 1, 20)
  }
}

# Optionally, remove rows where sgRNA is NA
pm8 <- pm8[!is.na(pm8$sgRNA), ]

pm9 <- pm8 %>%
  select(Chr, Pos, Ref, Alt, VAF, T_DP, N_DP, sgRNA)

# Step 12: Export results
# Write summary comment lines
cat(
  paste0(
    "# in VCF: ", nrow(combined_df), "\n",
    "# PASS: ", nrow(combined_df_pass), "\n",
    "# SBS: ", nrow(combined_df_filtered), "\n",
    "# Novel G/C: ", nrow(combined_df_GC), "\n",
    "# DP >= 18: ", nrow(combined_df_DP), "\n",
    "# VAF >= 0.3: ", nrow(combined_df_VAF), "\n",
    "# Novel PAM: ", nrow(pm9), "\n"
  ),
  file = "test_PamFinderoutput.txt"  ## Change!!!
)

# Append the actual table
write.table(
  pm9,
  file = "test_PamFinderoutput.txt",  ## Change!!!
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  append = TRUE
)