#============================================================#
# Archaic-derived alleles across original and imputed genomes
#
# Purpose:
# This script identifies SNPs that are derived in archaic
# humans and ancestral in Africans, and projects those sites
# onto original and imputed ancient genomes to assess allele
# state concordance across coverage and phasing.
#
# Scope:
# - Designed for a single genomic region or gene at a time
# - Input VCFs must already be restricted to the target region
# - Works with both unphased (original) and phased (imputed)
#   ancient genomes
#
# Output:
# - A sample × SNP matrix of allele states
# - A heatmap visualizing ancestral, derived, and heterozygous
#   states across individual and coverage levels
#
#============================================================#
#------------------------------------------------------------#
# Load required R packages
#------------------------------------------------------------#
# vcfR         : VCF parsing and genotype extraction
# tidyverse    : data manipulation and reshaping
# purrr        : iteration over multiple VCFs
# ggplot2      : visualization
# ComplexHeatmap / circlize are loaded for compatibility
# with alternative heatmap implementations
#------------------------------------------------------------#

library(vcfR)
library(tibble)
library(stringr)
library(dplyr)
library(tidyr)
library(purrr)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)

#------------------------------------------------------------#
# Load sample metadata
#------------------------------------------------------------#
# The metadata table links VCF sample identifiers to
# standardized sample names and population labels.
# Both primary (IID) and alternative IDs are supported.
#------------------------------------------------------------#

metadata_info <- as.data.frame(read.delim("../Data/MetadataInfo_ColorsPop.csv", header = TRUE, sep = ","))

#------------------------------------------------------------#
# Input sanity checks and safe VCF/BCF import
# Ensures all files exist and are readable before analysis
#------------------------------------------------------------#

check_exists <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  invisible(TRUE)
}

read_vcf_safe <- function(path, verbose = FALSE) {
  check_exists(path)
  vcf <- tryCatch(
    read.vcfR(path, verbose = verbose),
    error = function(e) {
      stop("Failed to read with vcfR::read.vcfR(): ", path, "\n", conditionMessage(e))
    }
  )
  return(vcf)
}


#------------------------------------------------------------#
# Define input VCFs for original and imputed genomes
#------------------------------------------------------------#
# Each VCF corresponds to a single coverage level and
# contains harmonized individual names of the form:
#   <SampleID>_<Coverage>
# e.g. BA64_Original, BA64_ImputedOC, BA64_2X, ...
#------------------------------------------------------------#

vcf_files <- c(
  Original   = "../Data/Haplotypes/Original_LEMD2_MLN.bcf.gz",
  ImputedOC = "../Data/Haplotypes/ImputedOC_LEMD2_MLN.bcf.gz",
  `2X`      = "../Data/Haplotypes/2X_LEMD2_MLN.bcf.gz",
  `1X`      = "../Data/Haplotypes/1X_LEMD2_MLN.bcf.gz",
  `0.5X`    = "../Data/Haplotypes/05X_LEMD2_MLN.bcf.gz",
  `0.0625X` = "../Data/Haplotypes/00625X_LEMD2_MLN.bcf.gz"
)

#------------------------------------------------------------#
# Load original (Original high-coverage) ancient genomes
#------------------------------------------------------------#
# These individuals represent non-imputed, high-coverage
# genotypes and serve as the reference for comparison
# with imputed datasets.
#------------------------------------------------------------#

vcf_original <- read_vcf_safe(
  vcf_files["Original"],
  verbose = FALSE
)

#------------------------------------------------------------#
# Step 11. Load imputed genomes across coverage levels
#------------------------------------------------------------#
# Imputed VCFs are read into a named list, where names
# correspond to sequencing coverage levels.
#------------------------------------------------------------#

vcf_imputed_list <- lapply(
  vcf_files[names(vcf_files) != "Original"],
  function(f) read_vcf_safe(f, verbose = FALSE)
)

names(vcf_imputed_list)


#------------------------------------------------------------#
# Load reference panels for ancestry inference
#------------------------------------------------------------#
# These VCFs are used to define ancestral and derived
# alleles at each site:
# - Chimpanzee defines the ancestral allele
# - Archaic humans define the derived allele
# - Africans act as a modern human control
#------------------------------------------------------------#

vcf_archaic <- read_vcf_safe("../Data/Haplotypes/4Archaic_LEMD2_MLN.vcf.gz", verbose = FALSE)
vcf_africans <- read_vcf_safe("../Data/Haplotypes/1000GP_YRI_ESN_MSL_LEMD2_MLN.vcf.gz", verbose = FALSE)
vcf_chimp    <- read_vcf_safe("../Data/Haplotypes/hg19_panTro6_LEMD2_MLN.vcf.gz", verbose = FALSE)

cat("Archaic samples:\n");  print(colnames(vcf_archaic@gt)[-1])
cat("African samples:\n");  print(colnames(vcf_africans@gt)[-1])
cat("Chimp sample:\n");     print(colnames(vcf_chimp@gt)[-1])


#------------------------------------------------------------#
# Define unique variant identifiers (CHR_POS)
#------------------------------------------------------------#
# To ensure consistent merging across VCFs originating
# from different datasets, each variant is assigned a
# unique identifier of the form:
#   <CHROM>_<POS>
# This identifier is used as the primary key throughout
# the remainder of the analysis.
#------------------------------------------------------------#

prepare_vcf <- function(vcf) {
  chr_pos <- paste0(vcf@fix[, "CHROM"], "_", vcf@fix[, "POS"])
  chr_pos_u <- make.unique(chr_pos)
  
  rownames(vcf@fix) <- chr_pos_u
  rownames(vcf@gt)  <- chr_pos_u
  vcf@fix[, "ID"]   <- chr_pos_u
  
  vcf
}

# Apply to all VCFs
vcf_original  <- prepare_vcf(vcf_original)
vcf_archaic   <- prepare_vcf(vcf_archaic)
vcf_africans  <- prepare_vcf(vcf_africans)
vcf_chimp     <- prepare_vcf(vcf_chimp)

vcf_imputed_list <- lapply(vcf_imputed_list, prepare_vcf)

#------------------------------------------------------------#
# Extract genotypes from reference panels
#------------------------------------------------------------#
# Genotypes from archaic humans, Africans, and chimpanzee
# are used to determine ancestral and derived alleles
# at each site.
#------------------------------------------------------------#

gt_archaic  <- extract.gt(vcf_archaic,  element = "GT")
gt_africans <- extract.gt(vcf_africans, element = "GT")

# Chimpanzee is haploid in practice; use the first column
gt_chimp <- extract.gt(vcf_chimp, element = "GT")[, 1]
names(gt_chimp) <- rownames(vcf_chimp@gt)

cat("Number of sites in reference panels:\n")
cat("Archaic:", nrow(gt_archaic), "\n")
cat("Africans:", nrow(gt_africans), "\n")
cat("Chimp:", length(gt_chimp), "\n")


anyDuplicated(rownames(vcf_archaic@gt))
cat("archaic dup:", anyDuplicated(rownames(vcf_archaic@gt)), "\n")
anyDuplicated(rownames(vcf_africans@gt))
cat("africans dup:", anyDuplicated(rownames(vcf_africans@gt)), "\n")
anyDuplicated(rownames(vcf_chimp@gt))
cat("chimp dup:", anyDuplicated(rownames(vcf_chimp@gt)), "\n")
sapply(vcf_imputed_list, function(v) anyDuplicated(rownames(v@gt)))
print(sapply(vcf_imputed_list, function(v) anyDuplicated(rownames(v@gt))))

#------------------------------------------------------------#
# Extract REF and ALT alleles for each site
#------------------------------------------------------------#
# REF and ALT alleles are taken from the archaic VCF,
# which is assumed to contain all sites shared across
# datasets in the target region.
#------------------------------------------------------------#

refalt_df <- tibble(
  CHR_POS = rownames(vcf_archaic@gt),
  REF     = vcf_archaic@fix[, "REF"],
  ALT     = vcf_archaic@fix[, "ALT"]
)

if (any(grepl(",", refalt_df$ALT))) {
  stop("Multi-allelic sites detected. This script assumes strictly biallelic SNPs.")
}

#------------------------------------------------------------#
# Infer ancestral and derived alleles using chimpanzee
#------------------------------------------------------------#
# For each site, chimpanzee defines the ancestral allele:
# - If chimp genotype is 0/0 → ancestral = REF
# - If chimp genotype is 1/1 → ancestral = ALT
# Sites where chimp is heterozygous or missing are excluded.
# The derived allele is then defined as the opposite allele.
#------------------------------------------------------------#

chimp_gt_df <- tibble(
  CHR_POS   = names(gt_chimp),
  Chimp_GT  = gt_chimp
)

infer_ancestral_allele <- function(ref, alt, chimp_gt) {
  chimp_alleles <- unique(unlist(strsplit(chimp_gt, "[|/]")))
  chimp_alleles <- chimp_alleles[!chimp_alleles %in% c(".", "", NA)]
  if (length(chimp_alleles) == 0) return(NA_character_)
  if (all(chimp_alleles == "0")) return(ref)
  if (all(chimp_alleles == "1")) return(alt)
  return(NA_character_)  # exclude heterozygous / ambiguous chimp calls
}

site_annot <- refalt_df %>%
  inner_join(chimp_gt_df, by = "CHR_POS") %>%
  mutate(
    Ancestral_Allele = mapply(infer_ancestral_allele, REF, ALT, Chimp_GT),
    Derived_Allele   = ifelse(Ancestral_Allele == REF, ALT,
                              ifelse(Ancestral_Allele == ALT, REF, NA_character_))
  ) %>%
  drop_na(Ancestral_Allele, Derived_Allele)

cat("Sites with unambiguous ancestral/derived alleles:", nrow(site_annot), "\n")


#------------------------------------------------------------#
# Identify archaic-derived SNPs (derived in archaics, ancestral in Africans)
#------------------------------------------------------------#
# We retain sites that satisfy BOTH conditions:
# (i) At least one archaic individual carries the derived allele
#     (including heterozygous genotypes).
# (ii) All African individuals are homozygous for the ancestral allele.
#------------------------------------------------------------#

archaic_ids <- c("Denisova", "AltaiNeandertal", "Vindija33.19", "Chagyrskaya-Phalanx")

missing_archaics <- setdiff(archaic_ids, colnames(gt_archaic))
if (length(missing_archaics) > 0) {
  stop("Missing archaic samples in archaic VCF: ",
       paste(missing_archaics, collapse = ", "))
}

# Helper: check if a genotype contains the derived allele (even as 0/1)
has_derived_allele <- function(gt, ref, alt, derived_allele) {
  if (is.na(gt) || gt %in% c("./.", ".|.", ".", "")) return(FALSE)
  alleles <- unlist(strsplit(gt, "[|/]"))
  alleles <- alleles[!alleles %in% c(".", "", NA)]
  if (length(alleles) == 0) return(FALSE)
  
  allele_bases <- vapply(alleles, function(a) {
    if (a == "0") return(ref)
    if (a == "1") return(alt)
    return(NA_character_)
  }, character(1))
  
  any(allele_bases == derived_allele, na.rm = TRUE)
}

# Helper: check if *all* African genotypes are homozygous ancestral at a site
# If STRICT_AFRICAN_FILTER = TRUE:
#   - all African individuals must be homozygous ancestral
#   - missing genotypes cause the site to be rejected
# If STRICT_AFRICAN_FILTER = FALSE:
#   - missing genotypes are allowed
#   - no African individual may carry the derived allele
STRICT_AFRICAN_FILTER <- TRUE

africans_all_hom_ancestral <- function(gt_vec, ancestral_allele, ref, alt) {
  
  # ancestral allele corresponds to GT code "0" if ancestral == REF else "1"
  ancestral_code <- ifelse(ancestral_allele == ref, "0", "1")
  
  ok <- vapply(gt_vec, function(gt) {
    
    # Missing genotype
    if (is.na(gt) || gt %in% c("./.", ".|.", ".", "")) {
      return(if (STRICT_AFRICAN_FILTER) FALSE else NA)
    }
    
    alleles <- unlist(strsplit(gt, "[|/]"))
    alleles <- alleles[!alleles %in% c(".", "", NA)]
    
    if (length(alleles) == 0) {
      return(if (STRICT_AFRICAN_FILTER) FALSE else NA)
    }
    
    # TRUE if homozygous ancestral
    all(alleles == ancestral_code)
    
  }, logical(1))
  
  if (STRICT_AFRICAN_FILTER) {
    # All Africans must be confidently homozygous ancestral
    all(ok)
  } else {
    # Allow missing, but forbid any derived allele
    all(ok | is.na(ok))
  }
}


filtered_sites <- site_annot %>%
  rowwise() %>%
  mutate(
    archaic_gt = list(gt_archaic[CHR_POS, archaic_ids]),
    african_gt = list(gt_africans[CHR_POS, , drop = TRUE]),
    
    archaic_has_derived = any(vapply(archaic_gt, has_derived_allele,
                                     ref = REF, alt = ALT, derived_allele = Derived_Allele,
                                     logical(1))),
    
    african_is_hom_ancestral = africans_all_hom_ancestral(
      gt_vec = african_gt,
      ancestral_allele = Ancestral_Allele,
      ref = REF,
      alt = ALT
    )
  ) %>%
  ungroup() %>%
  filter(archaic_has_derived, african_is_hom_ancestral)

target_chrpos <- filtered_sites$CHR_POS

cat("rchaic-derived SNPs retained:", length(target_chrpos), "\n")
print(head(filtered_sites %>% select(CHR_POS, REF, ALT, Ancestral_Allele, Derived_Allele), 10))

#------------------------------------------------------------#
# Restrict original and imputed genomes to archaic-derived SNPs
#------------------------------------------------------------#
# All downstream analyses are performed only on the set of
# SNPs identified as archaic-derived and African-ancestral.
#------------------------------------------------------------#

filter_vcf_by_chrpos <- function(vcf, chrpos_vec) {
  # Always recompute CHR_POS from FIX
  vcf_chrpos <- paste0(vcf@fix[, "CHROM"], "_", vcf@fix[, "POS"])
  
  keep <- which(vcf_chrpos %in% chrpos_vec)
  
  if (length(keep) == 0) {
    stop("No matching variants found after filtering.")
  }
  
  vcf_filt <- vcf[keep, ]
  
  # Enforce rownames AFTER filtering
  rownames(vcf_filt@gt)  <- paste0(vcf_filt@fix[, "CHROM"], "_", vcf_filt@fix[, "POS"])
  rownames(vcf_filt@fix) <- rownames(vcf_filt@gt)
  
  return(vcf_filt)
}

# Filter original genomes
vcf_original_filt <- filter_vcf_by_chrpos(vcf_original, target_chrpos)

# Filter imputed genomes (all coverages)
vcf_imputed_filt <- lapply(vcf_imputed_list, filter_vcf_by_chrpos,
                           chrpos_vec = target_chrpos)

# Sanity checks
cat("Original variants retained:", nrow(vcf_original_filt@gt), "\n")
cat("Imputed variants retained per coverage:\n")
print(sapply(vcf_imputed_filt, function(v) nrow(v@gt)))

#------------------------------------------------------------#
# Convert genotypes into allele states (Ancestral/Derived/Het)
#------------------------------------------------------------#
# We classify genotypes using the inferred ancestral/derived
# alleles. This makes genotypes comparable across datasets.
#
# Output labels:
# - "Ancestral"     : only ancestral allele(s)
# - "Derived"       : only derived allele(s)
# - "Heterozygous"  : one ancestral and one derived allele
# - "Different"     : allele(s) not matching ancestral/derived
# - "Missing"       : missing/unparseable genotype
#------------------------------------------------------------#

geno_to_state <- function(gt,
                          ref, alt,
                          ancestral_allele, derived_allele) {
  
  if (is.na(gt) || gt %in% c("./.", ".|.", ".", "")) {
    return("Missing")
  }
  
  # Normalize
  ref              <- toupper(trimws(ref))
  alt              <- toupper(trimws(alt))
  ancestral_allele <- toupper(trimws(ancestral_allele))
  derived_allele   <- toupper(trimws(derived_allele))
  
  # Split genotype
  alleles <- strsplit(gt, "[|/]")[[1]]
  alleles <- alleles[alleles != "."]
  
  if (length(alleles) == 0) {
    return("Missing")
  }
  
  # Translate GT → bases
  bases <- vapply(alleles, function(a) {
    if (a == "0") return(ref)
    if (a == "1") return(alt)
    NA_character_
  }, character(1))
  
  bases <- bases[!is.na(bases)]
  if (length(bases) == 0) {
    return("Missing")
  }
  
  unique_bases <- unique(bases)
  
  # Classify
  if (!all(unique_bases %in% c(ancestral_allele, derived_allele))) {
    return("Different")
  }
  
  if (length(unique_bases) == 1 && unique_bases == ancestral_allele) {
    return("Ancestral")
  }
  
  if (length(unique_bases) == 1 && unique_bases == derived_allele) {
    return("Derived")
  }
  
  if (length(unique_bases) == 2) {
    return("Heterozygous")
  }
  
  return("Different")
}

##------------------------------------------------------------#
# Extract African allele states and collapse to one row
#------------------------------------------------------------#
# Africans are used as a control group. We compute their
# allele states (for QC) and then collapse them to a single
# representative row for visualization.
#------------------------------------------------------------#
gt_africans <- gt_africans[target_chrpos, , drop = FALSE]
stopifnot(identical(rownames(gt_africans), target_chrpos))

allele_states_africans <- do.call(
  rbind,
  lapply(colnames(gt_africans), function(ind) {
    vapply(seq_along(target_chrpos), function(i) {
      geno_to_state(
        gt = gt_africans[i, ind],
        ref = filtered_sites$REF[i],
        alt = filtered_sites$ALT[i],
        ancestral_allele = filtered_sites$Ancestral_Allele[i],
        derived_allele   = filtered_sites$Derived_Allele[i]
      )
    }, character(1))
  })
)

rownames(allele_states_africans) <- colnames(gt_africans)
colnames(allele_states_africans) <- target_chrpos

# QC: Africans should be almost entirely ancestral
african_state_table <- table(unlist(allele_states_africans), useNA = "ifany")
print(african_state_table)

if (sum(allele_states_africans %in% c("Derived","Heterozygous","Different")) > 0) {
  warning("Some African genotypes are not strictly ancestral at retained sites.")
}

# Collapse Africans into a single row for plotting
allele_states_africans_collapsed <- matrix(
  "Ancestral",
  nrow = 1,
  ncol = length(target_chrpos),
  dimnames = list("Africans", target_chrpos)
)


#------------------------------------------------------------#
# Extract allele states for Chimp and Archaics
#------------------------------------------------------------#
# These reference groups are projected onto the same
# archaic-derived SNP set for visualization purposes.
#------------------------------------------------------------#

#--------------------#
# Chimp (ancestral)
#--------------------#

allele_states_chimp <- matrix(
  "Ancestral",
  nrow = 1,
  ncol = length(target_chrpos),
  dimnames = list("Chimp", target_chrpos)
)

#--------------------#
# Archaic individuals
#--------------------#
gt_archaic <- gt_archaic[target_chrpos, archaic_ids, drop = FALSE]
stopifnot(identical(rownames(gt_archaic), target_chrpos))

allele_states_archaic <- do.call(
  rbind,
  lapply(archaic_ids, function(ind) {
    vapply(seq_along(target_chrpos), function(i) {
      geno_to_state(
        gt = gt_archaic[i, ind],
        ref = filtered_sites$REF[i],
        alt = filtered_sites$ALT[i],
        ancestral_allele = filtered_sites$Ancestral_Allele[i],
        derived_allele   = filtered_sites$Derived_Allele[i]
      )
    }, character(1))
  })
)

rownames(allele_states_archaic) <- archaic_ids
colnames(allele_states_archaic) <- target_chrpos


#------------------------------------------------------------#
# Extract allele states for ORIGINAL (unphased diploid)
#------------------------------------------------------------#

gt_original <- extract.gt(vcf_original_filt, element = "GT")
gt_original <- gt_original[target_chrpos, , drop = FALSE]
stopifnot(identical(rownames(gt_original), target_chrpos))

individual_ids_original <- colnames(gt_original)

allele_states_original <- sapply(individual_ids_original, function(ind) {
  vapply(seq_along(target_chrpos), function(i) {
    geno_to_state(
      gt = gt_original[i, ind],
      ref = filtered_sites$REF[i],
      alt = filtered_sites$ALT[i],
      ancestral_allele = filtered_sites$Ancestral_Allele[i],
      derived_allele   = filtered_sites$Derived_Allele[i]
    )
  }, character(1))
})

# allele_states_original is sites × individuals; transpose to individuals × sites
allele_states_original <- t(allele_states_original)
colnames(allele_states_original) <- target_chrpos

# Sanity check: dimensions
cat("Original allele-state matrix:", dim(allele_states_original), "\n")

#------------------------------------------------------------#
# Extract allele states for IMPUTED datasets (phased; haplotypes A/B)
#------------------------------------------------------------#

haplotype_to_state <- function(gt, phase,
                               ref, alt,
                               ancestral_allele, derived_allele) {
  
  if (is.na(gt) || gt %in% c("./.", ".|.", ".", "")) {
    return("Missing")
  }
  
  sep <- ifelse(grepl("\\|", gt), "\\|", "/")
  alleles <- strsplit(gt, sep)[[1]]
  
  if (length(alleles) < phase) return("Missing")
  
  allele_code <- alleles[phase]
  if (allele_code == ".") return("Missing")
  
  ref              <- toupper(ref)
  alt              <- toupper(alt)
  ancestral_allele <- toupper(ancestral_allele)
  derived_allele   <- toupper(derived_allele)
  
  base <- if (allele_code == "0") ref else if (allele_code == "1") alt else NA
  if (is.na(base)) return("Missing")
  
  if (base == ancestral_allele) return("Ancestral")
  if (base == derived_allele)   return("Derived")
  
  "Different"
}

# If TRUE, append coverage label to haplotype names
# (required if coverage information is not encoded in sample IDs)
ADD_COVERAGE_TO_HAPLO_NAME <- FALSE

haplo_vecs <- list()

for (cov in names(vcf_imputed_filt)) {
  
  message("Processing imputed coverage: ", cov)
  
  vcf_obj <- vcf_imputed_filt[[cov]]
  gt_imp  <- extract.gt(vcf_obj, element = "GT")
  
  # Ensure variant order matches target SNP set
  gt_imp <- gt_imp[target_chrpos, , drop = FALSE]
  
  # Validate alignment before allele-state projection
  stopifnot(identical(rownames(gt_imp), target_chrpos))
  
  sample_ids <- colnames(gt_imp)
  
  # REF / ALT from THIS VCF
  refalt_cov <- tibble(
    CHR_POS = rownames(vcf_obj@gt),
    REF     = vcf_obj@fix[, "REF"],
    ALT     = vcf_obj@fix[, "ALT"]
  )
  
  refalt_cov <- refalt_cov[match(target_chrpos, refalt_cov$CHR_POS), ]
  stopifnot(all(refalt_cov$CHR_POS == target_chrpos))
  
  for (s in sample_ids) {
    
    # Optional coverage-aware haplotype naming
    pref <- if (ADD_COVERAGE_TO_HAPLO_NAME) {
      paste0(s, "_", cov)
    } else {
      s
    }
    
    # ---------- HAP_A ----------
    haplo_vecs[[paste0(pref, "_HAP_A")]] <- vapply(
      seq_along(target_chrpos),
      function(i) {
        haplotype_to_state(
          gt               = gt_imp[target_chrpos[i], s],
          phase            = 1,
          ref              = refalt_cov$REF[i],
          alt              = refalt_cov$ALT[i],
          ancestral_allele = filtered_sites$Ancestral_Allele[i],
          derived_allele   = filtered_sites$Derived_Allele[i]
        )
      },
      character(1)
    )
    
    # ---------- HAP_B ----------
    haplo_vecs[[paste0(pref, "_HAP_B")]] <- vapply(
      seq_along(target_chrpos),
      function(i) {
        haplotype_to_state(
          gt               = gt_imp[target_chrpos[i], s],
          phase            = 2,
          ref              = refalt_cov$REF[i],
          alt              = refalt_cov$ALT[i],
          ancestral_allele = filtered_sites$Ancestral_Allele[i],
          derived_allele   = filtered_sites$Derived_Allele[i]
        )
      },
      character(1)
    )
  }
}


allele_states_imputed <- do.call(rbind, haplo_vecs)
colnames(allele_states_imputed) <- target_chrpos

# Sanity check: dimensions
print(dim(allele_states_imputed))
print(table(allele_states_imputed, useNA = "ifany"))
cat("Imputed allele-state matrices:\n")
print(sapply(allele_states_imputed, dim))

#============================================================#
# Merge allele-state matrices for visualization
#============================================================#
# At this stage, we have already computed allele-state matrices
# for all relevant data types:
#
# - Chimp: synthetic ancestral reference (1 row)
# - Africans: collapsed control group (1 row)
# - Archaic individuals: Altai, Vindija, Denisova, Chagyrskaya
# - Original ancient genomes: unphased diploid genotypes
# - Imputed ancient genomes: phased haplotypes (A/B)
#
# All matrices:
#   - have rows = individuals / haplotypes
#   - have columns = target_chrpos
#   - are guaranteed to be aligned to target_chrpos
#
# We now concatenate them into a single matrix that will be
# used exclusively for visualization.
#============================================================#

allele_state_matrix <- do.call(
  rbind,
  list(
    allele_states_chimp,
    allele_states_africans_collapsed,
    allele_states_archaic,
    allele_states_original,
    allele_states_imputed
  )
)

# Final safety check: all columns must match the target SNP order
stopifnot(identical(colnames(allele_state_matrix), target_chrpos))


#============================================================#
# Convert allele-state matrix to long format (ggplot-ready)
#============================================================#
# The heatmap is produced using ggplot2, which expects data in
# long (tidy) format:
#
#   Individual | CHR_POS | Allele_State
#
# Each row represents the allele state of one individual (or
# haplotype) at one archaic-derived SNP.
#============================================================#

allele_df_long <- as.data.frame(allele_state_matrix) %>%
  mutate(Individual = rownames(.)) %>%
  pivot_longer(
    cols      = -Individual,
    names_to  = "CHR_POS",
    values_to = "Allele_State"
  )

# Explicit factor levels ensure consistent color mapping
# across plots and datasets
allele_df_long$Allele_State <- factor(
  allele_df_long$Allele_State,
  levels = c("Ancestral", "Derived", "Heterozygous", "Different", "Missing")
)


#============================================================#
# Define row ordering for visualization
#============================================================#
# We define an explicit ordering of rows (Individuals) so that
# the heatmap has a biologically meaningful structure.
#
# Conceptual order (top → bottom before reversal):
#   1) Chimp (ancestral reference)
#   2) Africans (collapsed control)
#   3) Archaic humans
#   4) Original (unphased) ancient genomes
#   5) Imputed (phased) haplotypes, grouped by coverage
#
# IMPORTANT:
# This block assumes that coverage labels (e.g. "_ImputedOC_",
# "_2X_") are present in Individual names. If this is not true
# (ADD_COVERAGE_TO_HAPLO_NAME = FALSE), this logic must be
# adapted to use an explicit metadata table instead.
#============================================================#

coverage_order <- c("Original", "ImputedOC", "2X", "1X", "0.5X", "0.0625X")

allele_df_long <- allele_df_long %>%
  mutate(
    Coverage = case_when(
      Individual == "Chimp"            ~ "Reference",
      Individual == "Africans"         ~ "Reference",
      Individual %in% archaic_ids      ~ "Reference",
      grepl("_Original$", Individual)  ~ "Original",
      grepl("_ImputedOC_", Individual) ~ "ImputedOC",
      grepl("_2X_", Individual)        ~ "2X",
      grepl("_1X_", Individual)        ~ "1X",
      grepl("_0.5X_", Individual)      ~ "0.5X",
      grepl("_0.0625X_", Individual)   ~ "0.0625X",
      TRUE                             ~ NA_character_
    )
  )

allele_df_long$Coverage <- factor(
  allele_df_long$Coverage,
  levels = c("Reference", coverage_order)
)

# Build final Individual ordering:
# - Reference groups first (Chimp, Africans, Archaics)
# - Then Original and Imputed, ordered by coverage
# - Within each group, keep Individuals together
individual_order <- allele_df_long %>%
  distinct(Individual, Coverage) %>%
  arrange(
    case_when(
      Individual == "Chimp"        ~ 1,
      Individual == "Africans"     ~ 2,
      Individual %in% archaic_ids  ~ 3,
      TRUE                         ~ 4
    ),
    Coverage,
    Individual
  ) %>%
  pull(Individual)

# Reverse order so the first entries appear at the top of the heatmap
allele_df_long$Individual <- factor(
  allele_df_long$Individual,
  levels = rev(unique(individual_order))
)

#------------------------------------------------------------#
# Define color scheme
#------------------------------------------------------------#

allele_colors <- c(
  "Ancestral"    = "white",
  "Derived"      = "darkgreen",
  "Heterozygous" = "orange",
  "Different"    = "red",
  "Missing"      = "grey80"
)

#------------------------------------------------------------#
# Draw heatmap
#------------------------------------------------------------#

heatmap_plot <- ggplot(
  allele_df_long,
  aes(x = CHR_POS, y = Individual, fill = Allele_State)
) +
  geom_tile(color = "grey70", linewidth = 0.05) +
  scale_fill_manual(values = allele_colors, name = NULL) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 8),
    legend.position = "right",
    panel.grid = element_blank()
  )

print(heatmap_plot)

#------------------------------------------------------------#
# Annotate Individuals as European vs non-European
#------------------------------------------------------------#
# Population assignment is derived from metadata and used
# only for visualization stratification.
#------------------------------------------------------------#

allele_df_long_Pop <- allele_df_long %>%
  mutate(
    # Remove haplotype suffix to match metadata IDs
    Individual_key = sub("_HAP_[AB]$", "", Individual)
  ) %>%
  left_join(
    metadata_info %>% select(Individual, Population),
    by = c("Individual_key" = "Individual")
  ) %>%
  mutate(
    Region_Group = case_when(
      Individual %in% c("Chimp", "Africans") ~ "Reference",
      Individual %in% archaic_ids            ~ "Reference",
      Population %in% c("European Mesolithic", "European Neolithic") ~ "European",
      TRUE                                   ~ "Non-European"
    )
  )

# Sanity checks
table(allele_df_long_Pop$Region_Group, useNA = "ifany")
table(is.na(allele_df_long_Pop$Population))

#------------------------------------------------------------#
# Subset to Europeans + references only
#------------------------------------------------------------#

allele_df_long_Pop_European <- allele_df_long_Pop %>%
  filter(Region_Group != "Non-European")

# Re-enforce Individual ordering after filtering
allele_df_long_Pop_European$Individual <- factor(
  allele_df_long_Pop_European$Individual,
  levels = levels(allele_df_long$Individual)
)

heatmap_plot_EU <- ggplot(
  allele_df_long_Pop_European,
  aes(x = CHR_POS, y = Individual, fill = Allele_State)
) +
  geom_tile(color = "grey70", linewidth = 0.05) +
  scale_fill_manual(values = allele_colors, name = NULL) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 8),
    legend.position = "right",
    panel.grid = element_blank()
  )

print(heatmap_plot_EU)

#------------------------------------------------------------#
# Session information (reproducibility)
#------------------------------------------------------------#
sessionInfo()
