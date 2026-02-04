# ================================================================
# Segment-level QC of imputed genotypes in archaic regions
# Publishable script header: packages, palettes, ordering, paths
# ================================================================

# ---- Package loading (fail fast with clear error) ------------------------
pkgs <- c(
  "dplyr", "tidyr", "ggplot2", "ggrain", "patchwork", "ggpubr",
  "purrr", "scales", "ggdist", "data.table"
)

missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them with install.packages(...) (or BiocManager::install(...) if needed)."
  )
}

# Attach packages (explicit order helps reproducibility)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrain)
library(patchwork)
library(ggpubr)
library(purrr)
library(scales)
library(ggdist)
library(data.table)

# ---- Aesthetics: palettes + plotting order -------------------------------
coverage_colors <- c(
  Original   = "#466F9D",
  ImputedOC  = "#00A087",
  `2X`       = "#D55E00",
  `1X`       = "#E69F00",
  `0.5X`     = "#56B4E9",
  `0.0625X`  = "#C77BA7"
)

coverage_order <- c("Original", "ImputedOC", "2X", "1X", "0.5X", "0.0625X")

individual_order <- c(
  "Ust_Ishim", "Yana1", "USR1", "Kolyma1", "KK1", "WC1", "sf12",
  "Loschbour", "Stuttgart", "NE1", "SRA62", "JP14", "BOT2016",
  "PB675", "2H10", "2H11", "atp016", "Yamnaya", "BA64", "Mota"
)

genotype_colors <- c(
  HomozygousReference    = "#677A8E",
  Heterozygous           = "#1B9E77",
  HomozygousAlternative  = "#D95F02"
)

# ---- Sanity checks (prevents silent mistakes in figures) -----------------
stopifnot(all(coverage_order %in% names(coverage_colors)))
if (length(names(coverage_colors)) != length(unique(names(coverage_colors)))) {
  stop("Duplicate names in coverage_colors.")
}
if (any(duplicated(individual_order))) {
  stop("Duplicate entries detected in individual_order.")
}

# # ---- Paths (publishable: avoid relying only on setwd) --------------------
# # Preferred: define a project root once and build paths from it.
# # If you keep setwd(), at least do it once here and keep all paths relative.
# 
# project_dir <- "C:/Users/capodifm/OneDrive - Trinity College Dublin/ArchaicInference_Manuscript/Figures/Figures2025/Script/Ready/CheckBam"
# 
# if (!dir.exists(project_dir)) {
#   stop("Project directory does not exist:\n", project_dir)
# }
# setwd(project_dir)
# 
# message("Working directory: ", getwd())

# -------------------------------------------------------------------
# Load per-individual published read-count tables
# -------------------------------------------------------------------
# Each file corresponds to one individual and contains per-site
# read support for imputed variants absent or missing in the Original VCF.

files_raw <- list.files(
  "../Data/BamCounts/",
  pattern = "_published\\.tsv\\.gz$",
  full.names = TRUE
)

dfs_raw <- lapply(
  files_raw,
  function(f) read.delim(gzfile(f), header = TRUE, stringsAsFactors = FALSE)
)

# Name list elements by file name (used later for provenance if needed)
names(dfs_raw) <- basename(files_raw)

# -------------------------------------------------------------------
# Per-individual processing function
# -------------------------------------------------------------------
# This function:
#   1) Interprets the imputed genotype (GT) into biological classes
#   2) Computes total read depth and allele proportions
#   3) Returns:
#        - a site-level ("wide") table
#        - an allele-level ("long") table
# -------------------------------------------------------------------

process_one <- function(df) {
  
  # ---------------------------------------------------------------
  # A. Interpret genotypes and compute read-based quantities
  # ---------------------------------------------------------------
  # - Normalise GT format (e.g. 0/1 -> 0|1)
  # - Map GT to HomozygousReference / Heterozygous / HomozygousAlternative
  # - Compute:
  #     dp_refalt   = ref_count + alt_count
  #     prop_ref_ra = ref_count / (ref + alt)
  #     prop_alt_ra = alt_count / (ref + alt)
  
  df_wide <- df %>%
    mutate(
      GT = as.character(GT),
      GT_norm = gsub("\\s*", "", GT),
      GT_norm = gsub("/", "|", GT_norm),
      
      Genotype = case_when(
        GT_norm == "0|0"           ~ "HomozygousReference",
        GT_norm == "1|1"           ~ "HomozygousAlternative",
        GT_norm %in% c("0|1","1|0") ~ "Heterozygous",
        TRUE                       ~ NA_character_
      ),
      
      dp_refalt   = ref_count + alt_count,
      prop_ref_ra = if_else(dp_refalt > 0, ref_count / dp_refalt, NA_real_),
      prop_alt_ra = if_else(dp_refalt > 0, alt_count  / dp_refalt, NA_real_)
    ) %>%
    select(-GT_norm)
  
  # ---------------------------------------------------------------
  # B. Convert to allele-level (long) format
  # ---------------------------------------------------------------
  # Each site becomes two rows:
  #   - one for REF allele
  #   - one for ALT allele
  # This is convenient for:
  #   - raincloud plots
  #   - distribution-based visualisations
  
  df_long <- df_wide %>%
    select(
      Individual,
      cov,
      Genotype,
      dp_refalt,
      prop_ref_ra,
      prop_alt_ra
    ) %>%
    pivot_longer(
      cols      = c(prop_ref_ra, prop_alt_ra),
      names_to  = "allele",
      values_to = "prop"
    ) %>%
    mutate(
      allele = recode(
        allele,
        prop_ref_ra = "REF",
        prop_alt_ra = "ALT"
      )
    )
  
  list(
    wide = df_wide,
    long = df_long
  )
}

# -------------------------------------------------------------------
# Apply processing to all individuals
# -------------------------------------------------------------------
# We keep results as a list so we can:
#   - combine later
#   - keep provenance if needed

res_list <- imap(dfs_raw, function(df, fname) {
  
  out <- process_one(df)
  
  # Keep track of source file (optional but useful)
  out$wide$source <- fname
  out$long$source <- fname
  
  out
})


# -------------------------------------------------------------------
# Combine processed results across individuals
# -------------------------------------------------------------------
# Each element of `res_list` contains:
#   - $wide : site-level information
#   - $long : allele-level information
# Here we bind them across individuals to obtain global tables.

wide_all <- dplyr::bind_rows(lapply(res_list, `[[`, "wide"))
long_all <- dplyr::bind_rows(lapply(res_list, `[[`, "long"))
long_all_NA <- long_all %>%
  filter(allele != "REF") %>%
  mutate(prop = if_else(is.na(prop), 0, prop))

rm(long_all)
rm(res_list)
gc()


# ================================================================
# Publishable segment-level false-positive check in archaic-new tracts
# Denominator: all sites missing from Original (wide_all)
# Numerator: ALT-genotype sites with <10% ALT support
#   - with reads: dp_refalt > 0
# Output: per-segment percentages + raincloud plot
# ================================================================
# ----------------------------------------------------------------
# Load and filter introgressed segments (Archaic New, Kept)
# ----------------------------------------------------------------

seg_path <- "C:/Users/capodifm/Desktop/ArchaicIntrogression/1.TestImputation/2.Tracks/Tracts/NewTractsNoLiftOver_Jun27/NewSimilarity_Jan2025/FileForSegmentAnalyses_March19/MergeFragmentsArchaic_SimilarityDistance_Sharing_Filtering_Masking_Sep04.csv"

segments_raw <- read.csv(seg_path, header = TRUE)

segments_ArchaicNew <- segments_raw %>%
  filter(Filter == "Kept", segment_type == "Archaic New")

str(segments_ArchaicNew)
seg <- as.data.table(segments_ArchaicNew)[
  , .(
    chr        = as.character(Chrom),
    start      = as.integer(start),
    end        = as.integer(end),
    Individual = as.character(individual),
    cov        = as.character(Coverage),
    Length     = as.integer(Length)
  )
]
seg[, segment_id := .I]
setkey(seg, chr, start, end)

if (nrow(seg) == 0) stop("No segments after filtering (Kept + Archaic New).")

# ----------------------------
# Build position tables (unique sites)
#    NOTE: wide_all is assumed to already be "missing-from-original" sites
# ----------------------------
# Helper: normalize chromosome naming
chr_norm <- function(x) {
  sub("^chr", "", as.character(x))
}

prep_pos <- function(dt) {
  x <- as.data.table(dt)[
    , .(
      Individual = as.character(Individual),
      cov        = as.character(cov),
      chr        = chr_norm(chr),
      pos        = as.integer(pos)
    )
  ]
  unique(x)
}

# Denominator: all missing sites
pos_all <- prep_pos(wide_all)

# Numerator (prop only): ALT genotypes with prop_alt_ra < 0.1
Prop_0.1_wide <- wide_all %>%
  filter(
    Genotype %in% c("HomozygousAlternative", "Heterozygous"),
    !is.na(prop_alt_ra),
    prop_alt_ra < 0.10
  )

# Numerator (prop + depth): ALT genotypes with reads and low ALT support
Prop_0.1_wide_Depth <- wide_all %>%
  filter(
    Genotype %in% c("HomozygousAlternative", "Heterozygous"),
    !is.na(prop_alt_ra),
    prop_alt_ra < 0.10,
    !is.na(dp_refalt),
    dp_refalt > 0
  )

pos_prop  <- prep_pos(Prop_0.1_wide)
pos_depth <- prep_pos(Prop_0.1_wide_Depth)

# Convert points to intervals for foverlaps
to_i <- function(p) setDT(copy(p))[, `:=`(start = pos, end = pos)]

i_all   <- to_i(pos_all);   setkey(i_all,   chr, start, end)
i_prop  <- to_i(pos_prop);  setkey(i_prop,  chr, start, end)
i_depth <- to_i(pos_depth); setkey(i_depth, chr, start, end)

# ----------------------------
# Efficient overlap: split by Individual x cov BEFORE foverlaps
#    (prevents huge intermediate joins)
# ----------------------------
join_and_count_split <- function(i_dt, seg_dt, count_name) {
  
  # split segments by Individual+cov
  seg_split <- split(seg_dt, by = c("Individual", "cov"), keep.by = TRUE)
  
  res <- lapply(seg_split, function(s) {
    ind <- s$Individual[1]
    cv  <- s$cov[1]
    
    i_sub <- i_dt[Individual == ind & cov == cv]
    if (nrow(i_sub) == 0) return(NULL)
    
    ov <- foverlaps(i_sub, s, nomatch = 0L)
    if (nrow(ov) == 0) return(NULL)
    
    ov[, .N, by = segment_id]
  })
  
  out <- rbindlist(res, fill = TRUE)
  if (nrow(out) == 0) out <- data.table(segment_id = integer(), N = integer())
  
  setnames(out, "N", count_name)
  out
}

cnt_total <- join_and_count_split(i_all,   seg, "n_total")
cnt_prop  <- join_and_count_split(i_prop,  seg, "n_lowALT")
cnt_depth <- join_and_count_split(i_depth, seg, "n_lowALT_depth")

# ----------------------------
# Merge + compute percentages
# ----------------------------
seg_counts <- copy(seg)[
  cnt_total, on = "segment_id"
][
  cnt_prop,  on = "segment_id"
][
  cnt_depth, on = "segment_id"
][
  , `:=`(
    n_total        = fifelse(is.na(n_total),        0L, n_total),
    n_lowALT       = fifelse(is.na(n_lowALT),       0L, n_lowALT),
    n_lowALT_depth = fifelse(is.na(n_lowALT_depth), 0L, n_lowALT_depth)
  )
][
  , `:=`(
    pct_lowALT       = fifelse(n_total > 0, 100 * n_lowALT       / n_total, NA_real_),
    pct_lowALT_depth = fifelse(n_total > 0, 100 * n_lowALT_depth / n_total, NA_real_)
  )
]

# enforce plotting order (and drop unused levels cleanly)
seg_counts <- seg_counts %>%
  mutate(
    Individual = factor(Individual, levels = individual_order),
    cov        = factor(cov,        levels = coverage_order)
  )

# Optional: remove segments with no denominator (rare but can happen)
# seg_counts <- seg_counts %>% filter(!is.na(pct_lowALT_depth))

# ----------------------------
# Raincloud 
# ----------------------------
raincloud_PrcSNPs_Depth_MultipleSamples <- ggplot(
  seg_counts,
  aes(x = 1, y = pct_lowALT_depth, fill = cov, color = cov)
) +
  geom_rain(
    alpha = 0.40,
    boxplot.args = list(color = "black", outlier.shape = NA)
  ) +
  facet_grid(rows = vars(cov), cols = vars(Individual), drop = TRUE) +
  scale_fill_manual(values = coverage_colors, name = "Coverage") +
  scale_color_manual(values = coverage_colors, name = "Coverage") +
  labs(
    x = NULL,
    y = "Per-segment % of missing sites with <10% ALT support (≥1 read)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    legend.position    = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background   = element_rect(fill = "grey95", color = NA),
    strip.text         = element_text(size = 9)
  )

raincloud_PrcSNPs_Depth_MultipleSamples

#Figure
# For imputed variants absent from the original VCF, the proportion of reads supporting the alternative allele increases with read depth
# and converges to the expected values for heterozygous and homozygous alternative genotypes,
# indicating that most imputed sites are supported by sequencing reads even at moderate coverage.
# -------------------------------------------------------------------
# Treat missing read support as zero
# -------------------------------------------------------------------
# For some imputed sites, no reads overlap the position in the BAM.
# We explicitly convert NA allele proportions to 0 so that:
#   - depth = 0 sites are retained
#   - they are visualised as zero support for the ALT allele

depth_bins <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 20, 50, 100, Inf)
depth_labels <- c(
  "0","1","2","3","4","5","6","7","8","9","10",
  "(11,20)","(20,50)","(50,100)","100+"
)

# -------------------------------------------------------------------
# Filter to ALT allele and non-reference genotypes
# -------------------------------------------------------------------
# We focus on:
#   - ALT allele rows only
#   - Sites where the imputed genotype is not HomozygousReference
# This removes trivial cases where ALT support is not expected.

df_all_NA <- long_all_NA %>%
  filter(allele == "ALT") %>%
  filter(Genotype != "HomozygousReference") %>%
  mutate(
    depth_bin = cut(
      dp_refalt,
      breaks = depth_bins,
      labels = depth_labels,
      right = FALSE,
      include.lowest = TRUE
    )
  )

# -------------------------------------------------------------------
# Summarise ALT allele proportion as a function of read depth
# -------------------------------------------------------------------
# For each combination of:
#   Individual × Coverage × Genotype × ReadDepth
# we compute:
#   - median ALT allele proportion
#   - interquartile range
#   - number of sites contributing

sum_by_dp_NA <- df_all_NA %>%
  filter(!is.na(prop), dp_refalt <= 40) %>%
  mutate(dp_refalt = as.integer(dp_refalt)) %>%
  group_by(Individual, cov, Genotype, dp_refalt) %>%
  summarise(
    n   = n(),
    q25 = quantile(prop, 0.25, na.rm = TRUE),
    med = median(prop, na.rm = TRUE),
    q75 = quantile(prop, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# -------------------------------------------------------------------
# Set factor levels for consistent plotting
# -------------------------------------------------------------------

sum_by_dp_NA <- sum_by_dp_NA %>%
  mutate(
    cov        = factor(cov, levels = coverage_order),
    Individual = factor(Individual, levels = individual_order),
    Genotype   = factor(
      Genotype,
      levels = c("Heterozygous", "HomozygousAlternative")
    )
  )

# -------------------------------------------------------------------
# Plot ALT allele proportion vs read depth
# -------------------------------------------------------------------

ProportionVSReads <- ggplot(
  sum_by_dp_NA,
  aes(
    x = dp_refalt,
    y = med,
    color = Genotype,
    fill  = Genotype,
    group = Genotype
  )
) +
  geom_ribbon(aes(ymin = q25, ymax = q75),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.6) +
  geom_point(aes(size = n),
             alpha = 0.6,
             show.legend = FALSE) +
  facet_grid(Individual ~ cov) +
  scale_color_manual(values = genotype_colors, name = "Genotype") +
  scale_fill_manual(values = genotype_colors, guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(limits = c(0, 40), breaks = seq(0, 40, 5)) +
  labs(
    x = "Number of reads",
    y = "Proportion of reads supporting the alternative allele"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    strip.text.y = element_text(angle = 0),
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7)
  )

# -------------------------------------------------------------------
# Filter to ALT allele and non-reference genotypes
# -------------------------------------------------------------------
# Raincloud plots show the distribution of alternative allele support for imputed variants absent from the original VCF.
# Across individuals and coverage levels, heterozygous and homozygous alternative genotypes
# display allele proportion distributions centred around their expected values,
# with increased dispersion at lower coverage.
# We restrict the analysis to:
#   - ALT allele rows only
#   - Sites with imputed genotypes that are not HomozygousReference
# This ensures that ALT allele support is biologically meaningful.

het_alt_NA <- long_all_NA %>%
  filter(allele == "ALT", Genotype != "HomozygousReference") %>%
  mutate(
    cov     = factor(cov, levels = coverage_order),
    Sample  = factor(Individual, levels = individual_order)
  )

# -------------------------------------------------------------------
# Plot ALT allele proportion distributions
# -------------------------------------------------------------------
# Raincloud plots summarise the distribution of ALT allele support
# across all imputed sites, stratified by genotype, coverage, and individual.

raincloud_MultipleSamples_NA <- ggplot(
  het_alt_NA,
  aes(1, prop, fill = Genotype, color = Genotype)
) + 
  geom_rain(
    alpha = 0.4,
    point.args   = list(size = 0, alpha = 0),
    boxplot.args = list(color = "black", outlier.shape = NA)
  ) +
  labs(
    x = "",
    y = "Proportion of reads supporting the alternative allele"
  ) +
  facet_grid(cov ~ Sample) +
  scale_fill_manual(values = genotype_colors, name = "Genotype") +
  scale_color_manual(values = genotype_colors, name = "Genotype") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"
  )

# For each individual and coverage, we quantified the fraction of imputed variants absent from the original VCF
# that show either no sequencing support or weak alternative allele support.
# Because sites with zero reads necessarily have an alternative allele proportion of zero, low allele support includes both low-depth and no-depth sites.
# Conditional probabilities further disentangle the contribution of read depth to low alternative allele support.
# -------------------------------------------------------------------
# INPUT
# - long_all_NA: long-format table with one row per allele per site
#
# ASSUMPTIONS
# - Only ALT rows are used (one row per site after filtering)
# - Sites with dp_refalt = 0 have prop = 0 by construction
# - Therefore, ALT < 10% includes both low-support and zero-depth sites
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Summarise ALT allele support per Individual × coverage
# -------------------------------------------------------------------

summary_by_ind_cov <- long_all_NA %>%
  filter(allele == "ALT", Genotype != "HomozygousReference") %>%
  group_by(Individual, cov) %>%
  summarise(
    # (A) Total number of imputed sites
    n_total = n(),
    
    # (B) Sites with no sequencing reads
    n_dp0    = sum(dp_refalt == 0),
    prop_dp0 = n_dp0 / n_total,
    
    # (C) Sites with low ALT allele support
    # NOTE: includes dp_refalt = 0
    n_prop10 = sum(prop < 0.10),
    prop_p10 = n_prop10 / n_total,
    
    # (D) Conditional on low ALT support
    prop_dp0_given_low  = if_else(
      n_prop10 > 0,
      sum(prop < 0.10 & dp_refalt == 0) / n_prop10,
      NA_real_
    ),
    prop_dp10_given_low = if_else(
      n_prop10 > 0,
      sum(prop < 0.10 & dp_refalt <= 10) / n_prop10,
      NA_real_
    ),
    
    # (E) Conditional on sites with reads
    n_gt0        = sum(dp_refalt > 0),
    prop_p10_gt0 = if_else(
      n_gt0 > 0,
      sum(dp_refalt > 0 & prop < 0.10) / n_gt0,
      NA_real_
    ),
    prop_alt0_gt0 = if_else(
      n_gt0 > 0,
      sum(dp_refalt > 0 & prop == 0) / n_gt0,
      NA_real_
    ),
    
    .groups = "drop"
  )
# -------------------------------------------------------------------
# Reshape summary table for plotting
# - Original coverage is excluded
# - Metrics include both unconditional and conditional probabilities
# -------------------------------------------------------------------

coverage_order2 <- setdiff(coverage_order, "Original")

plot_df <- summary_by_ind_cov %>%
  filter(cov != "Original") %>%
  mutate(
    cov        = factor(cov, levels = coverage_order2),
    Individual = factor(Individual, levels = individual_order)
  ) %>%
  select(
    Individual, cov,
    prop_dp0, prop_p10, prop_p10_gt0, prop_alt0_gt0,
    prop_dp0_given_low, prop_dp10_given_low
  ) %>%
  pivot_longer(
    c(prop_dp0, prop_p10, prop_p10_gt0, prop_alt0_gt0,
      prop_dp0_given_low, prop_dp10_given_low),
    names_to = "metric",
    values_to = "prop"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c(
        "prop_dp0",
        "prop_p10",
        "prop_p10_gt0",
        "prop_alt0_gt0",
        "prop_dp0_given_low",
        "prop_dp10_given_low"
      ),
      labels = c(
        "Depth = 0 (overall)",
        "ALT < 10% (overall)",
        "ALT < 10% (depth > 0)",
        "ALT reads = 0 (depth > 0)",
        "Depth = 0 | ALT < 10%",
        "Depth ≤ 10 | ALT < 10%"
      )
    )
  )
# -------------------------------------------------------------------
# Plot: Impact of read depth on ALT allele support
# -------------------------------------------------------------------

metric_colors <- c(
  "Depth = 0 (overall)"           = "#8DA0CB",
  "ALT < 10% (overall)"           = "#66C2A5",
  "ALT < 10% (depth > 0)"         = "#1B9E77",
  "ALT reads = 0 (depth > 0)"     = "#FFD92F",
  "Depth = 0 | ALT < 10%"         = "#FC8D62",
  "Depth ≤ 10 | ALT < 10%"        = "#E78AC3"
)

AltAllele_Impact <- ggplot(
  plot_df,
  aes(x = metric, y = prop, fill = metric)
) +
  geom_col(width = 0.75) +
  geom_text(
    aes(y = 0.5, label = scales::percent(prop, accuracy = 0.1)),
    size = 2.8
  ) +
  facet_grid(Individual ~ cov, drop = FALSE) +
  scale_fill_manual(values = metric_colors, guide = "none") +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(x = NULL, y = "Proportion of imputed sites") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x  = element_text(angle = 30, hjust = 1),
    strip.text.y = element_text(angle = 0)
  )

# ================================================================
# Summary statistics for raincloud plots
# Excludes Original; reports per-segment distributions
# ================================================================

# Use the same data frame used for plotting
plot_df_clean <- plot_df %>%
  filter(cov != "Original") %>%          # Original has no imputation
  filter(!is.na(prop)) %>%               # explicit NA handling
  mutate(
    prop = as.numeric(prop)              # defensive: ensure numeric
  )

# Sanity check: prop should be a proportion
stopifnot(all(plot_df_clean$prop >= 0 & plot_df_clean$prop <= 1))

# ----------------------------
# Overall statistics (all individuals & coverages pooled)
# ----------------------------
stats_overall <- plot_df_clean %>%
  group_by(metric) %>%
  summarise(
    n_segments = n(),
    min        = min(prop),
    q25        = quantile(prop, 0.25),
    median     = median(prop),
    q75        = quantile(prop, 0.75),
    max        = max(prop),
    .groups    = "drop"
  )

# ----------------------------
# Statistics stratified by coverage
# ----------------------------
stats_by_cov <- plot_df_clean %>%
  group_by(cov, metric) %>%
  summarise(
    n_segments = n(),
    min        = min(prop),
    q25        = quantile(prop, 0.25),
    median     = median(prop),
    q75        = quantile(prop, 0.75),
    max        = max(prop),
    .groups    = "drop"
  )

# ----------------------------
# Format for presentation (percent scale)
# ----------------------------
fmt_pct <- function(x) scales::percent(x, accuracy = 0.1)

stats_overall_print <- stats_overall %>%
  mutate(across(c(min, q25, median, q75, max), fmt_pct))

stats_by_cov_print <- stats_by_cov %>%
  mutate(across(c(min, q25, median, q75, max), fmt_pct))

# Display tables
stats_overall_print
stats_by_cov_print


