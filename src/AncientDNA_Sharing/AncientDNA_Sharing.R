# ============================================================
#  LIBRARIES
# ============================================================

# Data manipulation
library(dplyr)
library(tidyr)
library(tibble)
# Genomic interval handling
library(GenomicRanges)

# Visualization
library(ggplot2)
library(ggrepel)
library(patchwork)
library(pheatmap)
library(viridis)

# ================================================================
#  Multi-coverage analysis of archaic introgression sharing
# ================================================================
#
# This script performs a multi-intersection analysis of archaic
# introgressed genomic fragments across ancient individuals,
# stratified by sequencing coverage or imputation level.
#
# The analysis proceeds in three logically separated stages:
#
#   (1) DATA TRANSFORMATION
#       - Introgressed segments are split into minimal,
#         non-overlapping genomic intervals defined by shared
#         breakpoints across individuals.
#       - Presence or absence of introgression in each interval
#         is encoded as a binary matrix (regions × individuals).
#
#   (2) QUANTIFICATION
#       - Pairwise sharing between individuals is quantified as
#         the number of genomic intervals jointly introgressed.
#       - Sharing matrices are row-normalized to account for
#         individual-specific introgression burden.
#       - Multidimensional Scaling (MDS) coordinates are derived
#         from the normalized sharing matrix.
#
#   (3) VISUALIZATION
#       - MDS plots and heatmaps are generated exclusively from
#         precomputed tables.
#       - Numerical computation and plotting are fully decoupled
#         to ensure reproducibility and flexibility.
#
# Notes:
#   - Introgression is treated as binary (segment length is not
#     used as a weight).
#   - Coverage-specific results are always ordered according to
#     the user-defined 'coverage_values' vector.
#
# This script is designed to be publication-ready, modular,
# and easily extensible for additional coverage levels or
# visualization strategies.
# ================================================================

# ============================================================
#  HELPERS
# ============================================================

# Colors / orders
dataset_colors <- c(
  Original = "#466F9D",
  ImputedOC = "#00A087",
  `2X` = "#D55E00",
  `1X` = "#E69F00",
  `0.5X` = "#56B4E9",
  `0.0625X` = "#C77BA7"
)
dataset_order <- c("Original", "ImputedOC", "2X","1X","0.5X", "0.0625X")
coverage_values <- c("0.0625X", "0.5X", "1X", "2X", "Original", "ImputedOC")
imputed_coverages <- c("ImputedOC", "0.0625X", "0.5X", "1X", "2X")

individual_order <- c(
  "Ust_Ishim", "Yana1", "USR1", "Kolyma1", "KK1", "WC1", "sf12",
  "Loschbour", "Stuttgart", "NE1", "SRA62", "JP14", "BOT2016",
  "PB675", "2H10", "2H11", "atp016", "Yamnaya", "BA64", "Mota"
)

# Create a named vector of colors for each population
population_colors <- c("European Neolithic" =	"#8d2962ff",
                       "Steppe" = "#1b998bff"	,
                       "Caucasus Hunter-Gatherer" = "#6581e7"	,
                       "Ancient Siberian" = "#8c3821"	,
                       "European Mesolithic" =   "#7e5920ff",
                       "African Stone Age" = "#0c50b7" ,	
                       "Ancient Berinigian" = "#ffa737ff" ,
                       "Pre-LGM Eurasian" = "#ed217cff" ,	
                       "Iranian Neolithic" = "#008d00"
)

# ============================================================
#  INPUT DATA
# ============================================================
#
# metadata_info:
#   Sample-level metadata including population assignment
#   and color mapping used for visualization.
#
# IntrogressionMaps:
#   Genomic coordinates of introgressed segments inferred
#   for each individual and coverage level.

metadata_info <- as.data.frame(read.delim("../Data/MetadataInfo_ColorsPop.csv", header = TRUE, sep = ","))

IntrogressionMaps <- as.data.frame(read.delim("../Data/IntrogressionMaps.csv", header = TRUE, sep = ","))
IntrogressionMaps_Filtered <- IntrogressionMaps %>%
  filter(Filter == "Kept")

# ============================================================
#  PREPARATION OF INPUT DATA
# ============================================================

# Create a unique identifier for each individual–coverage pair
# (e.g. Ust_Ishim_Original, Ust_Ishim_1X, ...)
# This identifier will later be used as column names in the
# binary presence/absence matrix.
IntrogressionMaps_Filtered_ID <- 
  IntrogressionMaps_Filtered %>%
  mutate(Individual_selection = paste0(Individual, "_", Coverage))

# Remove the African reference individual (Mota) from the analysis
# to avoid influencing overlap frequencies and distances.
IntrogressionMaps_Filtered_ID_noMota <- 
  IntrogressionMaps_Filtered_ID %>%
  filter(Individual != "Mota")


# ============================================================
#  COVERAGE-SPECIFIC GENOMIC RANGES
# ============================================================
#
# For each coverage level, introgressed segments are converted
# into a GRanges object. The individual–coverage identifier
# is stored as metadata and used downstream to track sharing
# across individuals.

build_granges_for_coverage <- function(IntrogressionMaps, coverage_name) {
  
  # Ensure correct input type
  stopifnot(is.data.frame(IntrogressionMaps))
  
  # Required columns for downstream processing
  required_cols <- c("chr", "start", "end", "Coverage", "Individual_selection")
  missing_cols  <- setdiff(required_cols, colnames(IntrogressionMaps))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns in IntrogressionMaps: ",
         paste(missing_cols, collapse = ", "))
  }
  
  # Subset to the selected coverage
  df_cov <- IntrogressionMaps %>%
    dplyr::filter(Coverage == coverage_name)
  
  # If no segments are found, return NULL so the wrapper can skip
  if (nrow(df_cov) == 0) {
    message("⚠️  No segments found for coverage: ", coverage_name)
    return(NULL)
  }
  
  # Sanity check on genomic coordinates
  if (any(df_cov$start >= df_cov$end)) {
    stop("Invalid ranges detected (start >= end) for coverage: ", coverage_name)
  }
  
  # Convert to GRanges
  gr <- GRanges(
    seqnames = df_cov$chr,
    ranges   = IRanges(start = df_cov$start, end = df_cov$end),
    Individual_Coverage = df_cov$Individual_selection
  )
  
  message("✔️  Built GRanges for ", coverage_name,
          " (n = ", length(gr), " segments)")
  
  gr
}

# ============================================================
#  MULTI-INDIVIDUAL INTERSECTION ATOMIC BINS
# ============================================================
#
# This function decomposes introgressed segments into minimal
# non-overlapping genomic intervals ("atomic bins") defined by
# shared start and end breakpoints.
#
# For each atomic bin, it records:
#   - the number of individuals carrying introgression
#   - the list of individuals overlapping that bin
#
# This representation is the basis for constructing the
# binary presence/absence matrix.


process_chr_overlaps_generic <- function(gr, chr, label_col) {
  
  # Ensure the requested label exists in GRanges metadata
  stopifnot(label_col %in% names(mcols(gr)))
  
  # Subset GRanges to the current chromosome
  gr_chr <- gr[seqnames(gr) == chr]
  if (length(gr_chr) == 0) return(NULL)
  
  # Identify all unique breakpoints (segment starts and ends)
  breakpoints <- sort(unique(c(start(gr_chr), end(gr_chr))))
  if (length(breakpoints) < 2) return(NULL)
  
  # Define atomic bins between consecutive breakpoints
  split_ranges <- IRanges(
    start = head(breakpoints, -1),
    end   = tail(breakpoints, -1)
  )
  
  split_gr <- GRanges(seqnames = chr, ranges = split_ranges)
  
  # Find overlaps between bins and original segments
  overlaps <- findOverlaps(split_gr, gr_chr, minoverlap = 2)
  if (length(overlaps) == 0) return(NULL)
  
  # Initialise output data frame
  df <- data.frame(
    seqnames = as.character(seqnames(split_gr)),
    start    = start(split_gr),
    end      = end(split_gr),
    num_ind  = integer(length(split_gr)),   # number of individuals per bin
    ind_list = character(length(split_gr)), # list of individuals per bin
    stringsAsFactors = FALSE
  )
  
  # Extract individual labels
  labels <- mcols(gr_chr)[[label_col]]
  
  # Populate counts and individual lists
  for (i in seq_along(overlaps@from)) {
    region_idx <- overlaps@from[i]
    label_idx  <- overlaps@to[i]
    
    df$num_ind[region_idx] <- df$num_ind[region_idx] + 1
    df$ind_list[region_idx] <- paste(
      df$ind_list[region_idx],
      labels[label_idx],
      sep = ","
    )
  }
  
  # Clean leading/trailing commas
  df$ind_list <- gsub("^,|,$", "", df$ind_list)
  
  # Keep only bins overlapped by at least one individual
  df[df$num_ind > 0, ]
}


# ============================================================
#  COLLAPSED BINARY PRESENCE/ABSENCE MATRIX
# ============================================================
#
# Atomic bins are converted into a binary matrix indicating
# presence (1) or absence (0) of introgression for each
# individual.
#
# Adjacent bins with identical presence/absence patterns are
# collapsed to reduce redundancy while preserving sharing
# structure.

build_binary_matrix <- function(df_overlap,
                                collapse_adjacent = c("No", "Yes")) {
  
  collapse_adjacent <- match.arg(collapse_adjacent)
  
  stopifnot(is.data.frame(df_overlap))
  
  required_cols <- c("seqnames", "start", "end", "num_ind", "ind_list")
  missing_cols  <- setdiff(required_cols, colnames(df_overlap))
  
  if (length(missing_cols) > 0) {
    stop("df_overlap missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  
  if (nrow(df_overlap) == 0) {
    stop("df_overlap is empty – cannot build binary matrix")
  }
  
  # ------------------------------------------------------------
  # Build atomic-bin binary matrix
  # ------------------------------------------------------------
  
  all_ind <- unique(unlist(strsplit(
    paste(df_overlap$ind_list, collapse = ","), ","
  )))
  
  total_ind <- length(all_ind)
  
  df_bin <- df_overlap %>%
    mutate(freq = num_ind / total_ind) %>%
    mutate(ind_list = strsplit(ind_list, ",")) %>%
    tidyr::unnest(ind_list) %>%
    mutate(value = 1) %>%
    tidyr::pivot_wider(
      names_from  = ind_list,
      values_from = value,
      values_fill = list(value = 0)
    )
  
  # Standardise chromosome column name
  if ("seqnames" %in% colnames(df_bin)) { colnames(df_bin)[colnames(df_bin) == "seqnames"] <- "chrom" }
  
  ind_cols <- names(df_bin)[6:ncol(df_bin)]
  
  message("✔️  Built atomic-bin matrix: ",
          nrow(df_bin), " bins × ",
          length(ind_cols), " individuals")
  
  # ------------------------------------------------------------
  # Optional collapse of adjacent bins
  # ------------------------------------------------------------
  
  if (collapse_adjacent == "Yes") {
    
    message("🔧  Collapsing adjacent bins with identical sharing patterns")
    
    df_out <- df_bin %>%
      mutate(
        group = cumsum(
          (start - lag(end, default = start[1])) > 0 |
            (chrom != lag(chrom, default = chrom[1]))
        )
      ) %>%
      group_by(chrom, group) %>%
      summarise(
        start = min(start),
        end   = max(end),
        across(all_of(ind_cols), ~ as.integer(sum(.) > 0)),
        .groups = "drop"
      ) %>%
      dplyr::select(-group)
    
    message("✔️  Collapsed to ", nrow(df_out), " regions")
    
  } else {
    
    message("ℹ️  No collapse requested – keeping atomic bins")
    df_out <- df_bin
  }
  
  list(
    df_binary = df_out,
    ind_cols  = ind_cols,
    collapsed = collapse_adjacent
  )
}


# ============================================================
#  MDS COMPUTATION (TABLES ONLY)
# ============================================================
#
# This function computes:
#   - pairwise shared-region counts
#   - row-normalized sharing matrices
#   - two-dimensional MDS coordinates
#
# IMPORTANT:
#   - No plotting is performed here.
#   - The output consists exclusively of numerical tables.
#
# This separation ensures that visualization choices do not
# influence quantitative results.

compute_mds_from_collapsed <- function(df_collapsed, ind_cols, coverage_name) {
  
  mat <- as.matrix(df_collapsed[, ind_cols])
  storage.mode(mat) <- "numeric"
  
  mat_shared <- matrix(0, ncol = ncol(mat), nrow = ncol(mat))
  for (i in seq_len(ncol(mat))) {
    for (j in seq_len(ncol(mat))) {
      mat_shared[i, j] <- sum(pmin(mat[, i], mat[, j], na.rm = TRUE))
    }
  }
  
  rownames(mat_shared) <- colnames(mat)
  colnames(mat_shared) <- colnames(mat)
  diag(mat_shared) <- NA
  
  mat_norm <- sweep(mat_shared, 1, rowSums(mat_shared, na.rm = TRUE), FUN = "/")
  
  mds <- cmdscale(dist(mat_norm), k = 2)
  
  list(
    mds = data.frame(
      Dim1 = mds[,1],
      Dim2 = mds[,2],
      Individual = rownames(mds),
      Coverage = coverage_name
    ),
    mat_shared = mat_shared,
    mat_normalized = mat_norm
  )
}

# ============================================================
#  HEATMAP DATA COMPUTATION (TABLES ONLY)
# ============================================================
#
# This function computes sharing matrices suitable for heatmap
# visualization, including:
#   - raw shared-region counts
#   - normalized sharing proportions
#   - display-ready matrices
#   - population annotations
#
# No clustering or plotting is performed at this stage.


compute_heatmap_tables_from_collapsed <- function(
    df_collapsed,
    ind_cols,
    coverage_name,
    metadata_info,
    scale_percent = TRUE
) {
  
  stopifnot(is.data.frame(df_collapsed))
  stopifnot(all(ind_cols %in% colnames(df_collapsed)))
  
  # -------------------------------
  # 1. Binary matrix
  # -------------------------------
  mat <- as.matrix(df_collapsed[, ind_cols])
  storage.mode(mat) <- "numeric"
  
  if (ncol(mat) < 2) {
    stop("Need at least two individuals for heatmap: ", coverage_name)
  }
  
  # -------------------------------
  # 2. Shared-region matrix
  # -------------------------------
  mat_shared <- matrix(0, ncol = ncol(mat), nrow = ncol(mat))
  for (i in seq_len(ncol(mat))) {
    for (j in seq_len(ncol(mat))) {
      mat_shared[i, j] <- sum(pmin(mat[, i], mat[, j], na.rm = TRUE))
    }
  }
  
  rownames(mat_shared) <- colnames(mat)
  colnames(mat_shared) <- colnames(mat)
  diag(mat_shared) <- NA
  
  # -------------------------------
  # 3. Normalize
  # -------------------------------
  mat_norm <- sweep(
    mat_shared,
    1,
    rowSums(mat_shared, na.rm = TRUE),
    FUN = "/"
  )
  
  mat_display <- if (scale_percent) mat_norm * 100 else mat_norm
  
  # -------------------------------
  # 4. Population annotation
  # -------------------------------
  annotation <- metadata_info %>%
    dplyr::filter(Individual %in% rownames(mat_display)) %>%
    dplyr::select(Individual, Population) %>%
    tibble::column_to_rownames("Individual")
  
  annotation <- annotation[rownames(mat_display), , drop = FALSE]
  
  list(
    mat_shared     = mat_shared,
    mat_normalized = mat_norm,
    mat_display    = mat_display,
    annotation     = annotation,
    coverage       = coverage_name
  )
}

# ============================================================
#  WRAPPER: RUN THE FULL PIPELINE FOR ONE COVERAGE
# ============================================================

# This wrapper:
#   - builds GRanges for a given coverage
#   - computes overlaps chromosome by chromosome
#   - builds and collapses the binary matrix
#   - computes MDS coordinates
process_coverage_mds_and_heatmap <- function(
    coverage_name,
    IntrogressionMaps,
    metadata_info,
    population_colors,
    collapse_adjacent = "Yes"
) {
  
  message("🔹 Processing coverage: ", coverage_name)
  
  # -------------------------------
  # 1. Build GRanges
  # -------------------------------
  gr <- build_granges_for_coverage(IntrogressionMaps, coverage_name)
  if (is.null(gr)) return(NULL)
  
  # -------------------------------
  # 2. Atomic overlaps per chromosome
  # -------------------------------
  chromosomes <- unique(seqnames(gr))
  
  df_overlap <- dplyr::bind_rows(
    lapply(
      chromosomes,
      process_chr_overlaps_generic,
      gr = gr,
      label_col = "Individual_Coverage"
    )
  )
  
  if (nrow(df_overlap) == 0) return(NULL)
  
  # -------------------------------
  # 3. Binary matrix (collapse switch HERE)
  # -------------------------------
  bin_out <- build_binary_matrix(
    df_overlap,
    collapse_adjacent = collapse_adjacent
  )
  
  # -------------------------------
  # 4. MDS computation
  # -------------------------------
  mds_out <- compute_mds_from_collapsed(
    bin_out$df_binary,
    bin_out$ind_cols,
    coverage_name
  )
  
  # -------------------------------
  # 5. Heatmap tables
  # -------------------------------
  heatmap_tables <- compute_heatmap_tables_from_collapsed(
    bin_out$df_binary,
    bin_out$ind_cols,
    coverage_name,
    metadata_info
  )
  
  # -------------------------------
  # 6. Output
  # -------------------------------
  list(
    mds = mds_out$mds,
    mds_matrices = mds_out[c("mat_shared", "mat_normalized")],
    heatmap_tables = heatmap_tables,
    collapsed = bin_out$collapsed
  )
}

# ============================================================
#  VISUALIZATION FUNCTIONS
# ============================================================
#
# The functions below generate figures exclusively from
# precomputed tables.
#
# This design guarantees:
#   - consistent aesthetics across figures
#   - reproducible regeneration of plots
#   - flexible reuse of numerical results

# ============================================================
#  MDS VISUALIZATION
# ============================================================

plot_mds <- function(
    mds_df,
    metadata_info,
    population_colors,
    title = NULL
) {
  
  df_plot <- mds_df %>%
    left_join(metadata_info, by = "Individual")
  
  ggplot(df_plot, aes(Dim1, Dim2, color = Population)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey80") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(
      aes(label = Individual),
      size = 3,
      max.overlaps = 1000,
      show.legend = FALSE
    ) +
    scale_color_manual(values = population_colors) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      title = if (!is.null(title)) title else unique(mds_df$Coverage),
      x = "Dim 1",
      y = "Dim 2"
    )
}

# ============================================================
#  HEATMAP VISUALIZATION
# ============================================================

plot_heatmap <- function(
    mat_display,
    annotation,
    population_colors,
    title = NULL,
    cluster_rows = TRUE
) {
  
  # -------------------------------
  # 1. Extract row clustering
  # -------------------------------
  tmp <- pheatmap::pheatmap(
    mat_display,
    cluster_rows = cluster_rows,
    cluster_cols = FALSE,
    silent = TRUE
  )
  
  row_order <- tmp$tree_row$order
  sample_order <- rownames(mat_display)[row_order]
  
  mat_ordered <- mat_display[sample_order, sample_order]
  annotation  <- annotation[sample_order, , drop = FALSE]
  
  # -------------------------------
  # 2. Final heatmap
  # -------------------------------
  pheatmap::pheatmap(
    mat_ordered,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    treeheight_col = 0,
    display_numbers = round(mat_ordered, 1),
    main = title,
    color = viridis::viridis(100),
    annotation_row = annotation,
    annotation_col = annotation,
    annotation_colors = list(Population = population_colors),
    annotation_legend = TRUE,
    legend = TRUE
  )
}

# ============================================================
#  RUN EXAMPLE
# ============================================================

coverage_run <- "Original"
df_introgression_maps <- IntrogressionMaps_Filtered_ID

out_coverage_run <- process_coverage_mds_and_heatmap(
  coverage_run,
  df_introgression_maps,
  metadata_info,
  population_colors
)

# ---- tables ----
out_coverage_run$mds
out_coverage_run$heatmap_tables$mat_display

# ---- plots (optional!) ----
p_mds <- plot_mds(
  out_coverage_run$mds,
  metadata_info,
  population_colors,
  title = paste0("Dataset: ", coverage_run)
)

p_heat <- plot_heatmap(
  out_coverage_run$heatmap_tables$mat_display,
  out_coverage_run$heatmap_tables$annotation,
  population_colors,
  title = paste0("Dataset: ", coverage_run)
)

print(p_mds)
print(p_heat)

# ============================================================
#  MULTI-COVERAGE ANALYSIS
# ============================================================
#
# All coverage levels are processed in the order specified by
# 'coverage_values'. This ordering is enforced at the table
# level and propagated consistently to all visualizations.

df_IntrogressionMaps <- IntrogressionMaps_Filtered_ID


mds_results <- setNames(
  lapply(coverage_values, function(cov) {
    out <- process_coverage_mds_and_heatmap(
      coverage_name = cov,
      IntrogressionMaps = df_IntrogressionMaps,
      metadata_info = metadata_info,
      population_colors = population_colors
    )
    if (is.null(out)) return(NULL)
    out$mds
  }),
  coverage_values
)

# Drop coverages with no data but KEEP ORDER
mds_results <- mds_results[coverage_values[!sapply(mds_results, is.null)]]



# ============================================================
#  GENERATE ONE MDS PLOT PER COVERAGE
# ============================================================
# For each annotated MDS data frame, create a ggplot object.
# Each plot corresponds to a single coverage and shows:
#   - individuals positioned by MDS coordinates
#   - points coloured by population
#   - dashed reference lines at zero
#   - individual sample labels

mds_plots_named <- setNames(
  lapply(names(mds_results), function(cov) {
    plot_mds(
      mds_df = mds_results[[cov]],
      metadata_info = metadata_info,
      population_colors = population_colors,
      title = cov
    )
  }),
  names(mds_results)
)

# ============================================================
#  COMBINE ALL COVERAGE-SPECIFIC PLOTS INTO A SINGLE FIGURE
# ============================================================

# Arrange all coverage-specific MDS plots into a multi-panel figure.
# Panels are ordered according to the order of 'coverages' and
# automatically labelled (A, B, C, ...).
#
# The legend is collected and placed at the bottom for clarity.

combined_named_mds_plot <-
  wrap_plots(
    mds_plots_named,
    ncol = 3,
    guides = "collect"
  ) +
  plot_annotation(
    title = "MDS on Shared Archaic Segments",
    tag_levels = "A"
  ) &
  theme(legend.position = "bottom")


# ============================================================
#  SAVE THE COMBINED MDS PLOT TO FILE
# ============================================================

# Save the multi-panel MDS figure to disk.
# A large width/height is required because patchwork
# combines multiple panels with a shared legend.
# ggsave(
#   filename = "MDS_archaic_segments.png",
#   plot     = combined_named_mds_plot,
#   width    = 14,      # increase if you have many panels
#   height   = 10,
#   units    = "in",
#   dpi      = 300
# )


# ============================================================
#  RANDOMIZED COVERAGE CONTROL ANALYSIS
# ============================================================
#
# This section implements a control analysis in which sequencing
# coverage labels are randomly reassigned to individuals.
#
# The goal is to assess whether the structure observed in the
# coverage-specific MDS analyses could arise purely from the
# distribution of individuals across coverage levels, rather
# than from biologically meaningful sharing of introgressed
# segments.
#
# Coverage labels are reassigned in a balanced manner, such that
# the overall frequency of each coverage level is preserved.
#
# A fixed random seed is used to ensure full reproducibility.

seed <- 999999999

# ------------------------------------------------------------
# RANDOM, BALANCED ASSIGNMENT OF COVERAGE LABELS
# ------------------------------------------------------------
#
# Each individual is randomly assigned a coverage label, while
# preserving the total number of samples per coverage level.
#
# This ensures that downstream analyses are not biased by
# unequal coverage frequencies.
#
# The output is a two-column table mapping individuals to their
# randomly assigned coverage.

assign_random_coverages <- function(individuals, coverages, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  n_ind <- length(individuals)
  
  coverage_vec <- rep(coverages, length.out = n_ind)
  
  tibble(
    Individual = individuals,
    Coverage   = sample(coverage_vec)
  )
}

# Random balanced assignment
random_assignment <- assign_random_coverages(
  individuals = unique(metadata_info$Name),
  coverages   = unique(metadata_info$Coverage),
  seed        = seed
)

# ------------------------------------------------------------
# FILTER INTROGRESSION MAPS TO RANDOMIZED DATASET
# ------------------------------------------------------------
#
# The introgression maps are filtered to retain only those
# individual–coverage combinations present in the randomized
# assignment.
#
# Genomic coordinates are explicitly coerced to appropriate
# types and invalid ranges are removed as a safety check.


IntrogressionMaps_RandomBalanced <- IntrogressionMaps_Filtered_ID %>%
  inner_join(random_assignment, by = c("Individual", "Coverage")) %>%
  mutate(
    Individual_selection = paste0(Individual, "_", Coverage),
    chr   = as.character(chr),
    start = as.integer(start),
    end   = as.integer(end)
  ) %>%
  filter(!is.na(start), !is.na(end), start < end)

# ------------------------------------------------------------
# SANITY CHECKS
# ------------------------------------------------------------
#
# Verify that coverage frequencies are balanced and that each
# individual is represented exactly once in the randomized
# dataset.


table(random_assignment$Coverage)
length(unique(IntrogressionMaps_RandomBalanced$Individual))


# ============================================================
#  GRanges CONSTRUCTION FOR MIXED-COVERAGE DATASET
# ============================================================
#
# Introgressed segments from the randomized dataset are converted
# into a GRanges object. The individual–coverage identifier is
# stored as metadata and used for multi-individual intersection
# analyses.

gr_mixed <- GenomicRanges::GRanges(
  seqnames = IntrogressionMaps_RandomBalanced$chr,
  ranges   = IRanges::IRanges(
    start = IntrogressionMaps_RandomBalanced$start,
    end   = IntrogressionMaps_RandomBalanced$end
  ),
  Individual_Coverage = IntrogressionMaps_RandomBalanced$Individual_selection
)

# ============================================================
#  GENOME-WIDE MULTI-INDIVIDUAL INTERSECTION
# ============================================================
#
# The genome is decomposed into minimal non-overlapping intervals
# defined by shared segment boundaries across individuals.
#
# For each interval, the set of individuals carrying introgressed
# sequence is recorded.

chromosomes <- unique(seqnames(gr_mixed))

df_overlap <- dplyr::bind_rows(
  lapply(
    chromosomes,
    process_chr_overlaps_generic,
    gr = gr_mixed,
    label_col = "Individual_Coverage"
  )
)

# ============================================================
#  BINARY PRESENCE / ABSENCE MATRIX CONSTRUCTION
# ============================================================
#
# The overlap information is converted into a binary matrix
# indicating presence or absence of introgression for each
# individual across genomic intervals.
#
# Adjacent intervals with identical presence patterns are
# collapsed to reduce redundancy.

bin_out <- build_binary_matrix(
  df_overlap,
  collapse_adjacent = "Yes"   # or "No", but here "Yes" makes sense
)

# ============================================================
#  MDS ON RANDOMIZED MIXED-COVERAGE DATASET
# ============================================================
#
# Pairwise introgression sharing is quantified, normalized, and
# visualized using Multidimensional Scaling (MDS).
#
# Because coverage labels are randomized, any residual structure
# observed in this MDS reflects background similarities among
# individuals rather than coverage-driven effects.

mds_random_mixed_out <- compute_mds_from_collapsed(
  bin_out$df_binary,
  bin_out$ind_cols,
  coverage_name = "RandomMixedCoverages"
)

# ============================================================
#  VISUALIZATION OF RANDOMIZED CONTROL
# ============================================================
#
# MDS coordinates are annotated with population labels and
# visualized using the same color scheme and graphical settings
# as the main analysis, allowing direct qualitative comparison.
mds_random_mixed <- mds_random_mixed_out$mds %>%
  mutate(Individual_ID = sub("_(.*)$", "", Individual))

mds_random_mixed <- mds_random_mixed %>%
  left_join(metadata_info, by = "Individual")

p_mds_random_mixed <- ggplot(
  mds_random_mixed,
  aes(x = Dim1, y = Dim2, color = Population)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey90", linewidth = 0.01) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey90", linewidth = 0.01) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(
    values = population_colors,
    guide = guide_legend(override.aes = list(size = 4))
  ) +
  geom_text_repel(
    aes(label = Individual),
    size = 3.5,
    max.overlaps = 1000,
    box.padding = 0.5,
    show.legend = FALSE
  ) +
  labs(x = "Dim 1", y = "Dim 2") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "",
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.text = element_text(color = "black")
  )

p_mds_random_mixed


# ============================================================
#  SAVE THE MDS PLOT TO FILE
# ============================================================

# ggsave(
#   filename = "MDS_random_selection.png",
#   plot     = p_mds_random_mixed,
#   width    = 9,     
#   height   = 6,
#   units    = "in",
#   dpi      = 300
# )