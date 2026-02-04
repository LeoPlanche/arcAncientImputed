# ============================================================
#  LIBRARIES
# ============================================================

library(dplyr)
library(GenomicRanges)
library(biomaRt)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(future.apply)

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

# Create a unique identifier for each individual–coverage pair
# (e.g. Ust_Ishim_Original, Ust_Ishim_1X, ...)
IntrogressionMaps_ID <- IntrogressionMaps %>%
  mutate(Individual_selection = paste0(Individual, "_", Coverage))
IntrogressionMaps_Filtered_ID <- IntrogressionMaps_ID %>%
  filter(Filter == "Kept")

##### Region to test
# Load the data
gene_table <- read.delim("../Data/ListGenesToTest.txt", header = TRUE, na.strings = c(""))

#============================================#
#      Check 14 known candidate regions      #
#============================================#
#
# Goal
# ----
# Quantify introgression evidence at 14 previously proposed candidate loci by:
#   (i) building a transparent per-individual table WITHOUT applying any filtering,
#       keeping the Filter labels so users can inspect results themselves; and
#   (ii) producing a manuscript-ready gene × coverage summary using ONLY fragments
#       passing the filtering criterion ("Kept").
#
# Inputs
# ------
# gene_table:
#   Candidate regions (Chromosome, Start, End). Coordinates are interpreted on GRCh37/hg19.
#
# IntrogressionMaps_ID:
#   Fragment table containing genomic coordinates (Chrom/start/end) plus metadata:
#   Individual, Coverage, Filter, similarity/distance metrics, private SNP counts, etc.
#
# Output files
# ------------
# 1) SegmentsOverlapJan2026.txt
#    Per-gene × per-individual × per-coverage summary WITHOUT filtering.
#    The Filter column is retained for transparency.
#
# 2) FullGeneOverlap_Summary_SimilarityDistanceJan2026.txt
#    Gene × coverage table: number of individuals with evidence of introgression,
#    computed AFTER filtering to "Kept".
#
#============================================#
# Validate inputs
#============================================#
# Ensure the fragment table exists in the environment before proceeding

stopifnot(exists("IntrogressionMaps_ID"))

#============================================#
# Clean candidate region table
#============================================#
# Standardise candidate region coordinates:
#   - remove commas from Start/End
#   - convert Start/End to numeric
#   - remove any "chr" prefix from Chromosome
#   - remove incomplete rows
gene_table_clean <- gene_table %>%
  mutate(
    Start = as.numeric(gsub(",", "", Start)),
    End   = as.numeric(gsub(",", "", End)),
    Chromosome = as.character(Chromosome),
    Chromosome = gsub("^chr", "", Chromosome)
  ) %>%
  filter(!is.na(Start) & !is.na(End) & !is.na(Chromosome))

#============================================#
# STEP 3 — Retrieve protein-coding genes (GRCh37)
#============================================#
# Query Ensembl GRCh37 (hg19) for genes overlapping each candidate region,
# then retain only protein-coding genes with valid symbols.
ensembl_grch37 <- biomaRt::useMart(
  host    = "https://grch37.ensembl.org",
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl"
)


query_genes_one_region <- function(seqname, start, end, mart) {
  tryCatch({
    biomaRt::getBM(
      attributes = c(
        "chromosome_name", "start_position", "end_position",
        "ensembl_gene_id", "external_gene_name", "gene_biotype"
      ),
      filters = c("chromosome_name", "start", "end"),
      values  = list(seqname, start, end),
      mart    = mart
    )
  }, error = function(e) {
    data.frame()
  })
}

overlapping_genes_df <- dplyr::bind_rows(lapply(seq_len(nrow(gene_table_clean)), function(i) {
  query_genes_one_region(
    seqname = gene_table_clean$Chromosome[i],
    start   = gene_table_clean$Start[i],
    end     = gene_table_clean$End[i],
    mart    = ensembl_grch37
  )
}))

overlapping_genes_df_protein <- overlapping_genes_df %>%
  filter(gene_biotype == "protein_coding") %>%
  filter(!is.na(external_gene_name) & external_gene_name != "") %>%
  distinct(chromosome_name, start_position, end_position, external_gene_name, .keep_all = TRUE)

#============================================#
# Convert fragments and genes to GRanges
#============================================#
# Represent introgressed fragments as GRanges. All fragment metadata is stored
# in mcols so that it is retained after subsetting by overlap indices.
fragment_ranges <- GRanges(
  seqnames = gsub("^chr", "", IntrogressionMaps_ID$Chrom),
  ranges   = IRanges(
    start = IntrogressionMaps_ID$start,
    end   = IntrogressionMaps_ID$end
  ),
  mcols = IntrogressionMaps_ID %>%
    dplyr::select(-Chrom, -start, -end)
)

# ---------------------------------------------------
# Create GRanges for genes and find overlaps
# ---------------------------------------------------
gene_ranges <- GRanges(
  seqnames = overlapping_genes_df_protein$chromosome_name,
  ranges   = IRanges(
    overlapping_genes_df_protein$start_position,
    overlapping_genes_df_protein$end_position
  ),
  gene = overlapping_genes_df_protein$external_gene_name
)

#============================================#
# Compute gene–fragment overlaps (hit-level table)
#============================================#
# findOverlaps returns indices linking genes (queryHits) to fragments (subjectHits).
# We use these indices to build a "hit-level" table where each row is one overlap,
# retaining all fragment metadata.
hits <- findOverlaps(gene_ranges, fragment_ranges)
message("✔️  gene-fragment overlaps: ", length(hits))

if (length(hits) == 0) {
  
  overlap_hits_df <- tibble()
  
} else {
  
  # Subset GRanges using overlap indices
  gene_hit <- gene_ranges[queryHits(hits)]
  frag_hit <- fragment_ranges[subjectHits(hits)]
  
  # Build hit-level dataframe
  overlap_hits_df <- tibble(
    gene = mcols(gene_hit)$gene,
    
    gene_chr   = as.character(seqnames(gene_hit)),
    gene_start = start(gene_hit),
    gene_end   = end(gene_hit),
    
    fragment_chr   = as.character(seqnames(frag_hit)),
    fragment_start = start(frag_hit),
    fragment_end   = end(frag_hit),
    
    overlap_length = width(
      pintersect(ranges(gene_hit), ranges(frag_hit))
    ),
    
    gene_length = width(gene_hit)
  ) %>%
    bind_cols(as.data.frame(mcols(frag_hit)))
}

# Clean up any "mcols." prefix introduced when binding to tibble
overlap_hits_df <- overlap_hits_df %>%
  rename_with(~ gsub("^mcols\\.", "", .x))

#============================================#
# Define aggregation rules (SUM vs MEAN)
#============================================#
# Count-based evidence is summed across overlapping fragments, whereas
# similarity/distance metrics are averaged.
sum_cols <- c(
  "hmmPositionsOriginal",
  "hmmPositionsImp",
  grep("^PrivateSNPsOriginal_", colnames(overlap_hits_df), value = TRUE),
  grep("^PrivateSNPsImp_", colnames(overlap_hits_df), value = TRUE)
)

#============================================#
# Per-individual summary WITHOUT filtering (transparency table)
#============================================#
# For each gene × individual × coverage:
#   - compute % of gene overlapped by introgressed fragments
#   - retain Filter labels (not applied here)
#   - sum count-based metrics, average similarity/distance metrics
individuals_table <- overlap_hits_df %>%
  mutate(
    genomic_range = paste0(gene_chr, ":", gene_start, "-", gene_end)
  ) %>%
  group_by(gene, genomic_range, Individual, Filter, Coverage) %>%
  summarise(
    
    # gene coverage in this individual (no filtering)
    overlap_percentage =
      100 * sum(overlap_length) / unique(gene_length),
    
    # keep Filter info (do NOT apply it)
    Filter = paste(unique(Filter), collapse = ";"),
    
    # SUM evidence
    across(
      all_of(sum_cols),
      ~ sum(.x, na.rm = TRUE)
    ),
    
    # MEAN similarity / distance
    across(
      all_of(mean_cols),
      ~ mean(.x, na.rm = TRUE)
    ),
    
    n_fragments = n(),
    .groups = "drop"
  )


# # Save individual segments
write.table(
  individuals_table,
  "SegmentsOverlapJan2026.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

#============================================#
# Apply fragment filtering for manuscript summaries
#============================================#
# Restrict to fragments passing filtering criteria ("Kept"), then recompute
# gene × individual × coverage summaries (used for gene-level counts).

overlap_hits_df_Filtered <- overlap_hits_df %>%
  filter(Filter == "Kept")


gene_summary_df <- overlap_hits_df_Filtered %>%
  group_by(gene, Individual, Coverage) %>%
  summarise(
    
    # gene coverage
    overlap_percentage =
      100 * sum(overlap_length) / unique(gene_length),
    
    # SUM counts
    across(
      all_of(sum_cols),
      ~ sum(.x, na.rm = TRUE)
    ),
    
    # MEAN similarity / distance
    across(
      all_of(mean_cols),
      ~ mean(.x, na.rm = TRUE)
    ),
    
    # bookkeeping
    n_fragments = n(),
    
    .groups = "drop"
  )

#============================================#
# Ensure all candidate genes are represented (add zero-signal genes)
#============================================#
# Genes with no retained overlaps are included explicitly so they appear in the
# final gene × coverage table with zero counts.

genes_with_signal <- unique(overlap_hits_df_Filtered$gene)

missing_genes <- setdiff(
  unique(mcols(gene_ranges)$gene),
  genes_with_signal
)

missing_gene_ranges <- gene_ranges[
  mcols(gene_ranges)$gene %in% missing_genes
]

missing_genes_df <- tibble(
  gene = mcols(missing_gene_ranges)$gene,
  genomic_range = paste0(
    seqnames(missing_gene_ranges), ":",
    start(missing_gene_ranges), "-",
    end(missing_gene_ranges)
  )
)

gene_to_range <- overlap_hits_df_Filtered %>%
  distinct(gene, gene_chr, gene_start, gene_end) %>%
  mutate(
    genomic_range = paste0(gene_chr, ":", gene_start, "-", gene_end)
  ) %>%
  dplyr::select(gene, genomic_range) %>%
  bind_rows(missing_genes_df) %>%
  distinct()


#============================================#
# Gene × coverage summary (counts of individuals)
#============================================#
# Count how many distinct individuals show evidence of introgression per gene
# at each coverage level, applying an optional minimum overlap threshold.

presence_threshold <- 0   # % gene length

gene_coverage_counts <- gene_summary_df %>%
  filter(overlap_percentage >= presence_threshold) %>%
  group_by(gene, Coverage) %>%
  summarise(
    n_individuals = n_distinct(Individual),
    .groups = "drop"
  )
gene_coverage_wide <- gene_coverage_counts %>%
  tidyr::pivot_wider(
    names_from  = Coverage,
    values_from = n_individuals,
    values_fill = 0
  )
gene_info <- gene_to_range %>%
  distinct(gene, genomic_range)

gene_coverage_complete <- gene_info %>%
  left_join(gene_coverage_wide, by = "gene") %>%
  mutate(
    across(all_of(dataset_order), ~ replace_na(.x, 0))
  ) %>%
  dplyr::select(gene, genomic_range, all_of(dataset_order))

# # Save final summary table
write.table(
  gene_coverage_complete,
  "FullGeneOverlap_Summary_SimilarityDistanceJan2026.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

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
# NOT collapsed to reduce redundancy while preserving sharing
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
    # filter(!is.na(ind_list), ind_list != "") %>%
    mutate(value = 1) %>%
    tidyr::pivot_wider(
      names_from  = ind_list,
      values_from = value,
      values_fill = list(value = 0)
    )

  # Standardise chromosome column name
  if ("seqnames" %in% colnames(df_bin)) { colnames(df_bin)[colnames(df_bin) == "seqnames"] <- "chrom" }
  
  if (any(is.na(colnames(df_bin)))) {
    stop("NA column names detected after pivot_wider() — check ind_list cleaning")
  }
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

# ------------------------------------------------------------
# Run multi-intersection per coverage (no collapse)
# ------------------------------------------------------------

run_multiintersection_per_coverage <- function(coverage_name,
                                               IntrogressionMaps,
                                               collapse_adjacent = "No",
                                               workers = parallel::detectCores() - 1) {
  
  # Allow single coverage or vector of coverages
  if (is.null(coverage_name)) {
    coverages <- unique(IntrogressionMaps$Coverage)
  } else {
    coverages <- coverage_name
  }
  
  results <- vector("list", length(coverages))
  names(results) <- coverages
  
  for (cov in coverages) {
    
    message("🔬 Processing coverage: ", cov)
    
    # ------------------------------------------------------------
    # 1. Build GRanges
    # ------------------------------------------------------------
    gr <- build_granges_for_coverage(IntrogressionMaps, cov)
    if (is.null(gr)) {
      results[[cov]] <- NULL
      next
    }
    
    chromosomes <- unique(as.character(seqnames(gr)))
    
    # ------------------------------------------------------------
    # 2. Multi-intersection (parallel by chromosome)
    # ------------------------------------------------------------
    plan(multisession, workers = workers)
    
    chr_results <- future_lapply(chromosomes, function(chr) {
      process_chr_overlaps_generic(
        gr        = gr,
        chr       = chr,
        label_col = "Individual_Coverage"
      )
    })
    
    plan(sequential)
    
    df_overlap <- dplyr::bind_rows(chr_results)
    if (nrow(df_overlap) == 0) {
      results[[cov]] <- NULL
      next
    }
    
    # ------------------------------------------------------------
    # 3. Binary matrix (atomic bins unless specified)
    # ------------------------------------------------------------
    bin_out <- build_binary_matrix(
      df_overlap,
      collapse_adjacent = collapse_adjacent
    )
    
    # ------------------------------------------------------------
    # 4. Match original Multiinters_bedtools_format
    # ------------------------------------------------------------
    results[[cov]] <- bin_out$df_binary %>%
      dplyr::rename(
        seqnames = chrom
      )
  }
  
  results
}

mat_1x <- run_multiintersection_per_coverage(
  coverage_name = "1X",
  IntrogressionMaps = IntrogressionMaps_Filtered_ID,
  collapse_adjacent = "No"
)

# ============================================================
#  MULTI-COVERAGE ANALYSIS
# ============================================================
#
# All coverage levels are processed in the order specified by
# 'coverage_values'. This ordering is enforced at the table
# level and propagated consistently to all visualizations.

df_IntrogressionMaps <- IntrogressionMaps_Filtered_ID

mat_all <- run_multiintersection_per_coverage(
  coverage_name      = coverage_values,   # all coverages
  IntrogressionMaps  = df_IntrogressionMaps,
  collapse_adjacent  = "No"
)

# ============================================================
#  STANDARDISE MULTI-INTERSECTION OUTPUT FOR SCAN
# ============================================================
#
# The selection scan operates on introgressed atomic bins and
# evaluates how frequently introgression is shared across
# individuals within a population subset.
#
# In this pipeline, the fundamental unit is the *individual*
# (not technical samples). Therefore:
#   - num_ind already represents the correct count
#   - individual presence/absence columns are retained
#
# This step only harmonises column names and removes fields
# that are not required downstream.

standardise_for_scan <- function(mat_list) {
  stopifnot(is.list(mat_list))
  
  out <- vector("list", length(mat_list))
  names(out) <- names(mat_list)
  
  for (cov in names(mat_list)) {
    df <- mat_list[[cov]]
    
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
      out[[cov]] <- data.frame()
      next
    }
    
    # Required columns
    required_cols <- c("seqnames", "start", "end", "num_ind", "freq")
    missing_cols  <- setdiff(required_cols, colnames(df))
    if (length(missing_cols) > 0) {
      stop("Coverage ", cov, " is missing columns: ",
           paste(missing_cols, collapse = ", "))
    }
    
    # Identify individual columns (everything else)
    meta_cols <- intersect(
      colnames(df),
      c("seqnames", "start", "end", "num_ind", "ind_list", "freq")
    )
    individual_cols <- setdiff(colnames(df), meta_cols)
    
    # Clean individual names: remove coverage suffix if present
    clean_names <- gsub(paste0("_", cov, "$"), "", individual_cols)
    colnames(df)[match(individual_cols, colnames(df))] <- clean_names
    
    # Keep only scan-relevant fields
    df2 <- df %>%
      dplyr::select(
        seqnames, start, end,
        num_ind, freq,
        all_of(clean_names)
      )
    
    out[[cov]] <- df2
  }
  
  out
}

coverage_results_list <- standardise_for_scan(mat_all)


# Optional quick sanity check
lapply(coverage_results_list, function(x) {
  if (is.null(x) || nrow(x) == 0) return(c(nrow = 0, ncol = 0))
  c(
    nrow = nrow(x),
    ncol = ncol(x),
    has_meta = all(c("seqnames","start","end","num_ind","freq") %in% colnames(x))
  )
})

# ============================================================
#  DEFINE EUROPEAN INDIVIDUAL SET
# ============================================================
#
# Individual column names in the introgression matrices were
# cleaned upstream and correspond to the 'Name' field in
# metadata_info. Population assignment is therefore matched
# using 'Name', NOT Individual_selection.

european_individuals <- metadata_info %>%
  dplyr::filter(grepl("^European", Population)) %>%
  dplyr::pull(Name) %>%
  unique()

# ------------------------------------------------------------
# Sanity check: ensure metadata individuals are present
# in the introgression matrices
# ------------------------------------------------------------

example_cov <- names(coverage_results_list)[1]

present_individuals <- setdiff(
  colnames(coverage_results_list[[example_cov]]),
  c("seqnames", "start", "end", "num_ind", "freq")
)

message(
  "✔️  European individuals present in matrices: ",
  sum(european_individuals %in% present_individuals),
  " / ",
  length(european_individuals)
)

missing_eur <- setdiff(european_individuals, present_individuals)
if (length(missing_eur) > 0) {
  warning(
    "European individuals missing from matrices: ",
    paste(missing_eur, collapse = ", ")
  )
}

# ============================================================
#  SUBSET INTROGRESSION MATRICES TO EUROPEANS
# ============================================================
#
# For each coverage level, introgression sharing is recomputed
# using European individuals only. This avoids inflation of
# frequency estimates due to non-European carriers and ensures
# that downstream enrichment statistics reflect population-
# specific patterns.

subset_to_europeans <- function(coverage_results_list, european_individuals) {
  
  out <- vector("list", length(coverage_results_list))
  names(out) <- names(coverage_results_list)
  
  for (cov in names(coverage_results_list)) {
    
    df <- coverage_results_list[[cov]]
    
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
      out[[cov]] <- data.frame()
      next
    }
    
    # Identify individual columns
    meta_cols <- c("seqnames", "start", "end", "num_ind", "freq")
    individual_cols <- setdiff(colnames(df), meta_cols)
    
    # Keep only European individuals present in this coverage
    european_cols <- intersect(individual_cols, european_individuals)
    
    if (length(european_cols) == 0) {
      warning("No European individuals found for coverage: ", cov)
      out[[cov]] <- data.frame()
      next
    }
    
    # Subset and recompute statistics
    df_eur <- df %>%
      dplyr::select(seqnames, start, end, all_of(european_cols)) %>%
      mutate(
        num_ind = rowSums(across(all_of(european_cols)) > 0, na.rm = TRUE),
        freq    = num_ind / length(european_cols)
      ) %>%
      dplyr::filter(num_ind > 0)
    
    out[[cov]] <- df_eur
  }
  
  out
}

coverage_results_list_eur <- subset_to_europeans(
  coverage_results_list,
  european_individuals
)

# ============================================================
#  ONE-TAILED ENRICHMENT SCAN (EUROPEAN SUBSET)
# ============================================================
#
# For each coverage level, introgression frequency is converted
# into a Z-score relative to the genome-wide distribution.
# A one-tailed test is used to detect regions with unusually
# high sharing (overrepresentation) among Europeans.

compute_significant_segments <- function(df,
                                         alpha = 0.01) {
  
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
    return(list(
      df_with_zscores = data.frame(),
      significant_threshold = NA,
      significant_segments = data.frame()
    ))
  }
  
  # Genome-wide mean and standard deviation of introgression frequency
  mean_freq <- mean(df$freq, na.rm = TRUE)
  sd_freq   <- sd(df$freq, na.rm = TRUE)
  
  if (is.na(sd_freq) || sd_freq == 0) {
    warning("Zero or undefined SD in frequency distribution")
  }
  
  # Compute Z-scores
  df_z <- df %>%
    mutate(
      Z_score = (freq - mean_freq) / sd_freq
    )
  
  # One-tailed p-value for overrepresentation
  # (large positive Z)
  df_z <- df_z %>%
    mutate(
      p_value = pnorm(Z_score, lower.tail = FALSE),
      p_adj_bonferroni = pmin(1, p_value * nrow(df_z))
    )
  
  # Bonferroni-corrected Z threshold (for plotting / reference)
  significant_threshold <- qnorm(1 - (alpha / nrow(df_z)))
  
  # Core significant segments
  # Bonferroni correction across all tested bins within this coverage
  significant_segments <- df_z %>%
    filter(
      Z_score > 0,
      p_adj_bonferroni < alpha
    ) %>%
    arrange(p_adj_bonferroni)
  
  list(
    df_with_zscores        = df_z,
    significant_threshold = significant_threshold,
    significant_segments  = significant_segments
  )
}

# ------------------------------------------------------------
# Run enrichment scan for each coverage level
# ------------------------------------------------------------

coverage_significant_results_eur <- lapply(
  coverage_results_list_eur,
  compute_significant_segments
)

# ------------------------------------------------------------
# Extract core significant segments (Bonferroni < 0.01)
# ------------------------------------------------------------

significant_core_segments_per_coverage_eur <- lapply(
  coverage_significant_results_eur,
  function(res) res$significant_segments
)

# ============================================================
#  GENE ANNOTATION OF CORE SIGNIFICANT REGIONS
# ============================================================
#
# Protein-coding genes overlapping significantly enriched
# introgressed regions are retrieved using Ensembl (GRCh37),
# ensuring consistency with the reference genome used in
# introgression inference.

ensembl_grch37 <- useMart(
  host    = "https://grch37.ensembl.org",
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl"
)

# ------------------------------------------------------------
# Helper: query Ensembl genes for a genomic interval
# ------------------------------------------------------------

query_genes_in_region <- function(seqname, start, end, mart) {
  
  # Ensembl expects chromosomes without "chr"
  seqname <- gsub("^chr", "", seqname)
  
  tryCatch({
    getBM(
      attributes = c(
        "chromosome_name",
        "start_position",
        "end_position",
        "ensembl_gene_id",
        "external_gene_name",
        "gene_biotype"
      ),
      filters = c("chromosome_name", "start", "end"),
      values  = list(seqname, start, end),
      mart    = mart
    )
  }, error = function(e) {
    warning("Gene query failed for region: ",
            seqname, ":", start, "-", end)
    data.frame()
  })
}

# ------------------------------------------------------------
# Extract protein-coding genes from core significant regions
# ------------------------------------------------------------

protein_coding_genes_core_eur <- lapply(
  names(coverage_significant_results_eur),
  function(cov) {
    
    regions_df <- coverage_significant_results_eur[[cov]]$significant_segments
    
    if (is.null(regions_df) || nrow(regions_df) == 0) return(NULL)
    
    genes_list <- lapply(seq_len(nrow(regions_df)), function(i) {
      query_genes_in_region(
        seqname = regions_df$seqnames[i],
        start   = regions_df$start[i],
        end     = regions_df$end[i],
        mart    = ensembl_grch37
      )
    })
    
    genes_list <- genes_list[sapply(genes_list, nrow) > 0]
    if (length(genes_list) == 0) return(NULL)
    
    genes <- dplyr::bind_rows(genes_list) %>%
      dplyr::filter(gene_biotype == "protein_coding") %>%
      mutate(
        Coverage    = cov,
        region_type = "core"
      ) %>%
      distinct(
        chromosome_name,
        start_position,
        end_position,
        ensembl_gene_id,
        external_gene_name,
        .keep_all = TRUE
      )
    
    genes
  }
)

names(protein_coding_genes_core_eur) <- names(coverage_significant_results_eur)

# ============================================================
#  ANCHORED EXPANSION OF SIGNIFICANT REGIONS
# ============================================================
# Core significant bins identify loci with strong evidence of introgression enrichment.
# However, selective signals may extend beyond these narrow cores due to linkage and recombination.
# To capture the full genomic footprint of selection, we apply an anchored expansion strategy:
#   Anchor bins must pass a strict Bonferroni threshold
#   Adjacent bins are included if they show nominal evidence of enrichment# 
#   Expansion is restricted to contiguous bins on the same chromosome
# This yields extended candidate regions that preserve statistical rigor while improving biological interpretability.

extract_anchored_regions <- function(zscore_df,
                                     retain_threshold = 0.01,
                                     expand_threshold = 0.05) {
  
  if (is.null(zscore_df) || nrow(zscore_df) == 0) {
    return(data.frame())
  }
  
  # Ensure chromosome column
  if (!"chrom" %in% colnames(zscore_df)) {
    if ("seqnames" %in% colnames(zscore_df)) {
      zscore_df <- zscore_df %>%
        mutate(chrom = gsub("^chr", "", as.character(seqnames)))
    } else {
      stop("No chromosome information found")
    }
  }
  
  zscore_df <- zscore_df %>%
    arrange(as.numeric(chrom), start)
  
  results <- list()
  
  for (chr in unique(zscore_df$chrom)) {
    
    df_chr <- zscore_df %>% filter(chrom == chr)
    visited <- rep(FALSE, nrow(df_chr))
    
    for (i in seq_len(nrow(df_chr))) {
      
      if (visited[i]) next
      
      if (!is.na(df_chr$p_adj_bonferroni[i]) &&
          df_chr$p_adj_bonferroni[i] < retain_threshold) {
        
        start_idx <- i
        end_idx   <- i
        
        # Extend left
        while (start_idx > 1 &&
               df_chr$p_adj_bonferroni[start_idx - 1] < expand_threshold) {
          start_idx <- start_idx - 1
        }
        
        # Extend right
        while (end_idx < nrow(df_chr) &&
               df_chr$p_adj_bonferroni[end_idx + 1] < expand_threshold) {
          end_idx <- end_idx + 1
        }
        
        visited[start_idx:end_idx] <- TRUE
        
        region <- df_chr[start_idx:end_idx, ] %>%
          summarise(
            seqnames = dplyr::first(seqnames),
            start    = min(start),
            end      = max(end),
            num_ind  = max(num_ind, na.rm = TRUE),
            freq     = max(freq, na.rm = TRUE)
          )
        
        results[[length(results) + 1]] <- region
      }
    }
  }
  
  if (length(results) == 0) return(data.frame())
  bind_rows(results)
}

# ------------------------------------------------------------
# Apply anchored expansion to each coverage
# ------------------------------------------------------------

collapsed_segments_anchored_eur <- lapply(
  names(coverage_significant_results_eur),
  function(cov) {
    
    z_df <- coverage_significant_results_eur[[cov]]$df_with_zscores
    if (is.null(z_df) || nrow(z_df) == 0) return(data.frame())
    
    anchored <- extract_anchored_regions(
      zscore_df        = z_df,
      retain_threshold = 0.01,
      expand_threshold = 0.05
    )
    
    if (nrow(anchored) > 0) {
      anchored <- anchored %>%
        mutate(
          Coverage    = cov,
          region_type = "extended"
        )
    }
    
    anchored
  }
)

names(collapsed_segments_anchored_eur) <- names(coverage_significant_results_eur)

# ------------------------------------------------------------
# Extract protein-coding genes from anchored (extended) regions
# ------------------------------------------------------------

protein_coding_genes_extended_eur <- lapply(
  names(collapsed_segments_anchored_eur),
  function(cov) {
    
    regions_df <- collapsed_segments_anchored_eur[[cov]]
    if (is.null(regions_df) || nrow(regions_df) == 0) return(NULL)
    
    genes_list <- lapply(seq_len(nrow(regions_df)), function(i) {
      query_genes_in_region(
        seqname = regions_df$seqnames[i],
        start   = regions_df$start[i],
        end     = regions_df$end[i],
        mart    = ensembl_grch37
      )
    })
    
    genes_list <- genes_list[sapply(genes_list, nrow) > 0]
    if (length(genes_list) == 0) return(NULL)
    
    genes <- bind_rows(genes_list) %>%
      filter(gene_biotype == "protein_coding") %>%
      mutate(
        Coverage    = cov,
        region_type = "extended"
      ) %>%
      distinct(
        chromosome_name,
        start_position,
        end_position,
        ensembl_gene_id,
        external_gene_name,
        .keep_all = TRUE
      )
    
    genes
  }
)

names(protein_coding_genes_extended_eur) <- names(collapsed_segments_anchored_eur)

# ============================================================
#  MERGE CORE AND EXTENDED GENE SETS
# ============================================================
# Core regions represent loci with strong statistical evidence for introgression enrichment,
# while extended regions capture the broader genomic context around these signals.
# To avoid double-counting genes while preserving biological interpretation:
#   genes in core regions are prioritised
#   genes present only in extended regions are retained
#  coverage-specific information is preserved
# This yields a non-redundant, interpretable gene catalogue suitable for summary tables and downstream analyses.

genes_core_labeled <- lapply(protein_coding_genes_core_eur, function(df) {
  if (is.null(df)) return(NULL)
  df %>% mutate(region_type = "core")
})

genes_extended_labeled <- lapply(protein_coding_genes_extended_eur, function(df) {
  if (is.null(df)) return(NULL)
  df %>% mutate(region_type = "extended")
})

all_core_genes     <- bind_rows(genes_core_labeled)
all_extended_genes <- bind_rows(genes_extended_labeled)

# ------------------------------------------------------------
# Remove extended genes overlapping core genes
# ------------------------------------------------------------

extended_only_genes <- anti_join(
  all_extended_genes,
  all_core_genes,
  by = c(
    "ensembl_gene_id",
    "chromosome_name",
    "start_position",
    "end_position",
    "external_gene_name",
    "Coverage"
  )
)

# ------------------------------------------------------------
# Final non-redundant gene list
# ------------------------------------------------------------

combined_genes_eur <- bind_rows(
  all_core_genes,
  extended_only_genes
)

combined_genes_unique <- combined_genes_eur %>%
  distinct(
    chromosome_name,
    start_position,
    end_position,
    ensembl_gene_id,
    external_gene_name,
    .keep_all = TRUE
  )

# ------------------------------------------------------------
# Gene × Coverage presence table
# ------------------------------------------------------------

coverage_priority <- c("core", "extended")

gene_coverage_matrix <- combined_genes_eur %>%
  mutate(
    Gene = external_gene_name,
    region_rank = match(region_type, coverage_priority)
  ) %>%
  group_by(Gene, Coverage) %>%
  slice_min(region_rank, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(Gene, Coverage, region_type) %>%
  distinct() %>%
  mutate(Coverage = factor(Coverage, levels = dataset_order)) %>%
  pivot_wider(
    names_from  = Coverage,
    values_from = region_type,
    values_fill = ""
  ) %>%
  arrange(Gene)

# ============================================================
#  SAVE GENE × COVERAGE SUMMARY MATRIX
# ============================================================
# Enforce coverage column order
gene_coverage_matrix <- gene_coverage_matrix %>%
  dplyr::select(
    Gene,
    all_of(dataset_order)
  )


output_file <- "Table_GenePresence_Europeans_CoreExtended.csv"

write.csv(
  gene_coverage_matrix,
  file = output_file,
  row.names = FALSE
)


# ============================================================
#  MANHATTAN PLOT OF ENRICHMENT Z-SCORES
# ============================================================
# To visualise genome-wide patterns of introgression enrichment,
# we plot Z-scores across the genome for each coverage level using a Manhattan-style representation.
# Key features:
#   continuous genomic x-axis across autosomes
#   per-coverage facets to compare sequencing depth and imputation strategies
#   Bonferroni-corrected significance threshold
#   overlay of core and extended gene regions for biological interpretation
# This representation allows rapid identification of recurrent signals and coverage-specific effects.
plot_data <- bind_rows(lapply(
  names(coverage_significant_results_eur),
  function(cov) {
    
    df <- coverage_significant_results_eur[[cov]]$df_with_zscores
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    df %>%
      mutate(
        Coverage = cov,
        chr = gsub("^chr", "", as.character(seqnames))
      )
  }
)) %>%
  filter(chr %in% as.character(1:22)) %>%
  mutate(
    chr = factor(chr, levels = as.character(1:22)),
    seqnames = paste0("chr", chr)
  )

# ------------------------------------------------------------
# Build chromosome offsets for continuous x-axis
# ------------------------------------------------------------

chr_lengths <- plot_data %>%
  group_by(chr) %>%
  summarise(chr_len = max(end, na.rm = TRUE), .groups = "drop")

offset_table <- chr_lengths %>%
  arrange(as.numeric(as.character(chr))) %>%
  mutate(
    offset = lag(cumsum(as.numeric(chr_len)), default = 0),
    seqnames = paste0("chr", chr)
  ) %>%
  dplyr::select(seqnames, offset)

# ------------------------------------------------------------
# Add continuous genomic position
# ------------------------------------------------------------

plot_data <- plot_data %>%
  left_join(offset_table, by = "seqnames") %>%
  mutate(
    position_continuous =
      ((as.numeric(start) + as.numeric(end)) / 2) +
      as.numeric(offset)
  )


# ------------------------------------------------------------
# Axis labels and chromosome boundaries
# ------------------------------------------------------------
axis_df <- plot_data %>%
  group_by(chr, seqnames) %>%
  summarise(center = mean(position_continuous), .groups = "drop") %>%
  arrange(as.numeric(as.character(chr)))

chromosome_boundaries <- plot_data %>%
  group_by(chr) %>%
  summarise(chr_start = min(position_continuous), .groups = "drop") %>%
  arrange(as.numeric(as.character(chr))) %>%
  dplyr::slice(-1)


# ------------------------------------------------------------
# Gene annotation for Manhattan plot
# (keep genes per coverage — DO NOT deduplicate across coverages)
# ------------------------------------------------------------

gene_annotation <- combined_genes_eur %>%
  mutate(
    chr = factor(chromosome_name, levels = as.character(1:22)),
    seqnames = paste0("chr", chr),
    gene_mid = (start_position + end_position) / 2
  ) %>%
  left_join(offset_table, by = "seqnames") %>%
  mutate(
    position_continuous = as.numeric(gene_mid) + as.numeric(offset)
  )

# ------------------------------------------------------------
# Compute per-coverage label heights
# ------------------------------------------------------------

max_z <- plot_data %>%
  group_by(Coverage) %>%
  summarise(max_y = max(Z_score, na.rm = TRUE), .groups = "drop") %>%
  mutate(y_label = max_y + 2)

gene_annotation <- gene_annotation %>%
  left_join(max_z, by = "Coverage") %>%
  mutate(y_pos = y_label)

# ------------------------------------------------------------
# Bonferroni threshold
# ------------------------------------------------------------
bonferroni_threshold <- qnorm(
  1 - (0.01 / nrow(plot_data))
)

# ------------------------------------------------------------
# Draw the Manhattan plot
# ------------------------------------------------------------
# ------------------------------------------------------------
# Enforce coverage order for Manhattan plot
# ------------------------------------------------------------

plot_data$Coverage <- factor(
  plot_data$Coverage,
  levels = dataset_order
)

gene_annotation$Coverage <- factor(
  gene_annotation$Coverage,
  levels = dataset_order
)


manhattan_plot <- ggplot(
  plot_data,
  aes(x = position_continuous, y = Z_score)
) +
  geom_point(
    aes(color = Coverage),
    size = 0.8,
    alpha = 0.8
  ) +
  geom_vline(
    data = chromosome_boundaries,
    aes(xintercept = chr_start),
    linetype = "dotted",
    color = "grey80"
  ) +
  geom_rect(
    data = gene_annotation,
    aes(
      xmin = start_position + offset,
      xmax = end_position + offset,
      ymin = 0,
      ymax = y_pos,
      fill = region_type
    ),
    inherit.aes = FALSE,
    alpha = 0.25,
    color = NA
  ) +
  geom_text_repel(
    data = gene_annotation,
    aes(
      x = position_continuous,
      y = y_pos,
      label = external_gene_name,
      color = region_type
    ),
    inherit.aes = FALSE,
    size = 3.2,
    fontface = "italic",
    max.overlaps = Inf
  ) +
  geom_hline(
    yintercept = bonferroni_threshold,
    linetype = "dashed",
    color = "red"
  ) +
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$chr,
    expand = c(0.01, 0)
  ) +
  scale_color_manual(
    values = c(
      core = "firebrick",
      extended = "steelblue",
      dataset_colors
    )
  ) +
  scale_fill_manual(
    values = c(
      core = "firebrick",
      extended = "steelblue"
    )
  ) +
  facet_wrap(
    ~ Coverage,
    ncol = 1,
    strip.position = "right"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text.y.right = element_text(size = 11, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(y = "Z-score")

# ============================================================
#  SAVE MANHATTAN PLOT
# ============================================================

ggsave(
  "Manhattan_StrictFilters_OneTail_0.01_0.05.pdf",
  plot   = manhattan_plot,
  width  = 10,
  height = 8
)

ggsave(
  "Manhattan_StrictFilters_OneTail_0.01_0.05.svg",
  plot   = manhattan_plot,
  width  = 10,
  height = 8
)

ggsave(
  "Manhattan_StrictFilters_OneTail_0.01_0.05.png",
  plot   = manhattan_plot,
  width  = 10,
  height = 8,
  dpi    = 300
)

