suppressPackageStartupMessages({
  # This script uses the magrittr pipe (%>%)
  library(scales)
  library(hexbin)
  library(purrr)
  library(GenomicRanges)
  library(furrr)  
  library(future.apply)
  library(dplyr)
  library(ggplot2)
  library(biomaRt)
  library(readr)
  library(tibble)
  library(ggpubr)
  library(tidyr)
  library(cowplot)
})

#===================================#
#              SETUP                #
#===================================#
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

segment_colors <- c(
  "Archaic New" = "#009E73",
  "Archaic Shared" = "#E69F00",
  "Non Archaic" = "#0072B2"
)

#===================================#
#   Genome length (hg19 / GRCh37)   #
#===================================#

# Connect to Ensembl GRCh37 (hg19)
ensembl_hg19 <- useMart("ENSEMBL_MART_ENSEMBL", 
                        dataset = "hsapiens_gene_ensembl", 
                        host = "https://grch37.ensembl.org")

# Get gene-level coordinates and estimate chromosome lengths
gene_coords <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position"),
  filters = "chromosome_name",
  values = as.character(1:22),
  mart = ensembl_hg19
)

chr_lengths <- gene_coords %>%
  dplyr::filter(chromosome_name %in% as.character(1:22))  %>%
  dplyr::group_by(chromosome = chromosome_name)  %>%
  dplyr::summarise(chr_length = max(end_position), .groups = "drop") %>%
  dplyr::mutate(chromosome = paste0("chr", chromosome))

total_autosomal_length <- sum(chr_lengths$chr_length)

#===================================#
#              INPUT                #
#===================================#

metadata_info <- as.data.frame(read.delim("../Data/MetadataInfo_ColorsPop.csv", header = TRUE, sep = ","))

#====Archaic Segments===#
ArchaicFile <- as.data.frame(read.delim("../Data/IntrogressionMaps.csv", header = TRUE, sep = ","))
FragmentsArchaic = ArchaicFile
FragmentsArchaic_Filtered <- FragmentsArchaic %>%
  filter(Filter == "Kept")

#==== Not Archaic ====# 
NotArchaicFile <- as.data.frame(read.delim("../Data/NonIntrogressionMaps.csv", header = TRUE, sep = ","))
FragmentsNotArchaic = NotArchaicFile
FragmentsNotArchaic_Filtered <- FragmentsNotArchaic %>%
  filter(Filter == "Kept")

#===============================#
#      Compute  Proportion      #
#===============================#

# Filtered segments
introgressed_filtered <- FragmentsArchaic_Filtered %>%
  dplyr::group_by(Individual, Coverage)  %>%
  dplyr::summarise(total_introgressed_bp = sum(Length, na.rm = TRUE), .groups = "drop")  %>%
  dplyr::mutate(
    genome_length = total_autosomal_length,
    proportion_introgressed = total_introgressed_bp / (genome_length * 2),
    Filter = "After Filtering"
  )

# Not filtered segments
introgressed_unfiltered <- FragmentsArchaic %>%
  dplyr::group_by(Individual, Coverage)  %>%
  dplyr::summarise(total_introgressed_bp = sum(Length, na.rm = TRUE), .groups = "drop")  %>%
  dplyr::mutate(
    genome_length = total_autosomal_length,
    proportion_introgressed = total_introgressed_bp / (genome_length * 2),
    Filter = "Not Filtered"
  )

# Combine 
introgressed_proportions <- bind_rows(introgressed_unfiltered, introgressed_filtered) %>%
  dplyr::mutate(
    total_Mb = total_introgressed_bp / 1e6,
    Individual = factor(Individual, levels = individual_order),
    Coverage = factor(Coverage, levels = dataset_order),
    Filter = factor(Filter, levels = c("Not Filtered", "After Filtering"))
  )

# Compute scaling factor for secondary axis
scale_factor <- max(introgressed_proportions$total_Mb, na.rm = TRUE) /
  max(introgressed_proportions$proportion_introgressed * 100, na.rm = TRUE)

#===============================#
#    Plot Genome Proportion     #
#===============================#

introgressed_proportions_plot <- ggplot() +
  geom_bar(
    data = filter(introgressed_proportions, Filter == "Not Filtered"),
    aes(x = Individual, y = total_Mb, fill = Coverage),
    stat = "identity",
    position = position_dodge(width = 0.8),
    alpha = 0.6
  ) +
  geom_bar(
    data = filter(introgressed_proportions, Filter == "After Filtering"),
    aes(x = Individual, y = total_Mb, fill = Coverage, color = Filter),
    stat = "identity",
    position = position_dodge(width = 0.8),
    linewidth = 0.4,
    show.legend = TRUE
  ) +
  scale_y_continuous(
    name = "Total length introgressed (Mb)",
    sec.axis = sec_axis(~ . / scale_factor, name = "% of the genome detected as archaic")
  ) +
  scale_fill_manual(values = dataset_colors, name = "Coverage") +
  scale_color_manual(values = c("After Filtering" = "black")) +
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(title = "", override.aes = list(fill = NA))
  ) +
  theme_bw() +
  labs(title = "", x = "") +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    panel.grid.major.x = element_blank()
  )

#=============================#
# Compute Original Recovery   #
# i.e. Shared of the Original #
#=============================#
compute_original_recovery <- function(ind, imputed_coverages) {
  message("Processing individual: ", ind)
  
  MergeFragments <- FragmentsArchaic_Filtered %>%
    filter(Individual == ind)
  
  Original_sample <- MergeFragments %>%
    filter(Coverage == "Original")
  
  gr_Original <- GRanges(seqnames = as.character(Original_sample$Chrom),
                   ranges = IRanges(start = Original_sample$start, end = Original_sample$end))
  
  map_dfr(imputed_coverages, function(cov) {
    Imputed_sample <- MergeFragments %>%
      filter(Coverage == cov)
    
    gr_Imputed <- GRanges(seqnames = as.character(Imputed_sample$Chrom),
                          ranges = IRanges(start = Imputed_sample$start, end = Imputed_sample$end))
    
    overlaps <- findOverlaps(gr_Original, gr_Imputed)
    
    Original_sample %>%
      mutate(
        segment_type = ifelse(seq_along(start) %in% queryHits(overlaps), "Shared", "Not Shared"),
        Coverage = cov
      )
  }) %>% mutate(Individual = ind)
}

#========================================#
#     Compute Recovery - NotFiltered     #
#========================================#

compute_original_recovery_NotFiltered <- function(ind, imputed_coverages) {
  message("Processing individual: ", ind)
  
  MergeFragments <- FragmentsArchaic %>%
    filter(Individual == ind)
  
  Original_sample <- MergeFragments %>%
    filter(Coverage == "Original")
  
  gr_Original <- GRanges(seqnames = as.character(Original_sample$Chrom),
                   ranges = IRanges(start =  Original_sample$start, end =  Original_sample$end))
  
  map_dfr(imputed_coverages, function(cov) {
    Imputed_sample <- MergeFragments %>%
      filter(Coverage == cov)
    
    gr_Imputed <- GRanges(seqnames = as.character(Imputed_sample$Chrom),
                          ranges = IRanges(start = Imputed_sample$start, end = Imputed_sample$end))
    
    overlaps <- findOverlaps(gr_Original, gr_Imputed)
    
    Original_sample %>%
      mutate(
        segment_type = ifelse(seq_along(start) %in% queryHits(overlaps), "Shared", "Not Shared"),
        Coverage = cov
      )
  }) %>% mutate(Individual = ind)
}

#==============================#
#     Parallel Computation     #
#==============================#

plan(multisession)
# Filtered
distance_original_recovery <- future_map_dfr(
  unique(FragmentsArchaic_Filtered$Individual),
  ~ compute_original_recovery(.x, imputed_coverages),
  .progress = TRUE
)

# NotFiltered
distance_original_recovery_NotFiltered <- future_map_dfr(
  unique(FragmentsArchaic$Individual),
  ~ compute_original_recovery_NotFiltered(.x, imputed_coverages),
  .progress = TRUE
)

plan(sequential)

#=============================#
#      Summarize Recovery     #
#=============================#

summarize_recovery <- function(df) {
  df %>%
    group_by(Individual, Coverage, segment_type) %>%
    summarise(total_length = sum(Length, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = segment_type, values_from = total_length, values_fill = 0) %>%
    mutate(proportion_shared = Shared / (Shared + `Not Shared`)) %>%
    mutate(
      Individual = factor(Individual, levels = individual_order),
      Coverage = factor(Coverage, levels = dataset_order)
    )
}

# Summarize
recovery_summary_Match <- summarize_recovery(distance_original_recovery) %>%
  mutate(Filter = "After Filtering")

recovery_summary_NotFiltered_Match <- summarize_recovery(distance_original_recovery_NotFiltered) %>%
  mutate(Filter = "Not Filtered")

# Combine
recovery_combined_Match <- bind_rows(recovery_summary_Match, recovery_summary_NotFiltered_Match) %>%
  mutate(Filter = factor(Filter, levels = c("Not Filtered", "After Filtering")))

#=============================#
#       Plot with Legend      #
#=============================#

Proportion_RecoveryMatch <- ggplot() +
  # Background bars: Not Filtered (transparent fill)
  geom_bar(
    data = filter(recovery_combined_Match, Filter == "Not Filtered"),
    aes(x = Individual, y = proportion_shared * 100, fill = Coverage),
    stat = "identity",
    position = position_dodge(width = 0.8),
    alpha = 0.6
  ) +
  
  # Foreground bars: Filtered (black outline)
  geom_bar(
    data = filter(recovery_combined_Match, Filter == "After Filtering"),
    aes(x = Individual, y = proportion_shared * 100, fill = Coverage, color = Filter),
    stat = "identity",
    position = position_dodge(width = 0.8),
    linewidth = 0.4
  ) +
  # Manual fills and outlines
  scale_fill_manual(values = dataset_colors, name = "Coverage") +
  scale_color_manual(values = c("After Filtering" = "black")) +
  
  # Legend tweaks to show border meaning
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(title = "", override.aes = list(fill = NA))
  ) +
  
  theme_bw() +
  labs(title = "", x = "",
       y = "% of Original Introgressed Segments Recovered") +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    panel.grid.major.x = element_blank()
  )


#============================#
#     Original Similarity    #
#============================#
#============================#
#     Prepare plot data      #
#============================#

# Add "Not detected as Archaic" segments
df_not_archaic <- FragmentsNotArchaic %>%
  mutate(
    Segment_Type = "Non Archaic",
    Chrom = as.character(Chrom),
    Coverage = factor(Coverage, levels = dataset_order)
  )

# Add filter status annotation
df_archaic <- FragmentsArchaic %>%
  mutate(Filter_Status = ifelse(Filter == "Kept", "After Filtering", "Not Filtered"),
         Chrom = as.character(Chrom),
         Segment_Type = recode(Segment_Type,
                                "Original" = "Archaic Shared",
                                .default = Segment_Type),
         Coverage = factor(Coverage, levels = dataset_order)
  )


# Filtered segments (for outline only, new legend category)
df_filtered_archaic <- df_archaic %>%
  filter(Filter == "Kept") %>%
  mutate(Segment_Type = "After Filtering",
         Coverage   = factor(Coverage, levels = dataset_order)
         )  
binwidth <- 0.002 

#============================#
# STEP 3: Plotting
#============================#

similarity_filtered_plot <- ggplot() +
  
  # 1. Not detected as archaic (underlay)
  geom_histogram(
    data = df_not_archaic,
    aes(x = SimilarityOriginal_to_archaic, fill = Segment_Type),
    binwidth = binwidth,
    boundary = 0.86,
    position = "identity",
    alpha = 0.9
  ) +
  
  # 2. Archaic segments (stacked fill)
  geom_histogram(
    data = df_archaic ,
    aes(x = SimilarityOriginal_to_archaic, fill = Segment_Type),
    binwidth = binwidth,
    boundary = 0.86,
    position = "stack",
    alpha = 0.9
  ) +
  
  # 3. Filtered Archaic outline only (with legend via `color`)
  geom_histogram(
    data = df_filtered_archaic,
    aes(x = SimilarityOriginal_to_archaic, color = Segment_Type),
    binwidth = binwidth,
    boundary = 0.86,
    position = "identity",
    fill = NA,
    linewidth = 0.4
  ) +
  
  # 4. Styling
  scale_fill_manual(values = segment_colors) +
  scale_color_manual(
    values = c("After Filtering" = "black"),
    labels = c("After Filtering" = "After Filtering")
  ) +
  guides(
    fill = guide_legend(override.aes = list(alpha = 0.9)),
    color = guide_legend(override.aes = list(fill = NA, linewidth = 0.6))
  ) +
  # facet_wrap(~ Coverage, ncol = 3) 
  facet_wrap(~ Coverage, ncol = 2) + 
  coord_cartesian(xlim = c(0.90, 1.01)) +
  labs(
    x = "Similarity of high-coverage genomes to nearest archaic genome",
    y = "Count",
    fill = "",
    color = ""
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    legend.title = element_blank(),
    strip.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.margin = margin(1, 1, 1, 1)
  )

#==============================================#
#  Supp figures  Similarity vs Distance        #
#==============================================#
SimDist_Archaic <- ggplot(df_archaic, aes(x = SimilarityImp_to_archaic, y = log1p(DistanceImp_to_archaic))) +
  geom_hex(bins = 1000, alpha = 0.6) +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),
    trans = "log",
    breaks = c(1, 10, 100, 1000),
    labels = c("1", "10", "100", "1000")
  ) +
  geom_vline(xintercept = 0.99, linetype = "dashed", color = "red", linewidth = 0.5) +
  geom_hline(yintercept = log1p(0.2), linetype = "dashed", color = "red", linewidth = 0.5) +
  annotate("text", x = 0.90, y = log1p(0.25), label = "0.2", color = "red", hjust = -0.4, size = 4) +
  annotate("text", x = 0.95, y = 0, label = "0.99", color = "red", hjust = -0.4, size = 4) +
  facet_wrap(~ Coverage, nrow = 1) +
  theme_bw() +
  labs(
    x = "Similarity to closest archaic reference genome",
    y = "Distance to closest archaic reference genome (Log Scale)",
    title = "Archaic Segments"
  ) +
  theme(
    legend.position = "right",
    text = element_text(size = 14)
  )

#=======================================#
# Define Plot for Not Archaic Segments  #
#=======================================#

df_not_archaic_filter <- df_not_archaic %>%
  filter(Length > 40000)

SimDist_NotArchaic <- ggplot(df_not_archaic, aes(x = SimilarityImp_to_archaic, y = log1p(DistanceImp_to_archaic))) +
  geom_hex(bins = 1000, alpha = 0.6) +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),
    trans = "log",
    breaks = c(1, 10, 100, 1000),
    labels = c("1", "10", "100", "1000")
  ) +
  geom_vline(xintercept = 0.99, linetype = "dashed", color = "red", linewidth = 0.5) +
  geom_hline(yintercept = log1p(0.2), linetype = "dashed", color = "red", linewidth = 0.5) +
  annotate("text", x = 0.90, y = log1p(0.25), label = "0.2", color = "red", hjust = -0.4, size = 4) +
  annotate("text", x = 0.95, y = 0, label = "0.99", color = "red", hjust = -0.4, size = 4) +
  facet_wrap(~ Coverage, nrow = 1) +
  theme_bw() +
  labs(
    x = "Similarity to closest archaic reference genome",
    y = "Distance to closest archaic reference genome (Log Scale)",
    title = "Archaic Segments"
  ) +
  theme(
    legend.position = "right",
    text = element_text(size = 14)
  )

#==============================#
#  Compute Common Axes Limits  #
#==============================#

x_limits <- range(
  c(df_archaic$SimilarityImp_to_archaic, df_not_archaic$SimilarityImp_to_archaic),
  na.rm = TRUE
)

y_limits <- range(
  c(
    log1p(df_archaic$DistanceImp_to_archaic),
    log1p(df_not_archaic$DistanceImp_to_archaic)
  ),
  na.rm = TRUE
)

#==============================#
#   Apply Shared Axis Settings #
#==============================#

SimDist_Archaic <- SimDist_Archaic +
  coord_cartesian(xlim = c(0.90, 1), ylim = c(0, 1)) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(5, 5, 5, 5),
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) 

SimDist_NotArchaic <- SimDist_NotArchaic +
  coord_cartesian(xlim = c(0.90, 1), ylim = c(0, 1)) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(5, 5, 5, 5),
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) 
#===================#
#   Combine Plots   #
#===================#

combined_plot <- ggarrange(
  SimDist_Archaic, SimDist_NotArchaic,
  ncol = 1, nrow = 2,
  align = "v",
  common.legend = TRUE,
  legend = "top",
  labels = c("A", "B"),
  font.label = list(size = 14, face = "bold")
)

final_plot_with_shared_axes <- annotate_figure(
  combined_plot,
  left = text_grob(
    "Distance (Imp) to closest archaic reference genome (natural log scale)",
    rot = 90, size = 14, vjust = 0.5
  ),
  bottom = text_grob(
    "Similarity (Imp) to closest archaic reference genome",
    size = 14, vjust = 0.5
  )
)


#==========================#
#  Length vs Similarity    #
#     Supp figure          #
#==========================#
log_breaks <- c(1e4, 1e5, 1e6, 1e7)

df_plot <- FragmentsArchaic %>%
  filter(
    Coverage %in% c("ImputedOC", "2X", "1X", "0.5X", "0.0625X", "Original"),
    Length >= 40000
  ) %>%
  mutate(
    Coverage = factor(Coverage, levels = dataset_order),
    Plot_Group = case_when(
      Coverage == "Original" ~ "Original",
      TRUE ~ Segment_Type
    ),
    Plot_Group = factor(
      Plot_Group,
      levels = c("Archaic New", "Archaic Shared", "Original")
    )
  )

# Duplicate Original rows across imputed levels
df_Original_overlay <- expand.grid(
  row_id = seq_len(nrow(df_plot %>% filter(Coverage == "Original"))),
  Coverage_overlay = c("ImputedOC", "2X", "1X", "0.5X", "0.0625X")
) %>%
  left_join(
    df_plot %>% filter(Coverage == "Original") %>% mutate(row_id = row_number()),
    by = "row_id"
  ) %>%
  dplyr::select(-row_id) %>%
  mutate(
    Coverage = factor(Coverage_overlay, levels = levels(df_plot$Coverage)),
    Plot_Group = "Original"
  ) %>%
  dplyr::select(-Coverage_overlay)

df_imputed_combined <- df_plot %>%
  filter(Coverage != "Original") %>%
  bind_rows(df_Original_overlay) %>%
  mutate(Coverage = factor(Coverage, levels = dataset_order))

mean_lengths_annotated <- df_imputed_combined %>%
  group_by(Coverage, Plot_Group) %>%
  summarise(
    mean_length = mean(Length, na.rm = TRUE),
    y = max(table(cut(log10(Length), breaks = 60))),  # rough peak count estimate
    .groups = "drop"
  ) %>%
  mutate(
    x = mean_length,
    y = y * 1.05,  # Push a bit above histogram
    label = paste0(round(mean_length / 1000), " kb")
  )

PanelA_LengthOverlay <- ggplot() +
  geom_histogram(
    data = df_imputed_combined,
    aes(x = Length, fill = Plot_Group),
    bins = 60,
    alpha = 0.6,
    position = "identity"
  ) +
  geom_vline(
    data = mean_lengths_annotated,
    aes(xintercept = mean_length, color = Plot_Group),
    linetype = "dashed",
    linewidth = 0.6,
    show.legend = FALSE
  ) +
  geom_text(
    data = mean_lengths_annotated,
    aes(x = mean_length, y = y, label = label, color = Plot_Group),
    vjust = -0.5, size = 4, fontface = "bold",
    show.legend = FALSE
  ) +
  scale_x_log10(breaks = log_breaks, labels = scales::comma_format()) +
  scale_fill_manual(
    name = "Group",
    values = c(
      "Original" = "#D55E00",
      "Archaic New" = "#009E73",
      "Archaic Shared" = "#0072B2"
    )
  ) +
  scale_color_manual(
    name = "Group",
    values = c(
      "Original" = "#D55E00",
      "Archaic New" = "#009E73",
      "Archaic Shared" = "#0072B2"
    )
  ) +
  facet_wrap(~ Coverage, nrow = 1, labeller = labeller(
    ImputedOC = "Imputed OC"
  )) +
  theme_bw(base_size = 14) +
  labs(
    x = "Segment Length (bp, log scale)",
    y = "Count"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 14)
  )

PanelA_LengthOverlay

df_panelB <- FragmentsArchaic %>%
  filter(
    Coverage %in% c("ImputedOC", "2X", "1X", "0.5X", "0.0625X"),
    Segment_Type %in% c("Archaic New", "Archaic Shared"),
    !is.na(SimilarityImp_to_archaic),
    Length >= 40000
  ) %>%
  mutate(
    Coverage = factor(Coverage, levels = dataset_order),
    Segment_Type = factor(Segment_Type, levels = c("Archaic New", "Archaic Shared"))
  )

mean_similarity_panelB <- df_panelB %>%
  group_by(Coverage, Segment_Type) %>%
  summarise(
    mean_similarity = mean(SimilarityImp_to_archaic),
    max_length = max(Length),
    .groups = "drop"
  ) %>%
  mutate(
    x = max_length * 0.8,
    y = 0.902,
    label = paste0("mean = ", round(mean_similarity, 3))
  )

PanelB_LengthVsSimilarity <- ggplot(df_panelB, aes(x = Length, y = SimilarityImp_to_archaic)) +
  geom_hex(bins = 80, alpha = 0.8) +
  ylim(0.9, 1) +
  scale_x_log10(breaks = log_breaks, labels = scales::comma_format()) +
  geom_text(
    data = mean_similarity_panelB,
    aes(x = x, y = y, label = label, color = Segment_Type),
    hjust = 1, size = 4, fontface = "bold",
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c("Archaic New" = "#009E73", "Archaic Shared" = "#0072B2")
  ) +
  scale_fill_gradientn(
    colors = c("blue", "green", "yellow", "red"),
    trans = "log",
    breaks = c(1, 10, 100, 1000),
    labels = c("1", "10", "100", "1000"),
    name = "count"
  ) +
  facet_grid(Segment_Type ~ Coverage) +
  theme_bw(base_size = 14) +
  labs(
    x = "Segment Length (bp, log scale)",
    y = "Similarity to Archaic"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 14))


PanelA_LengthDistributions <- PanelA_LengthOverlay +
  theme(axis.title.x = element_blank())

# Align panels vertically (this aligns the x-axes)
Figure_Length_Similarity <- plot_grid(
  PanelA_LengthDistributions, PanelB_LengthVsSimilarity,
  ncol = 1,
  align = "v",      # align vertically => align x
  axis  = "lr",     # keep left/right margins consistent
  rel_heights = c(0.8, 1.2)
)
