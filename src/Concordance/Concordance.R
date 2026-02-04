## =======================
## Libraries
## =======================
library(ggplot2)   # all plots
library(dplyr)     # filter, mutate, group_by, summarise
library(readr)     # read_csv
library(ggpubr)    # ggarrange
library(ggh4x)     # facet_wrap2(), facetted_pos_scales()
library(grid)      # unit() used in legends
library(scales)    # alpha()

## =======================
## Metadata
## =======================
metadata_info <- read.delim(
  "../Data/InfoSamples.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

## =======================
## Color palettes
## =======================
dataset_order <- c("Original", "ImputedOC", "2X", "1X", "0.5X", "0.0625X")

dataset_colors <- c(
  Original   = "#466F9D",
  ImputedOC = "#00A087",
  `2X`       = "#D55E00",
  `1X`       = "#E69F00",
  `0.5X`     = "#56B4E9",
  `0.0625X`  = "#C77BA7"
)

category_colors <- c(
  "Archaic Shared"           = "#F4AB59",
  "Archaic New"              = "#6ABF72",
  "Not Detected As Archaic"  = "#5C9CC8"
)

individual_order <- c(
  "Ust_Ishim", "Yana1", "USR1", "Kolyma1", "KK1", "WC1", "sf12",
  "Loschbour", "Stuttgart", "NE1", "SRA62", "JP14", "BOT2016",
  "PB675", "2H10", "2H11", "atp016", "Yamnaya", "BA64", "Mota"
)

## =======================
## Load concordance data
## =======================
Concordance_Merged <- read_csv(
  "../Data/HeterozygousSensitivity.csv",
  show_col_types = FALSE
)

## =======================
## SUMMARY PLOT
## =======================
Concordance_Merged_Summary <- Concordance_Merged %>%
  filter(
    Individual != "Mota",
    Category %in% c(
      "Archaic Shared",
      "Archaic New",
      "Not Detected As Archaic"
    )
  ) %>%
  group_by(Coverage, Category) %>%
  summarise(
    mean_HET_SENSITIVITY = mean(HET_SENSITIVITY, na.rm = TRUE),
    sd_HET_SENSITIVITY   = sd(HET_SENSITIVITY, na.rm = TRUE),
    min_HET_SENSITIVITY  = min(HET_SENSITIVITY, na.rm = TRUE),
    max_HET_SENSITIVITY  = max(HET_SENSITIVITY, na.rm = TRUE),
    mean_SPECIFICITY     = mean(VAR_SPECIFICITY, na.rm = TRUE),
    sd_SPECIFICITY       = sd(VAR_SPECIFICITY, na.rm = TRUE),
    min_SPECIFICITY      = min(VAR_SPECIFICITY, na.rm = TRUE),
    max_SPECIFICITY      = max(VAR_SPECIFICITY, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Coverage = factor(Coverage, levels = dataset_order),
    Category = factor(Category, levels = names(category_colors))
  )

## ---- Heterozygous Sensitivity (summary) ----
hetsens_plot <- ggplot(
  Concordance_Merged_Summary,
  aes(
    x = Coverage,
    y = mean_HET_SENSITIVITY,
    color = Category,
    group = Category
  )
) +
  geom_ribbon(
    aes(
      ymin = min_HET_SENSITIVITY,
      ymax = max_HET_SENSITIVITY,
      fill = Category
    ),
    alpha = 0.2,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  coord_cartesian(ylim = c(0.65, 1)) +
  scale_color_manual(values = category_colors) +
  scale_fill_manual(values = category_colors) +
  labs(
    x = "",
    y = "Heterozygous Sensitivity"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = c(0.05, 0.02),
    legend.justification = c(0, 0),
    legend.background = element_rect(
      fill = alpha("white", 0.8),
      color = "black"
    ),
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

hetsens_plot

## =======================
## INDIVIDUAL-LEVEL PLOT
## =======================

dataset_order_ind <- c("ImputedOC", "2X", "1X", "0.5X", "0.0625X")

IndividualHetSens <- Concordance_Merged %>%
  filter(
    Category != "Archaic Original",
    Coverage %in% dataset_order_ind
  ) %>%
  mutate(
    Coverage   = factor(Coverage, levels = dataset_order_ind),
    Individual = factor(Individual, levels = individual_order),
    Category   = factor(
      Category,
      levels = c("Archaic New", "Archaic Shared", "Not Detected As Archaic")
    )
  )

IndividualHetSens_Plot <- ggplot(
  IndividualHetSens,
  aes(
    x = Individual,
    y = HET_SENSITIVITY,
    fill = Category
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  facet_wrap(
    ~ Coverage,
    ncol = 1,
    scales = "free_x",
    strip.position = "right"
  ) +
  scale_fill_manual(values = category_colors) +
  coord_cartesian(ylim = c(0.5, 1)) +
  labs(
    x = "",
    y = "Heterozygous Sensitivity",
    fill = ""
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey95"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

IndividualHetSens_Plot


#==============================
#  ABBA–BABA weighted GP analysis
#==============================

## ------------------------------
## Colors
## ------------------------------
category_colors_ABBA_BABA <- c(
  "ABBA" = "#F4AB59",
  "BABA" = "#5C9CC8"
)

## ------------------------------
## Load data
## ------------------------------
GP_ABBA_BABA <- read_csv(
  "../Data/GP_ABBA_BABA.csv",
  show_col_types = FALSE
)

## ------------------------------
## Weighted mean GP per individual
## ------------------------------
average_per_sample_pattern <- GP_ABBA_BABA %>%
  group_by(Individual, Coverage, Pattern) %>%
  summarise(
    Average_GP = weighted.mean(GP, w = Occurrences, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Coverage   = factor(Coverage, levels = c("ImputedOC", "2X", "1X", "0.5X", "0.0625X")),
    Individual = factor(Individual, levels = individual_order),
    Pattern    = factor(Pattern, levels = c("ABBA", "BABA"))
  )
#==============================
# Distribution across individuals
#==============================
average_per_sample_pattern_NoMota <- average_per_sample_pattern %>%
  filter(Individual != "Mota")

gp_sample_pattern_boxplot <- ggplot(
  average_per_sample_pattern_NoMota,
  aes(x = Coverage, y = Average_GP, fill = Pattern)
) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    size = 1.4,
    alpha = 0.7
  ) +
  scale_fill_manual(values = category_colors_ABBA_BABA) +
  labs(
    x = "",
    y = "Per-individual mean GP",
    fill = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = c(0.02, 0.02),
    legend.justification = c(0, 0),
    legend.background = element_rect(
      fill = alpha("white", 0.75),
      color = "black"
    ),
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

gp_sample_pattern_boxplot

#==============================
# Individual-level patterns
#==============================
IndividualABBA_BABA_Plot <- ggplot(
    average_per_sample_pattern,
    aes(x = Individual, y = Average_GP, fill = Pattern)
  ) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  facet_wrap2(
    ~ Coverage,
    ncol = 1,
    scales = "free_y"
  ) +
  scale_fill_manual(values = category_colors_ABBA_BABA) +
  facetted_pos_scales(
    y = list(
      Coverage == "ImputedOC" ~ scale_y_continuous(limits = c(0.992, 1)),
      Coverage == "2X"        ~ scale_y_continuous(limits = c(0.97, 1)),
      Coverage == "1X"        ~ scale_y_continuous(limits = c(0.96, 1)),
      Coverage == "0.5X"      ~ scale_y_continuous(limits = c(0.90, 1)),
      Coverage == "0.0625X"   ~ scale_y_continuous(limits = c(0.75, 1))
    )
  ) +
  labs(
    x = "",
    y = "GP",
    fill = ""
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey95"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

IndividualABBA_BABA_Plot


## ------------------------------------------------------------------
## Concordance analysis 
## ------------------------------------------------------------------

## Read accuracy summary (excluding Mota)
accuracy_nomota <- read.csv(
  "../Data/Accuracy_Summary_NoMota.csv",
  header = TRUE
)


plot_accuracy_nomota <- ggplot(
  accuracy_nomota,
  aes(
    x = MAF,
    y = r2,
    color = Coverage,
    linetype = SNP_Type
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  
  ## Coverage colours
  scale_color_manual(
    values = dataset_colors,
    labels = c("ImpHC" = "Imputed OC")
  ) +
  
  ## Archaic vs non-archaic SNPs
  scale_linetype_manual(
    values = c(
      "Archaic SNPs" = "dashed",
      "Non Archaic SNPs" = "solid"
    ),
    labels = c(
      "Archaic SNPs",
      "Non Archaic SNPs"
    )
  ) +
  
  ## Axes
  coord_cartesian(ylim = c(0.3, 1)) +
  labs(
    x = "1000 Genome Project MAF",
    y = "1 − MSE",
    color = NULL,
    linetype = NULL
  ) +
  
  ## Theme
  theme_bw(base_size = 14) +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    
    legend.position = c(0.40, 0.02),
    legend.justification = c(0, 0),
    legend.background = element_rect(
      fill = alpha("white", 0.7),
      color = "black"
    ),
    legend.text = element_text(size = 14),
    legend.key.width = unit(2.5, "lines")
  )

print(plot_accuracy_nomota)

## ------------------------------------------------------------------
## Multi-panel figure
## ------------------------------------------------------------------

figure3_ABC <- ggarrange(
  plot_accuracy_nomota,
  hetsens_plot,
  gp_sample_pattern_boxplot,
  ncol = 3,
  nrow = 1,
  labels = c("A", "B", "C"),
  font.label = list(size = 14, face = "bold"),
  align = "hv",
  widths = c(1, 1, 1),
  heights = c(1)
)

print(figure3_ABC)


