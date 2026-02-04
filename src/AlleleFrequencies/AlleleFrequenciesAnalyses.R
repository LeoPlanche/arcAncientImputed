################################################################################
# Allele-frequency based analyses
# D-statistics and f4-ratio plots across coverage and genotyping strategies
################################################################################

# ============================== #
#           Libraries            #
# ============================== #

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# ============================== #
#        Global settings         #
# ============================== #

dataset_order <- c(
  "Genotype Calls",
  "Original Coverage",
  "2X", "1X", "0.5X", "0.0625X"
)

individual_order <- c(
  "Ust_Ishim", "Yana1", "USR1", "Kolyma1", "KK1", "WC1", "sf12",
  "Loschbour", "Stuttgart", "NE1", "SRA62", "JP14", "BOT2016",
  "PB675", "2H10", "2H11", "atp016", "Yamnaya", "BA64", "Mota"
)

datatype_colors <- c(
  "Genotype Calls" = "#7A7A7A",
  "Imputed Raw" = "#4C78A8",
  "Imputed Filtered GP99" = "#E45756",
  "Imputed Filtered MAF05" = "lightblue",
  "Imputed GP Weighted" = "#9D5EA1",
  "Pseudohaploid" = "#59A14F",
  "Pseudohaploid Only Transversions" = "#F1CE63"
)

dataset_colors <- c(
  `Original Coverage` = "#466F9D",
  `Imputed Original Coverage` = "#00A087",
  `2X` = "#D55E00",
  `1X` = "#E69F00",
  `0.5X` = "#56B4E9",
  `0.0625X` = "#C77BA7"
)

# ============================== #
#        D STATISTICS            #
# ============================== #

Dstat_df <- read_tsv(
  "../Data/Dstatistics.txt",
  col_types = cols(
    Sample      = col_character(),
    Coverage    = col_character(),
    D           = col_double(),
    Z.score     = col_double(),
    Filter      = col_character(),
    `Data Type` = col_character(),
    Software    = col_character()
  )
)

Dstat_plot_df <- Dstat_df %>%
  mutate(
    dataset_order = case_when(
      `Data Type` == "Genotype Calls" & Coverage == "Original Coverage" ~ "Genotype Calls",
      TRUE ~ Coverage
    ),
    Highlight = if_else(Sample == "Mota", "Mota", "Other")
  )

# ============================== #
#        MAIN D FIGURES          #
# ============================== #

Fig1_df <- Dstat_plot_df %>%
  filter(
    Software == "Admixtools",
    !Filter %in% c(
      "Imputed GP Weighted",
      "Imputed Filtered MAF05",
      "Genotype Calls Only Transversions"
    )
  )

Fig1_df$dataset_order <- factor(Fig1_df$dataset_order, levels = dataset_order)

datasets <- levels(Fig1_df$dataset_order)
axis_colors <- unname(dataset_colors[datasets])
x_lines <- setdiff(seq(1.5, length(datasets) - 0.5, by = 1), 1.5)

Fig1_df$Filter <- factor(
  Fig1_df$Filter,
  levels = c(
    "Genotype Calls",
    "Imputed Raw",
    "Imputed Filtered GP99",
    "Pseudohaploid",
    "Pseudohaploid Only Transversions"
  )
)

# Figure 1C – D
Dstat_plot_WithDots <- ggplot(
  Fig1_df,
  aes(x = dataset_order, y = -D, fill = Filter)
) +
  geom_boxplot(position = position_dodge2(preserve = "single"), outliers = FALSE) +
  geom_point(
    color = "black", size = 1.5, alpha = 0.5,
    position = position_dodge(width = 0.75),
    show.legend = FALSE
  ) +
  geom_point(
    data = subset(Fig1_df, Highlight == "Mota"),
    color = "black", size = 4, shape = 22,
    position = position_dodge(width = 0.75),
    show.legend = FALSE
  ) +
  geom_vline(xintercept = x_lines, linetype = "dashed", color = "grey60") +
  scale_fill_manual(values = datatype_colors) +
  labs(x = "", y = expression(italic("D")), fill = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(color = axis_colors),
    text = element_text(size = 14),
    panel.grid.major.x = element_blank()
  )

# Figure S1A – |Z|
Dstat_Zscore_plot_WithDots <- ggplot(
  Fig1_df,
  aes(x = dataset_order, y = abs(Z.score), fill = Filter)
) +
  geom_boxplot(position = position_dodge2(preserve = "single"), outliers = FALSE) +
  geom_hline(yintercept = 3.3, linetype = "dashed", color = "red") +
  geom_point(
    color = "black", size = 1.5, alpha = 0.5,
    position = position_dodge(width = 0.75),
    show.legend = FALSE
  ) +
  geom_point(
    data = subset(Fig1_df, Highlight == "Mota"),
    color = "black", size = 4, shape = 22,
    position = position_dodge(width = 0.75),
    show.legend = FALSE
  ) +
  geom_vline(xintercept = x_lines, linetype = "dashed", color = "grey60") +
  scale_fill_manual(values = datatype_colors) +
  labs(x = "", y = expression(italic("|Z|")), fill = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(color = axis_colors),
    text = element_text(size = 14),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  )

# ============================== #
#        MAF05 FIGURES           #
# ============================== #

AdmixtoolsMAF05_selected <- Dstat_plot_df %>%
  filter(
    Software == "Admixtools",
    Filter == "Imputed Filtered MAF05"
  )

AdmixtoolsMAF05_selected$dataset_order <- factor(
  AdmixtoolsMAF05_selected$dataset_order,
  levels = dataset_order
)
AdmixtoolsMAF05_selected$dataset_order <-
  droplevels(AdmixtoolsMAF05_selected$dataset_order)

datasets <- levels(AdmixtoolsMAF05_selected$dataset_order)

# Special color remapping for MAF05
dataset_colors_MAF05 <- c(
  `Genotype Calls`    = unname(dataset_colors["Original Coverage"]),
  `Original Coverage` = unname(dataset_colors["Imputed Original Coverage"]),
  `2X`                = unname(dataset_colors["2X"]),
  `1X`                = unname(dataset_colors["1X"]),
  `0.5X`              = unname(dataset_colors["0.5X"]),
  `0.0625X`           = unname(dataset_colors["0.0625X"])
)

axis_colors_MAF05 <- dataset_colors_MAF05[datasets]

# Figure S1B – D (MAF05)
Dstat_plot_WithDots_MAF05 <- ggplot(
  AdmixtoolsMAF05_selected,
  aes(x = dataset_order, y = -D, fill = Filter)
) +
  geom_boxplot(
    position = position_dodge2(preserve = "single", padding = 0),
    outliers = FALSE
  ) +
  geom_point(
    color = "black",
    size = 1.5,
    shape = 16,
    alpha = 0.5,
    position = position_dodge(width = 0.75),
    show.legend = FALSE
  ) +
  geom_point(
    data = subset(AdmixtoolsMAF05_selected, Highlight == "Mota"),
    color = "black",
    size = 4,
    shape = 22,
    position = position_dodge(width = 0.75),
    show.legend = FALSE
  ) +
  geom_vline(
    xintercept = x_lines,
    linetype = "dashed",
    color = "grey60",
    linewidth = 0.5
  ) +
  scale_fill_manual(values = datatype_colors) +
  scale_x_discrete(
    limits = datasets,
    labels = c(
      "Genotype Calls"    = "Original Coverage\nGenotype Calls",
      "Original Coverage" = "Original Coverage\nImputed",
      "2X" = "2X",
      "1X" = "1X",
      "0.5X" = "0.5X",
      "0.0625X" = "0.0625X"
    )
  ) +
  labs(
    x = "",
    y = expression(italic("D")),
    fill = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(color = axis_colors_MAF05),
    text = element_text(size = 14),
    panel.grid.major.x = element_blank()
  ) +
  guides(fill = guide_legend(byrow = TRUE))

Zscore_plot_WithDots_MAF05 <- ggplot(
  AdmixtoolsMAF05_selected,
  aes(x = dataset_order, y = abs(Z.score), fill = Filter)
) +
  geom_boxplot(
    position = position_dodge2(preserve = "single", padding = 0),
    outliers = FALSE
  ) +
  
  # Dots for all individuals
  geom_point(
    aes(x = dataset_order, y = abs(Z.score), fill = Filter),
    color = "black",
    size = 1.5,
    shape = 16,
    alpha = 0.5,
    position = position_dodge(width = 0.75),
    show.legend = FALSE
  ) +
  
  # Highlight Mota
  geom_point(
    data = subset(AdmixtoolsMAF05_selected, Highlight == "Mota"),
    aes(x = dataset_order, y = abs(Z.score), fill = Filter),
    color = "black",
    size = 4,
    shape = 22,
    stroke = 1,
    position = position_dodge(width = 0.75),
    show.legend = FALSE
  ) +
  
  # Significance threshold
  geom_hline(
    yintercept = 3.3,
    color = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  
  scale_fill_manual(values = datatype_colors) +
  
  scale_x_discrete(
    limits = datasets,
    labels = c(
      "Genotype Calls"    = "Original Coverage\nGenotype Calls",
      "Original Coverage" = "Original Coverage\nImputed",
      "2X" = "2X",
      "1X" = "1X",
      "0.5X" = "0.5X",
      "0.0625X" = "0.0625X"
    )
  ) +
  
  labs(
    x = "",
    y = expression(italic("|Z|")),
    fill = ""
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(color = axis_colors_MAF05),
    text = element_text(size = 14),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(byrow = TRUE))


# ============================== #
#        GP-BASED FIGURE         #
# ============================== #

GP_df <- Dstat_plot_df %>%
  filter(
    Filter %in% c(
      "Genotype Calls",
      "Imputed Raw",
      "Imputed GP Weighted",
      "Imputed Filtered GP99"
    )
  )

GP_df$dataset_order <- factor(GP_df$dataset_order, levels = dataset_order)
GP_df$dataset_order <- droplevels(GP_df$dataset_order)

datasets <- levels(GP_df$dataset_order)
axis_colors <- unname(dataset_colors[datasets])

# Figure S3 – GP based
Dstat_GPBased_WithDots <- ggplot(
  GP_df,
  aes(x = dataset_order, y = -D, fill = Filter)
) +
  geom_boxplot(position = position_dodge2(preserve = "single"), outliers = FALSE) +
  geom_point(
    color = "black", size = 1.5, alpha = 0.5,
    position = position_dodge(width = 0.75),
    show.legend = FALSE
  ) +
  geom_point(
    data = subset(GP_df, Highlight == "Mota"),
    color = "black", size = 4, shape = 22,
    position = position_dodge(width = 0.75),
    show.legend = FALSE
  ) +
  scale_fill_manual(values = datatype_colors) +
  labs(x = "", y = expression(italic("D")), fill = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(color = axis_colors),
    text = element_text(size = 14),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  )

# ============================== #
#        PER-INDIVIDUAL D        #
# ============================== #

Dstat_plot_Ind_df <- Fig1_df %>%
  mutate(
    Sample = factor(Sample, levels = individual_order),
    Filter = factor(Filter, levels = names(datatype_colors))
  )

datasets <- levels(Dstat_plot_Ind_df$dataset_order)
axis_colors <- unname(dataset_colors[datasets])

# Figure S2
Dstat_plot_Ind <- ggplot(
  Dstat_plot_Ind_df,
  aes(x = dataset_order, y = -D, fill = Filter)
) +
  geom_col(position = position_dodge(width = 1)) +
  facet_wrap(~ Sample, ncol = 1, strip.position = "right") +
  scale_fill_manual(values = datatype_colors) +
  labs(x = "", y = expression(italic("D")), fill = "") +
  theme_bw() +
  theme(
    axis.text.x = element_text(color = axis_colors),
    text = element_text(size = 14),
    panel.grid.major.x = element_blank(),
    strip.text.y.right = element_text(angle = 0)
  )

# ============================== #
#           F4 RATIO             #
# ============================== #

f4_df <- read_tsv(
  "../Data/f4Ratio.txt",
  col_types = cols(
    Sample      = col_character(),
    Coverage    = col_character(),
    alpha       = col_double(),
    Filter      = col_character(),
    `Data Type` = col_character()
  )
)

f4_plot_df <- f4_df %>%
  mutate(
    dataset_order = case_when(
      `Data Type` == "Genotype Calls" & Coverage == "Original Coverage" ~ "Genotype Calls",
      TRUE ~ Coverage
    ),
    Highlight = if_else(Sample == "Mota", "Mota", "Other")
  )

Fig1D_df <- f4_plot_df %>%
  filter(
    !Filter %in% c(
      "Imputed GP Weighted",
      "Imputed Filtered MAF05",
      "Genotype Calls Only Transversions"
    )
  )

Fig1D_df$dataset_order <- factor(Fig1D_df$dataset_order, levels = dataset_order)

datasets <- levels(Fig1D_df$dataset_order)
axis_colors <- unname(dataset_colors[datasets])
x_lines <- setdiff(seq(1.5, length(datasets) - 0.5, by = 1), 1.5)

Fig1D_df$Filter <- factor(
  Fig1D_df$Filter,
  levels = c(
    "Genotype Calls",
    "Imputed Raw",
    "Imputed Filtered GP99",
    "Pseudohaploid",
    "Pseudohaploid Only Transversions"
  )
)

# Figure 1D – f4-ratio
Fig1D_f4Ratio_plot <- ggplot(
  Fig1D_df,
  aes(x = dataset_order, y = alpha, fill = Filter)
) +
  geom_boxplot(position = position_dodge2(preserve = "single"), outliers = FALSE) +
  geom_point(
    color = "black", size = 1.5, alpha = 0.5,
    position = position_dodge(width = 0.75),
    show.legend = FALSE
  ) +
  geom_point(
    data = subset(Fig1D_df, Highlight == "Mota"),
    color = "black", size = 4, shape = 22,
    position = position_dodge(width = 0.75),
    show.legend = FALSE
  ) +
  geom_vline(xintercept = x_lines, linetype = "dashed", color = "grey60") +
  scale_fill_manual(values = datatype_colors) +
  labs(x = "", y = expression(alpha), fill = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(color = axis_colors),
    text = element_text(size = 14),
    panel.grid.major.x = element_blank()
  )

# ============================== #
#      PER-INDIVIDUAL F4         #
# ============================== #

F4ratio_plot_Ind_df <- Fig1D_df %>%
  mutate(
    Sample = factor(Sample, levels = individual_order),
    Filter = factor(Filter, levels = names(datatype_colors))
  )

datasets <- levels(F4ratio_plot_Ind_df$dataset_order)
axis_colors <- unname(dataset_colors[datasets])

# Figure S4
F4ratio_plot_Ind <- ggplot(
  F4ratio_plot_Ind_df,
  aes(x = dataset_order, y = alpha, fill = Filter)
) +
  geom_col(position = position_dodge(width = 1)) +
  facet_wrap(~ Sample, ncol = 1, strip.position = "right") +
  scale_fill_manual(values = datatype_colors) +
  labs(
    x = "",
    y = expression(italic("f4-ratio (") * alpha * ")"),
    fill = ""
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(color = axis_colors),
    text = element_text(size = 14),
    panel.grid.major.x = element_blank(),
    strip.text.y.right = element_text(angle = 0)
  )

