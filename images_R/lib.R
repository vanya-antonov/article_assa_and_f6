
library(ggplot2)
# library(cowplot)   # To have ggplots side-by-side: plot_grid()

library(dplyr)
library(tidyr)     # separate(), gather() and spread() functions
library(tibble)    # for rownames_to_column() and column_to_rownames()

# library(ComplexHeatmap)
# library(circlize)   # colorRamp2()


###

DATA_DIR <- "../data/"
OUT_DIR <- "../images/"

LOG_COLORS <- colorRamp2(c(0, 2, 10), c("white", "gray", "black"))
LOG_COLORS_RED <- colorRamp2(c(0, 2, 10), c("white", "yellow", "red"))

# theme_set(theme_bw(base_size = 19))  # increase the font size: https://stackoverflow.com/a/11955412/310453

get_good_aso <- function(df)
{
  df %>%
    filter(hgd_pvalue < 0.01 | gsea_pvalue < 0.01) %>%
    pull(aso_id)
}

