
library(ggplot2)
# library(cowplot)   # To have ggplots side-by-side: plot_grid()

library(dplyr)
library(tidyr)     # separate(), gather() and spread() functions
library(tibble)    # for rownames_to_column() and column_to_rownames()

# library(ComplexHeatmap)
library(circlize)   # colorRamp2()


###

DATA_DIR <- "../data/"
OUT_DIR <- "../images/"

LOG_COLORS <- colorRamp2(c(0, 2, 10), c("white", "gray", "black"))
LOG_COLORS_RED <- colorRamp2(c(0, 2, 10), c("white", "yellow", "red"))

# We used 45480 as the total number of the human genes that was computed
# as the sum of the protein-coding (19954), long (17957) and small (7569) non-coding RNA genes
# as annotated in the GENCODE version 35.
NUM_HUMAN_GENES <- 45480

theme_set(theme_bw(base_size = 19))  # increase the font size: https://stackoverflow.com/a/11955412/310453

get_good_aso <- function(df)
{
  df %>%
    filter(hgd_pvalue < 0.01 | gsea_pvalue < 0.01) %>%
    pull(aso_id)
}


read_DE_summary_file <- function(fn) {
  read.delim(fn, as.is=TRUE, header=TRUE) %>%
    # p002@AC013394.2  =>  AC013394.2
    mutate(geneSymbol = gsub('^p\\d+\\@', '', geneSymbol)) %>%
    # SPNS1;p002@RP11-264B17   =>   SPNS1
    mutate(geneSymbol = gsub(';.*$', '', geneSymbol)) %>%
    # AC013394.2  =>  AC013394
    mutate(geneSymbol = gsub('\\.\\d+$', '', geneSymbol))
}

get_hgd_pvalue <- function(all_balls, red_balls, drawn_balls, drawn_red) {
  phyper(
    drawn_red, # number of red balls drawn without replacement from an urn
    red_balls, # the number of red balls in the urn
    all_balls - red_balls, # the number of white balls in the urn
    drawn_balls, # the number of balls drawn from the urn.
    lower.tail = FALSE
  )
}
