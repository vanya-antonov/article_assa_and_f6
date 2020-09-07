
source('lib.R')

# https://www.r-graph-gallery.com/14-venn-diagramm.html
# install.packages('VennDiagram')
library(VennDiagram)

###

co_df <- read.delim(paste0(DATA_DIR, '2_pvalues.co.txt'), as.is=TRUE, header=FALSE)
colnames(co_df) <- c('aso_id', 'all_trxs', 'de_trxs', 'assa_hits_all', 'assa_hits_de', 'hgd_pvalue', 'gsea_pvalue')

post_df <- read.delim(paste0(DATA_DIR, '2_pvalues.post.txt'), as.is=TRUE, header=FALSE)
colnames(post_df) <- c('aso_id', 'all_trxs', 'de_trxs', 'assa_hits_all', 'assa_hits_de', 'hgd_pvalue', 'gsea_pvalue')

get_good_aso <- function(df)
{
  df %>%
    filter(hgd_pvalue < 0.01 | gsea_pvalue < 0.01) %>%
    pull(aso_id)
}

good_co_v <- get_good_aso(co_df)
good_post_v <- get_good_aso(post_df)

# Do not create log file: https://stackoverflow.com/a/36812214/310453
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

venn.diagram(
  x = list(good_co_v, good_post_v),
  category.names = c(
    paste("Co-transcriptional\nASOs, n = ", length(good_co_v)),
    paste("Post-transcriptional\nASOs, n = ", length(good_post_v))),
  cat.fontface = "bold",
  cat.pos = c(-32, 27),
  cat.dist = c(0.055, 0.055),
  # height = 30,
  # width = 40,  
  filename = paste0(OUT_DIR, 'two_pvalues_venn.png'),
  output=TRUE
)

