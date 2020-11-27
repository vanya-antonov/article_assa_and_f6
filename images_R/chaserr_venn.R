
source('lib.R')

library(VennDiagram)

###

aso_07_df <- read_DE_summary_file(paste0(DATA_DIR, 'ASO_G0272888_AD_07.DE_Summary'))
aso_10_df <- read_DE_summary_file(paste0(DATA_DIR, 'ASO_G0272888_AD_10.DE_Summary'))

from_mouse_df <- read.delim(paste0(DATA_DIR, 'chaserr_orthologs.txt'), as.is=TRUE, header=TRUE) %>%
  # AC013394.2  =>  AC013394
  mutate(geneSymbol = gsub('\\.\\d+$', '', human_gene_name))
head(from_mouse_df)

aso_07_v <- unique(aso_07_df$geneSymbol)
aso_10_v <- unique(aso_10_df$geneSymbol)
mouse_v <- unique(from_mouse_df$geneSymbol)

common3_v <- Reduce(intersect, list(aso_07_v, aso_10_v, mouse_v))

# Do not create log file: https://stackoverflow.com/a/36812214/310453
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

?venn.diagram

venn.diagram(
  x = list(aso_07_v, aso_10_v, mouse_v),
  category.names = c(
    paste("ASO_G0272888_AD_07,\nn = ", length(aso_07_v)),
    paste("ASO_G0272888_AD_10,\nn = ", length(aso_10_v)),
    paste("Orthologs of the mEFs Chaserr-/- DEGs, n = ", length(mouse_v))),
  cat.fontface = "bold",
  cat.pos = c(-15, 15, 180),
  cat.dist = c(0.05, 0.05, 0.02),
  filename = paste0(OUT_DIR, 'other/chaserr_venn_3.png'),
  output=TRUE
)

venn.diagram(
  x = list(aso_07_v, aso_10_v),
  category.names = c(
    paste("ASO-07,\nn = ", length(aso_07_v)),
    paste("ASO-10,\nn = ", length(aso_10_v))),
  cat.fontface = "bold",
  cat.pos = c(-40, 40),
  cat.dist = c(0.04, 0.04),
  filename = paste0(OUT_DIR, 'other/chaserr_venn_aso.png'),
  output=TRUE
)



# Save common genes to file
aso_07_df %>%
  filter(geneSymbol %in% common3_v) %>%
  select(geneSymbol) %>%
  unique() %>%
  write.table(file = paste0(DATA_DIR, "chaserr_common3.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")



?phyper

n_07 = length(aso_07_v)
n_10 = length(aso_10_v)
both_aso = length(intersect(aso_07_v, aso_10_v))
all_genes = 20000
phyper(
  both_aso, # number of white balls drawn without replacement from an urn
  n_10, # the number of white balls in the urn
  all_genes - n_10, # the number of black balls in the urn
  n_07, # the number of balls drawn from the urn.
  lower.tail = FALSE
)




