
source('lib.R')

library(ComplexHeatmap)

###

common3_v <- read.delim(paste0(DATA_DIR, "chaserr_common3.txt"), as.is=TRUE, header=FALSE)$V1

aso_07_g_df <- read.delim(paste0(DATA_DIR, 'ASO_G0272888_AD_07.DE_Summary_gene'), as.is=TRUE, header=TRUE) %>%
  # AC013394.2  =>  AC013394
  mutate(geneSymbol = gsub('\\.\\d+$', '', geneSymbol)) %>%
  # https://stackoverflow.com/a/60861471/310453
  filter(!duplicated(geneSymbol)) %>%
  filter(geneSymbol %in% common3_v) %>%
  mutate(ASO_G0272888_AD_07 = ifelse(log2FC < 0, 'down', 'up')) %>%
  rename(log2FC_ASO_G0272888_AD_07 = log2FC) %>%
  select(geneSymbol, ASO_G0272888_AD_07, log2FC_ASO_G0272888_AD_07)
head(aso_07_g_df)

aso_10_g_df <- read.delim(paste0(DATA_DIR, 'ASO_G0272888_AD_10.DE_Summary_gene'), as.is=TRUE, header=TRUE) %>%
  # AC013394.2  =>  AC013394
  mutate(geneSymbol = gsub('\\.\\d+$', '', geneSymbol)) %>%
  # https://stackoverflow.com/a/60861471/310453
  filter(!duplicated(geneSymbol)) %>%
  filter(geneSymbol %in% common3_v) %>%
  mutate(ASO_G0272888_AD_10 = ifelse(log2FC < 0, 'down', 'up')) %>%
  rename(log2FC_ASO_G0272888_AD_10 = log2FC) %>%
  select(geneSymbol, ASO_G0272888_AD_10, log2FC_ASO_G0272888_AD_10)
head(aso_10_g_df)

from_mouse_df <- read.delim(paste0(DATA_DIR, 'chaserr_orthologs.txt'), as.is=TRUE, header=TRUE) %>%
  # AC013394.2  =>  AC013394
  mutate(geneSymbol = gsub('\\.\\d+$', '', human_gene_name)) %>%
  filter(!is.na(log2fc)) %>%
  filter(!duplicated(geneSymbol)) %>%
  rename(log2FC = log2fc) %>%
  mutate(mEF_knockout = ifelse(log2FC < 0, 'down', 'up')) %>%
  rename(log2FC_mEF_knockout = log2FC) %>%
  select(geneSymbol, mEF_knockout, log2FC_mEF_knockout)
head(from_mouse_df)



ht_df <- aso_07_g_df %>%
  left_join(aso_10_g_df, by='geneSymbol') %>%
  left_join(from_mouse_df, by='geneSymbol') %>%
  mutate(s_val = paste(ASO_G0272888_AD_07, ASO_G0272888_AD_10, mEF_knockout)) %>%
  arrange(s_val)

# ?Heatmap

ht <- Heatmap(ht_df[, c('ASO_G0272888_AD_07', 'ASO_G0272888_AD_10', 'mEF_knockout')],
        name = 'Expression\nChange',
        col = c('up' = 'darkgreen', 'down' = 'red'),
        show_row_names = FALSE,
        row_title_rot = 0)

pdf(paste0(OUT_DIR, 'other/chaserr_common_genes_ht.pdf'), width = 4, height = 8)
draw(ht, row_title = sprintf('%d common target genes', nrow(ht_df)))
dev.off()


# Number of 'down down down' genes = 126
n_ddd <- ht_df %>% filter(s_val == 'down down down') %>% nrow()
n_ddd / nrow(ht_df)
