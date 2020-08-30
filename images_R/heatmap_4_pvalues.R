
source('lib.R')

library(ComplexHeatmap)
library(circlize)   # colorRamp2()

###

LOG_COLORS <- colorRamp2(c(0, 2, 10), c("white", "gray", "black"))

###

genya_cols <- read.delim(paste0(DATA_DIR, 'genya.table_for_all_marks.tsv'), as.is=TRUE, header=FALSE, nrows = 2)
cols_v <- paste(genya_cols[1,], genya_cols[2,], sep='_')

genya_df <- read.delim(paste0(DATA_DIR, 'genya.table_for_all_marks.tsv'), as.is=TRUE, header=FALSE, skip = 2, col.names = cols_v) %>%
  separate(target_type, c('gene', 'aso_id'), sep='_ASO_') %>%
  mutate(aso_id = paste0('ASO_', aso_id)) %>%
  gather(key='type_pvalue', value='log_pvalue', -c(gene, aso_id)) %>%
  filter(log_pvalue >= 1.3) %>%
  # H3K27ac_wa  =>  H3K27ac
  mutate(mark = sub('_..$', '', type_pvalue) )
head(genya_df)

genya_aso <- genya_df %>%
  group_by(aso_id) %>%
  summarise(
    n_marks = n_distinct(mark),
    mark_str = paste(unique(mark), collapse=','))

co_df <- read.delim(paste0(DATA_DIR, '2_pvalues.co.txt'), as.is=TRUE, header=FALSE)
colnames(co_df) <- c('aso_id', 'all_trxs', 'de_trxs', 'assa_hits_all', 'assa_hits_de', 'hgd_pvalue', 'gsea_pvalue')
head(co_df)

post_df <- read.delim(paste0(DATA_DIR, '2_pvalues.post.txt'), as.is=TRUE, header=FALSE)
colnames(post_df) <- c('aso_id', 'all_trxs', 'de_trxs', 'assa_hits_all', 'assa_hits_de', 'hgd_pvalue', 'gsea_pvalue')
head(post_df)

aso_names <- read.delim(paste0(DATA_DIR, 'ASO_names.info'), as.is=TRUE, header=FALSE)
colnames(aso_names) <- c('aso_id', 'aso_id_short', 'prmtrID', 'prmtrName', 'trnscptID', 'trnscptName', 'trnscptType', 'geneID', 'geneName', 'geneType')
head(aso_names)


df1 <- co_df %>%
  filter(!is.na(all_trxs)) %>%
  mutate(co_hgd = ifelse(is.na(hgd_pvalue), 1, hgd_pvalue),
         co_gsea = ifelse(is.na(gsea_pvalue), 1, gsea_pvalue)) %>%
  select(aso_id, co_hgd, co_gsea)

df2 <- post_df %>%
  mutate(post_hgd = ifelse(is.na(hgd_pvalue), 1, hgd_pvalue),
         post_gsea = ifelse(is.na(gsea_pvalue), 1, gsea_pvalue)) %>%
  select(aso_id, post_hgd, post_gsea)


df <- full_join(df1, df2, by = 'aso_id') %>%
  replace_na(list(co_hgd = 1, co_gsea = 1, post_hgd = 1, post_gsea = 1)) %>%
  inner_join(aso_names, by = 'aso_id') %>%
  select(geneName, aso_id, co_hgd, co_gsea, post_hgd, post_gsea) %>%
  arrange(geneName, co_hgd)

good_lncrna <- df %>%
  filter(co_hgd < 0.01 | co_gsea < 0.01 | post_hgd < 0.01 | post_gsea < 0.01) %>%
  pull(geneName) %>%
  unique()


draw_ht <- function(good_lncrna, vanya_g_v=c(), show_row_names=FALSE)
{
  ht_df <- filter(df, geneName %in% good_lncrna) %>%
    mutate(co_hgd = ifelse(co_hgd < 1e-10, 10, -log10(co_hgd)),
           co_gsea = ifelse(co_gsea < 1e-10, 10, -log10(co_gsea)),
           post_hgd = ifelse(post_hgd < 1e-10, 10, -log10(post_hgd)),
           post_gsea = ifelse(post_gsea < 1e-10, 10, -log10(post_gsea))) %>%
    left_join(genya_aso, by='aso_id') %>%
    mutate(Epigenetics = factor(ifelse(is.na(n_marks), 'No', 'Yes'))) %>%
    mutate(Vanya = factor(ifelse(geneName %in% vanya_g_v, 'Yes', 'No')))
  rownames(ht_df)  <-  ht_df$aso_id

  log_ht <- Heatmap(ht_df[, c('co_hgd', 'co_gsea', 'post_hgd', 'post_gsea')],
                    name = '-log10(p-value)',
                    col = LOG_COLORS,
                    split = ht_df$geneName,
                    cluster_columns = FALSE,
                    cluster_rows = FALSE,
                    show_row_names = show_row_names,
                    #row_names_gp = gpar(fontsize = 1),
                    #gap = unit(1, "mm"),
                    row_title_rot = 0)
  
  genya_ht <- Heatmap(ht_df['Epigenetics'],
                      show_row_names = show_row_names,
                      name='Epigenetics',
                      col = c('No' = 'white', 'Yes' = 'red'))
  
  # vanya_ht <- Heatmap(ht_df['Vanya'],
  #                     show_row_names = show_row_names,
  #                     name='Vanya',
  #                     col = c('No' = 'white', 'Yes' = 'darkgreen'))
  # log_ht + marks_ht
  
  draw(log_ht + genya_ht)
}
#draw_ht(good_lncrna[1:half_good_n], vanya_genes)

vanya_genes <- c(
  'AC004980.7',
  'AC005592.2',
  'AC013394.2',
  'AC108488.3',
  'BOLA3−AS1',
  'CATG00000017883.1',
  'CD99P1',
  'DANCR',
  'FGD5-AS1',
  'MFI2-AS1',
  'RP11-660L16.2',
  'RP11-296O14.3',
  'RP11-33B1.1',
  'SERTAD4-AS1'
)

pdf(paste0(OUT_DIR, 'heatmap_4_pvalues_selected.pdf'), width = 8, height = 10)
draw_ht(vanya_genes, vanya_genes, show_row_names = TRUE)
dev.off()


half_good_n <- round(length(good_lncrna)/2, 0)

pdf(paste0(OUT_DIR, 'heatmap_4_pvalues_1.pdf'), width = 5, height = 10)
draw_ht(good_lncrna[1:half_good_n], vanya_genes)
dev.off()

pdf(paste0(OUT_DIR, 'heatmap_4_pvalues_2.pdf'), width = 5, height = 10)
draw_ht(good_lncrna[(half_good_n+1):length(good_lncrna)],  vanya_genes)
dev.off()




