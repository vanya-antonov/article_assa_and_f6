
source('lib.R')

library(ComplexHeatmap)
library(circlize)   # colorRamp2()

###

co_df <- read.delim(paste0(DATA_DIR, '2_pvalues.co.txt'), as.is=TRUE, header=FALSE)
colnames(co_df) <- c('aso_id', 'all_trxs', 'de_trxs', 'assa_hits_all', 'assa_hits_de', 'hgd_pvalue', 'gsea_pvalue')
head(co_df)

post_df <- read.delim(paste0(DATA_DIR, '2_pvalues.post.txt'), as.is=TRUE, header=FALSE)
colnames(post_df) <- c('aso_id', 'all_trxs', 'de_trxs', 'assa_hits_all', 'assa_hits_de', 'hgd_pvalue', 'gsea_pvalue')
head(post_df)

aso_names <- read.delim(paste0(DATA_DIR, 'ASO_names.info'), as.is=TRUE, header=FALSE)
colnames(aso_names) <- c('aso_id', 'aso_id_short', 'prmtrID', 'prmtrName', 'trnscptID', 'trnscptName', 'trnscptType', 'geneID', 'geneName', 'geneType')
head(aso_names)

gcorr_df <- read.delim(paste0(DATA_DIR, 'MARGI_GenometriCorrelation.txt'), as.is=TRUE, header=FALSE)
colnames(gcorr_df) <- c('aso_id', 'gcorr_pvalue')
head(gcorr_df)

good_co_v <- get_good_aso(co_df)
good_post_v <- get_good_aso(post_df)
good_two_v <- intersect(good_co_v, good_post_v)


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
head(df)

get_log_ht <- function(df, show_row_names=FALSE)
{
  ht_df <- df %>%
    mutate(co_hgd = ifelse(co_hgd < 1e-10, 10, -log10(co_hgd)),
           co_gsea = ifelse(co_gsea < 1e-10, 10, -log10(co_gsea)),
           post_hgd = ifelse(post_hgd < 1e-10, 10, -log10(post_hgd)),
           post_gsea = ifelse(post_gsea < 1e-10, 10, -log10(post_gsea)))
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
  return(log_ht)
}

get_margi_ht <- function(df, gcorr_df, show_row_names=FALSE)
{
  ht_df <- df %>%
    left_join(gcorr_df, by = 'aso_id') %>%
    mutate(gcorr_pvalue = ifelse(is.na(gcorr_pvalue), 1, gcorr_pvalue)) %>%
    mutate(MARGI = ifelse(gcorr_pvalue < 1e-10, 10, -log10(gcorr_pvalue)))
  rownames(ht_df)  <-  ht_df$aso_id
  
  log_ht <- Heatmap(ht_df['MARGI'],
                    name = '-log10(MARGI p-value)',
                    col = LOG_COLORS_RED,
                    # split = ht_df$geneName,
                    cluster_columns = FALSE,
                    cluster_rows = FALSE,
                    show_row_names = show_row_names,
                    #row_names_gp = gpar(fontsize = 1),
                    #gap = unit(1, "mm"),
                    row_title_rot = 0)
  
  return(log_ht)
}


good_df <- filter(df, aso_id %in% good_two_v)

pdf(paste0(OUT_DIR, 'two_pvalues_heatmap.pdf'), width = 8, height = 10)
get_log_ht(good_df, show_row_names=FALSE) +
  get_margi_ht(good_df, gcorr_df, show_row_names=TRUE)
dev.off()


# Number of unique genes = 31
good_df %>% pull(geneName) %>% unique() %>% length()


good_df %>%
  left_join(gcorr_df, by = 'aso_id') %>%
  mutate(gcorr_pvalue = ifelse(is.na(gcorr_pvalue), 1, gcorr_pvalue)) %>%
  filter(gcorr_pvalue < 0.01)
