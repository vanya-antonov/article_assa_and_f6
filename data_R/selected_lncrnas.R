
source('lib.R')

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
  select(geneName) %>%
  # Remove duplicates
  distinct(geneName, .keep_all = TRUE)

write.table(good_lncrna, file = paste0(DATA_DIR, "selected_lncrnas.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")
