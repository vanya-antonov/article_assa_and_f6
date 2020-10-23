
source('lib.R')

###
# "In total, we identified 104 cases (71 unique lncRNAs) where genes"
co_df <- read.delim(paste0(DATA_DIR, '2_pvalues.co.txt'), as.is=TRUE, header=FALSE)
colnames(co_df) <- c('aso_id', 'all_trxs', 'de_trxs', 'assa_hits_all', 'assa_hits_de', 'hgd_pvalue', 'gsea_pvalue')

post_df <- read.delim(paste0(DATA_DIR, '2_pvalues.post.txt'), as.is=TRUE, header=FALSE)
colnames(post_df) <- c('aso_id', 'all_trxs', 'de_trxs', 'assa_hits_all', 'assa_hits_de', 'hgd_pvalue', 'gsea_pvalue')

aso_names <- read.delim(paste0(DATA_DIR, 'ASO_names.info'), as.is=TRUE, header=FALSE)
colnames(aso_names) <- c('aso_id', 'aso_id_short', 'prmtrID', 'prmtrName', 'trnscptID', 'trnscptName', 'trnscptType', 'geneID', 'geneName', 'geneType')

good_co_v <- get_good_aso(co_df)
good_post_v <- get_good_aso(post_df)

unique_good_v <- c(good_co_v, good_post_v) %>% unique()

length(unique_good_v)
aso_names %>%
  filter(aso_id %in% unique_good_v) %>%
  pull(geneID) %>%
  unique() %>%
  length()


###
# 
rm(list = ls())
source('lib.R')

df <- read.delim(paste0(DATA_DIR, 'ASO_names.DE_Summary'), as.is=TRUE, header=TRUE)
head(df)

aso_names <- read.delim(paste0(DATA_DIR, 'ASO_names.info'), as.is=TRUE, header=FALSE)
colnames(aso_names) <- c('aso_id', 'aso_id_short', 'prmtrID', 'prmtrName', 'trnscptID', 'trnscptName', 'trnscptType', 'geneID', 'geneName', 'geneType')
head(aso_names)

aso_names %>%
  filter(aso_id %in% df$perturb_id) %>%
  group_by(geneType) %>%
  summarise(
    num_ASO = n(),
    num_genes = n_distinct(geneName))

aso_names %>%
  filter(aso_id %in% df$perturb_id) %>%
  pull(geneID) %>%
  unique() %>%
  length()
