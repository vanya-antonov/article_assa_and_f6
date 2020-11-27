
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


###
# Comparing with the RISE database  ----
rm(list = ls())
source('lib.R')

aso_07_df <- read_DE_summary_file(paste0(DATA_DIR, 'ASO_G0272888_AD_07.DE_Summary'))
aso_10_df <- read_DE_summary_file(paste0(DATA_DIR, 'ASO_G0272888_AD_10.DE_Summary'))
from_mouse_df <- read.delim(paste0(DATA_DIR, 'chaserr_orthologs.txt'), as.is=TRUE, header=TRUE) %>%
  # AC013394.2  =>  AC013394
  mutate(geneSymbol = gsub('\\.\\d+$', '', human_gene_name))
head(aso_07_df)
head(from_mouse_df)

rise_df <- read.csv(paste0(DATA_DIR, 'chaserr_RISE_db.csv'), as.is = TRUE)
head(rise_df)

aso_07_v <- unique(aso_07_df$geneSymbol)
aso_10_v <- unique(aso_10_df$geneSymbol)
mouse_v <- unique(from_mouse_df$geneSymbol)
rise_v <- unique(rise_df$Gene.name)

common3_v <- Reduce(intersect, list(aso_07_v, aso_10_v, mouse_v))

# "TNPO2", p-value = 0.1293645
rise_and_10_v <- intersect(rise_v, aso_10_v)
get_hgd_pvalue(all_balls=20000,
               red_balls=length(aso_10_v),
               drawn_balls=length(rise_v),
               drawn_red=length(rise_and_10_v))

# "NUP153", p-value = 0.2442752
rise_and_07_v <- intersect(rise_v, aso_07_v)
aso_07_df %>% filter(geneSymbol == 'NUP153')
get_hgd_pvalue(all_balls=20000,
               red_balls=length(aso_07_v),
               drawn_balls=length(rise_v),
               drawn_red=length(rise_and_07_v))


aso_07_and_10_v <- unique(c(aso_07_v, aso_10_v))
rise_and_07_10_v <- intersect(rise_v, aso_07_and_10_v)
get_hgd_pvalue(all_balls=NUM_HUMAN_GENES,
               red_balls=length(aso_07_and_10_v),
               drawn_balls=length(rise_v),
               drawn_red=length(rise_and_07_10_v))

