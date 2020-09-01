
source('lib.R')

###


df <- read.delim(paste0(DATA_DIR, '2_pvalues.co.txt'), as.is=TRUE, header=FALSE)
colnames(df) <- c('aso_id', 'all_trxs', 'de_trxs', 'assa_hits_all', 'assa_hits_de', 'hgd_pvalue', 'gsea_pvalue')
head(df)



good_df <- df %>%
  arrange(aso_id) %>%
  mutate(aso_id = factor(aso_id, levels = aso_id)) %>%
  mutate(type = ifelse(hgd_pvalue < 0.01 | gsea_pvalue < 0.01, 'good', 'bad')) %>%
  mutate(de_and_assa_hit = assa_hits_de,
         de_and_not_assa_hit = de_trxs - de_and_assa_hit) %>%
  mutate(non_de_trxs = all_trxs - de_trxs,
         non_de_and_assa_hit = assa_hits_all - assa_hits_de,
         non_de_and_not_assa_hit = non_de_trxs - non_de_and_assa_hit) %>%
  filter(type == 'good')



title <- sprintf('%s good ASOs', nrow(good_df))
good_df %>%
  select(aso_id, de_trxs, de_and_assa_hit, non_de_and_assa_hit) %>%
  gather('de_trxs', 'de_and_assa_hit', 'non_de_and_assa_hit',
         key = 'group', value = 'num_promoters') %>%
  mutate(group = factor(group, levels = c('de_trxs', 'de_and_assa_hit', 'non_de_and_assa_hit'))) %>%
  ggplot() +
  aes(x = aso_id, y = num_promoters, fill = group) +
  geom_bar(stat = 'identity', position = 'dodge') +
  coord_flip() +
  ggtitle(title, subtitle = 'ASSA pvalue < 0.01') +
  theme_bw() +
  theme(legend.position="bottom")
ggsave('assa_hits_among_de_and_non_de.CO.pdf', path = OUT_DIR,
       width = 6, height = 12)


###

rm(df)
rm(good_df)

df <- read.delim(paste0(DATA_DIR, '2_pvalues.post.txt'), as.is=TRUE, header=FALSE)
colnames(df) <- c('aso_id', 'all_trxs', 'de_trxs', 'assa_hits_all', 'assa_hits_de', 'hgd_pvalue', 'gsea_pvalue')
head(df)

good_df <- df %>%
  arrange(aso_id) %>%
  mutate(aso_id = factor(aso_id, levels = aso_id)) %>%
  mutate(type = ifelse(hgd_pvalue < 0.01 | gsea_pvalue < 0.01, 'good', 'bad')) %>%
  mutate(de_and_assa_hit = assa_hits_de,
         de_and_not_assa_hit = de_trxs - de_and_assa_hit) %>%
  mutate(non_de_trxs = all_trxs - de_trxs,
         non_de_and_assa_hit = assa_hits_all - assa_hits_de,
         non_de_and_not_assa_hit = non_de_trxs - non_de_and_assa_hit) %>%
  filter(type == 'good')



title <- sprintf('%s good ASOs', nrow(good_df))
good_df %>%
  select(aso_id, de_trxs, de_and_assa_hit, non_de_and_assa_hit) %>%
  gather('de_trxs', 'de_and_assa_hit', 'non_de_and_assa_hit',
         key = 'group', value = 'num_promoters') %>%
  mutate(group = factor(group, levels = c('de_trxs', 'de_and_assa_hit', 'non_de_and_assa_hit'))) %>%
  ggplot() +
  aes(x = aso_id, y = num_promoters, fill = group) +
  geom_bar(stat = 'identity', position = 'dodge') +
  coord_flip() +
  ggtitle(title, subtitle = 'ASSA pvalue < 0.01') +
  theme_bw() +
  theme(legend.position="bottom")
ggsave('assa_hits_among_de_and_non_de.POST.pdf', path = OUT_DIR, width = 6, height = 8)

