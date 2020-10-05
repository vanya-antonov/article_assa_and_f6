
source('lib.R')

###


aso_07_df <- read.delim(paste0(DATA_DIR, 'ASO_G0272888_AD_07.DE_Summary'), as.is=TRUE, header=TRUE) %>%
  # p002@AC013394.2  =>  AC013394.2
  mutate(geneSymbol = gsub('^p\\d+\\@', '', geneSymbol)) %>%
  # SPNS1;p002@RP11-264B17   =>   SPNS1
  mutate(geneSymbol = gsub(';.*$', '', geneSymbol)) %>%
  # AC013394.2  =>  AC013394
  mutate(geneSymbol = gsub('\\.\\d+$', '', geneSymbol)) %>%
  mutate(group = 'ASO_G0272888_AD_07') %>%
  select(geneSymbol, log2FC, group)
head(aso_07_df)

aso_10_df <- read.delim(paste0(DATA_DIR, 'ASO_G0272888_AD_10.DE_Summary'), as.is=TRUE, header=TRUE) %>%
  mutate(geneSymbol = gsub('p\\d+\\@', '', geneSymbol)) %>%
  # SPNS1;p002@RP11-264B17   =>   SPNS1
  mutate(geneSymbol = gsub(';.*$', '', geneSymbol)) %>%
  # AC013394.2  =>  AC013394
  mutate(geneSymbol = gsub('\\.\\d+$', '', geneSymbol)) %>%
  mutate(group = 'ASO_G0272888_AD_10') %>%
  select(geneSymbol, log2FC, group)
head(aso_10_df)

from_mouse_df <- read.delim(paste0(DATA_DIR, 'chaserr_orthologs.txt'), as.is=TRUE, header=TRUE) %>%
  mutate(group = 'Chaserr-/- mEFs') %>%
  # AC013394.2  =>  AC013394
  mutate(geneSymbol = gsub('\\.\\d+$', '', human_gene_name)) %>%
  rename(log2FC = log2fc) %>%
  select(geneSymbol, log2FC, group)
head(from_mouse_df)


df <- rbind(aso_07_df, aso_10_df, from_mouse_df) %>%
  filter(!is.na(log2FC)) %>%
  mutate(type = ifelse(log2FC > 0, 'up', 'down')) %>%
  group_by(group, type) %>%
  summarise(n=n())

UP_DOWN_COLS = c('up' = 'darkgreen', 'down' = 'red')

df %>%
  mutate(type = factor(type, levels=names(UP_DOWN_COLS))) %>%
  ggplot() +
  aes(x = group, y=n, fill = type) +
  geom_bar(stat = 'identity', col='black') +
#  theme_bw() +
  xlab('') +
  ylab('Number of differentially expressed genes') +
  coord_flip() +
  scale_fill_manual(
    name = 'Expression change',
    values = UP_DOWN_COLS,
    breaks = names(UP_DOWN_COLS)) +
  theme(legend.position = "top")
ggsave('chaserr_up_down.pdf', path = OUT_DIR, width = 9, height = 4)

