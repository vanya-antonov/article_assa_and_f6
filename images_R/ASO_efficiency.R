
source('lib.R')

###

df <- read.delim(paste0(DATA_DIR, 'ASO_names.DE_Summary'), as.is=TRUE, header=TRUE) %>%
 mutate(type = ifelse(log2FC < 0, 'log2FC < 0', 'log2FC > 0'))
head(df)

title <- sprintf('Total number of ASO = %s', nrow(df))
subT <- sprintf('Number of ASOs with log2FC < 0 = %s', filter(df, log2FC < 0) %>% nrow())
ggplot(df) +
  aes(x = log10(baseMean), y = log2FC, col = type) +
  geom_point() +
  ylim(-5, 5) +
  xlab('Target lncRNA expression, log10(baseMean)') +
  ylab('ASO knockdown efficiency, log2FC') +
  ggtitle(title, subtitle = subT) +
  theme_bw(base_size = 19)
ggsave('ASO_efficiency.pdf', path = OUT_DIR)

