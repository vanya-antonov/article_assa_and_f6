
source('lib.R')

###

draw_2_pvalues <- function(input_fn, title = '')
{
  df <- read.delim(paste0(DATA_DIR, input_fn), as.is=TRUE, header=FALSE)
  colnames(df) <- c('aso_id', 'all_trxs', 'de_trxs', 'assa_hits_all', 'assa_hits_de', 'hgd_pvalue', 'gsea_pvalue')
  head(df)
  
  good_name <- 'HGD p-value < 0.01  or  GSEA p-value < 0.01'
  bad_name  <- 'HGD p-value >= 0.01  and  GSEA p-value >= 0.01'
  df <- df %>%
    filter(!is.na(all_trxs)) %>%
    mutate(hgd_pvalue = ifelse(is.na(hgd_pvalue), 1, hgd_pvalue),
           gsea_pvalue = ifelse(is.na(gsea_pvalue), 1, gsea_pvalue)) %>%
    mutate(log10_hgd = ifelse(hgd_pvalue < 1e-6, 6, -log10(hgd_pvalue)),
           log10_gsea = ifelse(gsea_pvalue < 1e-6, 6, -log10(gsea_pvalue))) %>%
    mutate(type = ifelse(hgd_pvalue < 0.01 | gsea_pvalue < 0.01, good_name, bad_name))
  
  subT <- sprintf('Total number of ASOs = %s (%s good)',
                   nrow(df), filter(df, type == good_name) %>% nrow())
  ggplot(df) +
    aes(x = log10_hgd, y = log10_gsea, col = type) +
    geom_point() +
    xlim(0, 6) +
    ylim(0,6) +
    ylab('log10(GSEA p-value)') +
    xlab('log10(HGD p-value)') +
    ggtitle(title, subtitle = subT) +
    theme_bw() +
    theme(legend.position="bottom", legend.direction="vertical", legend.title=element_blank())
}

draw_2_pvalues(paste0(DATA_DIR, '2_pvalues.co.txt'), 'Co-transcriptional interactions')
ggsave('two_pvalues_co.pdf', path = OUT_DIR, width = 5, height = 5)

draw_2_pvalues(paste0(DATA_DIR, '2_pvalues.post.txt'), 'Post-transcriptional interactions')
ggsave('two_pvalues_post.pdf', path = OUT_DIR, width = 5, height = 5)

