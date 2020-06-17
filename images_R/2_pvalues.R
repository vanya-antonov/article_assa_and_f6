
source('lib.R')

###

draw_2_pvalues <- function(input_fn, out_pdf)
{
  df <- read.delim(paste0(DATA_DIR, input_fn), as.is=TRUE, header=FALSE)
  colnames(df) <- c('aso_id', 'all_trxs', 'de_trxs', 'assa_hits_all', 'assa_hits_de', 'hgd_pvalue', 'gsea_pvalue')
  head(df)
  
  df <- df %>%
    filter(!is.na(all_trxs)) %>%
    mutate(hgd_pvalue = ifelse(is.na(hgd_pvalue), 1, hgd_pvalue),
           gsea_pvalue = ifelse(is.na(gsea_pvalue), 1, gsea_pvalue)) %>%
    mutate(log10_hgd = ifelse(hgd_pvalue < 1e-6, 6, -log10(hgd_pvalue)),
           log10_gsea = ifelse(gsea_pvalue < 1e-6, 6, -log10(gsea_pvalue))) %>%
    mutate(type = ifelse(hgd_pvalue < 0.01 | gsea_pvalue < 0.01, 'good', 'bad'))
  
  subT <- 'Good ASO: GSEA_pvalue < 0.01 or HGD_pvalue < 0.01'
  title <- sprintf('Total number of ASOs = %s (%s good)',
                   nrow(df), filter(df, type == 'good') %>% nrow())
  ggplot(df) +
    aes(x = log10_hgd, y = log10_gsea, col = type) +
    geom_point() +
    xlim(0, 6) +
    ylim(0,6) +
    ggtitle(title, subtitle = subT) +
    theme_bw()
  ggsave(out_pdf, path = OUT_DIR, width = 5, height = 5)
}

draw_2_pvalues(paste0(DATA_DIR, '2_pvalues.co.txt'),   '2_pvalues_co.pdf')
draw_2_pvalues(paste0(DATA_DIR, '2_pvalues.post.txt'), '2_pvalues_post.pdf')

