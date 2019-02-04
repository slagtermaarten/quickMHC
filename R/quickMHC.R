format_query_table <- function(query_table) {
  setDT(query_table)
  query_table <- maartenutils::cond_setnames(query_table, 'hla', 'hla_allele')
  if ('hla_allele' %in% colnames(query_table)) {
    query_table[, hla_allele := shortenHLA(hla_allele)]
  }
  return(query_table)
}


#' Query a combination of peptides/MHCs
#'
#'
quickMHC <- function(
  query_table = data.frame('peptide' = c('SYFPEITHI'), 'hla_allele' = c('A0201')), 
  STS_percentile_rank = NULL, 
  verbose = F, batch_size = 1e5, ncores = 1) {

  query_table <- format_query_table(query_table)
  unique_hla_alleles <- setdiff(unique(query_table$hla_allele), c('', NA)) %>%
    shortenHLA

  plyr::ldply(unique_hla_alleles, function(hla_allele) {
    if (is.null(STS_percentile_rank)) {
      predictor <- BindingPredictor$new(
        hla_allele = hla_allele, 
        verbose = verbose
      )
    } else {
      predictor <- STSPredictor$new(
        hla_allele = hla_allele,
        STS_percentile_rank = STS_percentile_rank,
        verbose = verbose
      )
    }
    l_hla_allele <- hla_allele
    res <- cbind(
      'hla_allele' = l_hla_allele, 
      predictor$query(query_table[hla_allele == l_hla_allele, peptide],
        ncores = ncores, batch_size = batch_size)
    )
    rm(predictor)
    return(res)
  }, .parallel = F) %>% as.data.table
}
