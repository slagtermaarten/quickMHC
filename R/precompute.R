#' Precompute peptides
#'
#' As the loading of self lists is somewhat slow and heavy on RAM, we'd rather
#' not share this across multiple workers when RAM is not shared as is typically
#' the case for us. As an 'easy' fix, this function allows one to precompute all
#' relevant peptides (computation in parallel, lookups serially), such that they
#' can later be queried in parallel.
#'
#' @param peptides List of peptides to check. If NULL, all peptides from the
#' database are queried for their binding affinity and checked for STS coverage
#' in the associated STS database
#'
precompute_peps <- function(peptides = all_peps, hla_allele = 'B2705', 
  perform_STS = T, percentile_rank = 1.9, ncores = 10, batch_size = 1e6) {

  if (perform_STS) {
    job_description <- 'BA and STS'
  } else {
    job_description <- 'BA'
  }
  mymessage(sprintf('Precomputing all %s - %s - %.1f', 
      job_description, hla_allele, percentile_rank)) 
  BA_predictor <- BindingPredictor$new(
    hla_allele = hla_allele, verbose = T) 
  BA_predictor$index_SQL()

  # peptides <- c("AKARAAFSL", "AKARGHFQK", "AKARIIMLI", "AKARLVRYM", "AKARQLIGL",
  #   "AKARVVQLR", "AKASCILPV", "AKASCLPSK", "AKASFSLRL", "AKASIAFIF",
  #   "AKASMISKL")
  l_pr <- percentile_rank
  if (is.null(peptides) || length(peptides) == 0) {
    query_peps <- 
      dplyr::filter(BA_predictor$table_conn, percentile_rank <= l_pr) %>%
      select(peptide) %>%
      collect %>%
      unlist %>%
      setNames(NULL)
  } else {
    res <- BA_predictor$query(
      peptides, ncores = ncores, batch_size = batch_size)
    query_peps <- res[percentile_rank <= l_pr, peptide]
  }
  rm(BA_predictor)

  if (perform_STS) {
    STS_predictor <- STSPredictor$new(
      hla_allele = hla_allele,
      STS_percentile_rank = percentile_rank, verbose = T)
    # res <- RPostgreSQL::dbGetQuery(STS_predictor$conn,
    #   sprintf('SELECT EXISTS(SELECT * FROM "%s"  WHERE "peptide" IN (%s))',
    #     STS_predictor$table_name, paste(sprintf("'%s'", query_peps), collapse = ', ')))
    # res <- RPostgreSQL::dbGetQuery(STS_predictor$conn,
    #   sprintf('SELECT * FROM "%s"  WHERE "peptide" IN (%s)',
    #     STS_predictor$table_name, paste(sprintf("'%s'", query_peps), collapse = ', ')))
    # res <- STS_predictor$lookup(query_peps)
    # missing_peps <- setdiff(query_peps, res$peptide)
    res <- STS_predictor$query(query_peps, ncores = ncores, 
      batch_size = batch_size)
    rm(STS_predictor)
  }
  rm(res)
  rm(query_peps)
  gc()
  mymessage(sprintf('Finished precomputing all %s - %s - %.1f', 
      job_description, hla_allele, percentile_rank)) 
  invisible()
}
