#' Precompute peptides
#'
#' As the loading of self lists is somewhat slow and heavy on RAM,
#' we'd rather not share this across multiple workers (each loading it
#' for themselves...) during neo-antigen tallying. As an 'easy' fix,
#' this function allows one to precompute all candidate peptides such
#' that they can later be queried in parallel.
#'
#' @param peptides List of peptides to check. If NULL, all peptides
#' from the database are queried for their binding affinity and
#' checked for STS coverage in the associated STS database
#'
precompute_peps <- function(peptides = all_peps, 
                            hla_allele = 'B2705',
                            BP_class = 'BindingPredictor',
                            perform_STS = TRUE, 
                            percentile_rank = 1.9,
                            STS_percentile_rank = 1.9,
                            STS_threshold = 4,
                            ncores = 10, 
                            batch_size = 1e6,
                            STS_batch_size = 1e3,
                            STS_on_affinity_peps_only = TRUE,
                            verbose = FALSE,
                            index = FALSE) {

  if (perform_STS) {
    job_description <- 'BA and STS'
  } else {
    job_description <- 'BA'
  }
  if (STS_on_affinity_peps_only) {
    job_description <- glue::glue('{job_description}, STS ',
      'on peptides with percentile_rank <= {STS_threshold} only')
  }
  message(sprintf('Cmputing %s - %s - %.1f',
      job_description, hla_allele, percentile_rank))
  BA_predictor <- get(BP_class)$new(hla_allele = hla_allele,
                                    verbose = verbose)

  if (index) {
    BA_predictor$index_SQL()
  }

  # peptides <- c("AKARAAFSL", "AKARGHFQK", "AKARIIMLI", "AKARLVRYM",
  # "AKARQLIGL", "AKARVVQLR", "AKASCILPV", "AKASCLPSK", "AKASFSLRL",
  # "AKASIAFIF", "AKASMISKL")
  l_pr <- percentile_rank
  if (is.null(peptides) || length(peptides) == 0) {
    peptides <-
      dplyr::filter(
        BA_predictor$table_conn,
        percentile_rank <= l_pr
      ) %>%
      select(peptide) %>%
      collect %>%
      unlist %>%
      setNames(NULL)
  } else {
    res <- BA_predictor$query(
      peptides, 
      ncores = ncores,
      batch_size = batch_size,
      prefill_mode = FALSE
    )
    if (!null_dat(res)) {
      if (STS_on_affinity_peps_only) {
        STS_query_peps <- 
          res[percentile_rank <= STS_threshold, peptide]
      } else {
        STS_query_peps <- res[, peptide]
      }
    } else {
      STS_query_peps <- NULL
    }
    # STS_query_peps <- res[percentile_rank <= 4, peptide]
    # problem_peps <- c('TTFAHSSTV', 'IQNLQGLFV', 'RLVQARALI', 'LVLNCCRAV')
    # problem_peps %in% peptides
    # problem_peps %in% res$peptide
    # problem_peps %in% STS_query_peps
  }
  rm(BA_predictor)

  if (perform_STS) {
    STS_predictor <- STSPredictor$new(
      hla_allele = hla_allele,
      STS_percentile_rank = STS_percentile_rank,
      verbose = verbose
    )
    # res <- RPostgreSQL::dbGetQuery(STS_predictor$conn,
    #   sprintf('SELECT EXISTS(SELECT * FROM "%s"  WHERE "peptide" IN (%s))',
    #     STS_predictor$table_name, paste(sprintf("'%s'", STS_query_peps), collapse = ', ')))
    # res <- RPostgreSQL::dbGetQuery(STS_predictor$conn,
    #   sprintf('SELECT * FROM "%s"  WHERE "peptide" IN (%s)',
    #     STS_predictor$table_name, paste(sprintf("'%s'", STS_query_peps), collapse = ', ')))
    # res <- STS_predictor$lookup(STS_query_peps)
    # missing_peps <- setdiff(STS_query_peps, res$peptide)
    res <- STS_predictor$query(
      STS_query_peps,
      ncores = ncores,
      batch_size = STS_batch_size,
      prefill_mode = TRUE
    )
    rm(STS_predictor)
  }
  rm(res)
  rm(STS_query_peps)
  gc()
  message(sprintf('Finished precomputing all %s - %s - %.1f',
      job_description, hla_allele, percentile_rank))
  invisible()
}
