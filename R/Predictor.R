Predictor <- R6::R6Class('Predictor',
  public = list(
    hla_allele = NA,
    percentile_rank = NA,
    conn = NA,
    type = 'binding_affinity',
    table_name = NA,
    table_conn = NA,
    initialize = function(hla_allele, percentile_rank) {
      if (missing(hla_allele)) {
        self$hla_allele <- 'A0201'
      } else {
        self$hla_allele <- shortenHLA(hla_allele)
      }
      if (missing(percentile_rank)) {
        self$percentile_rank <- 1.9
      }
      self$conn <- DBI::dbConnect(RPostgreSQL::PostgreSQL(), 
        dbname = 'binding_affinity',
        user = 'm.slagter'
      )      
      self$table_name <- self$set_table_name(
        self$hla_allele, self$percentile_rank)
      if (self$table_name %in% RPostgreSQL::dbListTables(self$conn)) { 
        maartenutils::mymessage(sprintf('Connecting %s', self$table_name))
        self$table_conn <- tbl(self$conn, self$table_name)
      } else { 
        maartenutils::mymessage(sprintf('Table not found %s', self$table_name))
      }
    },
    finalize = function() {
      maartenutils::mymessage(sprintf('Disconnecting %s', self$table_name))
      on.exit(DBI::dbDisconnect(self$conn))
    },
    set_table_name = function(hla_allele, percentile_rank) {
      if (self$type == 'binding_affinity') {
        self$table_name <- hla_allele
      } else if (self$type == 'STS') {
        self$table_name <- sprintf('STS_%s_%.1f', hla_allele, percentile_rank)
      }
    },
    lookup = function(query_peps) {
      dplyr::filter(self$table_conn, peptide %in% query_peps) %>%
        dplyr::collect() %>%
        as.data.table
    },
    sanitize_query = function(peptides) {
      return(unique(setdiff(peptides, c('', NA))))
    },
    sanitize_res = function(dtf) {
      setDT(dtf)
      dtf <- dtf[, self$expected_res_columns, with = F]
      dtf <- self$set_col_types(dtf)
      return(unique(dtf))
    },
    set_col_types = function(dtf) {
      dtf <- dtf[!is.na(peptide)]
      if ('peptide' %in% colnames(dtf))
        dtf[, peptide := as.character(peptide)]

      if ('hla_allele' %in% colnames(dtf))
        dtf[, hla_allele := as.character(hla_allele)]

      if ('percentile_rank' %in% colnames(dtf))
        dtf[, percentile_rank := as.numeric(percentile_rank)]

      if ('affinity' %in% colnames(dtf))
        dtf[, affinity := as.numeric(affinity)]

      if ('different_from_self' %in% colnames(dtf))
        dtf[, different_from_self := as.logical(different_from_self)]
      dtf <- dtf[!is.na(peptide)]
      return(dtf)
    },
    store_res = function(dtf) {
      dtf <- self$sanitize_res(dtf)
      if (RPostgreSQL::dbExistsTable(self$conn, self$table_name)) {
        RPostgreSQL::dbWriteTable(self$conn, self$table_name, dtf, 
          row.names = F, append = T)
      } else {
        RPostgreSQL::dbWriteTable(self$conn, self$table_name, dtf, 
          row.names = F)
      }
    },
    compute = function(peptides, batch_size = 1e5, ncores = 1) {
      if (is.null(peptides) || length(peptides) == 0 || all(is.na(peptides))) 
        return(NULL)
      peptides <- self$sanitize_query(peptides)
      N_batches <- ceiling(length(peptides) / batch_size)
      query_peps_s <- split(peptides, 
        ceiling(seq(1, length(peptides)) / batch_size))
      if (ncores > 1) doParallel::registerDoParallel(ncores = ncores)
      query_res <- plyr::llply(seq_along(query_peps_s), function(idx) {
        qpl <- query_peps_s[[idx]]
        maartenutils::mymessage(
          msg = sprintf('Computing batch %d/%d (%d peptides)', idx, N_batches,
            length(qpl)))
        lres <- self$computer(qpl)
        self$store_res(lres)
        return(lres)
      }, .parallel = ncores > 1) %>% rbindlist(fill = T)
      return(self$sanitize_res(query_res))
    },
    query = function(peptides, batch_size = 1e5, ncores = 1) {
      peptides <- unique(peptides)
      already_computed <- self$lookup(peptides)
      already_computed <- self$sanitize_res(already_computed)
      to_be_computed <- setdiff(peptides, already_computed$peptide)
      extra_query_res <- rbindlist(
        list(already_computed, self$compute(to_be_computed)), fill = T)
      return(extra_query_res)
    },
    remove_peptides = function(peptides) {
      peptides <- unique(peptides)
      RPostgreSQL::dbSendQuery(self$conn, 
        sprintf('DELETE FROM "%s" WHERE "peptide" IN (%s)', self$table_name,
          paste(sprintf("'%s'", peptides), collapse = ', ')))
      invisible()
    }
  )
)