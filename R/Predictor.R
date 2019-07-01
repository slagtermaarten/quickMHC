Predictor <- R6::R6Class('Predictor',
  public = list(
    hla_allele = NA,
    STS_percentile_rank = NA,
    table_name = NA,
    conn = NA,
    table_conn = NA,
    verbose = F,
    initialize = function(hla_allele = 'A0201', verbose = F) {
      self$hla_allele <- shortenHLA(hla_allele)
      self$conn <- DBI::dbConnect(RPostgreSQL::PostgreSQL(max.con = 100),
        dbname = config$db_name,
        user = config$db_user
      )
      self$table_name <- self$set_table_name(
        self$hla_allele, self$STS_percentile_rank)
      self$verbose <- verbose
      if (self$table_name %in% RPostgreSQL::dbListTables(self$conn)) {
        if (self$verbose) {
          maartenutils::mymessage(instance = self$table_name,
            msg = sprintf('Connecting %s', self$table_name))
        }
        self$table_conn <- tbl(self$conn, self$table_name)
      } else {
        maartenutils::mymessage(instance = self$table_name,
          msg = sprintf('Table not found %s', self$table_name))
      }
    },
    finalize = function() {
      if (self$verbose) {
        maartenutils::mymessage(instance = self$table_name,
            msg = sprintf('Disconnecting %s', self$table_name))
      }
      DBI::dbDisconnect(self$conn)
    },
    set_table_name = function(hla_allele, STS_percentile_rank) {
      if (self$type == 'binding_affinity') {
        self$table_name <- hla_allele
      } else if (self$type == 'STS') {
        self$table_name <-
          sprintf('STS_%s_%.1f', hla_allele, STS_percentile_rank)
      }
    },
    lookup = function(query_peps) {
      query_peps <- self$sanitize_query(query_peps)
      quoted_peps <- paste(sprintf("'%s'", query_peps), collapse = ', ')
      sql <- sprintf('SELECT * FROM "%s" WHERE peptide IN (%s);', 
        self$table_name, quoted_peps)
      res <- RPostgreSQL::dbGetQuery(self$conn, sql) %>% as.data.table
      return(self$sanitize_res(res))
    },
    sanitize_query = function(peptides) {
      if (is.null(peptides) || length(peptides == 0)) return(NULL)
      return(unique(setdiff(peptides, c('', NA))))
    },
    sanitize_res = function(dtf) {
      if (maartenutils::null_dat(dtf)) return(NULL)
      setDT(dtf)
      dtf <- tryCatch(dtf[, self$expected_res_columns, with = F],
        error = function(e) {
          message('Not all expected columns present in compute res')
          return(NULL)
        })
      dtf <- self$set_col_types(dtf)
      setkey(dtf, peptide)
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
        tryCatch({
          RPostgreSQL::dbWriteTable(self$conn, self$table_name, dtf,
            row.names = F, append = T)
        }, error = function(e) { print(e) })
      } else {
        tryCatch({
          RPostgreSQL::dbWriteTable(self$conn, self$table_name, dtf,
            row.names = F)
        }, error = function(e) { print(e) })
        self$index_SQL()
      }
    },
    unique_SQL = function() {
      column_names <- paste(glue::glue('"{self$expected_res_columns}"'),
        collapse = ', ')
      tryCatch({
        RPostgreSQL::dbSendQuery(self$conn, glue::glue('
          CREATE TABLE "{self$table_name}_temp" (LIKE "{self$table_name}");
          INSERT INTO "{self$table_name}_temp"({column_names})
          SELECT DISTINCT ON ("peptide") {column_names} FROM "{self$table_name}";
          DROP TABLE "{self$table_name}" CASCADE;
          ALTER TABLE "{self$table_name}_temp" RENAME TO "{self$table_name}";
        '))
      }, error = function(e) { print(e) })
    },
    index_SQL = function() {
      self$unique_SQL()
      tryCatch({
        RPostgreSQL::dbSendQuery(self$conn,
          glue::glue('create unique index if not exists "{self$table_name}_idx"
            ON "{self$table_name}" (peptide);'))
      }, error = function(e) { print(e) })
    },
    compute = function(peptides, batch_size = 1e5, ncores = 1, overwrite = F) {
      if (is.null(peptides) || length(peptides) == 0 || all(is.na(peptides)))
        return(NULL)
      peptides <- self$sanitize_query(peptides)
      if (overwrite) {
        already_computed <- self$lookup(peptides)
        self$remove_peptides(already_computed$peptide)
      }

      if (length(peptides) == 0) {
        return(NULL) 
      } else {
        N_batches <- ceiling(length(peptides) / batch_size)
        query_peps_s <- split(peptides,
          ceiling(seq(1, length(peptides)) / batch_size))
        if (ncores > 1) doParallel::registerDoParallel(cores = ncores)
        compute_res <- plyr::llply(seq_along(query_peps_s), function(idx) {
          qpl <- query_peps_s[[idx]]
          maartenutils::mymessage(instance = self$table_name,
            msg = sprintf('computing batch %d/%d (%d peptides)', idx, N_batches,
              length(qpl)))
          if (self$verbose) {
            maartenutils::mymessage(instance = self$table_name,
              msg = sprintf('computing %s', paste(qpl, collapse = ', ')))
          }
          lres <- self$computer(qpl)
          return(lres)
        }, .parallel = ncores > 1) %>% rbindlist(fill = T)
        self$store_res(compute_res)
        return(compute_res)
      }
    },
    query = function(peptides, batch_size = 1e5, ncores = 1) {
      peptides <- self$sanitize_query(peptides)
      if (length(peptides) == 0) return(NULL)
      N_batches <- ceiling(length(peptides) / batch_size)
      query_peps_s <- split(peptides,
        ceiling(seq(1, length(peptides)) / batch_size))

      ## Lookups are performed in batches (probably faster but never really
      ## tested it) and serially in order to not multiply
      ## use the database connection associated with the R6 object across
      ## multiple threads; computations are done in parallel.
      # if (ncores > 1) doParallel::registerDoParallel(cores = ncores)
      query_res <- plyr::llply(seq_along(query_peps_s), function(idx) {
        qpl <- query_peps_s[[idx]]
        if (self$verbose) {
          maartenutils::mymessage(instance = self$table_name,
            msg = sprintf('querying batch %d/%d (%d peptides)', idx, N_batches,
              length(qpl)))
        }
        already_computed <- self$lookup(qpl)
        to_be_computed <- setdiff(qpl, already_computed$peptide)
        l_res <- rbindlist(
          list(already_computed,
            self$compute(to_be_computed, batch_size = 100, ncores = ncores)
          ), fill = T)
        return(l_res)
      }, .parallel = F) %>% rbindlist(fill = T)
      return(query_res)
    },
    remove_peptides = function(peptides) {
      peptides <- self$sanitize_query(peptides)
      if (is.null(peptides)) return(NULL)
      RPostgreSQL::dbSendQuery(self$conn,
        sprintf('DELETE FROM "%s" WHERE "peptide" IN (%s)', self$table_name,
          paste(sprintf("'%s'", peptides), collapse = ', ')))
      invisible()
    }
  )
)
