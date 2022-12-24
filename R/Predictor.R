Predictor <- R6::R6Class('Predictor',
  public = list(
    hla_allele = NA,
    STS_percentile_rank = NA,
    affinity_predictor = 'NetMHCpan3',
    table_name = NA,
    conn = NA,
    table_conn = NA,
    pred_wrapper = NA,
    verbose = F,
    initialize = function(hla_allele = 'A0201', 
                          verbose = F, max_conn = 1e3) {
      self$hla_allele <- shortenHLA(hla_allele)
      self$conn <- 
        DBI::dbConnect(
          RPostgreSQL::PostgreSQL(max.con = max_conn),
          dbname = config$db_name,
          user = config$db_user
        )
      ## Table to be created in SQL database
      self$table_name <- self$set_table_name(
        hla_allele = self$hla_allele, 
        STS_percentile_rank = self$STS_percentile_rank
      )
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
      if (self$affinity_predictor == 'NetMHCpan3') {
        affinity_pred_flag <- ''
      } else {
        affinity_pred_flag <- paste0('-', self$affinity_predictor)
      }
      if (self$type == 'binding_affinity') {
        self$table_name <- paste0(hla_allele, affinity_pred_flag)
      } else if (self$type == 'STS') {
        self$table_name <-
          sprintf('STS_%s_%.1f%s', 
                  hla_allele, STS_percentile_rank, affinity_pred_flag)
      }
    },
    #' Lookups are performed in serial batches 
    #' in order to not multiply use
    #' the database connection associated with the R6 object
    lookup = function(query_peps, batch_size = 1e5) {
      if (!self$check_table_existence()) {
        maartenutils::mymessage(instance = self$table_name,
                                msg = sprintf('Table not found %s', 
                                              self$table_name))
        return(NULL)
      }
      query_peps <- self$sanitize_query(query_peps)

      ## Lookup in batches such that the SQL string doesn't become too
      ## long
      N_batches <- ceiling(length(query_peps) / batch_size)
      query_peps_s <- split(query_peps,
        ceiling(seq(1, length(query_peps)) / batch_size))

      out <- plyr::llply(seq_along(query_peps_s), function(idx) {
        qpl <- query_peps_s[[idx]]
        quoted_peps <- paste(sprintf("'%s'", qpl), collapse = ', ')
        sql <- sprintf('SELECT * FROM "%s" WHERE peptide IN (%s);',
          self$table_name, quoted_peps)
        res <- RPostgreSQL::dbGetQuery(self$conn, sql)
        return(self$sanitize_res(res))
      }) %>% rbindlist(fill = TRUE)

      return(out)
    },
    sanitize_query = function(peptides) {
      if (is.null(peptides) || length(peptides) == 0) {
        return(NULL)
      } else {
        return(unique(setdiff(peptides, c('', NA))))
      }
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

      if ('score_el' %in% colnames(dtf))
        dtf[, score_el := as.numeric(score_el)]

      if ('rank_el' %in% colnames(dtf))
        dtf[, rank_el := as.numeric(rank_el)]

      if ('bind_level' %in% colnames(dtf))
        dtf[, bind_level := as.factor(bind_level)]

      dtf <- dtf[!is.na(peptide)]
      return(dtf)
    },
    store_res = function(dtf) {
      dtf <- self$sanitize_res(dtf)
      if (self$check_table_existence()) {
        tryCatch({
          RPostgreSQL::dbWriteTable(self$conn, self$table_name, dtf,
                                    row.names = F, append = T)
        }, error = function(e) { print(e) })
      } else {
        tryCatch({
          RPostgreSQL::dbWriteTable(self$conn, self$table_name, dtf,
            row.names = F)
        }, error = function(e) { print(e) })
      }
    },
    unique_SQL = function() {
      if (!self$check_table_existence()) {
        return(NULL)
      }
      dup_entries <- tryCatch({
        RPostgreSQL::dbGetQuery(self$conn, glue::glue('
          select count(*) from "{self$table_name}" group by \\
          "peptide" having count(*) > 1;'))
      }, error = function(e) { print(e) })
      if (maartenutils::null_dat(dup_entries) || 
          is.null(dup_entries) || is.na(dup_entries) || 
          length(dup_entries) == 0 || dup_entries == 0) {
        return(invisible(NULL))
      }
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
      if (!self$check_table_existence()) {
        maartenutils::mymessage(instance = sprintf('quickMHC: %s',
                                                   self$table_name),
                                msg = 'Table does not exist')
        return(NULL)
      }
      self$unique_SQL()
      tryCatch({
        RPostgreSQL::dbSendQuery(self$conn,
          glue::glue('create UNIQUE index if not exists "{self$table_name}_idx"
            ON "{self$table_name}" (peptide);'))
      }, error = function(e) { print(e) })
    },
    compute = function(peptides, batch_size = 1e5, ncores = 1, 
      verbose = F, overwrite = F) {
      if (is.null(peptides) || length(peptides) == 0 || 
          all(is.na(peptides)))
        return(NULL)
      peptides <- self$sanitize_query(peptides)
      already_computed <- self$lookup(peptides)
      if (overwrite) {
        self$remove_peptides(already_computed$peptide)
      }

      if (length(peptides) == 0) {
        return(NULL)
      } else {
        N_batches <- ceiling(length(peptides) / batch_size)
        query_peps_s <- split(peptides,
          ceiling(seq(1, length(peptides)) / batch_size))
        if (ncores > 1) 
          doParallel::registerDoParallel(cores = ncores)
        compute_res <- plyr::llply(seq_along(query_peps_s), 
          function(idx) {
          qpl <- query_peps_s[[idx]]
          maartenutils::mymessage(instance = self$table_name,
            msg = sprintf('computing batch %d/%d (%d peptides)', 
              idx, N_batches, length(qpl)))
          if (self$verbose) {
            maartenutils::mymessage(instance = self$table_name,
              msg = sprintf('computing %s', paste(qpl, collapse = ', ')))
          }
          lres <- self$computer(qpl)

          if (overwrite) {
            self$store_res(lres)
          } else {
            save_subs <- 
              lres[!peptide %in% already_computed$peptide]
            self$store_res(lres)
          }

          return(lres)
        }, .parallel = ncores > 1) %>% rbindlist(fill = T)

        self$index_SQL()

        return(compute_res)
      }
    },
    check_table_existence = function(appendix = '') {
      return(RPostgreSQL::dbExistsTable(self$conn, 
          sprintf('%s%s', self$table_name, appendix)))
    },
    #' Query the backend database, computing results that are missing
    #' from the query
    #'
    #'  across
    #' multiple threads; computations are done in parallel.
    #'
    #' @param prefill_mode Do not return results, just store them in
    #' the DB. If the latter is all we care about, we can save some
    #' RAM this way.
    #'
    query = function(peptides, 
      batch_size = 1e5, 
      prefill_mode = FALSE, 
      ncores = 1) {
      peptides <- self$sanitize_query(peptides)

      if (length(peptides) == 0) return(NULL)

      already_computed <- self$lookup(peptides)

      if (!maartenutils::null_dat(already_computed)) {
        to_be_computed <- setdiff(peptides, already_computed$peptide)
      } else {
        to_be_computed <- peptides
      }

      if (length(to_be_computed) > 0) {
        if (T) {
          computed_out <- self$compute(
            to_be_computed, batch_size = batch_size, ncores = ncores
          )
        } else {
          # N_batches <- ceiling(length(to_be_computed) / batch_size)
          # tbc_l <- split(to_be_computed,
          #   ceiling(seq(1, length(to_be_computed)) / batch_size))

          # if (ncores > 1) 
          #   doParallel::registerDoParallel(cores = ncores)
          # computed_out <- plyr::llply(seq_along(tbc_l), 
          #   function(idx) {
          #   qpl <- tbc_l[[idx]]
          #   # if (self$verbose) {
          #   #   maartenutils::mymessage(instance = self$table_name,
          #   #     msg = sprintf('querying batch %d/%d (%d peptides)', 
          #   #       idx, N_batches, length(qpl)))
          #   # }
          #   l_res <- self$compute(qpl, batch_size = length(qpl), 
          #         ncores = ncores)
          #   if (!prefill_mode) {
          #     return(l_res)
          #   } else {
          #     return(NULL)
          #   }
          # }, .parallel = F) %>% rbindlist(fill = T)
        }
      } else {
        computed_out <- NULL
      }

      if (prefill_mode) {
        out <- NULL
      } else {
        if (!maartenutils::null_dat(already_computed)) {
          out <- rbind(already_computed, computed_out)
        }
        out <- out[match(peptides, out$peptide, nomatch = 0L), ]
      }

      return(out)
    },
    remove_peptides = function(peptides) {
      peptides <- self$sanitize_query(peptides)
      if (is.null(peptides)) return(NULL)
      RPostgreSQL::dbSendQuery(self$conn,
        sprintf('DELETE FROM "%s" WHERE "peptide" IN (%s)', self$table_name,
          paste(sprintf("'%s'", peptides), collapse = ', ')))
      invisible()
    },
    drop_table = function() {
      RPostgreSQL::dbSendQuery(self$conn,
        sprintf('DROP TABLE "%s"', self$table_name))
      invisible()
    }
))
