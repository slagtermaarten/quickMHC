wrapper_netmhcpan4 <- function(
  peptides = c('SYFPEITHI', 'ELVISLIVE'), 
  hla_allele = 'A0201',
  peptidelength = 9,
  netMHCpan = '/DATA/users/m.slagter/libs/netMHCpan-4.1/netMHCpan') {

  temp_dir <- file.path(tempdir(), 
    sample(c(letters), size = 8, replace = T) %>% 
      paste(collapse = ''))
  dir.create(temp_dir, showWarnings = F)

  peptides_fn <- tempfile(tmpdir = temp_dir, fileext = '.fas')
  writeLines(peptides, con = peptides_fn)

  ## write peptides to disk in temp dir
  command <- paste(
    netMHCpan,
    '-a', expandHLA(hla_allele),
    '-l', peptidelength,
    # '-inptype 1', 
    # '-tdir', netmhc_temp,
    '-p',
    '-f', peptides_fn)
  output <- system(command = command, intern = TRUE)
  start_stop_idxs <- 
    which(grepl('----------', output))[c(2, 3)] + c(1, -1)
  # header <- 
  # gsub('%', '', tolower(strsplit(output[49], '\\s+')[[1]]))
  output_s <- output[start_stop_idxs[1]:start_stop_idxs[2]]
  dtf <- rbindlist(lapply(output_s, function(l) {
    fields <- strsplit(l, '\\s+')[[1]]
    c(as.list(fields[c(4, 13, 14)]), 
      list(paste0(fields[15], ' ', fields[16])))
  }))
  setnames(dtf, c('peptide', 'score_el', 'rank_el', 'bindlevel'))
  dtf[, score_el := as.numeric(score_el)]
  dtf[, rank_el := as.numeric(rank_el)]
  unlink(temp_dir)
  return(dtf)
}


wrapper_netmhcpan3 <- function(
  peptides = c('SYFPEITHI', 'ELVISLIVE'), 
  hla_allele = 'A0201',
  peptidelength = 9,
  netMHCpan = '/DATA/resources/predictors/netMHCpan-3.0/netMHCpan') {

  temp_dir <- file.path(tempdir(), 
    sample(c(letters), size = 8, replace = T) %>% paste(collapse = ''))
  dir.create(temp_dir, showWarnings = F)

  peptides_fn <- tempfile(tmpdir = temp_dir, fileext = '.fas')
  writeLines(peptides, con = peptides_fn)

  ## write peptides to disk in temp dir
  command <- paste(
    netMHCpan,
    '-a', expandHLA(hla_allele),
    '-l', peptidelength,
    # '-inptype 1', 
    # '-tdir', netmhc_temp,
    '-p',
    '-f', peptides_fn)
  output <- system(command = command, intern = TRUE)
  start_stop_idxs <- which(grepl('----------', output))[c(2, 3)] + c(1, -1)
  output <- output[start_stop_idxs[1]:start_stop_idxs[2]]
  dtf <- rbindlist(lapply(output, function(l) {
    fields <- strsplit(l, '\\s+')[[1]]
    as.list(fields[c(4, 13, 14, 15)])
  }))
  setnames(dtf, 
    c('peptide', 'peptide_score_log50k', 'affinity', 'percentile_rank'))
  dtf[, affinity := as.numeric(affinity)]
  dtf[, percentile_rank := as.numeric(percentile_rank)]
  dtf[, peptide_score_log50k := as.numeric(peptide_score_log50k)]
  unlink(temp_dir)
  return(dtf)
}


BindingPredictor <- BindingPredictor3 <- 
  R6::R6Class('BindingPredictor',
  inherit = Predictor,
  public = list(
    expected_res_columns = c('peptide', 'peptide_score_log50k', 
      'affinity', 'percentile_rank'),
    type = 'binding_affinity',
    affinity_predictor = 'NetMHCpan3',
    computer = function(peptides) {
      wrapper_netmhcpan3(peptides = peptides, 
        hla_allele = self$hla_allele) 
    }
  ),
  private = list()
)


BindingPredictor4 <- R6::R6Class('BindingPredictor',
  inherit = Predictor,
  public = list(
    expected_res_columns = c('peptide', 'score_el', 
      'rank_el', 'bindlevel'),
    type = 'binding_affinity',
    affinity_predictor = 'NetMHCpan4',
    computer = function(peptides) {
      wrapper_netmhcpan4(peptides = peptides, 
        hla_allele = self$hla_allele) 
    }
  ),
  private = list()
)
