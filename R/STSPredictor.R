load_self_epitope_list = function(
  hla_allele = 'A0201', STS_percentile_rank = 1.9, peptidelength = 9,
  processing_threshold = .5, 
  path = '/DATA/resources/neolution_selflists/pan_3.0') {

  availableSelfLists <- dir(path = path,
    pattern = paste0(hla_allele, '.*', paste0(peptidelength, 'mer')),
    include.dirs = FALSE,
    full.names = TRUE)

  if (length(availableSelfLists) == 1) {
    ## Prefilter using GAWK
    # awk_PF <- 'awk -F \",\" \'$9 >= %.1f && $8 <= %.1f {print $4, $9, $8} %s\''
    # maartenutils::systemf(awk_PF, processing_threshold, STS_percentile_rank,
    #   availableSelfLists)

    mymessage(msg = 'Loading self lists')
    selfEpitopes <- fread(availableSelfLists, header = TRUE,
      stringsAsFactors = FALSE, showProgress = F) %>%
      unique(by = 'peptide')
    pr_var <- sprintf('%spercentile_rank', shortenHLA(hla_allele))
    setkeyv(selfEpitopes, c(pr_var, 'processing_score'))
    selfEpitopes <- selfEpitopes[
      get(pr_var) <= STS_percentile_rank & 
      processing_score >= processing_threshold, peptide]
    mymessage(msg = 'Finished loading and filtering self list')
    return(selfEpitopes)
  } else {
    stop(paste('Zero or more than one self-epitope lists found, please make',
        'sure only 1 list per HLA & peptide length is present in ', path))
  }
}

load_selfsimilarity_matrix <- function() {
  scoreMatrix = matrix(c(1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,1,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,
                         0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,
                         0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,
                         0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,
                         0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,
                         0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,
                         0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,
                         0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,
                         0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,
                         0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1),
                       ncol = 21)
  rownames(scoreMatrix) = c('W', 'F', 'Y', 'I', 'V', 'L', 'M', 'C', 'D', 'E', 'G',
                            'A', 'P', 'H', 'K', 'R', 'S', 'T', 'N', 'Q', 'X')
  colnames(scoreMatrix) = c('W', 'F', 'Y', 'I', 'V', 'L', 'M', 'C', 'D', 'E', 'G',
                            'A', 'P', 'H', 'K', 'R', 'S', 'T', 'N', 'Q', 'X')
  return(scoreMatrix)
}

STSPredictor <- R6::R6Class('STSPredictor',
  inherit = Predictor,
  public = list(
    self_peptides = NULL,
    expected_res_columns = c('peptide', 'different_from_self'),
    type = 'STS',
    initialize = function(hla_allele=NULL, STS_percentile_rank=NULL) {
      if (!is.null(STS_percentile_rank)) {
        self$STS_percentile_rank <- STS_percentile_rank
      } else {
        self$STS_percentile_rank <- 1.9
      }
      super$initialize(hla_allele)
    },
    computer = function(peptides, inp_self_peptides = NULL, 
      score_matrix = load_selfsimilarity_matrix()) {
      ## Persist self_peptides throughout object's lifetime, loading it in can
      ## be costly timewise
      if (!is.null(inp_self_peptides) && is.null(self$self_peptides)) {
        self$self_peptides <- inp_self_peptides
      } else if (is.null(self$self_peptides)) {
        self$self_peptides <- load_self_epitope_list(
          hla_allele = self$hla_allele, STS_percentile_rank = self$STS_percentile_rank)
      }
      ## TODO properly integrate into R package so we don't have to recompile
      ## upon reloading
      Rcpp::sourceCpp(file.path('~/libs/quickMHC/cpp/match_ext_seq.cpp'))
      res <- data.table('peptide' = peptides)
      res$different_from_self <- match_seq_ext_Rcpp(query_peps = res$peptide, 
        ref_list = self$self_peptides, scorematrix = score_matrix)
      return(self$sanitize_res(res))
    }
  ),
  private = list(
  )
)
