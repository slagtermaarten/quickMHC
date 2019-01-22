shortenHLA <- function(hla_allele = 'HLA-A*02:01') {
  gsub('\\-|\\*|\\:|HLA', '', toupper(hla_allele))
}

expandHLA <- function(hla_allele = 'A*0201') {
  ## First get rid of all fluff, only to introduce some again
  hla_allele <- shortenHLA(hla_allele)
  sprintf('HLA-%s:%s', substring(hla_allele, 1, 3), substring(hla_allele, 4, 5)) 
}
