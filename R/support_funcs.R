shortenHLA <- function(hla_alleles = 'HLA-A*02:01') {
  vapply(hla_alleles, function(hla_allele) {
    toupper(hla_allele) %>%
      { gsub('\\-|\\*|\\:|HLA|_', '', .) } 
  }, character(1))
}

expandHLA <- function(hla_alleles = 'A*0201') {
  ## First get rid of all fluff, only to introduce some again in order for
  ## NetMHCpan to understand the HLA type
  hla_alleles <- shortenHLA(hla_alleles)
  vapply(hla_alleles, function(hla_allele)
    sprintf('HLA-%s:%s', substring(hla_allele, 1, 3),
      substring(hla_allele, 4, 5)),
    character(1))
}

ppHLA <- function(hla_alleles = 'A0201') {
  ## First get rid of all fluff, only to introduce some again in order for
  ## NetMHCpan to understand the HLA type
  hla_alleles <- shortenHLA(hla_alleles)
  vapply(hla_alleles, function(hla_allele)
    sprintf('HLA-%s*%s:%s', 
      substring(hla_allele, 1, 1),
      substring(hla_allele, 2, 3),
      substring(hla_allele, 4, 5)),
    character(1))
}
