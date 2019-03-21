rm(list = ls())
devtools::load_all(file.path('~/libs', 'maartenutils'))
devtools::load_all(file.path('~/libs', 'quickMHC'))

test <- BindingPredictor$new(hla_allele = 'A1101')
test$remove_peptides(c('SYFPEITHI', 'ELVISLIVE', 'ELVISLIVE'))
test$query(c('SYFPEITHI', 'ELVISLIVE', 'ELVISLIVE'))
test <- BindingPredictor$new(hla_allele = 'A0201')
test$remove_peptides(c('SYFPEITHI', 'ELVISLIVE', 'ELVISLIVE'))
test$query(c('SYFPEITHI', 'ELVISLIVE', 'ELVISLIVE'))
test <- BindingPredictor$new(hla_allele = 'B2705')
test$remove_peptides(c('SYFPEITHI', 'ELVISLIVE', 'ELVISLIVE'))
test$query(c('SYFPEITHI', 'ELVISLIVE', 'ELVISLIVE'))


# SP <- test$self_peptides
test <- STSPredictor$new(hla_allele = 'A1101', STS_percentile_rank = 1.9)
test$self_peptides
test$table_name
test$compute(c('SYFPEITHI', 'ELVISLIVE', 'ELVISLIVE'))
test$remove_peptides(c('SYFPEITHI', 'ELVISLIVE', 'ELVISLIVE'))
test$lookup(c('SYFPEITHI', 'ELVISLIVE', 'ELVISLIVE'))

peptides <- c('SYFPEITHI', 'ELVISLIVE', 'ELVISLIVE')
quickMHC(query_table = data.frame(
    'hla' = c('A0201', 'A0201', 'A0201'),
    'peptide' = c('SYFPEITHI', 'ELVISLIVE', 'ELVISLIVE')))
