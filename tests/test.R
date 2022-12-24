rm(list = ls())
devtools::load_all(file.path('~/libs', 'maartenutils'))
devtools::load_all(file.path('~/libs', 'quickMHC'))
get('BindingPredictor')

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
test <- STSPredictor$new(hla_allele = 'A1101', 
  STS_percentile_rank = 1.9)
test$self_peptides
test$table_name
test$compute(c('SYFPEITHI', 'ELVISLIVE', 'ELVISLIVE'))
test$remove_peptides(c('SYFPEITHI', 'ELVISLIVE', 'ELVISLIVE'))
test$lookup(c('SYFPEITHI', 'ELVISLIVE', 'ELVISLIVE'))

peptides <- c('SYFPEITHI', 'ELVISLIVE', 'ELVISLIVE')
quickMHC(query_table = data.frame(
    'hla' = c('A0201', 'A0201', 'A0201'),
    'peptide' = c('SYFPEITHI', 'ELVISLIVE', 'ELVISLIVE')))

devtools::load_all(file.path('~/libs', 'quickMHC'))
bp <- quickMHC::BindingPredictor$new(hla_allele = 'B7301')
bp <- quickMHC::BindingPredictor$new(hla_allele = 'A2401')
bp <- quickMHC::BindingPredictor$new(hla_allele = 'B1501')
bp <- quickMHC::BindingPredictor$new(hla_allele = 'A0201')
debugonce(bp$lookup)
debugonce(bp$query)
bp$lookup('AILEVCGUK')
bp$lookup('AILEVCGXK')
bp$lookup(c('AILEVCGXK', 'AILEVCGXK'))
bp$query(c('AILEVCGXK', 'AILEVCGUK'), ncores = 16)
bp$compute(c('AILEVCGXK', 'AILEVCGXK'), overwrite = F)
bp$remove_peptides('AILEVCGXK')

    rep_peps <- preds[, (.N / 4) %% 1 > 0, by = .(peptide)][V1 == T, peptide]
    rep_peps <- c('UKLGRFPQV', 'SPPPMAGGU', 'ILEVCGUKL', 'AILEVCGUK',
      'PPMAGGUGR', 'GAILEVCGU', 'PPPMAGGUG', 'EVCGUKLGR', 'VCGUKLGRF',
      'LEVCGUKLG', 'GUKLGRFPQ', 'CGUKLGRFP')
    preds[peptide %in% rep_peps]
    dput(preds[100:110, peptide])
    control_peps <- c('YLDNGVVFV', 'YLSGIAHFL', 'MMFDEPVLL', 'KLIEKNYFL',
      'FLQEYVANL', 'FMSEYLIEL', 'VLFEVAWEV', 'MLLELMTEV', 'YLLFGISKV',
      'SLMEVRFYV', 'FLYDHIQPV')
    devtools::load_all(file.path('~/libs', 'quickMHC'))
    bp <- quickMHC::BindingPredictor$new(hla_allele = 'A0201')
    bp$index_SQL()

    quoted_peps <- paste(sprintf("'%s'", rep_peps), collapse = ', ')
    sql <- sprintf('SELECT * FROM "%s" WHERE peptide IN (%s);', 
      'A0201', quoted_peps)
    res <- RPostgreSQL::dbGetQuery(bp$conn, sql) %>% as.data.table

    quoted_peps <- paste(sprintf("'%s'", control_peps), collapse = ', ')
    sql <- sprintf('SELECT * FROM "%s" WHERE peptide IN (%s);', 
      'A0201', quoted_peps)
    gsub('\\\\', '', sql, fixed = T)
    res <- RPostgreSQL::dbGetQuery(bp$conn, sql) %>% as.data.table

    bp$remove_peptides(rep_peps)
    bp$query(rep_peps)
    debugonce(bp$lookup)
    bp$lookup(rep_peps)
    bp$query(rep_peps)

    bp$query(control_peps)
    bp$lookup(control_peps)
    bp$remove_peptides(control_peps)
