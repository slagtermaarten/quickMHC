## List all potential config dirs
potential_config_locs <- 
  c(getwd(),
    '~/.config',
    file.path('~/libs', 'quickMHC'),
    sapply(.libPaths(), function(p_dir) file.path(p_dir, 'quickMHC')))

config_fn <-
  sapply(potential_config_locs, list.files, pattern = '.*quickMHC.*\\.yaml',
    full.names = T) %>%
  .[sapply(., length) == 1] %>%
  .[[1]]

if (is.null(config_fn) || is.na(config_fn) || length(config_fn) == 0) {
  stop('Could not find a neolution_config.yaml anywhere')
}

config <- yaml::read_yaml(config_fn)

if (is.null(config) || 
  config$db_user == 'db_user' || is.na(config$db_user ||
  config$db_name == 'db_name' || is.na(config$db_name))) {
  dir_string <- paste(potential_config_locs, collapse = ', ')
  stop(glue::glue('Define db_user and db_name in config files named\\
      quickMHC_config.yaml in one of these directories {dir_string}'))
}
