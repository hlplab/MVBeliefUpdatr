#' @importFrom stringr str_sub

write_NIW_ibbu_stanfit = function(fit, file) {
  assert_that(is.NIW_ibbu_stanfit(fit))
  if (str_sub(file, - 4, - 1) %nin% c(".rds", ".RDS"))
    saveRDS(fit, file, compress = T) else
      saveRDS(fit, paste0(file, ".rds", compress = T))
}


read_NIW_ibbu_stanfit = function(file) {
  assert_that(is_scalar_character(file))

  if (file.exists(file))
    fit = readRDS(file) else
      if (str_sub(file, - 4, - 1) %nin% c(".rds", ".RDS", ".Rds"))
        return(read_NIW_ibbu_stanfit(paste0(file, ".rds"))) else
          return(NULL)

  if (is.NIW_ibbu_stanfit(fit))
    return(fit) else stop(paste("File", file, "does not contain an NIW_ibbu_stanfit."))
}
