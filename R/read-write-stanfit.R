#' @importFrom stringr str_sub

write_NIW_ideal_adaptor_stanfit = function(fit, file) {
  assert_that(is.NIW_ideal_adaptor_stanfit(fit))
  assert_that(is_scalar_character(file))
  if (tolower(str_sub(file, - 4, - 1)) == ".rds")
    saveRDS(fit, file, compress = T) else
      saveRDS(fit, paste0(file, ".rds"), compress = T)
}


read_NIW_ideal_adaptor_stanfit = function(file) {
  assert_that(is_scalar_character(file))

  if (file.exists(file))
    fit = readRDS(file) else
      if (tolower(str_sub(file, - 4, - 1)) != ".rds")
        return(read_NIW_ideal_adaptor_stanfit(paste0(file, ".rds"))) else
          return(NULL)

  if (is.NIW_ideal_adaptor_stanfit(fit))
    return(fit) else stop(paste("File", file, "does not contain an NIW_ideal_adaptor_stanfit."))
}
