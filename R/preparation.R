if (0) {
  usethis::use_description(
    fields = list(
      Package     = "hintreg",
      Title       = "ML estimation of a heterogenous normal interval regression model",
      Description = "ML estimation of a heterogenous normal interval regression model,
        i.e., a normal interval regression model with conditional heteroskedasticity
        and heterogensou indifferece limens",
      `Authors@R` = person("Yasutomo", "Murasawa", , "yasutomo.murasawa@gmail.com",
        c("aut", "cre"), c(ORCID = "0000-0001-6716-0162"))
    )
  )
  usethis::use_gpl3_license()
  usethis::use_testthat(3)
  usethis::use_data_raw()
  usethis::use_vignette("hintreg-intro")
  usethis::create_github_token()
  gitcreds::gitcreds_set()
  usethis::use_github()
  usethis::use_package("generics", "Imports")
}
