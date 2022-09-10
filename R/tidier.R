#' @importFrom generics tidy
#' @export
generics::tidy
#' Tidy an intreg object
#'
#' @param lFit an object of class "intreg".
#' @param conf.int logical indicating whether or not to include a confidence interval in the tidied output.
#' @param conf.level the confidence level to use for the confidence interval if conf.int = TRUE.
#' @param ... additional arguments not used, but needed to match generic signature only.
#' @return a tidy [tibble::tibble()] summarizing component-level information about the model.
#' @examples
#' lout <- intreg(
#'   q2 ~ sex,
#'   data       = ias2009febmay,
#'   thresholds = c(-.5, .5, 1, 2, 3, 4, 5)
#' )
#' tidy(lout)
#' @export
tidy.intreg <- function(lFit, conf.int = FALSE, conf.level = .95, ...) {
  msummary <- summary(lFit)$summary |>
    tibble::as_tibble(rownames = "term") |>
    dplyr::rename(
      estimate  = Estimate,
      std.error = `Std. Error`,
      statistic = `t value`,
      p.value   = `Pr(>|t|)`
    )
  if (conf.int) {
    mci      <- confint.intreg(lFit, level = conf.level)
    msummary <- dplyr::left_join(msummary, mci, by = "term")
  }
  msummary
}
confint.intreg <- function(lFit, level = .95, ...) {
  vbeta     <- coef(lFit)
  dlsigma   <- lFit$logscale
  mhess     <- lFit$hess
  vpar      <- c(vbeta, dlsigma)
  mcov      <- solve(mhess)
  vse       <- sqrt(diag(mcov))
  vname     <- names(vpar)
  dalpha    <- (1 - level) / 2
  conf.low  <- vpar + qnorm(dalpha) * vse
  conf.high <- vpar + qnorm(1 - dalpha) * vse
  cbind(conf.low, conf.high) |> tibble::as_tibble(rownames = "term")
}
