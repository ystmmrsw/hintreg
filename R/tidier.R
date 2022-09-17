#' @importFrom generics tidy
#' @export
generics::tidy
#' Tidy an intreg object
#'
#' @param x an object of class "intreg".
#' @param conf.int logical indicating whether or not to include a confidence interval in the tidied output.
#' @param conf.level the confidence level to use for the confidence interval if conf.int = TRUE.
#' @param ... additional arguments not used, but needed to match generic signature only.
#' @return a tidy [tibble::tibble()] summarizing component-level information about the model.
#' @examples
#' lout <- intreg(
#'   q2 ~ sex,
#'   data       = ias2009febmay,
#'   thresholds = c(-.5, .5, 1:5)
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
  vbeta   <- coef(lFit)
  dlsigma <- lFit$logscale
  mhess   <- lFit$hess
  vpar    <- c(vbeta, dlsigma)
  mcov    <- solve(mhess)
  vse     <- sqrt(diag(mcov))
  dalpha  <- (1 - level) / 2
  tibble::tibble(
    term      = names(vpar),
    conf.low  = vpar + qnorm(dalpha) * vse,
    conf.high = vpar + qnorm(1 - dalpha) * vse
  )
}
#' Tidy an gintreg object
#'
#' @param x an object of class "gintreg".
#' @inheritParams tidy.intreg
#' @return a tidy [tibble::tibble()] summarizing component-level information about the model.
#' @examples
#' lout <- gintreg(
#'   q2 ~ sex,
#'      ~ sex,
#'   data       = ias2009febmay,
#'   thresholds = c(-.5, .5, 1:5)
#' )
#' tidy(lout)
#' @export
tidy.gintreg <- function(lFit, conf.int = FALSE, conf.level = .95, ...) {
  msummary <- summary(lFit)$summary |>
    tibble::as_tibble(rownames = "term") |>
    dplyr::rename(
      estimate  = Estimate,
      std.error = `Std. Error`,
      statistic = `t value`,
      p.value   = `Pr(>|t|)`
    )
  msummary <- tibble::tibble(
    coef.type = c(
      rep("location", length(lFit$coefL)),
      rep("scale", length(lFit$coefS))
    ),
    msummary
  )
  if (conf.int) {
    mci      <- confint.gintreg(lFit, level = conf.level)
    msummary <- dplyr::left_join(msummary, mci, by = c("coef.type", "term"))
  }
  msummary
}
confint.gintreg <- function(lFit, level = .95, ...) {
  vbeta  <- lFit$coefL
  vdelta <- lFit$coefS
  mhess  <- lFit$hess
  vpar   <- c(vbeta, vdelta)
  mcov   <- solve(mhess)
  vse    <- sqrt(diag(mcov))
  dalpha <- (1 - level) / 2
  tibble::tibble(
    coef.type = c(
      rep("location", length(lFit$coefL)),
      rep("scale", length(lFit$coefS))
    ),
    term      = names(vpar),
    conf.low  = vpar + qnorm(dalpha) * vse,
    conf.high = vpar + qnorm(1 - dalpha) * vse
  )
}
#' Tidy an hintreg object
#'
#' @param x an object of class "hintreg".
#' @inheritParams tidy.intreg
#' @return a tidy [tibble::tibble()] summarizing component-level information about the model.
#' @examples
#' lout <- hintreg(
#'   q2 ~ sex,
#'      ~ sex,
#'      ~ sex,
#'      ~ sex,
#'   data       = ias2009febmay,
#'   thresholds = c(-.5, .5, 1:5)
#' )
#' tidy(lout)
#' @export
tidy.hintreg <- function(lFit, conf.int = FALSE, conf.level = .95, ...) {
  msummary <- summary(lFit)$summary |>
    tibble::as_tibble(rownames = "term") |>
    dplyr::rename(
      estimate  = Estimate,
      std.error = `Std. Error`,
      statistic = `t value`,
      p.value   = `Pr(>|t|)`
    )
  msummary <- tibble::tibble(
    coef.type = c(
      rep("location", length(lFit$coefL)),
      rep("scale", length(lFit$coefS)),
      rep("lower", length(lFit$coefLower)),
      rep("upper", length(lFit$coefUpper))
    ),
    msummary
  )
  if (conf.int) {
    mci      <- confint.hintreg(lFit, level = conf.level)
    msummary <- dplyr::left_join(msummary, mci, by = c("coef.type", "term"))
  }
  msummary
}
confint.hintreg <- function(lFit, level = .95, ...) {
  vbeta    <- lFit$coefL
  vdelta   <- lFit$coefS
  vlambda  <- lFit$coefLower
  vupsilon <- lFit$coefUpper
  mhess    <- lFit$hess
  vpar     <- c(vbeta, vdelta, vlambda, vupsilon)
  mcov     <- solve(mhess)
  vse      <- sqrt(diag(mcov))
  dalpha   <- (1 - level) / 2
  tibble::tibble(
    coef.type = c(
      rep("location", length(lFit$coefL)),
      rep("scale", length(lFit$coefS)),
      rep("lower", length(lFit$coefLower)),
      rep("upper", length(lFit$coefUpper))
    ),
    term      = names(vpar),
    conf.low  = vpar + qnorm(dalpha) * vse,
    conf.high = vpar + qnorm(1 - dalpha) * vse
  )
}
