#' Fit a normal interval regression model to interval data
#'
#' \code{intreg} fits a normal interval regression model to interval data.
#' The dependent variable must be a factor representing J + 1 intervals,
#' and the user must provide J thresholds.
#'
#' @param formula an object of class "formula".
#' @param data a data.frame containing the variables in the model.
#' @param start an optional vector of starting values for the parameters.
#' @param weights an optional vector of weights.
#' @param thresholds a vector of J thresholds when there are J + 1 intervals.
#' @return an object of class "intreg"
#' @export
intreg <- function(formula, data, start, weights, thresholds) {
  mf      <- match.call(expand.dots = FALSE)
  m       <- match(c("formula", "data", "weights", "offset"), names(mf), 0)
  mf      <- mf[c(1, m)]
  mf[[1]] <- quote(model.frame)
  mf      <- eval.parent(mf)
  vy      <- model.response(mf)
  if (!is.factor(vy)) stop("response must be a factor")
  mx      <- model.matrix(formula, mf)
  vw      <- model.weights(mf)
  voffset <- model.offset(mf)
  vthresh <- thresholds
  cn      <- nrow(mx)
  ck      <- ncol(mx)
  if (is.null(voffset)) voffset <- rep(0, cn)
  if (length(vthresh) != nlevels(vy) - 1) stop("incorrect number of thresholds")
  if (missing(start)) {
    vstart_beta   <- numeric(ck)
    vstart_lsigma <- log(max(vthresh) - min(vthresh))
    vstart        <- c(vstart_beta, vstart_lsigma)
  } else {
    vstart <- start
  }
  lout            <- intreg.fit(vy, mx, vw, voffset, vthresh, vstart)
  vtheta          <- lout$par
  mhess           <- lout$hessian
  vbeta           <- vtheta[seq_len(ck)]
  dlsigma         <- vtheta[ck + 1]
  names(vbeta)    <- colnames(mx)
  names(dlsigma)  <- "log s.d."
  vparname        <- c(names(vbeta), names(dlsigma))
  dimnames(mhess) <- list(vparname, vparname)
  lfit            <- list(
    call         = match.call(),
    iteration    = lout$counts["function"],
    coefficients = vbeta,
    logscale     = dlsigma,
    logLik       = -lout$value,
    edf          = ck + 1,
    hess         = mhess
  )
  class(lfit) <- "intreg"
  lfit
}
intreg.fit <- function(vY, mX, vW, vOffset, vThresh, vStart) {
  optim(
    vStart,
    function(.) -intreg.loglikelihood(., vY, mX, vW, vOffset, vThresh),
    method  = "BFGS",
    hessian = TRUE
  )
}
intreg.loglikelihood <- function(vTheta, vY, mX, vW, vOffset, vThresh) {
  ck      <- ncol(mX)
  vbeta   <- vTheta[seq_len(ck)]
  dlsigma <- vTheta[ck + 1]
  vc      <- c(-Inf, vThresh, Inf)
  vy      <- unclass(vY)
  vmu     <- vOffset
  if (ck > 0) vmu <- vmu + drop(mX %*% vbeta)
  dsigma <- exp(dlsigma)
  vp     <- pnorm((vc[vy + 1] - vmu) / dsigma) - pnorm((vc[vy] - vmu) / dsigma)
  if (is.null(vW)) sum(log(vp)) else sum(vW * log(vp))
}
#' @export
print.intreg <- function(lFit, ...) {
  if (!is.null(lFit$call)) {
    cat("Call:\n")
    dput(lFit$call, control = NULL)
  }
  cat("\nlog-likelihood:", format(lFit$logLik, nsmall = 2))
  cat("\nAIC:", format(-2 * lFit$logLik + 2 * lFit$edf, nsmall = 2))
  cat("\n")
  if (length(coef(lFit)) > 0) {
    cat("\nCoefficients:\n")
    print(coef(lFit), ...)
  } else {
    cat("\nNo coefficients\n")
  }
  cat("\nScale:\n")
  print(lFit$logscale, ...)
  invisible(lFit)
}
#' @export
summary.intreg <- function(lFit, correlation = FALSE, ...) {
  vbeta    <- coef(lFit)
  dlsigma  <- lFit$logscale
  mhess    <- lFit$hess
  vpar     <- c(vbeta, dlsigma)
  mcov     <- solve(mhess)
  vse      <- sqrt(diag(mcov))
  vt       <- vpar / vse
  vp       <- 2 * pnorm(abs(vt), lower.tail = FALSE)
  msummary <- cbind(
    Estimate     = vpar,
    `Std. Error` = vse,
    `t value`    = vt,
    `Pr(>|t|)`   = vp
  )
  lFit$summary <- msummary
  if (correlation == TRUE) lFit$correlation <- mcov / outer(vse, vse)
  class(lFit) <- "summary.intreg"
  lFit
}
#' @export
print.summary.intreg <- function(lFit, digits = max(3, .Options$digits - 3), ...) {
  if (!is.null(lFit$call)) {
    cat("Call:\n")
    dput(lFit$call, control = NULL)
  }
  cat("\nNo. of iterations:", lFit$iteration)
  cat("\nlog-likelihood:", format(lFit$logLik, nsmall = 2))
  cat("\nAIC:", format(-2 * lFit$logLik + 2 * lFit$edf, nsmall = 2))
  cat("\n")
  cat("\nParameters:\n")
  printCoefmat(lFit$summary, digits = digits, signif.stars = TRUE, na.print = "NA", ...)
  if (!is.null(lFit$correlation)) {
    cat("\nCorrelation of parameters:\n")
    mcorr                    <- format(round(lFit$correlation, digits))
    mcorr[!lower.tri(mcorr)] <- ""
    print(mcorr[-1, -ncol(mcorr)], quote = FALSE, ...)
  }
  invisible(lFit)
}
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
