#' Fit a generalized normal interval regression model to interval data
#'
#' \code{gintreg} fits a generalized normal interval regression model,
#' i.e., a normal interval regression model with conditional heteroskedasticity,
#' to interval data.
#' The dependent variable must be a factor representing J + 1 intervals,
#' and the user must provide J thresholds.
#'
#' @param location an object of class "formula"
#' that specifies the equation for the conditional mean.
#' @param scale an object of class "formula"
#' that specifies the equation for the conditional standard deviation.
#' @inheritParams intreg
#' @return an object of class "gintreg"
#' @examples
#' gintreg(
#'   q2 ~ sex,
#'      ~ sex,
#'   data       = ias2009febmay,
#'   thresholds = c(-.5, .5, 1:5)
#' )
#' @export
gintreg <- function(location, scale, data, start, weights, thresholds) {
  cl <- match.call(expand.dots = FALSE)
  m  <- match(c("location", "data", "weights", "offset"), names(cl), 0)
  mf <- cl[c(1, m)]
  names(mf)[names(mf) == "location"] <- "formula"
  mf[[1]] <- quote(model.frame)
  mf      <- eval.parent(mf)
  vy      <- model.response(mf)
  if (!is.factor(vy)) stop("response must be a factor")
  mx       <- model.matrix(location, mf)
  vw       <- model.weights(mf)
  voffsetL <- model.offset(mf)
  cn       <- length(vy)
  if (is.null(voffsetL)) voffsetL <- rep(0, cn)
  if (missing(scale)) {
    mz           <- matrix(1, cn, 1)
    colnames(mz) <- "(Intercpet)"
    voffsetS     <- rep(0, cn)
  } else {
    m  <- match(c("scale", "data", "weights", "offset"), names(cl), 0)
    mf <- cl[c(1, m)]
    names(mf)[names(mf) == "scale"] <- "formula"
    mf[[1]]  <- quote(model.frame)
    mf       <- eval.parent(mf)
    mz       <- model.matrix(scale, mf)
    voffsetS <- model.offset(mf)
    if (is.null(voffsetS)) voffsetS <- rep(0, cn)
  }
  vthresh <- thresholds
  if (length(vthresh) != nlevels(vy) - 1) stop("incorrect number of thresholds")
  ck <- ncol(mx)
  cl <- ncol(mz)
  if (missing(start)) {
    vstart_beta  <- numeric(ck)
    vstart_delta <- c(log(max(vthresh) - min(vthresh)), numeric(cl - 1))
    vstart       <- c(vstart_beta, vstart_delta)
  } else {
    vstart <- start
  }
  lout            <- gintreg.fit(vy, mx, mz, vw, voffsetL, voffsetS, vthresh, vstart)
  vtheta          <- lout$par
  mhess           <- lout$hessian
  vbeta           <- vtheta[seq_len(ck)]
  vdelta          <- vtheta[ck + seq_len(cl)]
  names(vbeta)    <- colnames(mx)
  names(vdelta)   <- colnames(mz)
  vparname        <- c(names(vbeta), names(vdelta))
  dimnames(mhess) <- list(vparname, vparname)
  lfit            <- list(
    call      = match.call(),
    iteration = lout$counts["function"],
    coefL     = vbeta,
    coefS     = vdelta,
    logLik    = -lout$value,
    edf       = ck + cl,
    hess      = mhess
  )
  class(lfit) <- "gintreg"
  lfit
}
gintreg.fit <- function(vY, mX, mZ, vW, vOffsetL, vOffsetS, vThresh, vStart) {
  optim(
    vStart,
    function(vTheta) -gintreg.loglikelihood(vTheta, vY, mX, mZ, vW, vOffsetL, vOffsetS, vThresh),
    method  = "BFGS",
    hessian = TRUE
  )
}
gintreg.loglikelihood <- function(vTheta, vY, mX, mZ, vW, vOffsetL, vOffsetS, vThresh) {
  ck      <- ncol(mX)
  cl      <- ncol(mZ)
  vbeta   <- vTheta[seq_len(ck)]
  vdelta  <- vTheta[ck + seq_len(cl)]
  vmu     <- vOffsetL
  vlsigma <- vOffsetS
  if (ck > 0) vmu <- vmu + drop(mX %*% vbeta)
  if (cl > 0) vlsigma <- vlsigma + drop(mZ %*% vdelta)
  vsigma <- exp(vlsigma)
  vc     <- c(-Inf , vThresh, Inf)
  vy     <- unclass(vY)
  vp     <- pnorm((vc[vy + 1] - vmu) / vsigma) - pnorm((vc[vy] - vmu) / vsigma)
  if (is.null(vW)) sum(log(vp)) else sum(vW * log(vp))
}
#' @export
print.gintreg <- function(lFit, ...) {
  if (!is.null(lFit$call)) {
    cat("Call:\n")
    dput(lFit$call, control = NULL)
  }
  cat("\nlog-likelihood:", format(lFit$logLik, nsmall = 2))
  cat("\nAIC:", format(-2 * lFit$logLik + 2 * lFit$edf, nsmall = 2))
  cat("\n")
  if (length(lFit$coefL) > 0) {
    cat("\nLocation coefficients:\n")
    print(lFit$coefL, ...)
  } else {
    cat("\nNo location coefficients\n")
  }
  if (length(lFit$coefS) > 0) {
    cat("\nScale coefficients:\n")
    print(lFit$coefS, ...)
  } else {
    cat("\nNo scale coefficients\n")
  }
  invisible(lFit)
}
#' @export
summary.gintreg <- function(lFit, correlation = FALSE, ...) {
  vbeta    <- lFit$coefL
  vdelta   <- lFit$coefS
  mhess    <- lFit$hess
  vpar     <- c(vbeta, vdelta)
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
  lFit$summary  <- msummary
  if (correlation == TRUE) lFit$correlation <- mcov / outer(vse, vse)
  class(lFit) <- "summary.gintreg"
  lFit
}
#' @export
print.summary.gintreg <- function(lFit, digits = max(3, .Options$digits - 3), ...) {
  ck        <- length(lFit$coefL)
  cl        <- length(lFit$coefS)
  msummaryL <- lFit$summary[seq_len(ck), , drop = FALSE]
  msummaryS <- lFit$summary[ck + seq_len(cl), , drop = FALSE]
  if (!is.null(lFit$call)) {
    cat("Call:\n")
    dput(lFit$call, control = NULL)
  }
  cat("\nNo. of iterations:", lFit$iteration)
  cat("\nlog-likelihood:", format(lFit$logLik, nsmall = 2))
  cat("\nAIC:", format(-2 * lFit$logLik + 2 * lFit$edf, nsmall = 2))
  cat("\n")
  if (ck > 0) {
    cat("\nLocation coefficients:\n")
    printCoefmat(msummaryL, digits = digits, signif.stars = TRUE, signif.legend = FALSE, na.print = "NA", ...)
  } else {
    cat("\nNo location coefficients\n")
  }
  if (cl > 0) {
    cat("\nScale coefficients:\n")
    printCoefmat(msummaryS, digits = digits, signif.stars = TRUE, signif.legend = TRUE, na.print = "NA", ...)
  } else {
    cat("\nNo scale coefficients\n")
  }
  if (!is.null(lFit$correlation)) {
    cat("\nCorrelation of coefficients:\n")
    mcorr                    <- format(round(lFit$correlation, digits))
    mcorr[!lower.tri(mcorr)] <- ""
    print(mcorr[-1, -ncol(mcorr)], quote = FALSE, ...)
  }
  invisible(lFit)
}
