#' Fit a heterogeneous normal interval regression model to interval data
#'
#' \code{hintreg} fits a heterogeneous normal interval regression model
#' to interval data with an indifference limen,
#' i.e., a normal interval regression model with conditional heteroskedasticity
#' and heterogeneous indifference limen.
#' The dependent variable must be a factor representing J + 1 intervals,
#' and the user must provide thresholds below and above 0, respectively,
#' expect for an indifference limen (so there must be J - 2 known thresholds).
#'
#' @param lower an object of class "formula"
#' that specifies the equation for the lower bound of the indifference limen.
#' @param upper an object of class "formula"
#' that specifies the equation for the upper bound of the indifference limen.
#' @param threshbelow a vector of thresholds below 0
#' @param threshabove a vector of thresholds above 0
#' @inheritParams intreg
#' @inheritParams gintreg
#' @return an object of class "hintreg"
#' @export
hintreg <- function(location, scale, lower, upper, data, start, weights, threshbelow, threshabove) {
  cl <- match.call(expand.dots = FALSE)
  m  <- match(c("location", "data", "weights", "offset"), names(cl), 0)
  mf <- cl[c(1, m)]
  names(mf)[names(mf) == "location"] <- "formula"
  mf[[1]] <- quote(model.frame)
  mf      <- eval.parent(mf)
  vy      <- model.response(mf)
  if (!is.factor(vy)) stop("response is not a factor")
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
  if (missing(lower)) {
    mu           <- matrix(1, cn, 1)
    colnames(mu) <- "(Intercpet)"
    voffsetLB    <- rep(0, cn)
  } else {
    m  <- match(c("lower", "data", "weights", "offset"), names(cl), 0)
    mf <- cl[c(1, m)]
    names(mf)[names(mf) == "lower"] <- "formula"
    mf[[1]]   <- quote(model.frame)
    mf        <- eval.parent(mf)
    mu        <- model.matrix(lower, mf)
    voffsetLB <- model.offset(mf)
    if (is.null(voffsetLB)) voffsetLB <- rep(0, cn)
  }
  if (missing(upper)) {
    mv           <- matrix(1, cn, 1)
    colnames(mv) <- "(Intercpet)"
    voffsetUB    <- rep(0, cn)
  } else {
    m  <- match(c("upper", "data", "weights", "offset"), names(cl), 0)
    mf <- cl[c(1, m)]
    names(mf)[names(mf) == "upper"] <- "formula"
    mf[[1]]   <- quote(model.frame)
    mf        <- eval.parent(mf)
    mv        <- model.matrix(upper, mf)
    voffsetUB <- model.offset(mf)
    if (is.null(voffsetUB)) voffsetUB <- rep(0, cn)
  }
  vthresh <- c(threshbelow, NA, NA, threshabove)
  if (length(vthresh) != nlevels(vy) - 1) stop("incorrect number of thresholds")
  cbeta    <- ncol(mx)
  cdelta   <- ncol(mz)
  clambda  <- ncol(mu)
  cupsilon <- ncol(mv)
  if (missing(start)) {
    vstart_beta    <- numeric(cbeta)
    vstart_delta   <- c(log(max(vthresh, na.rm = TRUE) - min(vthresh, na.rm = TRUE)), numeric(cdelta - 1))
    vstart_lambda  <- numeric(clambda)
    vstart_upsilon <- numeric(cupsilon)
    vstart         <- c(vstart_beta, vstart_delta, vstart_lambda, vstart_upsilon)
  } else {
    vstart <- start
  }
  lout            <- hintreg.fit(vy, mx, mz, mu, mv, vw, voffsetL, voffsetS, voffsetLB, voffsetUB, vthresh, vstart)
  vtheta          <- lout$par
  mhess           <- lout$hessian
  vindex_beta     <- c(!logical(cbeta), logical(cdelta + clambda + cupsilon))
  vindex_delta    <- c(logical(cbeta), !logical(cdelta), logical(clambda + cupsilon))
  vindex_lambda   <- c(logical(cbeta + cdelta), !logical(clambda), logical(cupsilon))
  vindex_upsilon  <- c(logical(cbeta + cdelta + clambda), !logical(cupsilon))
  vbeta           <- vtheta[vindex_beta]
  vdelta          <- vtheta[vindex_delta]
  vlambda         <- vtheta[vindex_lambda]
  vupsilon        <- vtheta[vindex_upsilon]
  names(vbeta)    <- colnames(mx)
  names(vdelta)   <- colnames(mz)
  names(vlambda)  <- colnames(mu)
  names(vupsilon) <- colnames(mv)
  vparname        <- c(names(vbeta), names(vdelta), names(vlambda), names(vupsilon))
  dimnames(mhess) <- list(vparname, vparname)
  lfit            <- list(
    call      = match.call(),
    iteration = lout$counts["function"],
    coefL     = vbeta,
    coefS     = vdelta,
    coefLower = vlambda,
    coefUpper = vupsilon,
    logLik    = -lout$value,
    edf       = cbeta + cdelta + clambda + cupsilon,
    hess      = mhess
  )
  class(lfit) <- "hintreg"
  lfit
}
hintreg.fit <- function(vY, mX, mZ, mU, mV, vW, vOffsetL, vOffsetS, vOffsetLB, vOffsetUB, vThresh, vStart) {
  optim(
    vStart,
    function(.) -hintreg.loglikelihood(., vY, mX, mZ, mU, mV, vW, vOffsetL, vOffsetS, vOffsetLB, vOffsetUB, vThresh),
    method  = "BFGS",
    hessian = TRUE
  )
}
hintreg.loglikelihood <- function(vTheta, vY, mX, mZ, mU, mV, vW, vOffsetL, vOffsetS, vOffsetLB, vOffsetUB, vThresh) {
  cn             <- length(vY)
  cbeta          <- ncol(mX)
  cdelta         <- ncol(mZ)
  clambda        <- ncol(mU)
  cupsilon       <- ncol(mV)
  vindex_beta    <- c(!logical(cbeta), logical(cdelta + clambda + cupsilon))
  vindex_delta   <- c(logical(cbeta), !logical(cdelta), logical(clambda + cupsilon))
  vindex_lambda  <- c(logical(cbeta + cdelta), !logical(clambda), logical(cupsilon))
  vindex_upsilon <- c(logical(cbeta + cdelta + clambda), !logical(cupsilon))
  vbeta          <- vTheta[vindex_beta]
  vdelta         <- vTheta[vindex_delta]
  vlambda        <- vTheta[vindex_lambda]
  vupsilon       <- vTheta[vindex_upsilon]
  vmu            <- vOffsetL
  vlsigma        <- vOffsetS
  vlogit_l       <- vOffsetLB
  vlogit_u       <- vOffsetUB
  if (cbeta > 0) vmu <- vmu + drop(mX %*% vbeta)
  if (cdelta > 0) vlsigma <- vlsigma + drop(mZ %*% vdelta)
  if (clambda > 0) vlogit_l <- vlogit_l + drop(mU %*% vlambda)
  if (cupsilon > 0) vlogit_u <- vlogit_u + drop(mV %*% vupsilon)
  vsigma   <- exp(vlsigma)
  vlimen_l <- -plogis(vlogit_l)
  vlimen_u <- plogis(vlogit_u)
  vthresh  <- vThresh
  vy       <- unclass(vY)
  vlower   <- numeric(cn)
  vupper   <- numeric(cn)
  for (i in seq_len(cn)) {
    vthresh[is.na(vThresh)] <- c(vlimen_l[i], vlimen_u[i])
    vc                      <- c(-Inf, vthresh, Inf)
    vlower[i]               <- vc[vy[i]]
    vupper[i]               <- vc[vy[i] + 1]
  }
  vp <- pnorm((vupper - vmu) / vsigma) - pnorm((vlower - vmu) / vsigma)
  if (is.null(vW)) sum(log(vp)) else sum(vW * log(vp))
}
#' @export
print.hintreg <- function(lFit, ...) {
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
  if (length(lFit$coefLower) > 0) {
    cat("\nLower bound coefficients:\n")
    print(lFit$coefLower, ...)
  } else {
    cat("\nNo lower bound coefficients\n")
  }
  if (length(lFit$coefUpper) > 0) {
    cat("\nUpper bound coefficients:\n")
    print(lFit$coefUpper, ...)
  } else {
    cat("\nNo upper bound coefficients\n")
  }
  invisible(lFit)
}
#' @export
summary.hintreg <- function(lFit, correlation = FALSE, ...) {
  vbeta    <- lFit$coefL
  vdelta   <- lFit$coefS
  vlambda  <- lFit$coefLower
  vupsilon <- lFit$coefUpper
  mhess    <- lFit$hess
  vpar     <- c(vbeta, vdelta, vlambda, vupsilon)
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
  class(lFit) <- "summary.hintreg"
  lFit
}
#' @export
print.summary.hintreg <- function(lFit, digits = max(3, .Options$digits - 3), ...) {
  cbeta          <- length(lFit$coefL)
  cdelta         <- length(lFit$coefS)
  clambda        <- length(lFit$coefLower)
  cupsilon       <- length(lFit$coefUpper)
  vindex_beta    <- c(!logical(cbeta), logical(cdelta + clambda + cupsilon))
  vindex_delta   <- c(logical(cbeta), !logical(cdelta), logical(clambda + cupsilon))
  vindex_lambda  <- c(logical(cbeta + cdelta), !logical(clambda), logical(cupsilon))
  vindex_upsilon <- c(logical(cbeta + cdelta + clambda), !logical(cupsilon))
  msummaryL      <- lFit$summary[vindex_beta, , drop = FALSE]
  msummaryS      <- lFit$summary[vindex_delta, , drop = FALSE]
  msummaryLB     <- lFit$summary[vindex_lambda, , drop = FALSE]
  msummaryUB     <- lFit$summary[vindex_upsilon, , drop = FALSE]
  if (!is.null(lFit$call)) {
    cat("Call:\n")
    dput(lFit$call, control = NULL)
  }
  cat("\nNo. of iterations:", lFit$iteration)
  cat("\nlog-likelihood:", format(lFit$logLik, nsmall = 2))
  cat("\nAIC:", format(-2 * lFit$logLik + 2 * lFit$edf, nsmall = 2))
  cat("\n")
  if (cbeta > 0) {
    cat("\nLocation coefficients:\n")
    printCoefmat(msummaryL, digits = digits, signif.stars = TRUE, signif.legend = FALSE, na.print = "NA", ...)
  } else {
    cat("\nNo location coefficients\n")
  }
  if (cdelta > 0) {
    cat("\nScale coefficients:\n")
    printCoefmat(msummaryS, digits = digits, signif.stars = TRUE, signif.legend = FALSE, na.print = "NA", ...)
  } else {
    cat("\nNo scale coefficients\n")
  }
  if (clambda > 0) {
    cat("\nLower bound coefficients:\n")
    printCoefmat(msummaryLB, digits = digits, signif.stars = TRUE, signif.legend = FALSE, na.print = "NA", ...)
  } else {
    cat("\nNo lower bound coefficients\n")
  }
  if (cupsilon > 0) {
    cat("\nUpper bound coefficients:\n")
    printCoefmat(msummaryUB, digits = digits, signif.stars = TRUE, signif.legend = TRUE, na.print = "NA", ...)
  } else {
    cat("\nNo upper bound coefficients\n")
  }
  if (!is.null(lFit$correlation)) {
    cat("\nCorrelation of coefficients:\n")
    mcorr                    <- format(round(lFit$correlation, digits))
    mcorr[!lower.tri(mcorr)] <- ""
    print(mcorr[-1, -ncol(mcorr)], quote = FALSE, ...)
  }
  invisible(lFit)
}
