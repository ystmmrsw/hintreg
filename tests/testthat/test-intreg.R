# Set parameters

cN      <- 100
cL      <- 2
dBeta_0 <- 0
vBeta_1 <- rep(1, cL - 1)
dBeta_2 <- 2
vBeta   <- c(dBeta_0, vBeta_1, dBeta_2)
dSigma  <- 1
cCat    <- 8
vC_l    <- NULL
vC_u    <- 1:5
#cCat    <- 18
#vC_l    <- -5:-1
#vC_u    <- 1:10
vC_0    <- c(-.5, .5)
vC      <- c(vC_l, vC_0, vC_u)

# Data

vD      <- rbinom(cN, cL - 1, 1 / cL) |> factor()
mD      <- table(1:cN, vD)[, 2:cL, drop = FALSE]
vX      <- rnorm(cN)
vU      <- rnorm(cN, 0, dSigma)
vY_star <- cbind(1, mD, vX) %*% vBeta + vU
vY      <- cut(vY_star, c(-Inf, vC, Inf))
mData   <- data.frame(y = vY, d = vD, x = vX)

# Basic tests

logLik_oglmx <- function(formula) {
  oglmx::oglmx(
    formula,
    data        = mData,
    threshparam = vC
  )$loglikelihood[1] |> suppressWarnings()
}
logLik_intreg <- function(formula) {
  intreg(
    formula,
    data       = mData,
    thresholds = vC
  )$logLik
}
test_that("compare loglikelihoods of oglmx and intreg", {
  expect_equal(logLik_oglmx(y ~ d + x)    , logLik_intreg(y ~ d + x))
  expect_equal(logLik_oglmx(y ~ 0)        , logLik_intreg(y ~ 0))
  expect_equal(logLik_oglmx(y ~ 1)        , logLik_intreg(y ~ 1))
  expect_equal(logLik_oglmx(y ~ 0 + d + x), logLik_intreg(y ~ 0 + d + x))
})

# Offset

oglmx::oglmx(
  y ~ d + x,
  data        = mData,
  beta        = c(NA, NA ,2),
  threshparam = vC
) |> summary()
intreg(
  y ~ d + offset(2 * x),
  data       = mData,
  thresholds = vC
) |> summary()

# Use start value

vstart <- c(vBeta, log(dSigma))
oglmx::oglmx(
  y ~ d + x,
  data        = mData,
  start       = vstart,
  threshparam = vC
) |> summary()
intreg(
  y ~ d + x,
  data       = mData,
  start      = vstart,
  thresholds = vC
) |> summary()

# Use weights

vW <- rexp(cN)
oglmx::oglmx(
  y ~ d + x,
  data        = mData,
  weights     = vW,
  threshparam = vC
) |> summary()
intreg(
  y ~ d + x,
  data       = mData,
  weights    = vW,
  thresholds = vC
) |> summary()

# Comparison of alternative methods

list("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN") |>
  lapply(function(.) intreg(
    y ~ d + x,
    data       = mData,
    start      = vstart,
    thresholds = vC,
    method = .
  ) |> summary())
intreg(
  y ~ d + x,
  data       = mData,
  thresholds = vC
) |> broom::tidy()
