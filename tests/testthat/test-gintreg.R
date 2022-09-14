# Set parameters

cN       <- 100
cL       <- 2
dBeta_0  <- 0
vBeta_1  <- rep(1, cL - 1)
dBeta_2  <- 2
vBeta    <- c(dBeta_0, vBeta_1, dBeta_2)
dDelta_0 <- 0
vDelta_1 <- rep(1, cL - 1)
dDelta_2 <- 2
vDelta   <- c(dDelta_0, vDelta_1, dDelta_2)
cCat     <- 8
vC_l     <- NULL
vC_u     <- 1:5
#cCat     <- 18
#vC_l     <- -5:-1
#vC_u     <- 1:10
vC_0     <- c(-.5, .5)
vC       <- c(vC_l, vC_0, vC_u)

# Data

vD      <- rbinom(cN, cL - 1, 1 / cL) |> factor()
mD      <- table(1:cN, vD)[, 2:cL, drop = FALSE]
vX      <- rnorm(cN)
vZ      <- rnorm(cN)
vY_star <- cbind(1, mD, vX) %*% vBeta + exp(cbind(1, mD, vX) %*% vDelta) * vZ
vY      <- cut(vY_star, c(-Inf, vC, Inf))
mData   <- data.frame(y = vY, d = vD, x = vX)

# Basic tests

logLik_oglmx <- function(location, scale) {
  oglmx::oglmx(
    location,
    scale,
    data        = mData,
    threshparam = vC
  )$loglikelihood[1] |> suppressWarnings()
}
logLik_gintreg <- function(location, scale) {
  gintreg(
    location,
    scale,
    data       = mData,
    thresholds = vC
  )$logLik
}
test_that("compare loglikelihoods of oglmx and gintreg", {
  expect_equal(logLik_oglmx(y ~ d + x, ~ d + x)        , logLik_gintreg(y ~ d + x, ~ d + x))
  expect_equal(logLik_oglmx(y ~ 0, ~ 1)                , logLik_gintreg(y ~ 0, ~ 1))
  expect_equal(logLik_oglmx(y ~ 1, ~ 1)                , logLik_gintreg(y ~ 1, ~ 1))
  expect_equal(logLik_oglmx(y ~ 0 + d + x, ~ 0 + d + x), logLik_gintreg(y ~ 0 + d + x, ~ 0 + d + x))
})

# Offset

oglmx::oglmx(
  y ~ d + x,
    ~ d + x,
  data        = mData,
  beta        = c(NA, NA, 2),
  delta       = c(NA, NA, 2),
  threshparam = vC
) |> summary()
gintreg(
  y ~ d + offset(2 * x),
    ~ d + offset(2 * x),
  data       = mData,
  thresholds = vC
) |> summary()

# Use start value

lout <- intreg(
  y ~ d + x,
  data       = mData,
  thresholds = vC
)
vbeta   <- coef(lout)
dlsigma <- lout$logscale
#vstart  <- c(vbeta, dlsigma, numeric(length(vDelta) - 1))
vstart <- c(vBeta, vDelta)
oglmx::oglmx(
  y ~ d + x,
    ~ d + x,
  data        = mData,
  start       = vstart,
  threshparam = vC
) |> summary()
gintreg(
  y ~ d + x,
    ~ d + x,
  data       = mData,
  start       = vstart,
  thresholds = vC
) |> summary()

# Use weights

vW <- rexp(cN)
oglmx::oglmx(
  y ~ d + x,
    ~ d + x,
  data        = mData,
  weights     = vW,
  threshparam = vC
) |> summary()
lout <- gintreg(
  y ~ d + x,
    ~ d + x,
  data       = mData,
  weights    = vW,
  thresholds = vC
)
summary(lout)
tidy(lout, TRUE)
