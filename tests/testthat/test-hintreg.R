# Set parameters

cN         <- 1000
cL         <- 2
cCat       <- 8
vC_l       <- NULL
vC_u       <- 1:5
#cCat       <- 18
#vC_l       <- -5:-1
#vC_u       <- 1:10
dLimenlb   <- -1
dLimenub   <- 1
dBeta_0    <- 0
vBeta_1    <- rep(1, cL - 1)
dBeta_2    <- 2
vBeta      <- c(dBeta_0, vBeta_1, dBeta_2)
dDelta_0   <- 0
vDelta_1   <- rep(1, cL - 1)
dDelta_2   <- 2
vDelta     <- c(dDelta_0, vDelta_1, dDelta_2)
dLambda_0  <- 0
vLambda_1  <- rep(1, cL - 1)
dLambda_2  <- 2
vLambda    <- c(dLambda_0, vLambda_1, dLambda_2)
dUpsilon_0 <- 0
vUpsilon_1 <- rep(1, cL - 1)
dUpsilon_2 <- 2
vUpsilon   <- c(dUpsilon_0, vUpsilon_1, dUpsilon_2)

# Data

vD      <- rbinom(cN, cL - 1, 1 / cL) |> factor()
mD      <- table(1:cN, vD)[, 2:cL, drop = FALSE]
vX      <- rnorm(cN)
vZ      <- rnorm(cN)
vY_star <- cbind(1, mD, vX) %*% vBeta + exp(cbind(1, mD, vX) %*% vDelta) * vZ
vLower  <- -plogis(cbind(1, mD, vX) %*% vLambda)
vUpper  <- plogis(cbind(1, mD, vX) %*% vUpsilon)
vY      <- vapply(
  1:cN,
  function(i) findInterval(vY_star[i], c(-Inf, vC_l, vLower[i], vUpper[i], vC_u, Inf)),
  integer(1)
) |> factor(levels = 1:cCat)
mData <- data.frame(y = vY, d = vD, x = vX)

# Basic tests

logLik_oglmx <- function(location, scale) {
  oglmx::oglmx(
    location,
    scale,
    data        = mData,
    threshparam = c(vC_l, NA, NA, vC_u)
  )$loglikelihood[1] |> suppressWarnings()
}
logLik_hintreg <- function(location, scale) {
  hintreg(
    location,
    scale,
    data        = mData,
    threshbelow = vC_l,
    threshabove = vC_u,
    limenlb     = dLimenlb,
    limenub     = dLimenub
  )$logLik
}
test_that("compare loglikelihoods of oglmx and hintreg", {
  expect_equal(logLik_oglmx(y ~ d + x, ~ d + x)        , logLik_hintreg(y ~ d + x, ~ d + x))
  expect_equal(logLik_oglmx(y ~ 0, ~ 1)                , logLik_hintreg(y ~ 0, ~ 1))
  expect_equal(logLik_oglmx(y ~ 1, ~ 1)                , logLik_hintreg(y ~ 1, ~ 1))
  expect_equal(logLik_oglmx(y ~ 0 + d + x, ~ 0 + d + x), logLik_hintreg(y ~ 0 + d + x, ~ 0 + d + x))
})

# Offset

oglmx::oglmx(
  y ~ d + x,
  ~ d + x,
  data        = mData,
  beta        = c(NA, NA, 2),
  delta       = c(NA, NA, 2),
  threshparam = c(vC_l, NA, NA, vC_u)
) |> summary()
hintreg(
  y ~ d + offset(2 * x),
  ~ d + offset(2 * x),
  data        = mData,
  threshbelow = vC_l,
  threshabove = vC_u,
  limenlb     = dLimenlb,
  limenub     = dLimenub
) |> summary()

# Use start value

lout <- gintreg(
  y ~ d + x,
  ~ d + x,
  data       = mData,
  thresholds = c(vC_l, -.5, .5, vC_u)
)
vbeta        <- lout$coefL
vdelta       <- lout$coefS
vstart       <- c(vbeta, vdelta, -.5, .5)
#vstart       <- c(vBeta, vDelta, -.5, .5)
oglmx::oglmx(
  y ~ d + x,
  ~ d + x,
  data        = mData,
  start       = vstart,
  threshparam = c(vC_l, NA, NA, vC_u)
) |> summary()
hintreg(
  y ~ d + x,
  ~ d + x,
  data        = mData,
  start       = vstart,
  threshbelow = vC_l,
  threshabove = vC_u,
  limenlb     = dLimenlb,
  limenub     = dLimenub
) |> summary()

# Use weights

vW <- rexp(cN)
oglmx::oglmx(
  y ~ d + x,
  ~ d + x,
  data        = mData,
  weights     = vW,
  threshparam = c(vC_l, NA, NA, vC_u)
) |> summary()
hintreg(
  y ~ d + x,
  ~ d + x,
  data        = mData,
  weights     = vW,
  threshbelow = vC_l,
  threshabove = vC_u,
  limenlb     = dLimenlb,
  limenub     = dLimenub
) |> summary()

# Heterogenenous interval regression

lout <- hintreg(
  y ~ d + x,
  ~ d + x,
  ~ d + x,
  ~ d + x,
  data        = mData,
  threshbelow = vC_l,
  threshabove = vC_u,
  limenlb     = dLimenlb,
  limenub     = dLimenub
)
summary(lout)
tidy(lout, TRUE)
