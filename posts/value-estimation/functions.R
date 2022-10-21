#' Makes a suitable data drame out of pairwise agreement data.
#'
#' @param data Pairwise agreement data.
#' @param fixed Which question to fixed at 1.
#' @param keep_names If `TRUE`, keeps the names of the questions. This can
#'   make the data frame unwieldly. Defaults to `FALSE`.
#' @return An appropriate data frame.
make_frame <- function(data, fixed = 3, keep_names = FALSE) {
  levels <- levels(as.factor(c(data$source, data$target)))
  source <- as.numeric(factor(data$source, levels = levels))
  target <- as.numeric(factor(data$target, levels = levels))

  k <- length(levels)
  n <- nrow(data)
  
  d <- matrix(data = 0, nrow = n, ncol = k)
  for (i in seq(n)) {
    d[i, source[i]] <- -1
    d[i, target[i]] <- 1
  }

  selection <- setdiff(seq(k), fixed)
  y <- log(as.numeric(data$distance))
  data_frame <- data.frame(y, d[, selection])
  
  colnames(data_frame) <- if (keep_names) {
    c("distance", levels[selection])
  } else {
    c("distance", paste0("q", seq(k)[selection]))
  }
  
  data_frame
}

#' Estimate a pairwise model using linear regression
#'
#' @param data Pairwise agreement data.
#' @param fixed Which question to fixed at 1.
#' @param keep_names If `TRUE`, keeps the names of the questions. This can
#'   make the data frame unwieldly. Defaults to `FALSE`.
#' @return An `lm` object.
pairwise_model <- function(data, fixed = 3, keep_names = FALSE) {
  data_frame <- make_frame(data, fixed = fixed, keep_names = keep_names)
  lm(distance ~ . - 1, data = data_frame)
}

#' Estimate a pairwise mixed model using `lme4::lmer`
#'
#' @param data Pairwise agreement data.
#' @param fixed Which question to fixed at 1.
#' @param uncorrelated If `TRUE`, specifies multiple uncorrelated random 
#'   effects for the same grouping variable. If `FALSE`, does not assume
#'   uncorrelated effects.
#' @return An `lme4::lmer` object.
pairwise_mixed_model <- function(data, fixed = 3, uncorrelated = TRUE) {
  names <- names(data_list)
  model_list <- lapply(seq_along(data_list), \(i) {
    frame <- make_frame(data_list[[i]], fixed)
    cbind(frame, subject = names[i])
  })
  big_frame <- do.call(rbind, model_list)
  bar <- if (uncorrelated) " || " else " | "
  qs <- names(big_frame)[seq(ncol(big_frame) - 2) + 1]
  fixed_effects <- paste0(paste0(qs, collapse = " + "), " + 0")
  random_effects <- paste0("(", fixed_effects, bar, "subject)")
  formula <- paste0("distance ~ ", fixed_effects, " + ", random_effects)
  mod <- lme4::lmer(formula, data = big_frame, REML = TRUE)
}

#' Estimates a choice model with both binary and graded responses.
#' 
#' The parameters are identified if there is no perfectly separated point and
#'    the question graph is connected.
#'   
#' @data 
double_est <- function(data, maxit = 1000) {
  k <- ncol(data) - 2
  n <- nrow(data)
  
  indices <- !is.na(data$distance)
  y <- data$distance[indices]
  z <- data$binary[!indices]
  x <- as.matrix(data[, 3:ncol(data)])
  x_dist <- x[indices, ]
  x_bin <- x[!indices, ] * (2*z - 1)
  theta <- c(rep(0, k), 1)
  ui <- cbind(matrix(0, 1, k), 1)
  ci <- 0
  
  # fit <- constrOptim(
  #   theta = theta, 
  #   f = f,
  #   grad = \(x) numDeriv::grad(f, x),
  #   ui = ui,
  #   ci = ci,
  #   control = list(maxit = maxit))
  
  f <- function(theta) {
    betas <- theta[seq(k)]
    sigma <- theta[k + 1]
    
    mean <- x_dist %*% betas
    dist_contrib <- sum(dnorm(y / sigma, x_dist %*% betas, log = TRUE))
    bin_contrib <- sum(pnorm(x_bin %*% betas, log.p = TRUE))
    -(dist_contrib + bin_contrib)
  }
  
  # fit <- optim(
  #   par = theta, 
  #   fn = f,
  #   control = list(maxit = maxit),
  #   hessian = TRUE)
  
  fit <- nlm(
    p = theta, 
    f = f,
    iterlim = maxit,
    hessian = TRUE)
  
  names <- colnames(data[, 3:ncol(data)])
  beta <- stats::setNames(fit$estimate[seq(k)], names)
  sigma <- unname(fit$estimate[k + 1])

  list(
    beta_star = beta,
    beta = beta * sigma,
    sigma = sigma,
    likelihood = -fit$minimum,
    aic = 2 * fit$minimum + 2 * (k + 1),
    fit = fit,
    n = n,
    k = k
  )
}

cis <- function(x) {
  fit <- x$fit
  beta <- x$beta
  beta_star <- x$beta_star
  sigma <- x$beta
  n <- x$n
  k <- x$k

  sds <- diag(solve(x$fit$hessian))
  if(any(sds <= 0)) warning("Hessian contains non-positive elements.")
  sds <- pmax(sds, 0)
  ses <- sqrt(head(sds, -1) * sigma^2 + sds[k+1] * beta_star^2) / sqrt(n)
  
  matrix(rep(beta, 2), ncol = 2) +
    matrix(c(ses, ses), ncol = 2) *
    matrix(c(qt(0.025, n - k), qt(0.975, n - k)), ncol = 2, nrow = k, byrow = TRUE)
  
}
concordance = function(x) {
  n = nrow(x)
  r = ncol(x)
  sigma = cov(x) * (n - 1) / n
  mu = colMeans(x)
  trace = sum(diag(sigma))
  top = sum(sigma) - trace
  bottom = (r - 1) * trace + r ^ 2 * (mean(mu^2) - mean(mu)^2)
  top / bottom
}

