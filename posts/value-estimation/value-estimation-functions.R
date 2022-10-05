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
