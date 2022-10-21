#' Fits a density to quantiles.
#'
#' The method is based on monotone splines. We minimize the squared distance
#'    between the spline approximation and the identity function. Thus you
#'    minimize the distance between the CDF and the uniform CDF in a
#'    q-transformed space. The parameters `lambda` and `mu` are to tweak the
#'    resulting density into the desired shape.
#'
#' The function returns
#' @param x,y Vector of probabilities (`x`) and quantiles (`y`). The quantiles
#'    must have been transformed to the unit interval using your desired `q`.
#' @param m Number of internal knots in the spline.
#' @param degree Degree of the spline function. Defaults to `3`, which
#'   corresponds to cubic splines.
#' @param lambda The penalty term for the squared second derivative.
#' @param mu The penalty term for the squared derivative.
#' @return A list containing the fitted spline together with density functions,
#'   and so on.

fitter <- function(x, y, m = 20, degree = 3, lambda = 0.5, mu = 0) {
  knots <- (2:m) / (m + 1)
  boundary_knots <- c(y[1], y[length(y)])
  
  sf <- \(x, derivs = 0) splines2::bSpline(
    x,
    intercept = TRUE,
    knots = knots,
    degree = degree,
    Boundary.knots = boundary_knots,
    derivs = derivs
  )
  
  matrices <- get_penalties(knots, boundary_knots, degree)
  x_mat <- sf(x)
  f <- Vectorize(\(x, i) x * sf(x)[i])
  k <- ncol(x_mat)
  sx <- sapply(seq(k), \(i) integrate(f, lower = 0, upper = 1, i = i)$value)
  constraints <- increasing(ncol(x_mat) - 1)
  
  fit <- quadprog::solve.QP(
    Dmat = matrices$s2_mat + lambda * matrices$p_mat + mu * matrices$p1_mat,
    dvec = sx,
    Amat = t(rbind(x_mat, constraints$amat)),
    bvec = c(y, constraints$bvec),
    meq = length(y)
  )
  
  coefs <- fit$solution
  pdf_untrans = \(w) sf(w, derivs = 1) %*% coefs
  cdf_untrans = \(w) sf(w) %*% coefs
  
  list(
    x = x,
    y = y,
    m = m,
    degree = degree,
    lambda = lambda,
    mu = mu,
    fit = fit,
    pdf = \(w, d = dunif, p = punif) d(w) * pdf_untrans(p(w)),
    cdf = \(w, p = punif) cdf_untrans(p(w))
  )
}

#' Get penalty matrices for the fitter.
#' @param knots The internal knots.
#' @param boundary_knots The boundary knots.
#' @param degree Degree of the B-spline.
#' @return List of penalty matrices.
get_penalties <- function(knots, boundary_knots, degree) {
  basis_obj <- fda::create.bspline.basis(
    rangeval = boundary_knots,
    norder = degree + 1,
    breaks = c(boundary_knots[1], knots, boundary_knots[2])
  )
  inds <- seq(basis_obj$nbasis)
  list(
    s2_mat = fda::bsplinepen(basis_obj, Lfdobj = 0)[inds, inds],
    p_mat = fda::bsplinepen(basis_obj, Lfdobj = 2)[inds, inds],
    p1_mat = fda::bsplinepen(basis_obj, Lfdobj = 1)[inds, inds]
  )
}

#' Increasing constraints matrix and vector
#'
#' @keywords internal
#' @param m The number of knots in the spline.
#' @param intercept If `TRUE`, the model includes an intercept.
#' @return The increasing constraint matrix.

increasing <- function(m, intercept = TRUE) {
  amat <- cbind(0, -cbind(diag(m - 1), 0) + cbind(0, diag(m - 1)))
  amat <- rbind(0, amat)
  amat <- rbind(amat, 0)
  amat[m + 1, c(1, m + 1)] <- -1
  amat[1, 2] <- 1
  amat <- rbind(0, amat)
  amat[1, 1] <- 1
  
  if (!intercept) {
    amat <- amat[2:nrow(amat), 2:ncol(amat)]
  }
  
  bvec <- rep(0, nrow(amat))
  bvec[nrow(amat)] <- -1
  
  list(
    amat = amat,
    bvec = bvec
  )
}
