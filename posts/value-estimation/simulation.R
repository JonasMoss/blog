p <- 0.5
beta <- rnorm(15)
sigma <- rexp(15)
beta <- exp(beta - beta[3])
sigma[3] <- 0
sigma_ <- sigma[c(1:2, 4:15)]
beta_ <- beta[c(1:2, 4:15)]

sim_two_sigma <- function(p, beta, sigma) {
  stopifnot(length(beta) == length(sigma))
  k <- length(beta)
  is_connected = FALSE
  while(!is_connected)  {
    g <- igraph::erdos.renyi.game(
      k, p, directed = TRUE
    )
    is_connected = igraph::is.connected(g) 
  }
  
  x_ <- as.matrix(igraph::as_data_frame(g))
  colnames(x_) <- c("source", "target")
  y <- apply(x_, 1, \(x) 
             exp(beta[x[2]] - beta[x[1]] + sigma[x[1]] * rnorm(1) + sigma[x[2]] * rnorm(1)))
  data <- make_frame(data.frame(distances = y, x_))
  estimate(data, tau = 0)
}

mod <- sim_two_sigma(p, beta, sigma)
plot(mod$beta, log(beta_))
