data_list <- list(
  linch = jsonlite::fromJSON("posts/value-estimation/raw_data/linch-zhang.json"),
  finn = jsonlite::fromJSON("posts/value-estimation/raw_data/finn-moorhouse.json"),
  gavin = jsonlite::fromJSON("posts/value-estimation/raw_data/gavin-leech-complete.json"),
  jamie = jsonlite::fromJSON("posts/value-estimation/raw_data/jaime-sevilla.json"),
  misha = jsonlite::fromJSON("posts/value-estimation/raw_data/misha-yagudin.json"),
  ozzie = jsonlite::fromJSON("posts/value-estimation/raw_data/ozzie-gooen.json")
)

source("posts/value-estimation/functions.R")

data <- make_frame(data_list[[3]])
#data$x <- 1 * (data$distance >= 0)
lm(distance ~ . - 1, data = data)

estimate <- function(data, type = c("sqrt", "linear"), tau) {
  type <- match.arg(type)
  k <- ncol(data) - 1
  n <- nrow(data)
  initial <- lm(distance ~ . - 1, data = data)
  coefs <- coef(initial)
  sigma <- summary(initial)$sigma
  
  y <- data$distance
  x <- as.matrix(data[, 2:ncol(data)])
  indices <- apply(x, 1, \(x) which(x != 0))
  
  f <- function(theta) {
    betas <- theta[seq(k)]
    sigmas <- theta[seq(k) + k]
    if(type == "linear") {
      sds <- sapply(indices, \(x) sum(sigmas[x]))      
    } else {
      sds <- sqrt(sapply(indices, \(x) sum(sigmas[x]^2)))  
    }
    mean <- x %*% betas
    -sum(dnorm(y, mean, sds, log = TRUE))
  }
  
  theta <- c(coefs, rep(sigma, k))
  ui <- cbind(matrix(0, k, k), diag(k))
  ci <- rep(tau, k)
  
  fit <- constrOptim(
      theta = theta, 
      f = f,
      grad = NULL,
      ui = ui,
      ci = ci,
      control = list(maxit = 10000))
  
  # hessian <- numDeriv::hessian(f, fit$par)
  # 
  # cbind(-diag(1/hessian) / sqrt(n) * 1.96 + fit$par,
  #       diag(1/hessian) / sqrt(n) * 1.96 + fit$par)
  
  list(
    fit = fit,
    beta = fit$par[seq(k)],
    sigma = stats::setNames(fit$par[seq(k) + k], names(fit$par[seq(k)])),
    likelihood = -fit$value,
    aic = 2 * fit$value + 2 * 2 * k)
  
}


# new <- lapply(data_list, \(x) {
#   x = make_frame(x)
#   x$distance = x$distance + runif(n = nrow(x), 0, 0.1)
#   estimate(x, type = "sqrt")$aic
#   })
# old <- lapply(data_list, \(x) AIC(lm(distance ~ . - 1, data = make_frame(x))))
# cbind(new, old)

name <- "jamie"
data = make_frame(data_list[[name]])
mod <- estimate(data, tau = 0.00)
plot(mod$beta, type = "h")
points(coef(lm(distance ~ . - 1, data = data)), col = "blue")
summary(lm(distance ~ . - 1, data = data))$sigma
mod$sigma
mod$aic
AIC(lm(distance ~ . - 1, data = data))
cor(mod$sigma, mod$beta)


### -------------------------------------------------------------------------->
###  Graph
### -------------------------------------------------------------------------->
misha <- data_list[[name]]
levels <- levels(as.factor(c(misha$source, misha$target)))
source <- as.numeric(factor(misha$source, levels = levels))
target <- as.numeric(factor(misha$target, levels = levels))
graph <- igraph::graph_from_edgelist(cbind(source, target))
plot(graph)


x <- as.matrix(data[, 2:ncol(data)])
y <- data$distance


data = make_frame(data_list$gavin)
mod <- lm(distance ~ . - 1, data = data)
resid(mod)
x <- as.matrix(data[, 2:ncol(data)])
betas <- coef(mod)
x %*% betas
plot(resid(mod)[x[, 11] != 0], ylim = c(-3, 3), col = "blue")


z = (x != 0)

# pairwise_model_logit <- function(data, fixed = 3, keep_names = FALSE) {
#   data_frame <- make_frame(data, fixed = fixed, keep_names = keep_names)
#   data_frame$distance <- data_frame$distance >= 0
#   glm(distance ~ . - 1, data = data_frame, family = binomial(link = "logit"))
# }
# 
# 
# 
# 
# plotter <- \(i) {
#   x <- coef(pairwise_model(data_list[[i]]))
#   y <- coef(pairwise_model_logit(data_list[[i]]))
#   
#   plot(x, y, xlab = "Scale kept", ylab = "Scale not kept",
#        sub = paste0("Correlation: ", round(cor(x, y), 3)),
#        main = names(data_list)[i])
#   abline(lm(y ~ x))  
# }
# 
# plotter(6)


data = make_frame(data_list$gavin)
plot(lm(distance ~ . - 1, data = data), which = 1)


#robustbase::lmrob(distance ~ . - 1, data = data)
c1 <- coef(robust::lmRob(distance ~ . - 1, data = data))
c2 <- coef(lm(distance ~ . - 1, data = data))
plot(c1, c2)
abline(a = 0, b = 1)
cor(c1, c2)



data = make_frame(data_list$gavin)
#robustbase::lmrob(distance ~ . - 1, data = data)
f1 <- coef(robust::lmRob(distance ~ . - 1, data = data))
f2 <- coef(lm(distance ~ . - 1, data = data))




new <- lapply(data_list, \(x) {
  x = make_frame(x)
  summary(robust::lmRob(distance ~ . - 1, data = x))$sigma
})
new2 <- lapply(data_list, \(x) {
  x = make_frame(x)
  summary(MASS::rlm(distance ~ . - 1, data = x))$sigma
})
old <- lapply(data_list, \(x) summary(lm(distance ~ . - 1, data = make_frame(x)))$sigma)
cbind(new, old, new2)
