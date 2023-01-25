n <- 15
x <- mtcars$drat[1:n]
y <- mtcars$wt[1:n]
n <- length(x)

mat <- cbind(1, x, 0)
minimum = Inf

for(k in seq(0:n)) {
  icomb <- arrangements::icombinations(n, k)
  choose(n, k)
  for(i in seq(choose(n, k))) {
    indices <- icomb$getnext()
    mat[, 3] <- 0
    mat[indices, 3] <- 1
    mod <- lm.fit(mat, y)
    current <- sum(mod$residuals^2)
    if (current < minimum) {
      minimum = current
      z = indices
      coef <- mod$coefficients
    }
  }  
}

cor2 <- sqrt(1 - minimum/sum((y - mean(y))^2)) * sign(coef[2])
cor1 <- cor(x, y)

col <- rep("blue", n)
col[z] <- "red"
plot(x, y, col = col, xlab = "Rear axle ratio", ylab = "Weight",
     pch = 20, main = "Two or one line?")
abline(a = coef[1], b = coef[2], col = "blue")
abline(a = coef[1] + coef[3], b = coef[2], col = "red")
abline(lm(y ~ x))
