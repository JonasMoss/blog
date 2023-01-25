n = 1000

exp(lchoose(n, n/2)) / 2^n
exp(-log(0.5) * n)
2^n
2^n*(0.5)^(n)

n = 10 ** (1:60)
y = dbinom(n / 2, n, 0.5)
plot(n, y, log = "xy")
lm(I(log(y)) ~ log(n))

theta <- 1/10
n <- 100
x <- 0:n
plot(x, extraDistr::dbbinom(x, n, theta, theta))

theta <- 1
n = 10 ** (5:14)
y = extraDistr::dbbinom(n / 2, n, theta, theta + 7)
plot(n, y, log = "xy")
lines(n, exp(predict(lm(I(log(y)) ~ log(n) + log(2 + n)))))
lm(I(log(y)) ~ log(n))

exp(coef(lm(I(log(y)) ~ log(n)))[1])
sqrt(theta) * 2  / sqrt(pi)

## Using two
alpha <- 1
beta <- 1

calc3 = \(alpha, beta) {
  2^(2 - alpha - beta) / beta(alpha, beta)
}

calc4 = \(alpha, beta) {
  1/2 * dbeta(1/2, alpha, beta)
}

n = 10 ** (5:14)
y = extraDistr::dbbinom(n / 2, n, alpha, beta)
plot(n, y, log = "xy")
lines(n, exp(predict(lm(I(log(y)) ~ log(n) + log(2 + n)))))
lm(I(log(y)) ~ log(n))

exp(coef(lm(I(log(y)) ~ log(n)))[1])
calc3(alpha, beta)



sqrt(alpha) * 2  / sqrt(pi)





calc = \(alpha, beta) {
  pow = \(x) x^(x - 0.5)
  pow(alpha) * pow(beta) / pow(alpha + beta) / (
    0.5 ^ (alpha - 0.5) * 0.5 ^ (beta - 0.5) )
}

1 / calc(alpha, beta) / sqrt(pi) * sqrt(2)

calc2 = \(alpha, beta) {
  sqrt(2/pi) * sqrt(2*pi) * 2^(1 - alpha - beta) / beta(alpha, beta)
}
calc3 = \(alpha, beta) {
 2^(2 - alpha - beta) / beta(alpha, beta)
}
calc3(alpha, beta)

beta(0.1, 0.1 + 1)
beta(0.1, 0.1) 
2 * beta(0.1 + 1, 0.1)
24 * beta(0.1 + 1, 0.1 + 1)
beta(0.1, 0.1) / beta(0.1 + 2, 0.1 + 2)
beta(0.1, 0.1) / beta(0.1 + 3, 0.1 + 3)
beta(0.1, 0.1) / beta(0.1 + 4, 0.1 + 4)

n = 1:10
y = beta(theta + n, theta + n)
z = beta(n, n) * 0.5^(theta)
plot(n, y, log = "xy")
lines(n, z)


# Verification 1
f = \(m) sqrt(2*pi) * 2^(-2 * m + 1/2) * m^(-1/2) 
g = \(m) beta(m, m)/f(m)
g(100)

# Verification 2
f = \(m) sqrt(2*pi) * 2^(-2 * m + 1/2) * m^(-1/2) 
g = \(m, theta) beta(m + theta, m + theta) / f(m + theta)
g(300, 2)

# Verification 3
f1 = \(m, theta) 2^(-2*m) *(1 + m/theta)^(-1/2)
g1 = \(m, theta) beta(m + theta, m + theta) / beta(theta, theta)
h = \(m, theta) f(m + theta) / f(theta)
g1(400, 9) / f1(400, 9)





f = function(p, mean, lower) {
  (mean - p[1]/(p[1] + p[2]))^2 + (pbeta(lower, p[1], p[2]) - 0.1)^2
}

nlm(f, p = c(5, 5), mean = 0.516, lower = 0.49)
nlm(f, p = c(5, 5), mean = 0.52, lower = 0.51)

1/(dbeta(1/2, 81, 81) / 4)
1/(dbeta(1/2, 305, 285) / 4)
1/(dbeta(1/2, 2141.910, 1977.336) / 4)
