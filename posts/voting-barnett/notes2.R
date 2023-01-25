president <- readr::read_csv("posts/voting-barnett/1976-2020-senate.csv")
senate <- readr::read_csv("posts/voting-barnett/1976-2020-president.csv")

trans <- \(dat) {
  trans <- \(x) {
    groups <- dplyr::group_by(x, year, state)
    dplyr::filter(
      groups, 
      candidatevotes == max(candidatevotes)
    )    
  }

  democrat <- trans(dplyr::filter(dat, party_simplified == "DEMOCRAT"))
  republican <- trans(dplyr::filter(dat, party_simplified == "REPUBLICAN"))
  
  list(
    ns = democrat$candidatevotes + republican$candidatevotes, 
    votes = democrat$candidatevotes)
}

ns <- trans(president)$ns
votes <- trans(president)$votes
ps <- votes / ns

ns <- trans(senate)$ns
votes <- trans(senate)$votes
ps <- votes / ns

### Correlation.

plot(ns, ps, log = "x")






# ns <- c(trans(senate)$ns, trans(senate)$ns)
# votes <- c(trans(senate)$votes, trans(senate)$ns - trans(senate)$votes)
# ps <- votes / ns

ml <- function(x, n) {
  f <- \(p) {
    -sum(extraDistr::dbbinom(x, n, p[1], p[2], log = TRUE))
  }
  nlm(f, c(1, 1))$estimate
}



i <- ns < 10**6
y <- dbinom(floor(ns / 2), size = ns, p = ps, log = TRUE)
plot(ns, y)
abline(a = 0, b = -1)
summary(lm(y ~ ns))
par = univariateML::mlbeta(ps)
0.5 * (digamma(par[1]) + digamma(par[2]) - 2 * digamma(par[1] + par[2])) + log(2)

summary(lm(dbinom(floor(ns / 2), size = ns, p = ps) ~ I(1/ns)))

n <- length(ns)
params <- univariateML::mlbeta(ps)
ps_true <- univariateML::rml(n, univariateML::mlbeta(ps))
votes <- rbinom(n, ns, ps_true)
params <- univariateML::mlbeta(ps)
votes <- extraDistr::rbbinom(n, ns, params[1], params[2])
ps <- votes / ns

points(ns, dbinom(floor(ns / 2), size = ns, p = ps, log = TRUE), col = "red")
summary(lm(I(dbinom(floor(ns / 2), size = ns, p = ps, log = TRUE)) ~ ns))


N <- 10**6
ns <- sample(ns, N, replace = TRUE)
ps <- rbeta(N, 50, 50)
y <- dbinom(floor(ns / 2), size = ns, p = ps, log = TRUE)
#plot(log(ns), log(-y), col = "red")
summary(lm(y ~ ns + log(ns)))
abline(a = 0, b = -1)
digamma(50) - digamma(100) + log(2)

log(2 / pi) * 0.5
0.5
