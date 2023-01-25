jamie <- data_list$gavin
data <- make_frame(data_list$gavin)
lm(distance ~ . - 1, data = data)

levels <- levels(as.factor(c(gavin$source, gavin$target)))
source <- as.numeric(factor(gavin$source, levels = levels))
target <- as.numeric(factor(gavin$target, levels = levels))
graph <- igraph::graph_from_edgelist(cbind(source, target))
plot(graph)


k <- length(levels)
n <- nrow(data)

d <- matrix(data = 0, nrow = n, ncol = k)
for (i in seq(n)) {
  d[i, source[i]] <- -1
  d[i, target[i]] <- 1
}

y <- log(as.numeric(data$distance))

X = solve(t(x) %*% x)
D = MASS::ginv(t(d) %*% d)
X[1, 1]
D[1, 1] + D[3, 3] - 2*D[1, 3]

vol = sum(diag(t(d) %*% d))
sapply(1:nrow(D), \(i) D[i, i] + D[3, 3] - 2*D[i, 3]) * vol


k = ncol(d)
a <- MASS::ginv(d)
corr <- (diag(nrow(a)) - a %*% d)
a %*% y

con = (diag(nrow(a)) - a %*% d)[1, 1] 
con = (1 - (a %*% y)[3])

a %*% y - rep(1, k) * (a %*% y)[3]


ej = \(k, i) {
  out = rep(0, k)
  out[i] = 1
  out
}

(diag(k) - rep(1, k) %*% t(ej(k, 3))) %*% a %*% y


### File 2

data_list <- list(
  linch = jsonlite::fromJSON("posts/value-estimation/raw_data/linch-zhang.json"),
  finn = jsonlite::fromJSON("posts/value-estimation/raw_data/finn-moorhouse.json"),
  gavin = jsonlite::fromJSON("posts/value-estimation/raw_data/gavin-leech-complete.json"),
  jamie = jsonlite::fromJSON("posts/value-estimation/raw_data/jaime-sevilla.json"),
  misha = jsonlite::fromJSON("posts/value-estimation/raw_data/misha-yagudin.json"),
  ozzie = jsonlite::fromJSON("posts/value-estimation/raw_data/ozzie-gooen.json")
)
source("posts/value-estimation/functions.R")

gavin <- data_list$gavin
data <- make_frame(data_list$gavin)
levels <- levels(as.factor(c(gavin$source, gavin$target)))
source <- as.numeric(factor(gavin$source, levels = levels))
target <- as.numeric(factor(gavin$target, levels = levels))
k <- length(levels)
n <- nrow(data)
d <- matrix(data = 0, nrow = n, ncol = k)

for (i in seq(n)) {
  d[i, source[i]] <- -1
  d[i, target[i]] <- 1
}

y <- as.numeric(data$distance)
d_inv = MASS::ginv(d)

d_inv %*% y - rep(1, k) * (d_inv %*% y)[3]
coef(lm(distance ~ . - 1, data = data))


#diag(k) - d_inv %*% d


ej = \(k, i) {
  out = rep(0, k)
  out[i] = 1
  out
}

(diag(k) - rep(1, k) %*% t(ej(k, 3))) %*% d_inv %*% y


diag(k) - d_inv %*% d

c_mat = \(d) {
  ddd = MASS::ginv(t(d) %*% d)
  out = ddd * 0
  for (i in seq(nrow(ddd))) {
    for (j in seq(nrow(ddd))) {
      out[i, j] = ddd[i, i] + ddd[j, j] - 2*ddd[i, j]  
    }
  }
  out
}

### File 3

# Simulate.

m = 8
p = 20
k = 3
d = matrix(0, nrow = 3 * p + 5, ncol = 3 * m)

for (offset in c(0, 1, 2)) {
  for (i in (seq(p) + p * offset))  
    d[i, offset * m + sample(m, 2)] = c(-1, 1)
}

d[p * 3 + 1, c(m, m+1)] = c(-1, 1)
d[p * 3 + 2, c(m+1, 2 * m+1)] = c(-1, 1)
#d[p * 3 + 3, c(m+3, 2 * m+1)] = c(-1, 1)
#d[p * 3 + 4, c(m+4, 2 * m+2)] = c(-1, 1)
#d[p * 3 + 5, c(m+5, 2 * m+3)] = c(-1, 1)



c_mat = \(d) {
  ddd = MASS::ginv(t(d) %*% d)
  out = ddd * 0
  for (i in seq(nrow(ddd))) {
    for (j in seq(nrow(ddd))) {
      out[i, j] = ddd[i, i] + ddd[j, j] - 2*ddd[i, j]  
    }
  }
  out
}

dist = c_mat(d) * sum(diag(t(d) %*% d))
tsne_out <- Rtsne(dist, theta=0.0, perplexity = 3)
plot(tsne_out$Y,asp=1,type = "n")
text(tsne_out$Y,labels = seq(3*m))
eigen(MASS::ginv(t(d) %*% d))$value


y = rlnorm(nrow(d))
data <- data.frame(cbind(distance = y, d))

d_inv = MASS::ginv(d)
sigma2 <- sum((d %*% d_inv %*% y  - y)^2) / (nrow(d) - ncol(d) + 1)
#summary(lm(distance ~ . - 1, data = data[, -4]))$sigma^2
#summary(lm(distance ~ . - 1, data = data[, -3]))$sigma^2
sqrt(c_mat(d)[, 22] * sigma2)
#diag(solve(t(d[, -3]) %*% d[, -3]))


### Distances
eigs = eigen(t(d) %*% d)
lambda_dr = sqrt(diag(1/eigs$values))
lambda_dr[length(lambda_dr)] = 0
dist(t(lambda_dr %*% t(eigs$vectors)))^2
c_mat(d)

obs = lambda_dr %*% t(eigs$vectors)

eigs$vectors %*% lambda_dr^2 %*% t(eigs$vectors) - MASS::ginv(t(d) %*% d)

factoextra::fviz_dist(factoextra::get_dist(t(obs))^2)


### SVD

svd_ = svd(d)
singulars = 1/svd_$d
singulars[length(svd_$d)] = 0

sigma_ = matrix(0, nrow(svd_$u), length(singulars))
diag(sigma_) = svd_$d

svd_$u %*% diag(svd_$d) %*% t(svd_$v)
svd_$v %*% diag(singulars) %*% t(svd_$u) - MASS::ginv(d)
sqrt(diag(singulars)) %*% t(svd_$v)
obs

obs - diag(singulars) %*% t(svd_$v)

sum(abs(svd_$v %*% diag(eigs$values) %*% t(svd_$v) - t(d) %*% d))
sum(abs(eigs$vectors %*% diag(eigs$values) %*% t(eigs$vectors) - t(d) %*% d))

diag(lambda_dr)
singulars

## File 4

l = t(d) %*% d
d_ = d[, -3]
l_ = t(d_) %*% d_


r = c_mat(d)


l_inv = MASS::ginv(l)
r[3, ]
diag(solve(l_))

diag(l_inv) - l_inv[1:15, 3]^2/l_inv[3,3]



c(1, 1, 0.5, 1, 1) %*% t(c(0, 0, 1, 0, 0))
c(0, 0, 1, 0, 0) %*% t(c(1, 1, 0.5, 1, 1))

c(1, 1, -2, 1, 1) %*% t(c(1, 1, -2, 1, 1))

mat = matrix(1:16, 4, 4)
mat2 = diag(c(1, 1, 0, 1)) %*% mat %*% diag(c(1, 1, 0, 1))
mat3 = mat2[c(1, 2, 4), c(1, 2, 4)]

b = diag(nrow(l))
b[3, 3] = 0
(l_inv %*% l) %*% b - b %*% (l_inv %*% l)


b %*% l_inv %*% b
solve(l_)


t(d %*% b) %*% (d %*% b)




MASS::ginv(l) - MASS::ginv(d) %*% t(MASS::ginv(d))

j = matrix(1, 15, 15)
diag(15) - j / 15

MASS::ginv((diag(15) - j / 15) %*% b)



c_mat(d)[, 3]

prop = \(ld, j) {
  out = ld * 0
  for (i in seq(nrow(ddd))) {
    for (k in seq(nrow(ddd))) {
      out[i, k] = ld[i, k] + ld[j, j] - ld[i, j] - ld[k, j]  
    }
  }
  out
}

solve(l_) - prop(MASS::ginv(l), 3)[c(1, 2, 4:15), c(1, 2, 4:15)]


l[, -3] %*% 
  
  round(t(d) %*% d %*% MASS::ginv(b %*% t(d) %*% d %*% b) %*% t(d) %*% d)
round(t(d) %*% d %*% MASS::ginv(t(d) %*% d) %*% t(d) %*% d)

b = diag(4)
b[3, 3] = 0
a = matrix(1:16, 4, 4)
m = b %*% a %*% b
r = m[c(1, 2, 4), c(1, 2, 4)]
diags = diag(3)
diags = cbind(diags[, 1:2], rep(0, 3), diags[, 3])


a = 1/2 *(a + t(a))
b2 = matrix(0, 4, 4)
b2[3, ] = 1
a %*% b2
b3 = matrix(0, 4, 4)
b3[, 3] = 1
b3 %*% a 
-(b3 - diag(4)) %*% a %*% (b2 - diag(4))

l_inv = MASS::ginv(l)
l_iJ = matrix(rep(l_inv[, 3], 15), 15, byrow = TRUE)
l_Ji = t(l_iJ)
l_inv - l_iJ - l_Ji + l_inv[1, 1]

b = diag(nrow(l))
b[1, 1] = 0
MASS::ginv(b %*% l %*% b) - (l_inv - l_iJ - l_Ji + l_inv[1, 1])


#l_inv = MASS::ginv(l)
#l_iJ = matrix(rep(l[, 3], 15), 15, byrow = TRUE)
#l_Ji = t(l_iJ)
#MASS::ginv(b %*% l %*% b + l_iJ + l_Ji) - (l_inv + l_inv[3, 3])

j = 3
ind = setdiff(1:15, j)
l[j, j] - l[ind, j] %*% solve(l[ind, ind]) %*% l[ind, j]

l_inv

t(rep(1, 10)) %*% solve(diag(10) + rep(1, 10) %*% t(rep(1, 10)))
rep(1, 14) %*% solve(l[ind, ind]) / (k + 1)
solve(l[ind, ind]) %*% l[ind, j] * sum(solve(l[ind, ind])) / (15)^2

l[ind, j] %*% rep(1, 14)


l_inv = MASS::ginv(l)
l_iJ = matrix(rep(l_inv[ind, j], 14), 14, byrow = TRUE)
l_Ji = t(l_iJ)

l_star = (solve(l[ind, ind]))
mid = l[ind, j] %*% t(rep(1, 14)) + t(l[ind, j] %*% t(rep(1, 14)))
l_star %*% mid %*% l_star / (15) - l_iJ + l_Ji



ones = rep(1, 14) 
z = rep(1, 14) %*% 
  solve(diag(14) + l_star %*% l_star * c(t(l[ind, j]) %*% l[ind, j]))

-ones %*% solve(diag(14) + ones %*% t(ones))



1/15 * l_star %*% (l[ind, j] %*% t(ones) + ones %*% t(l[ind, j])) %*% l_star
(ones %*% t(ones) %*% l_star + l_star %*% ones %*% t(ones)) / 15

(ones %*% t(ones) %*% l_star + l_star %*% ones %*% t(ones)) / 15 - 
  (-l_iJ - l_Ji + 2 * l_inv[j, j])

l_iJ = matrix(rep(l_inv[ind, j], 14), 14, byrow = TRUE)
rep(1, 14) %*% t(rep(1, 14)) %*% 
  solve(diag(14) + l_star %*% l_star * c(t(l[ind, j]) %*% l[ind, j])) %*%
  l_star -
  l_iJ

(rep(1, 14) %*% t(rep(1, 14))) %*% l[ind, ind] * (1 + c(t(l[ind, j]) %*% l[ind, j]))



l[ind, ind] %*% l_inv[j, ind]
-l[ind, ind] %*% l_Ji - l[j, ind] %*% t(ones) * l_inv[j, j]


l_Ji - l_inv[j, ind] %*% t(ones)

(l_star %*% ones %*% t(ones))/15
-l_inv[ind, j] + l_inv[j, j]


l_star %*% ones  /15

(l_inv[ind, j] - l_inv[j, j] *ones)


MASS::ginv(t(d) %*% d) - MASS::ginv(d)%*%  t(MASS::ginv(d))

sing = 1/svd(d)$d
sing[15] = 0

### PCA

# Simulate.

m = 8
p = 20
k = 3
d = matrix(0, nrow = 3 * p + 5, ncol = 3 * m)

for (offset in c(0, 1, 2)) {
  for (i in (seq(p) + p * offset))  
    d[i, offset * m + sample(m, 2)] = c(-1, 1)
}

d[p * 3 + 1, c(m, m+1)] = c(-1, 1)
d[p * 3 + 2, c(m+1, 2 * m+1)] = c(-1, 1)
d[p * 3 + 3, c(m+2, 2 * m+2)] = c(-1, 1)
d[p * 3 + 4, c(m+3, 2 * m+3)] = c(-1, 1)
d[p * 3 + 5, c(m+4, 2 * m+4)] = c(-1, 1)


# dist = c_mat(d)
# tsne_out <- Rtsne(dist, theta=0.0, perplexity = 3)
# plot(tsne_out$Y,asp=1,type = "n")
# text(tsne_out$Y,labels = seq(3*m))
# eigen(MASS::ginv(t(d) %*% d))$value

svd_ = svd(d)
singulars = 1/svd_$d
singulars[length(svd_$d)] = 0

cols = c(rep("blue", 8), rep("red", 8), rep("black", 8))
res = umap::umap(t(diag(singulars) %*% t(svd_$v)))$layout
plot(res,asp=1,type = "n")
text(res,labels = seq(3*m), col = cols)
for (i in seq(nrow(d))) {
  indices = which(d[i, ] != 0)
  lines(res[indices, ])
}
