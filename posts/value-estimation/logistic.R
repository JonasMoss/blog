data_list <- list(
  linch = jsonlite::fromJSON("posts/value-estimation/raw_data/linch-zhang.json"),
  finn = jsonlite::fromJSON("posts/value-estimation/raw_data/finn-moorhouse.json"),
  gavin = jsonlite::fromJSON("posts/value-estimation/raw_data/gavin-leech-complete.json"),
  jamie = jsonlite::fromJSON("posts/value-estimation/raw_data/jaime-sevilla.json"),
  misha = jsonlite::fromJSON("posts/value-estimation/raw_data/misha-yagudin.json"),
  ozzie = jsonlite::fromJSON("posts/value-estimation/raw_data/ozzie-gooen.json")
)

source("posts/value-estimation/functions.R")

mod1 <- lm(distance ~ . - 1, data = make_frame(data_list$misha))

###============================================================================
### Find perfectly separated points.
###============================================================================

name <- "linch"
set.seed(11111) #9 does not work with p = 0.5
p <- 0.9
data = make_frame(data_list[[name]])
data$binary <- 1 * (data$distance >= 0)
indices <- rbinom(nrow(data), 1, p)
indices[1] <- 1
indices[21] <- 1
indices[10] <- 1
#indices <- rep(0, nrow(data))
data$distance[!indices] <- NA
data <- dplyr::relocate(data, binary, .after = 1)


mod2 <- double_est(data, maxit = 10000)
all = c(mod2$beta, coef(mod1))
range = c(min(all), max(all))
plot(coef(mod1), mod2$beta, xlim=range, ylim = range)
abline(a = 0, b = 1)

library("igraph")
misha <- data_list[[name]]
levels <- levels(as.factor(c(misha$source, misha$target)))
#misha <- misha[which(!!indices), ]
source <- as.numeric(factor(misha$source, levels = levels))
target <- as.numeric(factor(misha$target, levels = levels))
m <- length(source)
select <- which(!!indices)
source <- c(source, target[select])
target <- c(target, source[select])

graph <- igraph::graph_from_edgelist(cbind(source, target))
E(graph)$color <- "red"
E(graph)$color[which(!!indices)] <- "blue"
E(graph)$color[(m+1):length(target)] <- "blue"
plot(graph)

cis <- cbind(confint(mod1), cis(mod2))
apply(cis, 1, \(x) c(x[2] - x[1], x[4] - x[3]))

components(graph)
