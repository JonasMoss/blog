data_list <- list(
  linch = jsonlite::fromJSON("posts/value-estimation/raw_data/linch-zhang.json"),
  finn = jsonlite::fromJSON("posts/value-estimation/raw_data/finn-moorhouse.json"),
  gavin = jsonlite::fromJSON("posts/value-estimation/raw_data/gavin-leech-complete.json"),
  jamie = jsonlite::fromJSON("posts/value-estimation/raw_data/jaime-sevilla.json"),
  misha = jsonlite::fromJSON("posts/value-estimation/raw_data/misha-yagudin.json"),
  ozzie = jsonlite::fromJSON("posts/value-estimation/raw_data/ozzie-gooen.json")
)


data <- make_frame(data_list[[3]])
data$x <- data$distance >= 0

pairwise_model_logit <- function(data, fixed = 3, keep_names = FALSE) {
  data_frame <- make_frame(data, fixed = fixed, keep_names = keep_names)
  data_frame$distance <- data_frame$distance >= 0
  glm(distance ~ . - 1, data = data_frame, family = binomial(link = "logit"))
}




plotter <- \(i) {
  x <- coef(pairwise_model(data_list[[i]]))
  y <- coef(pairwise_model_logit(data_list[[i]]))
  
  plot(x, y, xlab = "Scale kept", ylab = "Scale not kept",
       sub = paste0("Correlation: ", round(cor(x, y), 3)),
       main = names(data_list)[i])
  abline(lm(y ~ x))  
}

plotter(6)
