library(data.table)
library(dplyr)
library(MASS)
library(tweedie)
library(COUNT)
library(modeldata)

# how offsets and weights work
set.seed(111)
n <- 1000
x1 <- sample(seq(3), n, replace = TRUE)
x1b <- .25
x2 <- sample(c(0, 1), n, replace = TRUE)
x2b <- 1.25
intb <- -1
u <- intb + x1 * x1b + x2 * x2b + rnorm(n)
y <- round(exp(u), 0)

dtab <- data.table(y, x1, x2)

m1 <- glm(y ~ x1 + x2, data = dtab, family = poisson(link = "log"))
m1 %>% summary()

dtab2 <- dtab %>%
  group_by(x1, x2) %>%
  summarise(
    y = sum(y),
    n = n()
  ) %>%
  as.data.table()

m2 <- glm(y ~ x1 + x2 + offset(log(n)), data = dtab2, family = poisson(link = "log"))
m2 %>% summary()

dtab2[, ye := y / n]
m3 <- glm(ye ~ x1 + x2, data = dtab2, weights = dtab2$n, family = poisson(link = "log"))
m3 %>% summary()

# model diagnostics
dtab <- modeldata::car_prices %>%
  as.data.table()

m1 <- glm(
  Price ~ I(Mileage / 1000)
    + Cylinder
    + Doors
    + convertible
    + hatchback
    + sedan,
  data = dtab,
  family = gaussian()
)

m2 <- glm(
  Price ~ I(Mileage / 1000),
  data = dtab,
  family = gaussian()
)

m3 <- glm(
  Price ~ 1,
  data = dtab,
  family = gaussian()
)

# univariate
univariate <- function(target, exposure, grp, plot = FALSE) {
  x <- data.frame(target, exposure, grp) %>%
    group_by(grp) %>%
    summarise(
      target = sum(target),
      exposure = sum(exposure),
      n = n()
    ) %>%
    transmute(
      group = grp,
      avgtarget = target / exposure,
      relativity = avgtarget / (sum(target) / sum(exposure)),
      target,
      exposure,
      n
    )
  if (plot) {
    # todo
  }
  return(x)
}
gbrk <- quantile(dtab$Price, seq(0, 1, 1 / 20)) %>% unique()
grp <- cut(dtab$Price, breaks = gbrk, include.lowest = TRUE, right = TRUE, dig.lab = 6, ordered_result = TRUE)
univariate(dtab$Price, rep(1, nrow(dtab)), grp)

# or wtd.quantile

# lift table
lift <- function(target, exposure, score, bins = 10, plot = FALSE) {
  require(dplyr)
  x <- data.frame(target, exposure, score) %>%
    arrange(score) %>%
    mutate(
      g = pmin(floor(1 + bins * cumsum(exposure) / sum(exposure)), bins)
    )
  x <- x %>%
    group_by(group = g) %>%
    summarise(
      target = sum(target),
      exposure = sum(exposure),
      n = n()
    ) %>%
    transmute(
      group,
      avgtarget = target / exposure,
      relativity = avgtarget / (sum(target) / sum(exposure)),
      target,
      exposure,
      n
    )
  if (plot) {
    barplot(x$avgtarget, main = "Lift Chat", col = "black")
  }
  return(list(
    data = x,
    lift = as.vector(unlist(x[which.max(x$group), "relativity"] / x[which.min(x$group), "relativity"]))
  ))
       
}
lift(dtab$Price, rep(1, nrow(dtab)), predict(m1), bins = 10, plot = TRUE)

# double lift
doublelift <- function(target, exposure, score1, score2, bins = 10, plot = FALSE) {
  require(dplyr)
  x <- data.frame(target, exposure, score1, score2) %>%
    arrange(score1 / score2) %>%
    mutate(
      g = pmin(floor(1 + bins * cumsum(exposure) / sum(exposure)), bins)
    ) %>%
    group_by(group = g) %>%
    summarise(
      target = sum(target),
      exposure = sum(exposure),
      n = n()
    ) %>%
    transmute(
      group,
      avgtarget = target / exposure,
      relativity = avgtarget / (sum(target) / sum(exposure)),
      target,
      exposure,
      n
    )
  if (plot) {
    barplot(x$avgtarget, main = "Double Lift Chart", col = "black")
  }
  return(list(
    data = x,
    lift = as.vector(unlist(x[which.max(x$group), "relativity"] / x[which.min(x$group), "relativity"]))
  ))
}
doublelift(dtab$Price, rep(1, nrow(dtab)), predict(m1), predict(m2), bins = 10, plot = TRUE)
doublelift(dtab$Price, rep(1, nrow(dtab)), predict(m1), predict(m3), bins = 10, plot = TRUE)
doublelift(dtab$Price, rep(1, nrow(dtab)), predict(m1), predict(m1), bins = 10, plot = TRUE)

# lorenz and gini
gini <- function(target, exposure, score) {
  require(dplyr)
  x <- data.frame(target, exposure, score) %>%
    arrange(score) %>%
    mutate(
      ct = cumsum(target) / sum(target),
      ce = cumsum(exposure) / sum(exposure)
    )
  ginif <- approxfun(c(0, x$ce), c(0, x$ct), na.rm = TRUE)
  ginii <- integrate(ginif, min(x$ce), max(x$ce))$value
  return(list(gini_func = ginif, gini_index = 1 - 2 * ginii))
}

lorenz <- function(target, exposure, score, plot = FALSE) {
  g <- gini(target, exposure, score)
  x <- seq(0, 1, 1 / 100)
  y <- g$gini_func(x)
  ng <- gini(target, exposure, target)
  ny <- ng$gini_func(x)
  if (plot) {
    plot(x, y, type = "l", main = "Lorenz Curve", col = "red")
    lines(y, y, col = "black")
    lines(x, ny, col = "blue")
  }
  return(list(
    data = data.table(x, y),
    gini_index = g$gini_index,
    gini_index_normalized = g$gini_index / ng$gini_index
  ))
}

gini(dtab$Price, rep(1, nrow(dtab)), predict(m1))
lorenz(dtab$Price, rep(1, nrow(dtab)), predict(m1), plot = TRUE)

# time saving things for GLMs and similar models
# list methods for getting statistics
m3_ws <- m3
m3_ws$coefficients[1] <- 0
m3_ws %>% coef()
predict(m3_ws, newdata=dtab[seq(10)]) %>% head()
predict(m3_ws) %>% head()
predict(m3, newdata=dtab[seq(10)]) %>% head()
predict(m3) %>% head()

m1_coefs <- m1 %>% coef()
m1_mm <- model.matrix(m1)
m1_coefs %*% t(head(m1_mm))
predict(m1) %>% head()

# model plots
# resid vs fitted
x <- fitted(m1)
y <- residuals(m1)
l1 <- 1
l2 <- 0
plot(x,y)
lines(loess.smooth(x,y), col = "red")
abline(0, 0) # horizontal line

# qqplot, see below

# sq std pearson vs. fitted
x <- fitted(m1)
y <- sqrt(residuals(m1, type="pearson"))

# prediction vs. confidence in predict interval arg

# loess smoothing
# you can fit a spline and this isnt needed
#  or you can fit binned categories and smooth using loess

# rounding in R vs excel

# qq plot
set.seed(556)
x <- rpois(1000, 3)
xx <- qpois(ppoints(length(x)), 3)[order(order(x))]
plot(xx, x, main = "QQ Plot for Poisson dist.", pch = 20, cex = .5)
xxq <- quantile(xx, probs = c(.25, .75))
xq <- quantile(x, probs = c(.25, .75))
slope <- diff(xq) / diff(xxq)
int <- xq[1] - slope * xxq[1]
abline(int, slope)
