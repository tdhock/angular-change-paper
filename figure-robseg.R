if(FALSE){
  remotes::install_github("tdhock/robust-fpop@interval-count")
}
##double lBorder = X_- lthrs_;
##double rBorder = X_ + rthrs_;
library(data.table)
set.seed(1)
param <- function(lambda, lthreshold, lslope){
  data.table(lambda, lthreshold, lslope)
}
param.dt <- rbind(
  param(100, 1, 0),
  param(0.1, 1, 0),
  param(100, 0, -1),
  param(0.1, 0, -1))
N <- 1000
data.vec <- rnorm(N)
for(param.i in 1:nrow(param.dt)){
  param.row <- param.dt[param.i]
  param.list <- c(list(x=data.vec), as.list(param.row))
  L <- do.call(robseg::Rob_seg, param.list)
  robseg.dt <- print(with(L, data.table(param.row, intervals, path)))
}
count.vec <- rpois(N, 50)
for(lambda in c(100, 0.1)){
  fit <- PeakSegOptimal::PeakSegFPOP(count.vec, penalty=lambda)
  fit.dt <- with(fit, data.table(
    lambda, path=rev(ends.vec), intervals=t(intervals.mat)))
  print(fit.dt)
}

N=10
count.vec <- as.integer(2^seq(1,N))
pen.vec <- c(1000000, 0.1)
for(lambda in pen.vec){
  fit <- PeakSegOptimal::PeakSegFPOP(count.vec, penalty=lambda)
  fit.dt <- with(fit, data.table(
    lambda, path=rev(ends.vec), intervals=t(intervals.mat)))
  print(fit.dt)
  rob.fit <- robseg::Rob_seg(count.vec, lambda=lambda, lthreshold=1)
  print(with(rob.fit, data.table(lambda, t.est, intervals)))
}

