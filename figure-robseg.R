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
data.vec <- rnorm(1000)
for(param.i in 1:nrow(param.dt)){
  param.row <- param.dt[param.i]
  param.list <- c(list(x=data.vec), as.list(param.row))
  L <- do.call(robseg::Rob_seg, param.list)
  print(with(L, data.table(param.row, intervals, path)))
}
