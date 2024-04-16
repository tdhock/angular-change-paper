rob_intervals <- function(data.vec,penalty,lthreshold,lslope){
  L <- robseg::Rob_seg(data.vec,penalty,lthreshold=lthreshold,lslope=lslope)
  with(L, data.table(
    segments=as.numeric(length(L$t.est)),
    mean.intervals=mean(intervals),
    max.intervals=as.numeric(max(intervals))))
}
expr.list <- atime::atime_grid(
  list(penalty=c(0.01,1000)),
  Poisson={
    fit <- PeakSegOptimal::PeakSegFPOP(count.vec, penalty=penalty)
    with(fit, data.table(
      segments=as.numeric(sum(ends.vec >= 0)),
      mean.intervals=mean(intervals.mat),
      max.intervals=as.numeric(max(intervals.mat))))
  },
  geodesic={
    fit <- geodesichange::geodesicFPOP_vec(angle.vec, penalty)
    with(fit$loss, data.table(
      segments=as.numeric(segments),
      mean.intervals,
      max.intervals=as.numeric(max.intervals)))
  },
  L2=rob_intervals(data.vec,penalty,lthreshold=1000,lslope=0),
  L1=rob_intervals(data.vec,penalty,lthreshold=0,lslope=-1),
  L2rob=rob_intervals(data.vec,penalty,lthreshold=1,lslope=0))
atime.result <- atime::atime(
  N=2^seq(2, 100),
  expr.list=expr.list,
  setup={
    data.vec <- rnorm(N)
    count.vec <- rpois(N, 50)
    angle.vec <- ifelse(data.vec<0, data.vec+2*pi, data.vec)
  },
  seconds.limit=10,
  result=TRUE)
saveRDS(atime.result, "figure-geodesicFPOP-robseg-data.rds")
