library(ggplot2)
library(data.table)
data.per.seg <- 10
set.seed(1)
angle.vec <- c(
  runif(data.per.seg,-10,10),
  runif(data.per.seg,180,190),
  runif(data.per.seg,-10,10))
max.angle <- 360
data.dt <- data.table(
  data.i=seq_along(angle.vec),
  angle=ifelse(angle.vec<0, angle.vec+max.angle, angle.vec))
param.dt <- data.table(
  param=as.numeric(seq(0, max.angle))
)[-.N]
gg <- ggplot()+
  theme_bw()+
  geom_hline(aes(
    yintercept=param),
    color="grey",
    data=param.dt)+
  geom_point(aes(
    data.i, angle),
    data=data.dt)
gg

check_angle <- function(x){
  stopifnot(0 <= x & x < max.angle)
}
geodesic_dist <- function(data.vec, param.vec){
  check_angle(data.vec)
  check_angle(param.vec)
  data.table(data.val=data.vec, param.val=param.vec)[, ifelse({
    data.val==0
  }, ifelse(
    param.val<max.angle/2, param.val,
    max.angle-param.val
  ),
  ifelse({
    data.val<max.angle/2
  }, ifelse(
    param.val<data.val, data.val-param.val,
    ifelse(
      param.val<data.val+max.angle/2,
      param.val-data.val,
      max.angle+data.val-param.val)
  ),
  ifelse({
    data.val==max.angle/2
  }, ifelse(
      param.val<max.angle/2, max.angle/2-param.val,
      param.val-max.angle/2
  ), ifelse(
    param.val<data.val-max.angle/2, param.val+max.angle-data.val,
    ifelse(
      param.val<data.val,
      data.val-param.val,
      param.val-data.val))))
  )]
}
geodesic_dist <- function(data.vec, param.vec){
  check_angle(data.vec)
  check_angle(param.vec)
  d <- abs(data.vec-param.vec)
  ifelse(d>max.angle/2, max.angle-d, d)
}
plot(geodesic_dist(300, param) ~ param, param.dt)
cum.mat <- sapply(param.dt$param, function(param.val){
  c(0,cumsum(geodesic_dist(data.dt$angle, param.val)))
})

penalty <- 100
n.out <- nrow(data.dt)+1
best.change <- rep(NA,n.out)
best.cost <- rep(NA,n.out)
best.cost[1] <- 0
best.param <- rep(NA,n.out)
for(data.t in 1:nrow(data.dt)){
  last.seg.start <- seq(1, data.t)
  out.t <- data.t+1
  total.up.to.t <- cum.mat[out.t,]
  total.before.start <- t(cum.mat[last.seg.start,,drop=FALSE])
  last.seg.cost.mat <- total.up.to.t-total.before.start
  best.param.i <- apply(last.seg.cost.mat, 2, which.min)
  (best.params <- param.dt$param[best.param.i])
  last.seg.cost <- apply(last.seg.cost.mat, 2, min)
  prev.cost <- best.cost[last.seg.start]+penalty
  total.cost <- prev.cost+last.seg.cost
  best.start <- which.min(total.cost)
  best.change[out.t] <- best.start
  best.cost[out.t] <- total.cost[best.start]
  best.param[out.t] <- best.params[best.start]
}
cost.dt <- data.table(best.change, best.cost, best.param)[-1]

##decoding.
seg.dt.list <- list()
last.i <- cost.dt[, .N]
while({
  cost.row <- cost.dt[last.i]
  seg.dt.list[[paste(first.i,last.i)]] <- cost.row[, data.table(
    first.i=best.change,
    last.i,
    param=best.param
  )]
  last.i <- cost.row$best.change-1L
}>0)TRUE
(seg.dt <- setkey(rbindlist(seg.dt.list),first.i))
seg.color <- "red"
dist.dt <- seg.dt[
  data.dt,
  data.table(data.i, angle, param),
  on=.(first.i<=data.i,last.i>=data.i)
][
, dist := abs(angle-param)
][
, segs := ifelse(dist<max.angle/2, 1, 2)
][]
residual.dt <- dist.dt[, if(segs==1){
  data.table(y=angle,yend=param)
}else{
  yend <- c(angle,param)
  data.table(y=ifelse(yend<max.angle/2,0,max.angle),yend)
}, by=.(segs,data.i)]
gg <- ggplot()+
  ggtitle(paste0(
    "Approximately optimal partitioning (APART) with penalty=",
    penalty))+
  theme_bw()+
  theme(panel.grid.minor=element_blank())+
  geom_segment(aes(
    first.i-0.5, param,
    xend=last.i+0.5, yend=param),
    color=seg.color,
    linewidth=1,
    data=seg.dt)+
  geom_vline(aes(
    xintercept=first.i-0.5),
    color=seg.color,
    data=seg.dt[-1])+
  geom_segment(aes(
    data.i, y,
    xend=data.i, yend=yend),
    data=residual.dt)+
  geom_point(aes(
    data.i, angle),
    shape=1,
    data=data.dt)+
  scale_x_continuous(
    "Data sequence index",
    breaks=data.dt$data.i)+
  scale_y_continuous(
    "Angle (degrees)")
png("figure-approx-algo-1d.png", 6, 2, units="in", res=200)
print(gg)
dev.off()

