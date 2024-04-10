library(ggplot2)
library(data.table)
data.per.seg <- 10
set.seed(1)
angle.vec <- c(
  runif(data.per.seg,-10,10),
  runif(data.per.seg,170,190),
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
  d <- abs(data.vec-param.vec)
  ifelse(d>max.angle/2, max.angle-d, d)
}
plot(geodesic_dist(300, param) ~ param, param.dt)

penalty <- 100
n.out <- nrow(data.dt)
best.change <- rep(n.out)
best.cost <- rep(n.out)
best.cost[1] <- 0
best.param <- rep(n.out)
cost.model <- rep(0,nrow(param.dt))
change.model <- rep(1L,nrow(param.dt))
for(data.t in 1:nrow(data.dt)){
  loss.vec <- geodesic_dist(data.dt[data.t, angle], param.dt$param)
  if(data.t>1){
    cost.of.change <- best.cost[data.t-1]+penalty
    change.better <- cost.of.change<cost.model
    change.model[change.better] <- data.t
    cost.model[change.better] <- cost.of.change
  }
  cost.model <- cost.model+loss.vec
  best.param.i <- which.min(cost.model)
  best.change[data.t] <- change.model[best.param.i]
  best.param[data.t] <- param.dt$param[best.param.i]
  best.cost[data.t] <- cost.model[best.param.i]
}
(cost.dt <- data.table(best.change, best.cost, best.param))
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
rbind(
  check=penalty*(nrow(seg.dt)-1)+residual.dt[, sum(abs(y-yend))],
  recursion=cost.dt[.N, best.cost])
