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

APART <- function(angle.vec, penalty, param.vec){
  n.data <- length(angle.vec)
  best.change <- rep(NA,n.data)
  best.cost <- rep(NA,n.data)
  best.param <- rep(NA,n.data)
  n.param <- length(param.vec)
  cost.model <- rep(0,n.param)
  change.model <- rep(1L,n.param)
  for(data.t in 1:n.data){
    loss.vec <- geodesic_dist(angle.vec[data.t], param.vec)
    if(data.t>1){
      cost.of.change <- best.cost[data.t-1]+penalty
      change.better <- cost.of.change<cost.model
      change.model[change.better] <- data.t
      cost.model[change.better] <- cost.of.change
    }
    cost.model <- cost.model+loss.vec
    best.param.i <- which.min(cost.model)
    best.change[data.t] <- change.model[best.param.i]
    best.param[data.t] <- param.vec[best.param.i]
    best.cost[data.t] <- cost.model[best.param.i]
  }
  data.table(best.change, best.cost, best.param)
}
APART_decode <- function(DT){
  seg.dt.list <- list()
  last.i <- DT[, .N]
  while({
    cost.row <- DT[last.i]
    seg.dt.list[[paste(last.i)]] <- cost.row[, data.table(
      first.i=best.change,
      last.i,
      param=best.param
    )]
    last.i <- cost.row$best.change-1L
  }>0)TRUE
  setkey(rbindlist(seg.dt.list),first.i)
}

angle.breaks <- seq(0, 360, by=90)
show.seg.list <- list()
show.resid.list <- list()
for(penalty in c(1, 100, 10000)){
  (cost.dt <- APART(data.dt$angle, penalty, param.dt$param))
  (seg.dt <- APART_decode(cost.dt))
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
  print(rbind(
    check=penalty*(nrow(seg.dt)-1)+residual.dt[, sum(abs(y-yend))],
    recursion=cost.dt[.N, best.cost]))
  show.seg.list[[paste(penalty)]] <- data.table(penalty, seg.dt)
  show.resid.list[[paste(penalty)]] <- data.table(penalty, residual.dt)
  gg <- ggplot()+
    ggtitle(paste0(
      "Approximate partitioning (APART), geodesic loss, penalty=",
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
      "Angle (degrees)",
      breaks=angle.breaks)
  out.png <- paste0("figure-approx-algo-1d-", penalty, ".png")
  png(out.png, 6, 2, units="in", res=200)
  print(gg)
  dev.off()
}

show.resid <- rbindlist(show.resid.list)
show.seg <- rbindlist(show.seg.list)
gg <- ggplot()+
  ggtitle(
    "Approximate partitioning (APART), geodesic loss")+
  theme_bw()+
  theme(panel.grid.minor=element_blank())+
  geom_segment(aes(
    first.i-0.5, param,
    xend=last.i+0.5, yend=param),
    color=seg.color,
    linewidth=1,
    data=show.seg)+
  geom_vline(aes(
    xintercept=first.i-0.5),
    color=seg.color,
    data=show.seg[, .SD[-1], by=penalty])+
  geom_segment(aes(
    data.i, y,
    xend=data.i, yend=yend),
    data=show.resid)+
  geom_point(aes(
    data.i, angle),
    shape=1,
    data=data.dt)+
  scale_x_continuous(
    "Data sequence index",
    breaks=data.dt$data.i)+
  scale_y_continuous(
    "Angle (degrees)",
    breaks=angle.breaks)+
  facet_grid(penalty ~ ., labeller=label_both)
png("figure-approx-algo-1d.png", 6, 3.5, units="in", res=200)
print(gg)
dev.off()
