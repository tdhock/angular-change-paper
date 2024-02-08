library(data.table)
library(ggplot2)
data.vec <- c(3, 4)
param.grid <- seq(0, 2*pi, l=201)
loss.dt.list <- list()
zero.dt.list <- list()
get.zeros_ <- function(constant, sin.coef, cos.coef){
  under.sqrt <- -constant^2+sin.coef^2+cos.coef^2
  if(under.sqrt<1e-6) under.sqrt <- 0
  numerator <- c(1,-1)*sqrt(under.sqrt)-sin.coef
  denominator <- constant-cos.coef
  atan.term <- atan(numerator/denominator)
  int <- 0
  out <- (atan.term+pi*int)*2
  ifelse(out<0, out+2*pi, out)
}
get.zeros <- function(loss.coefs){
  do.call(get.zeros_, as.list(loss.coefs))
}
get.coefs <- function(data.val){
  dx <- sin(data.val)
  dy <- cos(data.val)
  c(2, -2*dx, -2*dy)
}
get.loss_ <- function(param.grid, constant, sin.coef, cos.coef){
  constant+sin.coef*sin(param.grid)+cos.coef*cos(param.grid)
}
get.loss <- function(param.grid, loss.coefs){
  do.call(get.loss_, c(list(param.grid), as.list(loss.coefs)))
}
coef.list <- list()
for(data.i in seq_along(data.vec)){
  data.val <- data.vec[[data.i]]
  loss.coefs <- get.coefs(data.val)
  coef.list[[data.i]] <- loss.coefs
  zero <- get.zeros(loss.coefs)
  loss <- get.loss(param.grid, loss.coefs)
  loss.dt.list[[data.i]] <- data.table(
    data.val,
    param.grid,
    loss)
  zero.dt.list[[data.i]] <- data.table(
    data.val,
    zero)
}
loss.dt <- rbindlist(loss.dt.list)
zero.dt <- rbindlist(zero.dt.list)
ggplot()+
  geom_hline(
    yintercept=0,
    color="grey50")+
  geom_point(aes(
    param.grid, loss),
    data=loss.dt)+
  geom_vline(aes(
    xintercept=zero),
    data=zero.dt,
    color="red")+
  coord_cartesian(xlim=c(0,2*pi))+
  facet_grid(. ~ data.val, labeller=label_both)+
  theme_bw()

diff.coefs <- coef.list[[1]]-coef.list[[2]]
diff.dt <- data.table(
  param.grid,
  loss=get.loss(param.grid, diff.coefs))
ggplot()+
  geom_hline(
    yintercept=0,
    color="grey50")+
  geom_point(aes(
    param.grid, loss),
    data=diff.dt)+
  coord_cartesian(xlim=c(0,2*pi))+
  theme_bw()

both.dt <- rbind(
  diff.dt[, .(y="diff", data.val="diff", param.grid, loss)],
  data.table(y="loss", loss.dt))
diff.zeros <- data.table(
  zero=get.zeros(diff.coefs))
gg <- ggplot()+
  geom_hline(aes(
    yintercept=yintercept),
    data=data.table(yintercept=0, y="diff"),
    color="grey50")+
  geom_point(aes(
    param.grid, loss, color=data.val),
    data=both.dt)+
  geom_vline(aes(
    xintercept=zero),
    data=diff.zeros)+
  coord_cartesian(xlim=c(0,2*pi))+
  scale_x_continuous(breaks=seq(0, 10))+
  facet_grid(y ~ ., scales="free")+
  theme_bw()
png("figure-roots-sin-cos.png", width=6, height=6, units="in", res=100)
print(gg)
dev.off()
