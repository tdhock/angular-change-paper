library(data.table)
library(ggplot2)
data.vec <- c(1, 4)
N.grid <- 101
param.grid <- seq(0, 2*pi, l=N.grid)
loss.dt.list <- list()
sin_cos_mat <- function(x){
  cbind(sin=sin(x), cos=cos(x))
}
norm.list <- list(
  L1=identity,
  L2=function(x)x^2)
for(data.i in seq_along(data.vec)){
  data.val <- data.vec[[data.i]]
  d <- abs(data.val-param.grid)
  angle.dist <- ifelse(d>pi, 2*pi-d, d)
  for(norm.name in names(norm.list)){
    norm.fun <- norm.list[[norm.name]]
    sin.cos.loss <- rowSums(norm.fun(abs(
      sin_cos_mat(param.grid)-sin_cos_mat(rep(data.val,N.grid)))))
    base.dt <- data.table(data.val, param.grid)
    loss.dt.list[[paste(data.i, norm.name)]] <- rbind(
      data.table(base.dt, norm.name, loss.name="angle", loss.value=norm.fun(angle.dist)),
      data.table(base.dt, norm.name, loss.name="sin.cos", loss.value=sin.cos.loss))
  }
}
loss.dt <- rbindlist(loss.dt.list)
min.dt <- loss.dt[
, .SD[which.max(loss.value)]
, by=.(data.val, norm.name, loss.name)]
gg <- ggplot()+
  facet_grid(norm.name ~ data.val, labeller=label_both, scales="free")+
  geom_vline(aes(
    xintercept=param.grid, color=loss.name),
    data=min.dt)+
  geom_point(aes(
    param.grid, loss.value, color=loss.name),
    data=loss.dt)
png("figure-loss.png", width=6, height=6, res=100, units="in")
print(gg)
dev.off()
