library(data.table)
library(ggplot2)
max.angle <- 360
seg <- function(state.i, n.data=30){
  data.table(state.i, n.data)
}
seg.dt <- rbind(
  seg(1),
  seg(2),
  seg(1),
  seg(2))
params <- function(p1, p2){
  data.table(p1, p2)
}
param.dt <- rbind(
  params(10, 100),
  params(200, 350))
set.seed(1)
data.dt.list <- list()
for(seg.i in 1:nrow(seg.dt)){
  seg.row <- seg.dt[seg.i]
  for(dimension in 1:2){
    param.row <- param.dt[seg.row$state.i]
    param <- param.row[[dimension]]
    data.vec <- rnorm(seg.row$n.data, param, 20)
    data.dt.list[[paste(seg.i, dimension)]] <- data.table(
      seg.i, dimension, param.name=names(param.row)[[dimension]],
      angle=data.vec %% max.angle)
  }
}
data.dt <- rbindlist(data.dt.list)[
, data.i := 1:.N, by=dimension
][]
ggplot()+
  geom_point(aes(
    data.i, angle),
    shape=1,
    data=data.dt)+
  facet_grid(dimension ~ ., labeller=label_both)

data.wide <- dcast(data.dt, data.i + seg.i ~ param.name, value.var="angle")
ggplot()+
  geom_point(aes(
    p1, p2),
    shape=1,
    data=data.wide)+
  facet_grid(. ~ seg.i, labeller=label_both)+
  coord_equal()

n.states <- 2
log.A.mat <- log(matrix(1/n.states, n.states, n.states))
y.mat <- as.matrix(data.wide[, .(p1,p2)])
mean.mat <- y.mat[1:n.states,]
N.data <- nrow(y.mat)
## Need to compute emission whenever parameters are updated.
check_angle <- function(x){
  stopifnot(0 <= x & x < max.angle)
}
geodesic_dist <- function(data.vec, param.vec){
  check_angle(data.vec)
  check_angle(param.vec)
  d <- abs(data.vec-param.vec)
  ifelse(d>max.angle/2, max.angle-d, d)
}
emission <- function(){
  e.mat <- matrix(NA, nrow(y.mat), n.states)
  for(state.i in 1:nrow(mean.mat)){
    mean.vec <- mean.mat[state.i,]
    mean.rep <- matrix(mean.vec, N.data, ncol(mean.mat), byrow=TRUE)
    stopifnot(identical(dim(y.mat),dim(mean.rep)))
    dist.mat <- geodesic_dist(y.mat, mean.rep)
    e.mat[,state.i] <- apply(dist.mat, 1, sum)
  }
  -e.mat
}
log.emission.mat <- emission()
log.pi.vec <- log(rep(1/n.states, n.states))
prob.dt.list <- list()
lik.dt.list <- list()
trans.dt.list <- list()
mean.dt.list <- list()
viterbi.dt.list <- list()
prior.dt.list <- list()
done <- FALSE
prev.log.lik <- -Inf
iteration <- 1
while(!done){
  fwd.list <- plotHMM::forward_interface(
    log.emission.mat, log.A.mat, log.pi.vec)
  log.alpha.mat <- fwd.list[["log_alpha"]]
  log.beta.mat <- plotHMM::backward_interface(log.emission.mat, log.A.mat)
  log.gamma.mat <- plotHMM::multiply_interface(log.alpha.mat, log.beta.mat)
  prob.mat <- exp(log.gamma.mat)
  log.xi.array <- plotHMM::pairwise_interface(
    log.emission.mat, log.A.mat, log.alpha.mat, log.beta.mat)
  ## update rules. could use https://en.wikipedia.org/wiki/Von_Mises_distribution#Estimation_of_parameters
  ## mean below comes from https://en.wikipedia.org/wiki/Circular_mean#Definition
  (log.pi.vec <- log.gamma.mat[1,])
  for(dimension in 1:ncol(mean.mat)){
    angle.vec <- y.mat[,dimension]
    atan2.args <- list()
    trig.list <- list(y=sin,x=cos)
    for(fun.name in names(trig.list)){
      fun <- trig.list[[fun.name]]
      radian.vec <- 2*pi*angle.vec/max.angle
      atan2.args[[fun.name]] <- colSums(
        fun(radian.vec)*prob.mat
      )/colSums(prob.mat)
    }
    new.angles <- max.angle*do.call(atan2, atan2.args)/(2*pi)
    mean.mat[,dimension] <- new.angles %% max.angle
  }
  log.A.mat <- plotHMM::transition_interface(
    log.gamma.mat[-N.data,], log.xi.array)
  ## M step done, now store stuff to plot.
  log.emission.mat <- emission()
  print(log.lik <- fwd.list[["log_lik"]])
  viterbi.result <- plotHMM::viterbi_interface(
    log.emission.mat, log.A.mat, log.pi.vec)
  viterbi.dt.list[[iteration]] <- data.table(
    iteration,
    data.i=1:N.data,
    state=factor(viterbi.result[["state_seq"]]))
  if(FALSE){
    trans.dt.list[[iteration]] <- data.table(
      iteration,
      prob=as.numeric(exp(log.A.mat)),
      from.state=as.integer(row(log.A.mat)),
      to.state=as.integer(col(log.A.mat)))
    prior.dt.list[[iteration]] <- data.table(
      iteration,
      state=1:n.states,
      prob=exp(log.pi.vec))
    lik.dt.list[[iteration]] <- data.table(
      iteration,
      log.lik)
    mean.dt.list[[iteration]] <- data.table(
      iteration,
      state=factor(seq_along(mean.vec)),
      sd=sd.param,
      mean=mean.vec)
    prob.dt.list[[iteration]] <- data.table(
      iteration,
      prob=as.numeric(prob.mat),
      data.i=as.integer(row(prob.mat)),
      state=factor(as.integer(col(prob.mat))))
  }
  iteration <- iteration+1
  if(log.lik <= prev.log.lik){
    done <- TRUE
  }
  prev.log.lik <- log.lik
}
viterbi.dt <- do.call(rbind, viterbi.dt.list)
prob.dt <- do.call(rbind, prob.dt.list)
mean.dt <- do.call(rbind, mean.dt.list)
trans.dt <- do.call(rbind, trans.dt.list)
prior.dt <- do.call(rbind, prior.dt.list)
lik.dt <- do.call(rbind, lik.dt.list)

best.states <- viterbi.dt[iteration==max(iteration)]
rlist <- rle(as.integer(best.states$state))
infer.seg.dt <- with(rlist, data.table(
  state.i=values, n.data=lengths
))[
, end := cumsum(n.data)+0.5
][
, start := end-n.data
][]
infer.param.dt <- data.table(
  param=as.numeric(mean.mat),
  state.i=as.integer(row(mean.mat)),
  dimension=as.integer(col(mean.mat)))
seg.param.dt <- infer.seg.dt[
  infer.param.dt, on="state.i"
][
, state := factor(state.i)
][]
gg <- ggplot()+
  ggtitle("Changes detected by 2 state HMM for 2d angular data")+
  theme_bw()+
  geom_vline(aes(
    xintercept=start),
    color="grey",
    data=infer.seg.dt[-1])+
  geom_segment(aes(
    start, param,
    xend=end, yend=param,
    color=state),
    size=2,
    data=seg.param.dt)+
  geom_point(aes(
    data.i, angle),
    shape=1,
    data=data.dt)+
  facet_grid(dimension ~ ., labeller=label_both)+
  scale_x_continuous(
    "Data sequence index")+
  scale_y_continuous(
    "Angular data value")
png('figure-2d-hmm.png', width=5, height=3, units="in", res=200)
print(gg)
dev.off()
