library(ggplot2)
library(data.table)
N.per.seg <- 100
set.seed(1)
N.segs <- 3
N.dim <- 2
max.angle <- 2*pi
(mean.mat <- matrix(
  runif(N.segs*N.dim, 0, max.angle),
  N.segs, N.dim))
sim.dt.list <- list()
true.mean.dt.list <- list()
kappa.vec <- c(0.1,1,10)
for(seg.i in 1:N.segs){
  for(dim.i in 1:N.dim){
    mean.angle <- mean.mat[seg.i, dim.i]
    for(kappa in kappa.vec){
      circ.vec <- suppressWarnings(circular::rvonmises(N.per.seg, mean.angle, kappa))
      num.vec <- as.numeric(circ.vec)
      true.mean.dt.list[[paste(seg.i, dim.i)]] <- data.table(
        seg.i, dim.i, mean.angle,
        first.i=(seg.i-1)*N.per.seg+1, 
        last.i=seg.i*N.per.seg, 
        parameter="true"
      )[, `:=`(
        start=first.i-0.5,
        end=last.i+0.5
      )]
      sim.dt.list[[paste(seg.i, dim.i, kappa)]] <- data.table(
        seg.i, dim.i, kappa, i.in.seg=1:N.per.seg, angle=num.vec)
    }
  }
}
true.mean.dt <- rbindlist(true.mean.dt.list)
(sim.dt <- rbindlist(
  sim.dt.list
)[
, data.i := seq_along(seg.i)
, by=.(kappa,dim.i)][])
Ngrid <- 10
(grid.1d <- seq(0, max.angle, l=Ngrid+1)[-(Ngrid+1)])
grid.1d.dt <- data.table(param=grid.1d, parameter="grid")
ggplot()+
  geom_point(aes(
    data.i, angle),
    shape=1,
    data=sim.dt)+
  facet_grid(dim.i ~ kappa, labeller=label_both)+
  geom_segment(aes(
    first.i, mean.angle,
    color=parameter,
    xend=last.i, yend=mean.angle),
    size=2,
    data=true.mean.dt)+
  geom_hline(aes(
    yintercept=param,
    color=parameter),
    size=1,
    data=grid.1d.dt)

grid.2d <- as.matrix(CJ(dim1=grid.1d, dim2=grid.1d))
check_angle <- function(x){
  stopifnot(0 <= x & x < max.angle)
}
geodesic_dist <- function(data.vec, param.vec){
  check_angle(data.vec)
  check_angle(param.vec)
  d <- abs(data.vec-param.vec)
  ifelse(d>max.angle/2, max.angle-d, d)
}
geodesic_L2 <- function(data.vec, param.mat){
  if(is.matrix(data.vec))stop("data.vec should be vector not matrix")
  stopifnot(identical(length(data.vec), ncol(param.mat)))
  rep.data.mat <- matrix(
    data.vec, nrow(param.mat), ncol(param.mat), byrow=TRUE)
  l1.mat <- geodesic_dist(rep.data.mat, param.mat)
  rowSums(l1.mat^2)
}
OPART <- function(x.mat, penalty){
  n.data <- nrow(x.mat)
  best.change <- rep(NA,n.data)
  best.cost <- rep(NA,n.data)
  n.dim <- ncol(x.mat)
  best.param <- matrix(NA,n.data,n.dim)
  cum.mat <- rbind(0,apply(sin.cos.mat,2,cumsum))
  sq.mat <- rbind(0,apply(sin.cos.mat^2,2,cumsum))
  best.cost <- best.change <- rep(NA,nrow(x.mat))
  cost <- function(start,end){
    sum.mat <- matrix(cum.mat[end+1,],length(start),ncol(cum.mat),byrow=TRUE)-cum.mat[start,]
    N.vec <- end-start+1
    mean.mat <- sum.mat/N.vec
    squares <- matrix(sq.mat[end+1,],length(start),ncol(cum.mat),byrow=TRUE)-sq.mat[start,]
    rowSums(squares-sum.mat^2/N.vec)
  }
  best.cost[1] <- cost(1,1)
  for(up.to in 2:length(best.cost)){
    last.start <- seq(2, up.to)
    prev.end <- last.start-1L
    prev.cost <- best.cost[prev.end]
    last.cost <- cost(last.start, up.to)
    total.cost <- prev.cost+last.cost+penalty
    best.i <- which.min(total.cost) 
    best.cost[up.to] <- total.cost[best.i]
    best.change[up.to] <- best.i
  }
  best.param <- NA#TODO
  data.table(best.change, best.cost, best.param)
}
APART <- function(angle.mat, penalty, param.mat){
  n.data <- nrow(angle.mat)
  best.change <- rep(NA,n.data)
  best.cost <- rep(NA,n.data)
  n.param <- nrow(param.mat)
  n.dim <- ncol(angle.mat)
  stopifnot(identical(n.dim, ncol(param.mat)))
  best.param <- matrix(NA,n.data,n.dim)
  cost.model <- rep(0,n.param)
  change.model <- rep(1L,n.param)
  for(data.t in 1:n.data){
    loss.vec <- geodesic_L2(angle.mat[data.t,], param.mat)
    if(data.t>1){
      cost.of.change <- best.cost[data.t-1]+penalty
      change.better <- cost.of.change<cost.model
      change.model[change.better] <- data.t
      cost.model[change.better] <- cost.of.change
    }
    cost.model <- cost.model+loss.vec
    best.param.i <- which.min(cost.model)
    best.change[data.t] <- change.model[best.param.i]
    best.param[data.t,] <- param.mat[best.param.i,]
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
      data.frame(.SD)[grepl("^V", names(.SD))]
    )]
    last.i <- cost.row$best.change-1L
  }>0)TRUE
  setkey(rbindlist(seg.dt.list),first.i)
}

APART.fit.list <- list()
both.seg.list <- list()
for(kappa in kappa.vec){
  print(kappa)
  select.dt <- data.table(kappa)
  select.sim <- sim.dt[select.dt,on="kappa"]
  data.2d <- as.matrix(
    dcast(select.sim, data.i ~ dim.i, value.var="angle")
  )[,-1]
  sin.cos.df <- data.frame(sin=sin(data.2d),cos=cos(data.2d))
  sin.cos.mat <- as.matrix(sin.cos.df)
  fpop::multiBinSeg(sin.cos.mat, 2)
  ## should try OPART.
  OPART(sin.cos.mat)
  unlink("figure-2d-hmm-sim-noise-data*.csv")
  fwrite(data.table(data.2d),"figure-2d-hmm-sim-noise-data.csv")
  py.status <- system("conda activate angular-change-paper && python figure_msmbuilder_data.py figure-2d-hmm-sim-noise-data.csv")
  if(py.status != 0)stop(py.status)
  data.list <- list()
  for(data.type in c("means", "states")){
    data.csv <- sprintf("figure-2d-hmm-sim-noise-data_%s.csv", data.type)
    data.list[[data.type]] <- fread(data.csv,header=TRUE)
  }
  setnames(data.list$means,c("V1","V2"))
  state.rle <- rle(data.list$states[[1]])
  last.i <- cumsum(state.rle$lengths)
  first.i <- 1L+c(0L,last.i[-length(last.i)])
  (pen.vec <- c(10,15,50,60,70,80,90,100,120,1000))
  hmm.mean.mat <- data.list$means[state.rle$values+1,]
  both.seg.list[[paste(kappa,"HMM")]] <- data.table(
    kappa,
    algorithm="VonMisesHMM",
    parameter.name="num.states",
    parameter.value=3,
    first.i,
    last.i,
    hmm.mean.mat,
    seg.i=seq_along(last.i),
    start=first.i-0.5,
    end=last.i+0.5)
  pen.vec <- 70
  for(penalty in pen.vec){
    afit <- APART(data.2d, penalty, grid.2d)
    asegs <- APART_decode(afit)[, `:=`(
      algorithm = "APART",
      seg.i = .I,
      start=first.i-0.5,
      end=last.i+0.5
    )][]
    APART.fit.list[[paste(kappa,penalty)]] <- data.table(
      kappa, penalty, afit)
    both.seg.list[[paste(kappa,penalty)]] <- data.table(
      kappa,
      parameter.name="penalty",
      parameter.value=penalty,
      asegs)
  }    
  ## TODO - change 2d data to 4d sin/cos run L2.
  ## - L2 on 2d data? obviously bad.
}
APART.fit <- rbindlist(APART.fit.list)
(both.seg <- rbindlist(both.seg.list, use.names=TRUE))

save(
  APART.fit, both.seg, sim.dt,
  file="figure-2d-hmm-sim-noise-data.RData")
