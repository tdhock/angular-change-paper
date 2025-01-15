library(ggplot2)
library(data.table)
N.per.seg <- 30
set.seed(1)
N.means <- 3
N.dim <- 2
max.angle <- 2*pi
(mean.mat <- matrix(
  runif(N.means*N.dim, 0, max.angle),
  N.means, N.dim))
sim.dt.list <- list()
true.mean.dt.list <- list()
N.segs <- 12
seg.ord.vec <- sample(rep(1:N.means, l=N.segs))
is.true.change <- diff(seg.ord.vec) != 0
between.segs <- seq(N.per.seg, N.per.seg*(N.segs-1), N.per.seg)+0.5
true.change.vec <- between.segs[is.true.change]
N.true.changes <- sum(is.true.change)
N.true.segs <- N.true.changes+1
kappa.vec <- c(
  ##0.5,
  1, 2, 5)
for(seg.i in 1:N.segs){
  seg.mean.index <- seg.ord.vec[seg.i]
  seg.mean.vec <- mean.mat[seg.mean.index,]
  for(dim.i in 1:N.dim){
    mean.angle <- seg.mean.vec[dim.i]
    for(kappa in kappa.vec){
      circ.vec <- suppressWarnings(
        circular::rvonmises(N.per.seg, mean.angle, kappa))
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
APART.mean.mat.list <- list()
APART.mean.dt.list <- list()
for(kappa in kappa.vec){
  print(kappa)
  select.dt <- data.table(kappa)
  select.sim <- sim.dt[select.dt,on="kappa"]
  data.2d <- as.matrix(
    dcast(select.sim, data.i ~ dim.i, value.var="angle")
  )[,-1]
  sin.cos.dt <- data.table(sin=sin(data.2d), cos=cos(data.2d))
  sin.cos.mat <- as.matrix(sin.cos.dt)
  km.fit <- kmeans(sin.cos.mat, N.means)
  km.angles <- melt(
    data.table(mean.i=1:nrow(km.fit$centers),km.fit$centers),
    measure.vars=measure(value.name, dim.i, pattern="(sin|cos)[.]([12])")
  )[
  , mean.radians := atan2(sin,cos)
  ][
  , param := ifelse(mean.radians < 0, mean.radians+2*pi, mean.radians)
  ][]
  APART.mean.mat.list[[paste(kappa)]] <- as.matrix(dcast(km.angles, mean.i ~ dim.i, value.var="param")[,-1])
  APART.mean.dt.list[[paste(kappa)]] <- data.table(
    kappa,
    km.angles,
    parameter="kmeans")
}
APART.mean.dt <- rbindlist(APART.mean.dt.list)
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
    data=APART.mean.dt)

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
  cum.mat <- rbind(0,apply(x.mat,2,cumsum))
  sq.mat <- rbind(0,apply(x.mat^2,2,cumsum))
  best.cost <- best.change <- best.Nsegs <- rep(NA,nrow(x.mat))
  for(last.seg.end in seq_along(best.cost)){
    last.seg.start <- seq(1, last.seg.end)
    get_sum <- function(m){
      matrix(
        m[last.seg.end+1,], last.seg.end, ncol(m), byrow=TRUE
      )-m[last.seg.start,]
    }
    N.vec <- last.seg.end-last.seg.start+1
    candidate.dt <- data.table(
      last.seg.start,
      last.seg.cost=rowSums(get_sum(sq.mat)-get_sum(cum.mat)^2/N.vec)
    )[
      last.seg.start==1, `:=`(
        prev.cost=0,
        Nsegs=1
      )
    ][
      last.seg.start>1, `:=`(
        prev.cost=best.cost[last.seg.start-1]+penalty,
        Nsegs=best.Nsegs[last.seg.start-1]+1L
      )
    ][
    , total.cost := prev.cost + last.seg.cost
    ][]
    best.row <- candidate.dt[which.min(total.cost)]
    best.cost[last.seg.end] <- best.row$total.cost
    best.Nsegs[last.seg.end] <- best.row$Nsegs
    best.change[last.seg.end] <- best.row$last.seg.start
  }
  data.table(best.change, best.cost, best.Nsegs)
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
  setkey(rbindlist(seg.dt.list),first.i)[]
}

APART.fit.list <- list()
both.seg.list <- list()
for(kappa in kappa.vec){
  print(kappa)
  APART.mean.mat <- APART.mean.mat.list[[paste(kappa)]]
  select.dt <- data.table(kappa)
  select.sim <- sim.dt[select.dt,on="kappa"]
  data.2d <- as.matrix(
    dcast(select.sim, data.i ~ dim.i, value.var="angle")
  )[,-1]
  sin.cos.dt <- data.table(sin=sin(data.2d),cos=cos(data.2d))
  sin.cos.mat <- as.matrix(sin.cos.dt)
  bs.fit <- fpop::multiBinSeg(sin.cos.mat,N.true.changes)
  after.change <- sort(bs.fit$t.est)
  first.i <- c(1,after.change+1)
  last.i <- c(after.change,nrow(sin.cos.mat))
  cum.mat <- rbind(0,apply(sin.cos.mat,2,cumsum))
  sum.mat <- cum.mat[last.i+1,]-cum.mat[first.i,]
  N.vec <- last.i-first.i+1
  bs.mean <- melt(
    data.table(first.i,last.i,sum.mat/N.vec),
    measure.vars=measure(value.name, dim.i, pattern="(sin|cos)[.]([12])")
  )[
  , mean.radians := atan2(sin,cos)
  ][]
  (bs.wide <- dcast(
    bs.mean[, V := paste0("V",dim.i)],
    first.i+last.i ~ V,
    value.var="mean.radians"))
  both.seg.list[[paste(kappa,"SinCos")]] <- data.table(
    kappa,
    algorithm="SinCosL2BinSeg",
    parameter.name="num.segs",
    parameter.value=N.true.segs,
    bs.wide,
    seg.i=seq_along(last.i),
    start=first.i-0.5,
    end=last.i+0.5)
  if(FALSE){
    ## should try OPART.
    (ofit <- OPART(sin.cos.mat, 5.5))
    APART_decode(ofit)
    solver_OPART <- function(data, penalty){
      result <- OPART(data, penalty)
      changes <- result[.N, best.Nsegs]-1
      list(
        penmap_loss=result[.N, best.cost]-changes*penalty,
        penmap_size=changes,
        details=result)
    }
    pmap <- penmap::penmap_create(solver_OPART, sin.cos.mat)
    penmap::penmap_solve_plot(pmap)
    pmap$penmap$df()
    stop(1)
  }
  unlink("figure-2d-hmm-sim-noise-data*.csv")
  fwrite(data.table(data.2d),"figure-2d-hmm-sim-noise-data.csv")
  py.status <- system("/home/tdhock/miniconda3/bin/conda run -n angular-change-paper python figure_msmbuilder_data.py figure-2d-hmm-sim-noise-data.csv")
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
  solver_APART <- function(data, penalty){
    result <- APART(data, penalty, APART.mean.mat)
    asegs <- APART_decode(result)[, `:=`(
      algorithm = "APART",
      seg.i = .I,
      start=first.i-0.5,
      end=last.i+0.5
    )][]
    changes <- nrow(asegs)-1
    pen.term <- if(penalty<Inf)changes*penalty else 0
    list(
      penmap_loss=result[.N, best.cost]-pen.term,
      penmap_size=changes,
      details=result)
  }
  pmap <- penmap::penmap_create(solver_APART, data.2d)
  penmap::penmap_solve(pmap)
  while({
    pmap.dt <- data.table(pmap$penmap$df())
    good.pen.row <- pmap.dt[size_on==N.true.changes]
    more.changes <- pmap.dt[size_on>N.true.changes][.N]
    less.changes <- pmap.dt[size_on<N.true.changes][1]
    print(h.vec <- pmap$penmap$helpful())
    new.pen <- h.vec[more.changes$penalty < h.vec & less.changes$penalty > h.vec]
    nrow(good.pen.row)==0 && length(new.pen)
  }){
    penmap::penmap_solve(pmap, new.pen)
  }
  closest <- pmap.dt[which.min(abs(size_on-N.true.changes))]
  afit <- pmap$results[[paste(closest$size_on)]]$details
  asegs <- APART_decode(afit)[, `:=`(
    algorithm = "APART",
    seg.i = .I,
    start=first.i-0.5,
    end=last.i+0.5
  )][]
  penalty <- closest$penalty
  APART.fit.list[[paste(kappa)]] <- data.table(
    kappa, penalty, afit)
  both.seg.list[[paste(kappa)]] <- data.table(
    kappa,
    parameter.name="penalty",
    parameter.value=penalty,
    asegs)
}
APART.fit <- rbindlist(APART.fit.list)
(both.seg <- rbindlist(both.seg.list, use.names=TRUE))
save(
  APART.fit, both.seg, sim.dt, true.change.vec,
  file="figure-2d-hmm-sim-noise-many-data.RData")
