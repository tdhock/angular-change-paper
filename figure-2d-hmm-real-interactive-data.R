library(data.table)
col.i.show <- 1010
col.i.show <- 0
real.dt.list <- list()
data.per.window <- 500
col.i.vec <- 0:1 + col.i.show
for(col.i in col.i.vec){
  col.i.csv.gz <- file.path(
    "avocados",
    "trajectory1458",
    paste0(col.i,".csv.gz"))
  col.i.dt <- fread(file=col.i.csv.gz, col.names="radians")
  row.i <- 1:nrow(col.i.dt)
  in.window <- function(x){
    ((x-1) %% data.per.window)+1
  }
  row.in.window <- in.window(row.i)
  window <- ((row.i-1) %/% data.per.window) + 1
  real.dt.list[[paste(col.i)]] <- data.table(
    col.i, row.i, row.in.window, window, col.i.dt
  )[, degrees := 360*radians/(2*pi)]
}
(real.dt <- rbindlist(real.dt.list)[row.i<=1e5])
range(real.dt$degrees)
(summary.dt <- real.dt[, {
  L <- as.list(quantile(degrees, c(0,0.01,0.25,0.5,0.75,1,0.99)))
  L$min.row <- min(row.i)
  L$max.row <- max(row.i)
  L
}, by=.(col.i, window)
][
, mid.row := (min.row+max.row)/2
][])

prob.dt.list <- list()
lik.dt.list <- list()
trans.dt.list <- list()
mean.dt.list <- list()
viterbi.dt.list <- list()
prior.dt.list <- list()
max.states <- 5
for(n.states in 2:max.states){
  cat(sprintf("%4d / %4d states\n", n.states, max.states))
  log.A.mat <- log(matrix(1/n.states, n.states, n.states))
  data.wide <- dcast(real.dt, row.i ~ col.i, value.var="degrees")
  max.angle <- 360
  y.mat <- as.matrix(data.wide[,-1])
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
  done <- FALSE
  prev.log.lik <- -Inf
  iteration <- 1
  while(!done){
    it.time.vec <- system.time({
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
    })#time
    viterbi.dt.list[[paste(n.states, iteration)]] <- data.table(
      n.states, iteration,
      data.i=1:N.data,
      state=factor(viterbi.result[["state_seq"]]))
    mean.dt.list[[paste(n.states, iteration)]] <- data.table(
      n.states, iteration,
      state.i=as.integer(row(mean.mat)),
      col.i=col.i.vec[as.integer(col(mean.mat))],
      mean=as.numeric(mean.mat))
    lik.dt.list[[paste(n.states, iteration)]] <- data.table(
      n.states, iteration,
      seconds=it.time.vec[["elapsed"]],
      log.lik)
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
}
viterbi.dt <- do.call(rbind, viterbi.dt.list)
mean.dt <- do.call(rbind, mean.dt.list)
lik.dt <- do.call(rbind, lik.dt.list)

save(
  viterbi.dt, mean.dt, lik.dt, real.dt, summary.dt,
  file="figure-2d-hmm-real-interactive-data.RData")
