library(data.table)
if(FALSE){
  install.packages(c("circular","moveHMM"))
}
max.angle <- 2*pi
n.grid <- 360
(param.grid <- seq(0,2*pi,l=n.grid)[-n.grid])
check_angle <- function(x){
  stopifnot(0 <= x & x < max.angle)
}
geodesic_dist <- function(data.vec, param.vec){
  check_angle(data.vec)
  check_angle(param.vec)
  d <- abs(data.vec-param.vec)
  ifelse(d>max.angle/2, max.angle-d, d)
}
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

##from ?moveHMM::fitHMM
mu<-c(15,50)
sigma<-c(10,20)
angleMean <- c(pi,0)
kappa <- c(0.7,1.5)
mu0 <- c(20,70)
sigma0 <- c(10,30)
atime.result <- atime::atime(
  verbose=TRUE,
  N=as.integer(10^seq(1, 6, by=0.5)),
  setup={
    set.seed(1)
    sim.data <- moveHMM::simData(
      nbAnimals=1,
      nbStates=2,
      stepDist="gamma",
      angleDist="vm",
      stepPar=c(mu,sigma),
      anglePar=c(angleMean,kappa),
      nbCovs=2,
      zeroInflation=FALSE,
      obsPerAnimal=N)
    no.na <- sim.data[-c(1,N),]#rm missing first/last
    c.vec <- circular::circular(no.na$angle)
    pos.vec <- with(no.na, ifelse(angle<0,angle+2*pi,angle))
  },
  seconds.limit=100,
  ##seconds.limit=0.01,
  "moveHMM"={
    hmm <- moveHMM::fitHMM(
      data=sim.data,
      nbStates=2,
      stepPar0=c(mu0,sigma0),
      anglePar0=c(1,1),
      formula=~cov1+cos(cov2),
      stepDist="gamma",
      angleDist="vm",
      angleMean=angleMean)
  },
  geodesicFPOP={
    geodesichange::geodesicFPOP_vec(pos.vec, 1)
  },
  APART_grid=APART(pos.vec, log(N), param.grid),
  APART_kmeans={
    km <- kmeans(pos.vec, 2)
    two.params <- as.numeric(km$centers)
    APART(pos.vec, log(N), two.params)
  },
  result=FALSE,
  "circular"=circular::change.point(c.vec))
plot(atime.result)
saveRDS(atime.result, "figure-compare-time-sim-data.rds")
