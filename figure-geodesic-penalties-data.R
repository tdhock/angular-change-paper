library(data.table)
if(FALSE){
  install.packages(c("circular","moveHMM"))
}
##from ?moveHMM::fitHMM
mu<-c(15,50)
sigma<-c(10,20)
angleMean <- c(pi,0)
kappa <- c(0.7,1.5)
mu0 <- c(20,70)
sigma0 <- c(10,30)
expr.list <- atime::atime_grid(
  list(penalty=10^seq(-1,2)),
  geodesic={
    fit <- geodesichange::geodesicFPOP_vec(pos.vec, penalty)
    with(fit$loss, data.table(
      segments=as.numeric(segments),
      mean.intervals,
      max.intervals=as.numeric(max.intervals)))
  })
atime.result <- atime::atime(
  verbose=TRUE,
  N=2^seq(2,100),
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
  seconds.limit=1,
  expr.list=expr.list,
  result=TRUE)
saveRDS(atime.result, "figure-geodesic-penalties-data.rds")
