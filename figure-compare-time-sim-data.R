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
atime.result <- atime::atime(
  N=10^seq(1, 1e6, by=0.5),
  setup={
    sim.data <- moveHMM::simData(
      nbAnimals=1,
      nbStates=2,
      stepDist="gamma",
      angleDist="vm",
      stepPar=c(mu,sigma),
      anglePar=c(angleMean,kappa),
      nbCovs=2,
      zeroInflation=FALSE,
      obsPerAnimal=N+2)[-c(1,N+2),]#rm missing first/last
    c.vec <- circular::circular(sim.data$angle)
  },
  seconds.limit=100,
  "moveHMM::fitHMM"={
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
  "circular::change.point"=circular::change.point(c.vec))
saveRDS(atime.result, "figure-compare-time-sim-data.rds")
