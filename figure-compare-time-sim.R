atime.result <- readRDS("figure-compare-time-sim-data.rds")
plot(atime.result)+
  ggplot2::ggtitle("Tests on data from moveHMM::simData")
atime.refs <- atime::references_best(atime.result)
png("figure-compare-time-sim.png", width=6, height=6, units="in", res=200)
plot(atime.refs)+
  ggplot2::ggtitle("Tests on data from moveHMM::simData")
dev.off()

