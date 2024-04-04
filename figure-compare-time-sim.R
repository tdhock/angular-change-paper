atime.result <- readRDS("figure-compare-time-sim-data.rds")
plot(atime.result)+
  ggplot2::ggtitle("Tests on data from moveHMM::simData")
atime.refs <- atime::references_best(atime.result)
atime.pred <- predict(atime.refs, seconds=atime.result$seconds.limit, kilobytes=8e6)
png("figure-compare-time-sim-pred.png", width=6, height=6, units="in", res=200)
plot(atime.pred)+
  ggplot2::ggtitle("Tests on data from moveHMM::simData")+
  ggplot2::scale_x_log10("N", limits=c(NA, 5e5))
dev.off()
png("figure-compare-time-sim.png", width=6, height=6, units="in", res=200)
plot(atime.refs)+
  ggplot2::ggtitle("Tests on data from moveHMM::simData")
dev.off()

