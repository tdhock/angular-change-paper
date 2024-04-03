atime.result <- readRDS("figure-compare-time-data.rds")
atime.refs <- atime::references_best(atime.result)
png("figure-compare-time.png", width=6, height=6, units="in", res=200)
plot(atime.refs)+
  ggplot2::ggtitle("circular::change.point on two columns of trajectory 1458")
dev.off()
