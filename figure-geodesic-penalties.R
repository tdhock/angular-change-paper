atime.result <- readRDS("figure-geodesic-penalties-data.rds")
library(data.table)
more.dt <- rbindlist(atime.result$meas$result)
atime.result$measurements <- data.table(
  atime.result$measurements,
  more.dt)
atime.refs <- atime::references_best(
  atime.result,
  more.units=names(more.dt))
png("figure-geodesic-penalties.png",width=8,height=8,units="in",res=200)
plot(atime.refs)
dev.off()
