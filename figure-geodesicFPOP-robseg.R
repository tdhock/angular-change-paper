atime.result <- readRDS("figure-geodesicFPOP-robseg-data.rds")
library(data.table)
more.dt <- rbindlist(atime.result$meas$result)
atime.result$measurements <- data.table(
  atime.result$measurements,
  more.dt)
atime.refs <- atime::references_best(
  atime.result,
  more.units=names(more.dt))
png("figure-geodesicFPOP-robseg.png",width=18,height=8,units="in",res=200)
plot(atime.refs)
dev.off()
