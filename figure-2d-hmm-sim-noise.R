library(ggplot2)
library(data.table)
(objs=load("figure-2d-hmm-sim-noise-data.RData"))
both.long <- melt(both.seg, measure.vars=measure(dim.i=as.integer, pattern="V([12])"), value.name="mean.angle")

gg <- ggplot()+
  geom_point(aes(
    data.i, angle),
    shape=1,
    data=sim.dt)+
  facet_grid(dim.i ~ kappa, labeller=label_both)+
  geom_segment(aes(
    start, mean.angle,
    color=algorithm,
    size=algorithm,
    xend=end, yend=mean.angle),
    size=2,
    data=both.long)
png("figure-2d-hmm-sim-noise.png", width=6, height=3, units="in", res=200)
print(gg)
dev.off()
