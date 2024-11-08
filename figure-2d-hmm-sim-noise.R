library(ggplot2)
library(data.table)
(objs=load("figure-2d-hmm-sim-noise-data.RData"))
both.long <- melt(both.seg, measure.vars=measure(dim.i=as.integer, pattern="V([12])"), value.name="mean.angle")[, mean.radians := ifelse(mean.angle<0, mean.angle+2*pi, mean.angle)][, feature := dim.i]
sim.dt[, feature := dim.i]
true.dt <- data.table(index=c(100,200)+0.5)
gg <- ggplot()+
  theme_bw()+
  geom_point(aes(
    data.i, angle),
    shape=1,
    color="grey50",
    data=sim.dt)+
  geom_vline(aes(
    xintercept=index),
    color="black",
    data=true.dt)+
  facet_grid(feature ~ kappa, labeller=label_both)+
  geom_segment(aes(
    start, mean.radians,
    color=algorithm,
    size=algorithm,
    xend=end, yend=mean.radians),
    data=both.long)+
  scale_size_manual(values=c(
    APART=1,
    SinCosL2BinSeg=3,
    VonMisesHMM=2))+
  scale_x_continuous(
    "Index in data sequence")+
  scale_y_continuous(
    "Angle (radians)")
png("figure-2d-hmm-sim-noise.png", width=8, height=3, units="in", res=200)
print(gg)
dev.off()
