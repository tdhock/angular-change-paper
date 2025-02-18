library(ggplot2)
library(data.table)
(objs=load("figure-2d-hmm-sim-noise-many-data.RData"))
show.kappa <- c(1,5)
both.long <- melt(
  both.seg,
  measure.vars=measure(
    dim.i=as.integer, pattern="V([12])"),
  value.name="mean.angle"
)[
, mean.radians := ifelse(mean.angle<0, mean.angle+2*pi, mean.angle)
][
, feature := dim.i
][
  kappa %in% show.kappa
]
sim.show <- sim.dt[
, feature := dim.i
][
  kappa %in% show.kappa
]
true.dt <- data.table(index=true.change.vec)
algo <- function(algorithm, Algorithm, color, size){
  data.table(algorithm, Algorithm, color, size)
}
algo.dt <- rbind(
  ##algo("LR","LR","00b50d",1),
  algo("APART", "APART", "4562b8", 1),
  algo("SinCosL2BinSeg", "BINSEG", "fa9919", 3),
  algo("VonMisesHMM", "HMM", "d6001f", 2))
algo.colors <- algo.dt[, structure(paste0("#",color), names=Algorithm)]
algo.sizes <- algo.dt[, structure(size, names=Algorithm)]
both.disp <- both.long[algo.dt, on="algorithm"][order(-size)]
gg <- ggplot()+
  theme_bw()+
  geom_point(aes(
    data.i, angle),
    shape=1,
    color="grey50",
    data=sim.show)+
  geom_vline(aes(
    xintercept=index),
    color="black",
    data=true.dt)+
  facet_grid(feature ~ kappa, labeller=label_both)+
  geom_segment(aes(
    start, mean.radians,
    color=Algorithm,
    size=Algorithm,
    xend=end,
    yend=mean.radians),
    data=both.disp)+
  scale_color_manual(values=algo.colors)+
  scale_size_manual(values=algo.sizes)+
  scale_x_continuous(
    "Index in data sequence")+
  scale_y_continuous(
    "Angle (radians)")
png("figure-2d-hmm-sim-noise-many.png", width=8, height=3, units="in", res=200)
print(gg)
dev.off()
