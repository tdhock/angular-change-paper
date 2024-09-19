library(data.table)
library(ggplot2)
(objs <- load("figure-2d-hmm-sim-data.RData"))
selected.states <- 3
best.mean <- mean.dt[n.states==selected.states][iteration==max(iteration)]
best.states <- viterbi.dt[n.states==selected.states][iteration==max(iteration)]
rlist <- rle(as.integer(best.states$state))
infer.seg.dt <- with(rlist, data.table(
  state.i=values, n.data=lengths, parameter="HMM"
))[
, end := cumsum(n.data)+0.5
][
, start := end-n.data
][]
seg.param.dt <- infer.seg.dt[
  best.mean, on="state.i", allow.cartesian=TRUE
][
, state := factor(state.i)
][]
tall.segs <- melt(
  APART.seg[penalty==100],
  measure.vars=measure(dim.i=as.integer, pattern="V(.*)"),
  value.name="mean.angle")
mean.all <- rbind(
  tall.segs[, .(parameter, start, end, dim.i, mean.angle)],
  true.mean.dt[, .(parameter, start, end, dim.i, mean.angle)],
  seg.param.dt[, .(parameter, start, end, dim.i=col.i, mean.angle=mean)])
change.all <- mean.all[start>1]

gg <- ggplot()+
  theme_bw()+
  geom_point(aes(
    data.i, angle),
    shape=1,
    data=sim.dt)+
  facet_grid(dim.i ~ ., labeller=label_both)+
  geom_segment(aes(
    start, mean.angle,
    color=parameter,
    size=parameter,
    xend=end, yend=mean.angle),
    data=mean.all)+
  scale_size_manual(values=c(
    true=0.5,
    grid=1,
    HMM=2,
    APART=2))
png("figure-2d-hmm-sim.png", width=8, height=3, units="in", res=200)
print(gg)
dev.off()
