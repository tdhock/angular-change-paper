library(animint2)
library(data.table)
(objs <- load("figure-2d-hmm-real-interactive-data.RData"))
in.window <- function(x){
  ((x-1) %% data.per.window)+1
}
data.per.window <- 500

ggplot()+
  geom_point(aes(
    iteration, log.lik),
    data=lik.dt)+
  facet_grid(n.states ~ ., labeller=label_both)

ggplot()+
  geom_point(aes(
    n.states, log.lik),
    data=lik.dt[, .SD[.N], by=n.states])

show.states <- 2
best.states <- viterbi.dt[, .SD[iteration==max(iteration)], by=n.states]
best.states <- viterbi.dt[n.states==show.states][iteration==max(iteration)]
show.mean <- mean.dt[n.states==show.states][iteration==max(iteration)]
rlist <- rle(as.integer(best.states$state))
infer.seg.dt <- with(rlist, data.table(
  state.i=values, n.data=lengths
))[
, end := cumsum(n.data)+0.5
][
, start := end-n.data
][]
seg.param.dt <- infer.seg.dt[
  show.mean, on="state.i", allow.cartesian=TRUE
][
, state := factor(state.i)
][]
setkey(seg.param.dt, col.i, start, end)
setkey(summary.dt, col.i, min.row, max.row)
seg.overlaps <- foverlaps(
  summary.dt, seg.param.dt, nomatch=0L
)[, `:=`(
  start.in.window=in.window(start+0.5)-0.5,
  end.in.window=in.window(end-0.5)+0.5
)][, `:=`(
  start.in.window=c(0.5, start.in.window[-1]),
  end.in.window=c(end.in.window[-.N], data.per.window+0.5)
), by=.(window, col.i)
][]

win.id <- 59
one.window <- real.dt[window==win.id]
summary.dt[window==win.id]
ggplot()+
  geom_point(aes(
    row.i, degrees),
    shape=1,
    data=one.window)+
  facet_grid(col.i ~ ., labeller=label_both)

ggplot()+
  facet_grid(col.i ~ ., labeller=label_both)+
  geom_point(aes(
    window, `50%`),
    fill="blue",
    color=NA,
    shape=21,
    data=summary.dt)+
  geom_segment(aes(
    window, `25%`,
    xend=window, yend=`75%`),
    data=summary.dt)+
  geom_segment(aes(
    window, `99%`,
    xend=window, yend=`100%`),
    data=summary.dt)+
  geom_segment(aes(
    window, `0%`,
    xend=window, yend=`1%`),
    data=summary.dt)

viz <- animint(
  title="Hidden Markov Model for 2d angular data",
  source="https://github.com/tdhock/angular-change-paper/blob/main/figure-2d-hmm-real-interactive.R",
  out.dir="figure-2d-hmm-real-interactive",
  summary=ggplot()+
    ggtitle("Summary of all windows, click to select window")+
    theme_bw()+
    theme_animint(width=1000, height=300)+
    theme(panel.margin=grid::unit(1,"lines"))+
    coord_cartesian(expand=FALSE)+
    scale_y_continuous(
      "degrees (quantiles in window)",
      limits=c(0,360),
      breaks=seq(0,360,by=90))+
    scale_x_continuous(
      "Index/row in real data sequence")+
    facet_grid(col.i ~ ., labeller=label_both)+
    geom_segment(aes(
      mid.row, `25%`,
      xend=mid.row, yend=`75%`),
      data=summary.dt)+
    geom_point(aes(
      mid.row, `50%`),
      fill=NA,
      color="blue",
      data=summary.dt)+
    geom_point(aes(
      mid.row, `0%`),
      fill="white",
      data=summary.dt)+
    geom_point(aes(
      mid.row, `100%`),
      fill="white",
      data=summary.dt)+
    geom_segment(aes(
      mid.row, `75%`,
      xend=mid.row, yend=`99%`),
      color="grey50",
      size=1,
      data=summary.dt)+
    geom_segment(aes(
      mid.row, `25%`,
      xend=mid.row, yend=`1%`),
      color="grey50",
      size=1,
      data=summary.dt)+
    geom_tallrect(aes(
      xmin=min.row, xmax=max.row),
      fill="black",
      alpha=0.5,
      clickSelects="window",
      data=summary.dt),
  zoom=ggplot()+
    ggtitle("Zoom to data in selected window")+
    theme_bw()+
    theme_animint(width=1000, height=300)+
    theme(panel.margin=grid::unit(1,"lines"))+
    coord_cartesian(expand=FALSE)+
    scale_x_continuous(
      "Index/row relative to selected window")+
    scale_y_continuous(
      "degrees (individual data points)",
      limits=c(0,360),
      breaks=seq(0,360,by=90))+
    geom_segment(aes(
      start.in.window, mean,
      xend=end.in.window, yend=mean),
      data=seg.overlaps,
      showSelected="window",
      color="red")+
    geom_point(aes(
      row.in.window, degrees),
      fill=NA,
      showSelected="window",
      data=real.dt)+
    facet_grid(col.i ~ ., labeller=label_both)
)
viz
if(FALSE){
  animint2pages(viz, "2024-07-04-HMM-angular-data")
}
