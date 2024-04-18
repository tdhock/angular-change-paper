atime.result <- readRDS("figure-compare-time-sim-data.rds")
library(atime)
library(data.table)
library(ggplot2)
atime.refs <- atime::references_best(atime.result)
atime.pred <- predict(atime.refs, seconds=atime.result$seconds.limit, kilobytes=8e6)
png("figure-compare-time-sim-pred.png", width=6, height=6, units="in", res=200)
plot(atime.pred)+
  ggplot2::ggtitle("Tests on data from moveHMM::simData")+
  ggplot2::scale_x_log10("N", limits=c(NA, 5e5))
dev.off()
png("figure-compare-time-sim-nopred.png", width=4, height=3, units="in", res=300)
ref.dt <- data.table(
  unit="seconds",
  seconds=c(1,60),
  label=c("1 second","1 minute"))
gg <- ggplot2::ggplot() +
  ggplot2::theme_bw() +
  ggplot2::geom_hline(ggplot2::aes(
    yintercept=seconds),
    color="grey",
    data=ref.dt)+
  ggplot2::geom_text(ggplot2::aes(
    0, seconds, label=label),
    vjust=-0.5,
    size=3,
    hjust=0,
    color="grey50",
    data=ref.dt)+
  ggplot2::geom_ribbon(ggplot2::aes(
    N, 
    ymin = min, ymax = max, fill = expr.name),
    data = data.table(atime.result$meas, unit = "seconds"),
    alpha = 0.5) +
  ggplot2::geom_line(ggplot2::aes(
    N, 
    median, color = expr.name),
    data = atime.refs$meas) +
  ggplot2::facet_grid(
    unit ~ ., scales = "free") +
  ggplot2::ggtitle("Tests on data from moveHMM::simData")+
  ggplot2::scale_x_log10(
    "Number of data to segment",
    breaks=10^seq(2,6,by=2))+
  coord_cartesian(
    ylim=c(1e-3,1e4),
    xlim=10^c(2, 8.5))+
  ggplot2::scale_y_log10("median line, min/max band")
dl <- directlabels::direct.label(gg, list(cex=0.5,"right.polygons"))
print(dl)
dev.off()
png("figure-compare-time-sim-nopred-nomem.png", width=5, height=3, units="in", res=300)
ref.dt <- data.table(
  unit="seconds",
  seconds=c(1,60),
  label=c("1 second","1 minute"))
gg <- ggplot2::ggplot() +
  ggplot2::theme_bw() +
  ggplot2::geom_hline(ggplot2::aes(
    yintercept=seconds),
    color="grey",
    data=ref.dt)+
  ggplot2::geom_text(ggplot2::aes(
    0, seconds, label=paste0(" ",label)),
    vjust=-0.5,
    size=3,
    hjust=0,
    color="grey50",
    data=ref.dt)+
  ggplot2::geom_ribbon(ggplot2::aes(
    N, 
    ymin = min, ymax = max, fill = expr.name),
    data = data.table(atime.result$meas, unit = "seconds"),
    alpha = 0.5) +
  ggplot2::geom_line(ggplot2::aes(
    N, 
    median, color = expr.name),
    data = atime.result$meas) +
  ggplot2::scale_x_log10(
    "N = Number of data to segment",
    breaks=10^seq(2,6,by=2))+
  coord_cartesian(
    ylim=c(1e-3,1e3),
    xlim=10^c(2, 8.4))+
  ggplot2::scale_y_log10(
    "Computation time (seconds)\nmin/median/max over 10 timings",
    breaks=10^seq(-3,3))
dl <- directlabels::direct.label(gg, list(cex=0.7,"right.polygons"))
print(dl)
dev.off()
png("figure-compare-time-sim.png", width=12, height=6, units="in", res=200)
plot(atime.refs)+
  ggplot2::ggtitle("Tests on data from moveHMM::simData")
dev.off()


