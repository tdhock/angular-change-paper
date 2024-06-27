atime.result <- readRDS("figure-geodesicFPOP-robseg-data.rds")
library(data.table)
more.dt <- rbindlist(atime.result$meas$result)
atime.result$unit.col.vec <- c(seconds="median", "kilobytes", names(more.dt))
atime.result$measurements <- data.table(
  atime.result$measurements,
  more.dt)
atime.refs <- atime::references_best(atime.result)
png("figure-geodesicFPOP-robseg.png",width=18,height=8,units="in",res=200)
plot(atime.refs)
dev.off()

p <- function(base_intervals,loss){
  data.table(base_intervals,loss)
}
piece.meta <- rbind(
  p("2-3","geodesic"),
  p("2","L1"),
  p("3","L2rob"),
  p("1","L2"),
  p("1","Poisson")
)[
, loss.pieces := sprintf("%s(%s)",loss,base_intervals)
][
  order(base_intervals)
][
, `Loss(Pieces)` := factor(loss.pieces,loss.pieces)
][]
some.meas <- piece.meta[
  nc::capture_first_df(
    atime.refs$meas[unit%in%c("mean.intervals","seconds")],
    expr.name=list(
      loss=".*?",
      " ",
      nc::field("penalty","=","[0-9.]+",as.numeric))),
  on="loss"
][, `:=`(
  Penalty = paste0(
    ifelse(penalty==1000,"Worst case (slow) for losses with >1 piece\nAverage case (fast) for losses with 1 piece","Best case (all losses fast)"),
    "\nPenalty=",penalty,
    " Changes=",ifelse(penalty==1000,"none","O(N)")
  ),
  Unit = ifelse(unit=="mean.intervals","intervals",unit)
)][]
library(ggplot2)
lab <- function(x,y,label,hjust,vjust,Unit,penalty){
  data.table(x,y,label,hjust,vjust,Unit,penalty)
}
pen.meta <- unique(some.meas[,.(penalty,Penalty)])
lab.dt <- pen.meta[rbind(
  lab(1e4, 10, "O(1) intervals", 0.5, 0,"intervals",0.01),
  lab(4e4, 1e4, "O(N) intervals\ngeodesic/L2rob/L1",0,1,"intervals",1000),
  lab(1e7, 30, "O(log N) intervals\nPoisson/L2",1,0,"intervals",1000),
  lab(1e6, 20, "O(N) time", 0.5,0,"seconds",0.01),
  lab(5e4, 20, "O(N^2) time\ngeodesic/L2rob/L1",1,0,"seconds",1000),
  lab(2e5, 20, "O(N log N) time\nPoisson/L2", 0, 0, "seconds",1000)),
  on="penalty"]
blank.dt <- data.table(
  x=1e6,
  y=4e2,
  Unit="seconds")
ref.dt <- pen.meta[data.table(
  seconds=c(1),
  penalty=0.01,
  Unit="seconds",
  label=c("1 second")),
  on="penalty"]
RColorBrewer::display.brewer.all()
dput(RColorBrewer::brewer.pal(Inf,"Dark2"))
loss.colors <- c(
  geodesic="#1B9E77",
  L1="#D95F02",
  L2="#7570B3",
  L2rob="#E7298A",
  Poisson="#66A61E", "#E6AB02", 
  "#A6761D", "#666666"
)[piece.meta$loss]
names(loss.colors) <- piece.meta[["Loss(Pieces)"]]
gg <- ggplot()+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  geom_hline(aes(
    yintercept=seconds),
    data=ref.dt[,.(seconds,Unit)],
    color="grey")+
  geom_text(aes(
    0,seconds,label=paste0(" ",label)),
    hjust=0,
    color="grey50",
    vjust=-0.5,
    size=2.5,
    data=ref.dt)+
  geom_text(aes(
    x,y,label=label,hjust=hjust,vjust=vjust),
    size=2.5,
    data=lab.dt)+
  geom_blank(aes(x,y),data=blank.dt)+
  geom_ribbon(aes(
    N,
    ymin=min,
    ymax=max,
    fill=`Loss(Pieces)`),
    alpha=0.5,
    data=some.meas[Unit=="seconds"])+
  scale_color_manual(values=loss.colors)+
  scale_fill_manual(values=loss.colors)+
  geom_line(aes(
    N, empirical, color=`Loss(Pieces)`),
    data=some.meas)+
  facet_grid(Unit ~ Penalty, scales="free")+
  scale_x_log10(
    "T = Number of simulated data")+
  scale_y_log10("")
png("figure-geodesicFPOP-robseg-simple.png",width=7,height=3,units="in",res=300)
print(gg)
dev.off()

some.meas[, `:=`(
  Penalty = factor(penalty),
  Base_intervals=paste0("Base intervals=",base_intervals),
  Loss=paste0("Loss=",loss)
)]
ggplot()+
  theme_bw()+
  geom_ribbon(aes(
    N,
    ymin=min,
    ymax=max,
    fill=Penalty),
    alpha=0.5,
    data=some.meas[Unit=="seconds"])+
  geom_line(aes(
    N, empirical, color=Penalty),
    data=some.meas)+
  facet_grid(Unit ~ Base_intervals + Loss, scales="free")+
  scale_x_log10()+
  scale_y_log10()
