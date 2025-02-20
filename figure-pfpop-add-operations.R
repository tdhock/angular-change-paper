library(ggplot2)
library(data.table)

pfpop_map_verbose <- function(degrees_vec, penalty=Inf, weight_vec = rep(1, length(degrees_vec)), verbose_file=tempfile()){
  result <- pfpop::pfpop_map(degrees_vec, penalty, weight_vec, verbose_file)
  result$clusters <- fread(verbose_file)
  result
}

mean_cost <- function(result)melt(
  data.table(result$iterations)[, data_i := .I-1],
  measure.vars=measure(limit, value.name, pattern="(min|max)_(.*)")
)[
, N := data_i+1
][, let(
  Lmean=Linear/N,
  Cmean=Constant/N,
  cost_mean=cost/N
)][]
cldt <- function(data_i, opt, start, end){
  sedt <- data.table(
    data_i,
    start=as.numeric(start),
    end=as.numeric(end)
  )[
    ##start != end
  ][
  , noInf := start <= end
  ]
  data.table(opt, rbind(
    sedt[noInf==TRUE],
    sedt[noInf==FALSE][, let(start=-Inf)],
    sedt[noInf==FALSE][, let(end=Inf)]))
}
plot_check <- function(tit, suffix, ...){
  data_vec <- round(c(...),1)
  (result <- pfpop_map_verbose(data_vec))
  gres <- geodesichange::geodesicFPOP_vec(data_vec, Inf, verbose=1)
  gres$data[, let(
    Value=radians,
    data_i=start)]
  map_dt <- mean_cost(result)
  cluster.dt <- result$clusters[, rbind(
    cldt(data_i, "before", first_param, opt_param),
    cldt(data_i, "after", opt_param, last_param))]
  lab.dt <- melt(
    result$clusters,
    measure.vars=measure(
      pointer, value.name, pattern="(first|opt|last)_(param|diff)"))
  add_meta <- function(DT){
    value <- gres$data[DT, radians, on="data_i"]
    DT[, Value := value]
  }
  add_meta(lab.dt)
  add_meta(cluster.dt)
  add_meta(map_dt)
  add_meta(gres$model)
  gg <- ggplot2::ggplot() +
    ggtitle(tit)+
    ggplot2::theme_bw() +
    ggplot2::geom_vline(ggplot2::aes(
      xintercept = x), 
      color = "grey", data = gres$model[
      , .SD[, .(x = unique(c(min_param, max_param)))]
      , by = .(data_i,Value)]) +
    ggplot2::geom_segment(ggplot2::aes(
      min_param, 
      min_param * Linear + Constant, xend = max_param, 
      yend = max_param * Linear + Constant),
      data = gres$model) + 
    theme_bw()+
    facet_grid(data_i + Value ~ ., scales="free", labeller=label_both)+
    geom_rect(aes(
      xmin=start, xmax=end,
      fill=opt,
      ymin=-Inf, ymax=Inf),
      data=cluster.dt,
      alpha=0.5,
      color="black")+
    geom_label(aes(
      param, ifelse(diff<0, -Inf, Inf),
      vjust=ifelse(diff<0, 0, 1),
      label=diff),
      data=lab.dt,
      alpha=0.5)+
    scale_color_discrete(
      "global")+
    scale_fill_manual(
      "Position\nrelative\nto cluster\noptimum", values=c(
      before="grey50",
      after="violet"))+
    geom_point(aes(
      param, cost_mean, color=limit),
      size=4,
      shape=21,
      fill=NA,
      data=map_dt)+
    geom_abline(aes(
      slope=Lmean, intercept=Cmean, color=limit),
      data=map_dt)+
    geom_point(aes(
      Value, -Inf),
      data=gres$data)+
    scale_x_continuous(
      "Geodesic center parameter",
      breaks=seq(0,360,by=90))+
    scale_y_continuous(
      "Mean loss")
  f.png <- sprintf(
    "figure-pfpop-add-operations-%s.png",
    suffix)
  png(f.png, width=6, height=10, units="in", res=200)
  print(gg)
  dev.off()
}

data_vec <- c(30, 20, 335, 10, 325, 340, 330, 320, 310, 300)
data_vec <- c(10, 200, 40, 50, 60, 70, 205, 210, 215, 220)
data_vec <- c(10, 20,    30,    40,50,    60, 70, 80, 90, 100)#best case
data_vec <- c(10, 20+180,40,50+180,80,90+180)#worst case
N <- 10
data_vec <- seq(0,90,l=N)+rep(c(0,180),l=N)

plot_check(
  "Base case: constant number of clusters of breakpoints (2)",
  "best",
  10, 20,    30,    40,50,    60, 70, 80, 90, 100)

plot_check(
  "Worst case: linear number of clusters of breakpoints",
  "worst",
  10, 20+180,40,50+180,80,90+180)

set.seed(1)
plot_check(
  "Typical data: 3 clusters of data values, 6 clusters of breakpoints",
  "typical",
  rnorm(3,60,10),
  rnorm(4, 180,10),
  rnorm(5, 300,10),
  ##rnorm(3,60,10),
  NULL)

plot_check(
  "Cluster which is not a local optimum",
  "crit",
  50, 190, 60, 200, 300, 310, 320, 330, 340,
  ##55, 45, 65, 70, 40,
  NULL)
