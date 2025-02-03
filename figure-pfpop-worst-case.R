library(ggplot2)
library(data.table)

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
plot_check <- function(data_value){
  (result <- pfpop::pfpop_map(data_value, Inf))
  data_dt <- data.table(data_value, data_i=seq_along(data_value)-1)
  cum_dt <- data.table(i=seq_along(data_value))[, .(
    data_value=data_value[1:i],
    data=c(rep("previous",i-1),"current")
  ), by=.(data_i=i-1)]
  gres <- geodesichange::geodesicFPOP_vec(data_value, Inf, verbose=1)
  map_dt <- mean_cost(result)
  ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_vline(ggplot2::aes(
      xintercept = x),
      size=0.5,
      color = "grey",
      data = gres$model[, .SD[, .(
        x = unique(c(min_param, max_param))
      )], by = data_i]) +
    scale_size_manual(values=c(
      previous=1,
      current=2))+
    geom_vline(aes(
      xintercept=data_value, size=data),
      data=cum_dt,
      color="grey50")+
    geom_text(aes(
      data_value, Inf, label=paste0(" ",data_value)),
      data=data_dt,
      vjust=1.1,
      hjust=0)+
    ggplot2::geom_segment(ggplot2::aes(
      min_param, min_param * Linear + Constant,
      xend = max_param, 
      yend = max_param * Linear + Constant),
      size=1.5,
      data = gres$model) + 
    ggplot2::facet_grid(data_i ~ .) +
    theme_bw()+
    facet_grid(data_i ~ ., scales="free")+
    geom_point(aes(
      param, cost_mean, color=limit),
      size=4,
      shape=21,
      fill=NA,
      data=map_dt)+
    geom_abline(aes(
      slope=Lmean, intercept=Cmean, color=limit),
      size=0.5,
      data=map_dt)+
    scale_x_continuous("Geodesic center parameter",breaks = seq(0, 360, by = 90))+
    scale_y_continuous("Mean L1 geodesic loss")+
    ggtitle("Iterations of FPOP with L1 geodesic loss and penalty=Inf")
}

if(FALSE){
atime::atime(
  setup={
    data_vec <- (seq(0, 1, l=N)^2+1)*90/2+rep(c(0,180),l=N)
  },
  foo=(result <- pfpop::pfpop_map(data_vec, Inf)),
  )
}

png("figure-pfpop-worst-case-typical.png", width=6, height=6, units="in", res=200)
gg <- plot_check(c(40, 50, 20, 80, 60, 70))
print(gg)
dev.off()

png("figure-pfpop-worst-case.png", width=6, height=6, units="in", res=200)
gg <- plot_check(c(10, 200, 40, 230, 80, 270))
print(gg)
dev.off()


