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
all_funs <- function(data_value){
  (result <- pfpop::pfpop_map(data_value, Inf))
  data_dt <- data.table(data_value, data_i=seq_along(data_value)-1)
  cum_dt <- data.table(i=seq_along(data_value))[, .(
    data_value=data_value[1:i],
    data=c(rep("previous",i-1),"current")
  ), by=.(data_i=i-1)]
  gres <- geodesichange::geodesicFPOP_vec(data_value, Inf, verbose=1)
  gres$data[, cum.weight := cumsum(end-start)]
  gres$model[
  , weight := gres$data$cum.weight[data_i+1]
  ][, let(
    Linear_total=Linear*weight,
    Constant_total=Constant*weight
  )][]
  model_long <- melt(
    melt(
      gres$model,
      measure.vars=measure(param_end, pattern="(min|max)_param"),
      value.name="param_value"),
    measure.vars=measure(value.name, cost.type, pattern="(Linear|Constant)(.*)")
  )[
  , cost_value := param_value*Linear+Constant
  ][]
  map_dt <- mean_cost(result)
  min_max_dt <- model_long[cost.type=="_total"][
  , data.table(
    .SD[c(which.min(cost_value), which.max(cost_value))],
    opt_name=c("min","max"),
    vjust = c(-0.5, 1.5))
  , by=data_i]
  range_dt <- min_max_dt[, .(diff=diff(cost_value)), by=data_i]
  list(result=result, data_dt=data_dt, cum_dt=cum_dt, gres=gres, map_dt=map_dt, model_long=model_long, min_max_dt=min_max_dt, range_dt=range_dt)
}
plot_mean <- function(data_value)with(all_funs(data_value), {
  ggplot2::ggplot() +
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
      data_value, Inf, label=paste0(" data=",data_value)),
      data=data_dt,
      vjust=1.1,
      hjust=0)+
    ggplot2::geom_segment(ggplot2::aes(
      min_param, min_param * Linear + Constant,
      xend = max_param, 
      yend = max_param * Linear + Constant),
      size=1.5,
      data = gres$model) + 
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
})

if(FALSE){
atime::atime(
  setup={
    data_vec <- (seq(0, 1, l=N)^2+1)*90/2+rep(c(0,180),l=N)
  },
  foo=(result <- pfpop::pfpop_map(data_vec, Inf)),
  )
}

png("figure-pfpop-worst-case-typical.png", width=6, height=6, units="in", res=200)
gg <- plot_mean(c(40, 50, 20, 80, 60, 70))
print(gg)
dev.off()

png("figure-pfpop-worst-case.png", width=6, height=6, units="in", res=200)
gg <- plot_mean(c(10, 200, 40, 230, 80, 270))
print(gg)
dev.off()

plot_total <- function(data_value)with(all_funs(data_value), {
  ggplot2::ggplot() +
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
      data_value, Inf, label=paste0("data=",data_value)),
      data=data_dt,
      vjust=1.1)+
    ggplot2::geom_segment(ggplot2::aes(
      min_param, min_param * Linear_total + Constant_total,
      xend = max_param,
      yend = max_param * Linear_total + Constant_total),
      size=1.5,
      data = gres$model) +
    theme_bw()+
    facet_grid(data_i ~ ., scales="free")+
    geom_point(aes(
      param, cost, color=limit),
      size=4,
      shape=21,
      fill=NA,
      data=map_dt)+
    geom_abline(aes(
      slope=Linear, intercept=Constant, color=limit),
      size=0.5,
      data=map_dt)+
    geom_label(aes(
      param_value, cost_value,
      vjust=vjust,
      label=paste0(opt_name,"=",round(cost_value,2))),
      alpha=0.5,
      data=min_max_dt,
      color="red")+
    geom_point(aes(
      param_value, cost_value),
      data=min_max_dt,
      color="red")+
    geom_text(aes(
      -Inf, Inf, label=paste0("max-min=",round(diff,2)," ")),
      data=range_dt,
      vjust=1.1,
      hjust=1,
      angle=90,
      color="red")+
    scale_x_continuous("Geodesic center parameter",breaks = seq(0, 360, by = 90))+
    scale_y_continuous("Total L1 geodesic loss")+
    ggtitle("Iterations of FPOP with L1 geodesic loss and penalty=Inf")
})

plot_new <- function(data_values, new_values){
  fun5.dt.list <- list()
  min.dt.list <- list()
  for(last_value in new_values){
    ares <- all_funs(c(data_values, last_value))
    fun5.dt.list[[paste(last_value)]] <- data.table(
      last_value, ares$gres$model[data_i==max(data_i)])
    min.dt.list[[paste(last_value)]] <- data.table(
      last_value, ares$min_max_dt[data_i==max(data_i) & opt_name=="min"])
  }
  (fun5.dt <- rbindlist(fun5.dt.list))
  (min.dt <- rbindlist(min.dt.list))
  some.dt <- ares$min_max_dt[data_i==max(data_i)-1 & opt_name=="min"]
  ggplot2::ggplot() +
    geom_vline(aes(
      xintercept=param_value),
      data=some.dt,
      color="grey")+
    ggplot2::geom_segment(ggplot2::aes(
      min_param, min_param * Linear_total + Constant_total,
      xend = max_param,
      yend = max_param * Linear_total + Constant_total),
      data = fun5.dt) +
    theme_bw()+
    facet_wrap("last_value")+
    scale_x_continuous("Geodesic center parameter",breaks = seq(0, 360, by = 90))+
    scale_y_continuous("Total L1 geodesic loss")+
    ggtitle("Iterations of FPOP with L1 geodesic loss and penalty=Inf")+
    geom_point(aes(
      param_value, cost_value),
      data=min.dt,
      color="red")
}

initial.vec <- c(270)
plot_total(initial.vec)
plot_new(initial.vec, seq(0,350, by=50))

## diff 200 one local min
initial.vec <- c(270, 290, 80)
plot_total(initial.vec)
plot_new(initial.vec, seq(0,270, by=45))

## diff=170 ???
initial.vec <- c(270, 290, 105)
plot_total(initial.vec)
plot_new(initial.vec, c(75,95,115))

## no new local min
initial.vec <- c(270, 290, 300, 280)
plot_total(initial.vec)
plot_new(initial.vec, seq(0,315, by=45))

initial.vec <- c(200, 310)
plot_total(initial.vec)
plot_new(initial.vec, c(90,180,270,340))

initial.vec <- c(200, 310, 90)
plot_total(initial.vec)
plot_new(initial.vec, c(80, 255, 300))

initial.vec <- c(200, 310, 90, 255, 240)
plot_total(initial.vec)

initial.vec <- c(270, 310, 50)
plot_total(initial.vec)
plot_new(initial.vec, c(90,95))

## something like the worst case?
plot_total(c(10, 200, 40, 230, 80, 270))

plot_total(c(10, 200, 40, 230, 35, 270))

plot_total(c(10, 100, 200, 250))

