fun_list <- list(
  best=function(N)rnorm(N,100),
  worst=function(N)seq(0,90,l=N)+rep(c(0,180),l=N))
(expr.list <- atime::atime_grid(
  list(Data=names(fun_list)),
  map={
    mres <- pfpop::pfpop_map(data_list[[Data]], Inf)
    with(mres$iterations, data.frame(
      cost=min_cost[N]/N,
      ##space=max(map_size),
      pointers=max(list_size),
      moves=sum(num_moves)))
  },
  list={
    lres <- pfpop::pfpop_list(data_list[[Data]], Inf)
    data.frame(
      cost=lres$iterations$best_cost[N],
      ##space=max(lres$iterations$num_pieces),
      pointers=NA,
      moves=NA)
  }))
map_list_atime <- atime::atime(
  setup={
    set.seed(1)
    data_list <- list()
    for(Data in names(fun_list)){
      fun <- fun_list[[Data]]
      data_list[[Data]] <- fun(N)
    }
  },
  seconds.limit=0.1,
  result=TRUE,
  expr.list=expr.list)

png("figure-pfpop-atime-compare.png", width=5, height=4, units="in", res=200)
gg <- plot(map_list_atime)
print(gg)
dev.off()

library(data.table)
dcast(
  map_list_atime$measurements,
  N + Data ~ expr.grid,
  value.var="cost"
)[
, cost.diff := list - map
][order(-abs(cost.diff))]
## biggest difference is ~1e-14

map_list_refs <- atime::references_best(map_list_atime)

png("figure-pfpop-atime.png", width=5, height=6, units="in", res=200)
gg <- plot(map_list_refs)
print(gg)
dev.off()
