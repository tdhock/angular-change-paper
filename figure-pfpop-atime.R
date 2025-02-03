map_list_atime <- atime::atime(
  setup={
    set.seed(1)
    data_vec <- runif(N,0,360)
  },
  map={
    mres <- pfpop::pfpop_map(data_vec, Inf)
    with(mres$iterations, data.frame(
      cost=min_cost[N]/N,
      space=max(map_size),
      moves=sum(pointer_moves)))
  },
  list={
    lres <- pfpop::pfpop_list(data_vec, Inf)
    data.frame(
      cost=lres$iterations$best_cost[N],
      space=max(lres$iterations$num_pieces),
      moves=NA)
  },
  seconds.limit=0.1,
  result=TRUE
)

png("figure-pfpop-atime-compare.png", width=4, height=4, units="in", res=200)
gg <- plot(map_list_atime)
print(gg)
dev.off()


library(data.table)
dcast(
  map_list_atime$measurements,
  N ~ expr.name,
  value.var="cost"
)[list<map]

map_list_refs <- atime::references_best(map_list_atime)

png("figure-pfpop-atime.png", width=5, height=6, units="in", res=200)
gg <- plot(map_list_refs)
print(gg)
dev.off()
