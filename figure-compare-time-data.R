library(data.table)
if(FALSE){
  install.packages(c("circular","moveHMM"))
}
(traj.csv.vec <- Sys.glob("avocados/*/*.csv.gz"))
csv.i.vec <- seq_along(traj.csv.vec)
csv.i.vec <- 1:2
traj.list <- list()
expr.list <- list()
for(csv.i in csv.i.vec){
  traj.csv <- traj.csv.vec[[csv.i]]
  traj.dt <- fread(traj.csv)
  traj.vec <- traj.dt[[1]]
  traj.list[[csv.i]] <- circular::circular(traj.vec)
  expr.list[[paste0("col",csv.i)]] <- substitute(
    circular::change.point(sub.list[[I]]),
    list(I=csv.i))
}

(N.seq <- as.integer(10^seq(1, log10(length(traj.vec)), by=0.5)))
atime.result <- atime::atime(
  N=N.seq,
  setup={
    sub.list <- lapply(traj.list, "[", 1:N)
  },
  seconds.limit=1,
  expr.list=expr.list)
saveRDS(atime.result, "figure-compare-time-data.rds")

