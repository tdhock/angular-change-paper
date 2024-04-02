library(data.table)
if(FALSE){
  install.packages(c("circular","moveHMM"))
}
traj.csv.vec <- Sys.glob("avocados/*/*.csv.gz")
csv.i.vec <- seq_along(traj.csv.vec)
csv.i.vec <- 1:2
for(csv.i in csv.i.vec){
  traj.csv <- traj.csv.vec[[csv.i]]
  traj.dt <- fread(traj.csv)
  c.vec <- circular::circular(traj.dt[[1]])
  circular::change.point(c.vec)
}
