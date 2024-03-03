library(data.table)
library(ggplot2)
N.grid <- 51
param.grid <- seq(0, 2*pi, l=N.grid)
get.loss <- function(data.val, param.val){
  d <- abs(data.val-param.val)
  ifelse(d>pi, 2*pi-d, d)
}
set.seed(1)
(data.mat <- matrix(1:4, 2))
mean1.mat <- matrix(data.mat[1,],2,2,byrow=TRUE)
mean2.mat <- matrix(data.mat[2,],2,2,byrow=TRUE)
loss.dt.list <- list()
loss.sf.list <- list()
for(data.i in 1:nrow(data.mat)){
  data.row <- data.mat[data.i,]
  data.vals <- paste(data.row, collapse=",")
  iseq <- 1:N.grid
  pair.dt <- CJ(i1=iseq, i2=iseq)
  piece.dt.list <- list()
  for(dim.i in seq_along(data.row)){
    iname <- paste0("i",dim.i)
    data.val <- data.row[[dim.i]]
    piece <- function(Linear, Constant, min_param, max_param){
      data.table(iname, Linear, Constant, min_param, max_param)
    }
    piece.dt.list[[iname]] <- if(data.val == 0)rbind(
      piece(1, 0, 0, pi),
      piece(-1, 2*pi, pi, 2*pi)
    )else if(data.val < pi)rbind(
      piece(-1, data.val, 0, data.val),
      piece(1, -data.val, data.val, data.val+pi),
      piece(-1, (2*pi+data.val), data.val+pi, 2*pi)
    )else if(data.val == pi)rbind(
      piece(-1, pi, 0, pi),
      piece(1, -pi, pi, 2*pi)
    )else rbind(
      piece(1, 2*pi-data.val, 0, data.val-pi),
      piece(-1, data.val, data.val-pi, data.val),
      piece(1, -data.val, data.val, 2*pi)
    )
    dim.dt <- data.table(
    )[
    , paste0("loss",dim.i) := get.loss(data.val, param.grid)
    ][
    , paste0("param",dim.i) := param.grid
    ][
    , (iname) := iseq
    ]
    pair.dt <- pair.dt[dim.dt, on=iname]
  }
  dim.list <- lapply(piece.dt.list, function(DT)1:nrow(DT))
  index.dt <- do.call(CJ, dim.list)
  poly.sf.list <- list()
  for(polygon.i in 1:nrow(index.dt)){
    index.row <- index.dt[polygon.i]
    irow.list <- list()
    for(iname in names(index.row)){
      irow <- index.row[[iname]]
      irow.list[[iname]] <- piece.dt.list[[iname]][irow]
    }
    poly.sf.list[[polygon.i]] <- with(irow.list, data.table(
      data.i,
      Constant=i1$Constant+i2$Constant,
      Linear1=i1$Linear,
      Linear2=i2$Linear,
      geometry=sf::st_sfc(sf::st_polygon(list(rbind(
        c(i1$min_param, i2$min_param),
        c(i1$max_param, i2$min_param),
        c(i1$max_param, i2$max_param),
        c(i1$min_param, i2$max_param),
        c(i1$min_param, i2$min_param)
      ))))
    ))
  }
  (poly.sf <- do.call(rbind, poly.sf.list))
  loss.sf.list[[data.i]] <- poly.sf
  loss.dt.list[[data.i]] <- data.table(
    data.i, data.vals, pair.dt[, loss := loss1+loss2])
}
(loss.dt <- rbindlist(loss.dt.list))
loss.wide <- dcast(
  loss.dt[, loss.i := paste0("loss.data",data.i)],
  param1+param2 ~ loss.i,
  value.var="loss")

s <- function(step, title){
  data.table(step, title)
}
step.info <- rbind(
  s(1, "loss.data1"),
  s(2, "min.data1"),
  s(3, "loss.data2"),
  s(4, "cost.data2"))
last.cost.dt.list <- list()
last.best.dt.list <- list()
last.tile.dt.list <- list()
for(penalty in c(0.5, 5)){
  Penalty <- paste0("Penalty=",penalty)
  loss.sf.current <- data.table(
    loss.sf.list[[1]]
  )[, `:=`(
    loss = Constant+Linear1*1+Linear2*3,
    polygon.i = .I
  )]
  ggplot()+
    geom_tile(aes(
      param1, param2, fill=loss),
      data=loss.dt[data.i==2])+
    scale_fill_gradient(low="white", high="red")+
    geom_sf(aes(
      geometry=geometry),
      fill=NA,
      data=loss.sf.current)+
    geom_sf_text(aes(
      geometry=geometry, label=polygon.i),
      data=loss.sf.current)
  min1.viz.list <- list()
  for(polygon.i in 1:nrow(loss.sf.current)){
    loss.sf.piece <- loss.sf.current[polygon.i,]
    intercept <- (penalty-loss.sf.piece$Constant)/loss.sf.piece$Linear2
    slope <- -loss.sf.piece$Linear1/loss.sf.piece$Linear2
    roots.line <- sf::st_linestring(rbind(
      c(0,intercept),
      c(2*pi, intercept+slope*2*pi)
    ))
    loss.sf.geom <- loss.sf.piece$geometry[[1]]
    split.collection <- lwgeom::st_split(loss.sf.geom, roots.line)
    cmat <- sapply(split.collection, sf::st_centroid)
    change.is.better <- loss.sf.piece[
    , Constant+cmat[1,]*Linear1+cmat[2,]*Linear2>penalty]
    min1.viz.list[[polygon.i]] <- loss.sf.piece[, data.table(
      Constant, Linear1, Linear2,
      geometry=sf::st_sfc(split.collection[]),
      polygon.i, change=FALSE
    )][
      change.is.better,
      `:=`(Constant=penalty, Linear1=0, Linear2=0, change=TRUE)
    ]
    ## TODO visualize these lines on top of heat map.
    ##https://stackoverflow.com/questions/63815365/how-can-i-split-clip-a-polygon-by-lines-in-r
    ##https://r-spatial.github.io/lwgeom/
  }
  (min1.viz <- rbindlist(min1.viz.list))
  gg <- ggplot()+
    geom_sf(aes(
      fill=change,
      geometry=geometry),
      data=min1.viz)+
    geom_sf_text(aes(
      label=polygon.i,
      geometry=geometry),
      data=min1.viz)+
    coord_sf(xlim=c(0,2*pi),ylim=c(0,2*pi))    
  mySimplify <- function(DT){
    DT[, .(
      geometry=sf::st_sfc(sf::st_simplify(sf::st_multipolygon(geometry)))
    ), by=.(Constant, Linear1, Linear2, change)][, geometry := sf::st_sfc(geometry)]
  }
  min1.simple <- mySimplify(min1.viz)
  ggplot()+
    geom_sf(aes(
      fill=change,
      geometry=geometry),
      data=min1.simple)+
    coord_sf(xlim=c(0,2*pi),ylim=c(0,2*pi))
  last.cost.dt.list[[paste(Penalty,"Min1")]] <- data.table(
    Penalty, step=2, min1.simple)[, polygon.i := .I]
  ggplot()+
    geom_tile(aes(
      param1, param2, fill=loss),
      data=loss.dt[data.i==2])+
    geom_sf(data=loss.sf.list[[2]]$geometry, fill=NA)+
    coord_sf(xlim=c(0,2*pi),ylim=c(0,2*pi))+
    scale_fill_gradient(low="white", high="red")
  #?st_intersections says spatial index is built for first arg, TODO
  #should first arg be smaller or larger? you build an index on a table
  #that you want to query, so that is the large set, should be the
  #current cost as first arg, new data point loss as second arg.
  cost2.dt <- data.table(sf::st_intersection(
    sf::st_sf(min1.simple[, .(
      min.Constant=Constant,
      min.Linear1=Linear1,
      min.Linear2=Linear2,
      geometry,
      change
    )]),
    sf::st_sf(loss.sf.list[[2]][, .(
      data.Constant=Constant,
      data.Linear1=Linear1,
      data.Linear2=Linear2,
      geometry
    )])
  ))[, `:=`(
    Constant=min.Constant+data.Constant,
    Linear1=min.Linear1+data.Linear1,
    Linear2=min.Linear2+data.Linear2,
    polygon.i = .I
  )][]
  gg <- ggplot()+
    geom_sf(aes(
      geometry=geometry,
      fill=change),
      data=cost2.dt)+
    geom_sf_text(aes(
      geometry=geometry,
      label=polygon.i),
      data=cost2.dt)
  cost2.simple <- mySimplify(cost2.dt)
  ggplot()+
    geom_sf(aes(
      geometry=geometry,
      fill=change),
      data=cost2.simple)+
    coord_sf(xlim=c(0,2*pi),ylim=c(0,2*pi))  
  ## Final optimization (and optimization that we should do at each
  ## time step.
  vertex.dt <- cost2.simple[, {
    vertex.mat <- sf::st_coordinates(sf::st_boundary(geometry[[1]]))
    data.table(
      vertex.mat[,1:2],
      cost=Constant+vertex.mat[,1]*Linear1+vertex.mat[,2]*Linear2)
  }, by=.(Constant, Linear1, Linear2, change)]
  best.vertices <- vertex.dt[min(cost)==cost]
  gg+
    geom_point(aes(
      X,Y),
      data=best.vertices)
  last.cost.dt.list[[Penalty]] <- data.table(
    Penalty, step=4, cost2.simple)[, polygon.i := .I]
  last.best.dt.list[[Penalty]] <- data.table(
    Penalty, step=4, best.vertices)
  last.tile.dt.list[[Penalty]] <- data.table(Penalty, melt(loss.wide[
  , min.data1 := pmin(loss.data1, penalty)
  ][
  , cost.data2 := min.data1+loss.data2
  ],
  measure.vars=step.info$title,
  variable.name="title",
  value.name="cost"
  )[step.info, on="title"])
}
addStep <- function(DF){
  out <- data.table(DF)
  if(!"title" %in% names(out)){
    out <- out[step.info, on="step", nomatch=0L]
  }
  out[, Step := paste("Step", step)]
}
last.cost.dt <- addStep(do.call(
  rbind, lapply(last.cost.dt.list, sf::st_sf)
))[, `:=`(
  vertices = sapply(last.cost.dt$geometry, function(g)length(unlist(g))/2),
  polygons = sapply(last.cost.dt$geometry, length)
)][]
total.dt <- last.cost.dt[, .(
  unique.parameters=.N,
  polygons=sum(polygons),
  vertices=sum(vertices)
), by=.(Step, title, Penalty)]
last.best.dt <- addStep(rbindlist(last.best.dt.list))
last.tile.dt <- addStep(rbindlist(last.tile.dt.list))[
, relative.cost := (cost-min(cost))/(max(cost)-min(cost)), by=.(Step, Penalty)]
last.loss.dt <- addStep(rbindlist(
  loss.sf.list
)[, step := ifelse(data.i==1, 1, 3)])
gg <- ggplot()+
  ggtitle("Change-point detection, geodesic loss, data1=(1,3), data2=(2,4)")+
  theme_bw()+
  scale_color_manual(values=c(
    "TRUE"="black",
    "FALSE"="deepskyblue"))+
  geom_tile(aes(
    param1, param2, fill=relative.cost),
    data=last.tile.dt)+
  geom_sf(aes(
    geometry=geometry,
    color=change),
    fill=NA,
    ##linewidth=1,
    data=last.cost.dt)+
  ## geom_sf_text(aes(
  ##   geometry=geometry,
  ##   label=polygon.i),
  ##   size=2.5,
  ##   data=last.cost.dt)+
  geom_text(aes(
    0,-0.1,label=sprintf(
      "%d unique parameters,\n%d polygons, %d vertices",
      unique.parameters, polygons, vertices)),
    hjust=0,
    size=2,
    vjust=1,
    data=total.dt)+
  scale_shape_manual(values=c(min=21))+
  geom_point(aes(
    X,Y, shape=vertex),
    fill="white",
    color="grey50",
    data=last.best.dt[, vertex := "min"])+
  geom_sf(aes(
    geometry=geometry),
    fill=NA,
    data=last.loss.dt)+
  scale_fill_gradient(low="white", high="red")+
  facet_grid(Penalty ~ Step + title)+
  scale_x_continuous(
    "Angle parameter 1")+
  scale_y_continuous(
    "Angle parameter 2",
    breaks=seq(0,6,by=1),
    limits=c(-1, 6.5))
png("figure-min-l1-2d.png", width=6.5, height=4, units="in", res=300)
print(gg)
dev.off()

