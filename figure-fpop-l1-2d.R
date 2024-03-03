library(ggplot2)
library(data.table)
N <- 5
set.seed(1)
data.mat <- rbind(
  cbind(runif(N, -pi/4, pi/4), runif(N, 2, 3)),
  cbind(runif(N, 2, 3), runif(N, -1,1)))
data.mat <- rbind(
  cbind(runif(N, 5, 5.1), runif(N, 2, 2.1)),
  cbind(runif(N, 3, 3.1), runif(N, 4, 4.1)))
## data.mat <- rbind(
##   c(5, 2),
##   c(5.1, 2.1),
##   c(3, 4),
##   c(3.1, 4.1))
print(data.mat)
data.mat[data.mat<0] <- data.mat[data.mat<0]+2*pi
data.wide <- data.table(data.mat)[, data.i := .I]
data.long <- melt(
  data.wide,
  id.vars="data.i")
ggplot()+
  geom_point(aes(
    data.i, value),
    data=data.long)+
  facet_grid(variable ~ .)
    
piece <- function(Linear, Constant, min_param, max_param){
  data.table(Linear, Constant, min_param, max_param)
}
loss_for_data <- function(data.row){
  piece.dt.list <- list()
  for(dim.i in seq_along(data.row)){
    iname <- paste0("i",dim.i)
    data.val <- data.row[[dim.i]]
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
      Constant=i1$Constant+i2$Constant,
      Linear1=i1$Linear,
      Linear2=i2$Linear,
      change=NA_integer_,
      geometry=sf::st_sfc(sf::st_polygon(list(rbind(
        c(i1$min_param, i2$min_param),
        c(i1$max_param, i2$min_param),
        c(i1$max_param, i2$max_param),
        c(i1$min_param, i2$max_param),
        c(i1$min_param, i2$min_param)
      ))))
    ))
  }
  do.call(rbind, poly.sf.list)
}
loss_for_data(c(1,3))

mySimplify <- function(DT){
  out.dt <- DT[, {
    ##print(geometry)
    g1 <- geometry[[1]]
    multi <- if(inherits(g1, "MULTIPOLYGON")){
      g1
    }else if(inherits(g1, "POLYGON")){
      sf::st_multipolygon(geometry)
    }
    if(is.null(multi)){
      data.table()
    }else .(
      geometry=sf::st_sfc(sf::st_simplify(multi))
    )
  }, by=.(Constant, Linear1, Linear2, change)
  ]
  g.out <- sf::st_sfc(out.dt$geometry)
  out.dt[
  , geometry := if(nrow(out.dt)==1)list(g.out) else g.out
  ][]
}

sum_of_two <- function(min.dt, data.dt){
  min.sf <- sf::st_sf(min.dt[, .(
    min.Constant=Constant,
    min.Linear1=Linear1,
    min.Linear2=Linear2,
    geometry,
    change
  )])
  data.sf <- sf::st_sf(data.dt[, .(
    data.Constant=Constant,
    data.Linear1=Linear1,
    data.Linear2=Linear2,
    geometry
  )])
  int.sf <- suppressWarnings(sf::st_intersection(min.sf, data.sf))
  data.table(int.sf)[, `:=`(
    Constant=min.Constant+data.Constant,
    Linear1=min.Linear1+data.Linear1,
    Linear2=min.Linear2+data.Linear2
  )][]
}  

compare_penalty <- function(loss.sf.current, data.i, cost.of.change){
  min1.viz.list <- list()
  for(polygon.i in 1:nrow(loss.sf.current)){
    loss.sf.piece <- loss.sf.current[polygon.i,]
    numerator <- cost.of.change-loss.sf.piece$Constant
    roots.line <- sf::st_linestring(if(loss.sf.piece$Linear2 != 0){
      intercept <- numerator/loss.sf.piece$Linear2
      slope <- -loss.sf.piece$Linear1/loss.sf.piece$Linear2
      rbind(
        c(0, intercept),
        c(2*pi, intercept+slope*2*pi))
    }else{
      xint <- numerator/loss.sf.piece$Linear1
      rbind(
        c(xint, 0),
        c(xint, 2*pi))
    })
    loss.sf.geom <- loss.sf.piece$geometry[[1]]
    split.collection <- lwgeom::st_split(loss.sf.geom, roots.line)
    cmat <- sapply(split.collection, sf::st_centroid)
    if(length(cmat)){
      change.is.better <- loss.sf.piece[
      , Constant+cmat[1,]*Linear1+cmat[2,]*Linear2>cost.of.change]
      min1.viz.list[[polygon.i]] <- loss.sf.piece[, data.table(
        Constant, Linear1, Linear2,
        change,
        geometry=sf::st_sfc(split.collection[])
      )][
        change.is.better,
        `:=`(Constant=cost.of.change, Linear1=0, Linear2=0, change=data.i)
      ]
    }
  }
  rbindlist(min1.viz.list)
}  

min_cost_vertices <- function(cost.dt){
  vertex.dt <- cost.dt[, {
    vertex.mat <- sf::st_coordinates(sf::st_boundary(geometry[[1]]))
    data.table(
      vertex.mat[,1:2],
      cost=Constant+vertex.mat[,1]*Linear1+vertex.mat[,2]*Linear2)
  }, by=.(Constant, Linear1, Linear2, change)]
  vertex.dt[min(cost)==cost]
}  

penalty <- 1
best.param.dt.list <- list()
cost.dt.list <- list()
for(data.i in 1:nrow(data.mat)){
  print(data.i)
  data.vec <- data.mat[data.i,]
  data.loss.dt <- loss_for_data(data.vec)
  cost.i.dt <- if(data.i==1){
    data.loss.dt
  }else{
    ## TODO special cases for penalty=0 or Inf.
    change.cost <- penalty+best.vertices$cost[1]
    min.dt <- compare_penalty(cost.prev.dt, data.i-1, change.cost)
    min.simple <- mySimplify(min.dt)
    sum.dt <- sum_of_two(min.simple, data.loss.dt)
    mySimplify(sum.dt)
  }
  cost.dt.list[[data.i]] <- cost.i.dt
  best.vertices <- min_cost_vertices(cost.i.dt)
  best.param.dt.list[[data.i]] <- print(dcast(
    best.vertices,
    cost + change ~ .,
    list(min, max),
    value.var=c("X","Y")
  )[, `:=`(
    vertices = sum(sapply(cost.i.dt$geometry, function(g)length(unlist(g))/2)),
    polygons = sum(sapply(cost.i.dt$geometry, length))
  )][])
  cost.prev.dt <- cost.i.dt
}
best.param.dt <- rbindlist(best.param.dt.list)[, N.data := .I]

ggplot()+
  geom_point(aes(
    N.data, vertices),
    data=best.param.dt)
