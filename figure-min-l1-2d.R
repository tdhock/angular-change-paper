library(data.table)
library(ggplot2)
N.grid <- 41
param.grid <- seq(0, 2*pi, l=N.grid)
get.loss <- function(data.val, param.val){
  d <- abs(data.val-param.val)
  ifelse(d>pi, 2*pi-d, d)
}
set.seed(1)
(data.mat <- matrix(1:4, 2))
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
      data.table(Linear, Constant, min_param, max_param)
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
    poly.sf.list[[polygon.i]] <- with(irow.list, sf::st_sf(
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
  plot(poly.sf)
  loss.sf.list[[data.i]] <- poly.sf
  loss.dt.list[[data.i]] <- data.table(
    data.i, data.vals, pair.dt[, loss := loss1+loss2])
}
(loss.dt <- rbindlist(loss.dt.list))

sf::st_intersects(loss.sf.list[[1]], loss.sf.list[[2]])

ggplot()+
  geom_sf(data=loss.sf.list[[1]])

ggplot()+
  geom_sf(data=loss.sf.list[[2]])

penalty <- 1.5
loss.sf.current <- loss.sf.list[[1]]
for(polygon.i in 1:nrow(loss.sf.current)){
  loss.sf.piece <- loss.sf.current[polygon.i,]
  intercept <- (penalty-loss.sf.piece$Constant)/loss.sf.piece$Linear2
  slope <- -loss.sf.piece$Linear1-loss.sf.piece$Linear2
  roots.line <- sf::st_linestring(rbind(
    c(0,intercept),
    c(2*pi, intercept+slope*2*pi)
  ))
  int.feat <- sf::st_intersection(loss.sf.piece$geometry, roots.line)
  if(length(int.feat)){
    print(int.feat)
    collection <- sf::st_geometrycollection(list(loss.sf.piece$geometry[[1]], int.feat[[1]]))
    gg <- ggplot()+
      geom_sf(data=collection)+
      coord_sf(xlim=c(0,2*pi),ylim=c(0,2*pi))    
    print(gg)
    browser()
    ## TODO visualize these lines on top of heat map.
  }
}

color.scale <- scale_color_manual(values=c(
  "1,3"="black",
  "2,4"="grey50"))
gg <- ggplot()+
  geom_tile(aes(
    param1, param2, fill=loss),
    data=loss.dt)+
  color.scale+
  geom_rect(aes(
    xmin=min1, xmax=max1,
    ymin=min2, ymax=max2,
    color=data.vals),
    fill=NA,
    data=draw.rect.dt)+
  geom_text(aes(
    (min1+max1)/2,
    (min2+max2)/2,
    color=data.vals,
    label=rect.i),
    data=draw.rect.dt)+
  facet_grid(. ~ data.vals, labeller=label_both)+
  scale_fill_gradient(low="white", high="red")+
  coord_equal()
png("figure-min-l1-2d-loss.png", width=8, height=4, units="in", res=100)
print(gg)
dev.off()

(loss.wide <- dcast(
  loss.dt[, lossi := paste0("loss", data.i)],
  param1+param2 ~ lossi,
  value.var="loss"))
gg <- ggplot()+
  geom_tile(aes(
    param1, param2, fill=loss1-loss2),
    data=loss.dt)+
  color.scale+
  geom_rect(aes(
    xmin=min1, xmax=max1,
    color=data.vals,
    ymin=min2, ymax=max2),
    fill=NA,
    data=draw.rect.dt)+
  geom_text(aes(
    (min1+max1)/2,
    (min2+max2)/2,
    color=data.vals,
    label=rect.i),
    data=draw.rect.dt)+
  scale_fill_gradient2()+
  coord_equal()
png("figure-roots-l1-2d.png", width=4, height=4, units="in", res=100)
print(gg)
dev.off()

coll.dt <- draw.rect.dt[
, poly.list := list(list(sf::st_polygon(list(rbind(
  c(min1,min2),
  c(min1,max2),
  c(max1,max2),
  c(max1,min2),
  c(min1,min2)
)))))
, by=.(data.i, i1, i2)
][
, .(poly.coll=list(sf::st_sfc(poly.list)))
, by=.(data.i, data.vals)]
plot(coll.dt$poly.coll[[1]])
plot(coll.dt$poly.coll[[2]])
draw.rect.dt[, .(data.i, min1, max1, min2, max2)]
(ilist <- coll.dt[, sf::st_intersects(poly.coll[[1]], poly.coll[[2]])])
## this list has an element for each item in second set of
## polygons. The element is an integer vector of overlapping polygons
## from the first set. How many of these have zeros?
sf::st_point#see ex on https://r-spatial.github.io/sf/reference/geos_binary_pred.html
