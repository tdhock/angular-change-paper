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
draw.rect.dt.list <- list()
for(data.i in 1:nrow(data.mat)){
  data.row <- data.mat[data.i,]
  data.vals <- paste(data.row, collapse=",")
  iseq <- 1:N.grid
  mseq <- 1:3
  Mseq <- mseq+1L
  rect.dt <- CJ(i1=mseq, i2=mseq)[, rect.i := .I]
  pair.dt <- CJ(i1=iseq, i2=iseq)
  for(dim.i in seq_along(data.row)){
    iname <- paste0("i",dim.i)
    data.val <- data.row[[dim.i]]
    too.many.edges <- c(0, data.val-pi, data.val, data.val+pi, 2*pi)
    edge.vec <- too.many.edges[too.many.edges %between% c(0,2*pi)]
    min.max.dt <- data.table()[
    , paste0("min",dim.i) := edge.vec[mseq]
    ][
    , paste0("max",dim.i) := edge.vec[Mseq]
    ][
    , (iname) := mseq
    ]
    rect.dt <- rect.dt[min.max.dt, on=iname]
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
  pair.dt[, loss := loss1+loss2]
  loss.dt.list[[data.i]] <- data.table(
    data.i, data.vals, pair.dt)
  draw.rect.dt.list[[data.i]] <- data.table(
    data.i, data.vals, rect.dt)
}
(loss.dt <- rbindlist(loss.dt.list))
(draw.rect.dt <- rbindlist(draw.rect.dt.list))
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
png("figure-roots-l1-2d-loss.png", width=8, height=4, units="in", res=100)
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
