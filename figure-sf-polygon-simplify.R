m.triangle <- rbind(
  c(0,0),
  c(0,0.5),
  c(0,1),
  c(1,1),
  c(0,0)
)
m.rect <- rbind(
  c(0,0.4),
  c(0,2),
  c(-1,2),
  c(-1,0.4),
  c(0,0.4)
)
## from https://r-spatial.github.io/sf/articles/sf1.html
##POLYGON	geometry with a positive area (two-dimensional); sequence of points form a closed, non-self intersecting ring; the first ring denotes the exterior ring, zero or more subsequent rings denote holes in this exterior ring
p <- sf::st_polygon(list(m.triangle, m.rect)) #m.rect is the hole
##The way this is printed is called well-known text, and is part of the standards. The word MULTIPOLYGON is followed by three parentheses, because it can consist of multiple polygons, in the form of MULTIPOLYGON(POL1,POL2), where POL1 might consist of an exterior ring and zero or more interior rings, as of (EXT1,HOLE1,HOLE2). Sets of coordinates belonging to a single polygon are held together with parentheses, so we get ((crds_ext)(crds_hole1)(crds_hole2)) where crds_ is a comma-separated set of coordinates of a ring. This leads to the case above, where MULTIPOLYGON(((crds_ext))) refers to the exterior ring (1), without holes (2), of the first polygon (3) - hence three parentheses.
plot(p)
plot(sf::st_simplify(p))
two.polys <- sf::st_multipolygon(list(list(m.triangle), list(m.rect)))
plot(two.polys)
(two.polys.joined <- sf::st_simplify(two.polys))
plot(two.polys.joined)
sf.obj <- sf::st_sf(Linear=1.0, Constant=0, params=sf::st_sfc(two.polys.joined))
sf.obj$params[[1]]#is the polygon
library(ggplot2)
ggplot()+
  theme_bw()+
  geom_sf(data=sf.obj)
