# These are used for dist_from_poly()
library(sp)
library(rgeos)

dist_from_poly <- function(x, y, xp, yp) {
  stopifnot(length(xp) >= 2)

  the_point_sp <- SpatialPoints(cbind(x, y))
  one_poly <- cbind(xp, yp)

  if (length(xp) > 2) {
    # this avoids the warning "less than 4 coordinates in polygon".
    # (does not change the end result, though)
    one_poly <- rbind(one_poly, one_poly[1, ])
    pts_sp <- SpatialPoints(one_poly)
    # does not work with a line (that's why we condition on length(xp) > 2)
    psk <- Polygons(list(Polygon(pts_sp)), ID = "not_used")
    def <- SpatialPolygons(list(psk))
  } else {
    pts_sp <- SpatialPoints(one_poly)
    psk <- Lines(list(Line(pts_sp)), ID = "not_used")
    def <- SpatialLines(list(psk))
  }

  dist_ <- gDistance(the_point_sp, def)

  return(dist_)
}


do_qconvexify <- function(x, i = i, id = id) {
  # TODO change name (?)
  idi <- id[i]

  l <- 1
  u <- length(id)
  while (u > l){
    p <- (u + l)/2
    I <- id[1:floor(p)]
    if (p < 3){
      if ((idi == I[1])||(idi == I[length(I)])) {
        u <- floor(p)
      } else {
        l <- ceiling(p)
      }
    }
    if (p >= 3){
      q <- chull(x[I,])
      d_ <- dist_from_poly(x[id[i],1], x[id[i],2],x[I[q],1], x[I[q],2])
      # TODO ideally we would just use inpolygon(). However, there seems to
      # be a numeric precision issue somewhere.
      # TODO this tolerance is arbitrary. We should either figure out
      # the source of the numerical precision issue, or make the tolerance
      # relative to the grid in some way.
      tol_ <- 1e-7
      if (d_ < tol_){
        u <- floor(p)
      } else {
        l <- ceiling(p)
      }
    }
  }
  return(u)
}

# compute qconvexified f(x) d=2
qconvexify_2d <- function(x, y, parallel = FALSE) {
  # Because of a bug in R, R versions 3.4.0 and below could see a memory
  # leak from this function if JIT is enabled (which it is by default on
  # 3.4.0). This R bug was fixed by Luke Tierney in in R-devel at r72788)
  # and R-patched at (r72789) and should be included in 3.4.1.

  id <- order(y)
  fqc <- rep(0, length(y))
  the_us <- rep(0, length(y))
  if (parallel) {
    # TODO move library call comewhere else
    library(parallel)
    the_us <- simplify2array(mclapply(X = 1:length(y), FUN = function(i) do_qconvexify(i = i, x = x, id = id)))
  } else {
    the_us <- sapply(X = 1:length(y), FUN = function(i) do_qconvexify(i = i, x = x, id = id))
  }
  for (i in 1:length(y)) {
    fqc[id[i]] <- y[id[the_us[i]]]
  }
  fqc
}
