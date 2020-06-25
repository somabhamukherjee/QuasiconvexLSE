library(Rearrangement)

vec2gridmat <- function(xygrid, z) {
  # xygrid is something created from e.g. expand.grid()
  #
  # only the *number* of unique values in each column of xygrid
  # is used. The actual values are not used.
  stopifnot(ncol(xygrid) == 2)
  x_unique <- unique(xygrid[, 1])
  y_unique <- unique(xygrid[, 2])

  # "nu" = "number of unique values"
  x_nu <- length(x_unique)
  y_nu <- length(y_unique)
  zlen <- length(z)

  if (zlen != x_nu * y_nu)
    stop("length of z not consistent with xygrid")

  z_mat <- matrix(z, nrow = x_nu, ncol = y_nu)
}


# This function provides a wrapper around Rearrangement::rearrangement()
# so that we can call this operator in the same way as the others

# same arguments as qconvexify_2d
# ... is passed to rearrangement()
monotonize <- function(x, y, ...) {
  stopifnot(ncol(x) == 2)
  x1_unique <- unique(x[, 1])
  x2_unique <- unique(x[, 2])

  if (length(x1_unique) * length(x2_unique) != length(y))
    stop("monotonize() only works if x is a balanced grid")

  # alternative:
  # y_mat <- matrix(y, nrow = length(x1_unique), ncol = length(x2_unique))
  y_mat <- vec2gridmat(xygrid = x, z = y)

  # c() unwinds the matrix in the same way that matrix() winds it (by default).
  mono <- c(rearrangement(x = list(x1_unique, x2_unique), y = y_mat))
  return(mono)
}
