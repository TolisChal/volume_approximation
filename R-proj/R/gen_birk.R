#' Generator function for hypercubes
#' 
#' This function can be used to generate the \eqn{d}-dimensional unit hypercube \eqn{[-1,1]^d} in H- or V-representation.
#' 
#' @param n The dimension of the hypercube
#' 
#' @return A polytope class representing the unit \eqn{d}-dimensional hypercube in H- or V-representation.
#' @examples 
#' # generate a 10-dimensional hypercube in H-representation
#' P = gen_birk(5)
#' 
#' # generate a 15-dimension hypercube in V-representation
#' P = gen_birk(3)
#' @export
gen_birk <- function(n) {
  
  kind_gen = 7
  m_gen = 0
  
  
  Mat = poly_gen(kind_gen, FALSE, FALSE, n, m_gen)
  
  # first column is the vector b
  b = Mat[,1]
  Mat = Mat[,-c(1)]
  
  P = Hpolytope$new(Mat, b)
  
  
  return(P)
  
}
