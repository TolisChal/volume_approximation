#' Generator function for cross polytopes
#' 
#' This function can be used to generate the \eqn{d}-dimensional cross polytope in H- or V-representation.
#' 
#' @param n The dimension of the cross polytope.
#' @param s The dimension of the cross polytope.
#' 
#' @return A polytope class representing a cross polytope in H- or V-representation.
#' @examples 
#' # generate a 10-dimensional cross polytope in H-representation
#' P = gen_everest(10)
#' 
#' @export
gen_everest <- function(n, s) {
  
  
  Mat = poly_gen(7, TRUE, FALSE, n, s)
  
  # first column is the vector b
  b = Mat[,1]
  Mat = Mat[,-c(1)]
  
  
  P = Vpolytope$new(Mat, (prod(1:((n+1)*s)) / prod(1:(n*s))) / prod(1:s)^(n+1) )
  
  return(P)
  
}