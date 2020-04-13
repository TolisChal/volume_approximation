#' Generator function for cross polytopes
#' 
#' This function can be used to generate the \eqn{d}-dimensional cross polytope in H- or V-representation.
#' 
#' @param dimension The dimension of the cross polytope.
#' 
#' @return A polytope class representing a cross polytope in H- or V-representation.
#' @examples 
#' # generate a 10-dimensional cross polytope in H-representation
#' P = gen_dual_knapsack(10)
#' 
#' @export
gen_dual_knapsack <- function(dimension) {
  
  
  Mat = poly_gen(6, TRUE, FALSE, dimension, 0)
  
  # first column is the vector b
  b = Mat[,1]
  Mat = Mat[,-c(1)]
  

  P = Vpolytope$new(Mat)
  
  return(P)
  
}