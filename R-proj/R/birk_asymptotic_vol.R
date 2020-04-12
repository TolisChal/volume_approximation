birk_asymptotic_vol <- function(n){
  
  return(exp(-(n-1)^2*log(n) + n^2 - (n-0.5)*log(2*pi) + 1/3))
  
  
}