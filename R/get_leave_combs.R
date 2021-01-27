#' For a given number of leaves, compute all possible leave configurations over a number of components.
#'
#' @param c Number of components.
#' @param l Number which the sum of leaves over all components must equal.
#
#' @return A matrix with c columns, where each row is a possible combination of leaves over the c components.

get_leave_combs <- function(c, l)
{
  # Preliminary checks
  if(!is.numeric(c) | !is.numeric(l)) stop("Number of components c and l must be numeric.")
  if(!(c == floor(c)) | !(l == floor(l))) stop("Number of components c and l must be integer.")
  if(c < 0 | l < 0) stop("Number of components c and l must be nonnegative.")

  # Check for special cases
  if(l == 0) return(rep(0, c)) # we treat the empty graph as a tree with no leaf
  else if(l == 1) return(diag(c)) # we treat a point as a tree with one leaf
  if(c == 1) return(l) # easy solution for one component

  # Obtain the possible leave combinations through recursion
  leavecombs <- t(sapply(c(0:l), function(k){
    if(k < l) return(cbind(k, get_leave_combs(c - 1, l - k)))
    else return(c(k, rep(0, c - 1)))
  }))

  # Concatenate results
  if(c >  2) leavecombs <- do.call("rbind", leavecombs)
  colnames(leavecombs) <- NULL

  # Return output
  return(leavecombs)
}
