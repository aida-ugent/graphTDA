#' Compute the elbows for a set of CLOF cost functions (one per component) by minimizing the second-order finite differentiates.
#'
#' @param C CLOF cost function, as stored in a backbone object.
#'
#' @return data.frame where each row contains the elbow (leaves + cost) for the corresponding component.

get_elbow <- function(C)
{
  # Preliminary checks
  if(!ncol(C) == 3) stop("Cost must have three columns (two curve coordinates + component per point).")

  # Compute elbow for each curve
  elbows <- sapply(levels(C[,"component"]), function(c){
    I <- which(C[,"component"] == c)
    if(length(I) <= 3) return(c(NA, NA)) # no elbow for defined for a curve with at most three points
    elbow_idx <- which.min(diff(diff(C[I, "cost"]))[-1]) + 1
    return(c(elbow_idx, C[I[elbow_idx], "cost"]))
  })

  # Convert to a data.frame
  elbows <- data.frame(leaves=as.integer(elbows[1,]), cost=elbows[2,])

  # Output the elbows
  return(elbows)
}
