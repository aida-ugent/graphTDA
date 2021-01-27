#' Compute the boundary coefficients of a weighted graph.
#'
#' @param G             The graph of which boundary coefficients are to be computed.
#' @param constrain_mem (Optional) If TRUE and pairD is missing, the entire distance matrix on G as is computed first, from which the hop-k-approximation is extracted afterwards.
#'                      If FALSE, entries that are not needed in pairD are never computed, requiring less storage.
#'                      Standard is FALSE, which leads to faster results with our current implementation.
#'
#' @return A vector of |V(G)| boundary coefficients.

boundary_coefficient <- function(G, constrain_mem=FALSE)
{
  # Preliminary checks
  if(!is_igraph(G)) stop("G must be an igraph object.")

  # Compute hop-2-approximation
  hop2 <- bounded_hop_pairG(G, constrain_mem=constrain_mem)

  # Compute boundary coefficients
  hop1 <- as.matrix(summary(as_adj(G)))[,c(1, 2)]
  A <- spam(x=list(indices=hop1, values=hop2[hop1]))
  oneOvA <- 1 / A
  BC <- (apply.spam(A, 1, sum) * apply.spam(oneOvA, 1, sum) - diag(oneOvA %*% hop2^2%*% oneOvA) / 2) / (degree(G)^2)

  # Return boundary coefficients
  return(BC)
}
