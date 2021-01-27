#' Compute the hop-k-approximation of a pairwise distance matrix of a given graph G.
#' We suggest a future implementation directly in (and called from) C++.
#'
#' @param G             The graph from which the distance matrix is to be approximated.
#' @param k             (Optional) Maximum allowed hops between pairs of nodes for which the weighted distances are to be computed.
#'                      Standard is k = 2.
#' @param constrain_mem (Optional) If TRUE and pairD is missing, the entire distance matrix on G as is computed first, from which the hop-k-approximation is extracted afterwards.
#'                      If FALSE, entries that are not needed in pairD are never computed, requiring less storage.
#'                      Standard is FALSE, which leads to faster results with our current implementation.
#'
#' @return A sparse matrix storing the weighted distance D_uv between nodes u and v if the unweighted distance between u and v is less than k.

bounded_hop_pairG <- function(G, k=2, constrain_mem=FALSE)
{
  # Preliminary checks
  if(!is_igraph(G)) stop("G must be an igraph object.")

  # Compute hop-k-neighborhoods
  to <- ego(G, order=k)
  i <- as.integer(unlist(sapply(1:length(to), function(n) rep(n, length(to[[n]])))))

  # Case 1: compute hop-k-approximation (fast)
  if(!constrain_mem){
    to <- unlist(to)
    hopk <- spam(x=list(i=i, j=to, values=distances(G)[cbind(i, to)]))
  }
  # Case 2:Compute hop-k-approximation (memory efficient)
  else{
    values <- unlist(lapply(1:length(V(G)), function(v) distances(G, v=v, to=to[[v]])))
    to <- unlist(to)
    hopk <- spam(x=list(i=i, j=to, values=values))
  }

  # Return hop-k-approximation
  return(hopk)
}
