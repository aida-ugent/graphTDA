#' Compute softmax weights for a graph using a predefined function
#'
#' @param G  Graph containing the edges for which the softmax values have to be computed
#' @param f  vector of function values of length |V(G)|
#'
#' @return a vector containing the sum weight for each edge in G, ordered as E(G)

sum_weights <- function(G, f)
{
  # Preliminary checks
  if(!is_igraph(G)) stop("G must be an igraph object.")
  if(!is.numeric(f) | !length(V(G)) == length(f)) stop("f must be a numeric vector of length |V(G)|.")

  # Get matrix representation of edges
  edges <- get.edges(G, E(G))

  # Compute sum weights
  W <- f[edges[,1]] + f[edges[,2]]

  # Return sum weights
  return(W)
}
