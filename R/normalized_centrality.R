#' Compute the normalized centrality for each node in a given graph G.
#'
#' @param G The input graph.
#'
#' @return A numeric vector of length(V(G)) containing the degrees of centrality, computed per component.

normalized_centrality <- function(G)
{
  # Preliminary checks
  if(!is_igraph(G)) stop("G must be of class igraph.")

  # Compute vertex eccentricities
  vertex_eccentricities <- apply(distances(G), 1, function(r) max(r[r < Inf]))

  # Negate and scale
  m <- min(vertex_eccentricities )
  M <- max(vertex_eccentricities )
  normalized_centralities  <- (M - vertex_eccentricities) / (M - m)

  # Return the normalized centralities
  return(normalized_centralities)
}
