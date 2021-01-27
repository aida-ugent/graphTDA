#' Convert an igraph object to a simplicial complex.
#'
#' @param G The input graph.
#'
#' @return A simplicial complex S, compatible with the TDA library in R.

graph_complex <- function(G)
{
  # Preliminary checks
  if(!is_igraph(G)) stop("G must be of class igraph.")

  # Construct simplical complex object
  S <- c(lapply(V(G), function(v) as.integer(v)),
         split(t(ends(G, E(G), names=FALSE)), rep(1:length(E(G)), each=2)))
  class(S) <- "cmplx"

  # Return simplical complex object
  return(S)
}
