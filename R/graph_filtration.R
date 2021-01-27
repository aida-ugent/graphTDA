#' Compute a filtration from a vertex-valued defined over a graph G.
#'
#' @param S        A graph represented by either an igraph or cmplx object.
#' @param f        A vertex-valued function defined over the nodes of G.
#' @param sublevel (Optional) If FALSE, superlevel filtration is computed instead.
#'                 Standard is TRUE.
#'
#' @return A filtration defined by f on G, compatible with the TDA library in R.

graph_filtration <- function(S, f, sublevel=TRUE)
{
  # Preliminary checks
  if(!is_igraph(S) & class(S) != "cmplx") stop("G must be of a graph or complex object.")

  # Convert to complex if necessary
  if(is_igraph(S)) S <- graph_complex(S)

  # Compute the filtration
  FltFun <- funFiltration(FUNvalues=f, cmplx=S, sublevel=sublevel)

  # Return the filtration
  return(FltFun)
}
