#' Compute zero-dimensional persistence from a vertex-valued defined over a graph G.
#'
#' @param FltFun   (Optional) Filtration of which to compute persistence.
#'                 If missing, both the graph and function values must be provided.
#' @param S        (Optional) A graph represented by either an igraph or cmplx object.
#' @param f        (Optional) A vertex-valued function defined over the nodes of G.
#' @param sublevel (Optional) Only used if FltFun is missing.
#'                 If FALSE, superlevel filtration is computed instead.
#'                 Standard is TRUE.
#'
#' @return The result of the persistent homology computation for the specified filtration through the TDA library in R.

graph_persistence <- function(FltFun=NULL, S=NULL, f=NULL, sublevel=TRUE)
{
  # Preliminary checks
  if(is.null(FltFun) & (is.null(S) | is.null(f))) stop("Either the filtration, or the graph and vertex-valued function must be provided.")

  # Construct filtration if necessary
  if(is.null(FltFun)) FltFun <- graph_filtration(S, f, sublevel)

  # Perform persistent homology computation through the TDA library in R
  diag <- filtrationDiag(filtration=FltFun, maxdimension=0, library="Dionysus")$diagram

  # Return the persitent homology computation results
  return(diag)
}
