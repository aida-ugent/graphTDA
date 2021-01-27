#' Identify the endpoints of a proximity graph G in the 2D Euclidean plane through the data from which it was constructed.
#'
#' @param df Two-dimensional data.frame where the rownames correspond to the vertex names in G (not necessarily in order).
#' @param G  A proximity graph constructed from the points (represented) in df.
#'
#' @return   A four-dimensional data.frame containing the coordinates of the endpoints of in the form (x1, y1, x2, y2).

get_edges2D <- function(df, G)
{
  # Preliminary checks
  if(ncol(df) != 2) stop("df must be two-dimensional.")
  if(!is.igraph(G)) stop("G must be an igraph object.")
  if(!all(names(V(G)) %in% rownames(df))) stop("all vertex names in G bust me row names in df.")
  if(length(E(G)) == 0) return(NULL) # no coordinates to be returned

  # Construct coordinate matrix
  df.E <- ends(G, E(G))
  df.E <- data.frame(df[df.E[,1],], df[df.E[,2],])
  colnames(df.E) <- c("x1", "y1", "x2", "y2")
  rownames(df.E) <- 1:nrow(df.E)

  # Return coordinate matrix
  return(df.E)
}
