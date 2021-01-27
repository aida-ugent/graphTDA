#' Construct a Vietoris-Rips graph from by connecting each point to all points within a certain distance eps.
#'
#' @param X   A numeric data.frame or matrix containing observations, or a distance object containing distances between pairs of observations.
#' @param eps Upper bound on edge weights/lengths.
#'
#' @return Undirected Vietoris-Rips graph G (igraph object).
#'         Note that the order of nodes in X and G may differ.
#'         Rownames in X correspond to vertex names in G.

VRgraph <- function(X, eps)
{
  # Preliminary checks
  if(!is.numeric(eps)) stop("eps must be a given numeric value.")
  if(!is.data.frame(X) & !is.matrix(X) & !class(X) == "dist") stop("X must be a data.frame, matrix, or dist object.")

  # Case 1: X is a point cloud data set
  if(is.data.frame(X) | is.matrix(X)){
    adj <- bind_rows(lapply(1:(nrow(X) - 1), function(n){
      d <- pdist(X[n,], X[(n + 1):nrow(X),])@dist
      to <- which(d < eps)
      return(data.frame(from=rep(n, length(to)), to=((n + 1):nrow(X))[to], weight=d[to]))
    }))
    G <- graph.data.frame(adj, directed=FALSE)
    V(G)$name <- rownames(X)[as.integer(V(G)$name)]
  }

  # Case 2: X is a distance object
  else if(class(X) == "dist"){
    adj <- which(as.matrix(X) < eps, arr.ind=TRUE)
    adj <- adj[adj[,1] < adj[,2],]
    G <- graph.data.frame(adj, directed=FALSE)
    V(G)$name <- attr(X, "Labels")[as.integer(V(G)$name)]
    E(G)$weight <- as.matrix(X)[adj]
  }

  # Return Rips graph
  return(G)
}
