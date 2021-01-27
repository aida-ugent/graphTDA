#' Construct an undirected kNN graph by connecting each point to its k closest neighbors.
#'
#' @param X A numeric data.frame or matrix containing observations, or a distance object containing distances between pairs of observations.
#' @param k Number of neighbors in the kNN graph.
#'          Standard is 10.
#'
#' @return  Undirected k-nearest neighbor graph G (igraph object).
#'          Note that the order of nodes in X and G may differ.
#'          Rownames in X correspond to vertex names in G.

kNNgraph <- function(X, k=10){
  # Preliminary checks
  if(!is.numeric(k) | !floor(k) == k | k < 0) stop("k must be a nonnegative integer.")
  if(!is.data.frame(X) & !is.matrix(X) & !class(X) == "dist") stop("X must be a data.frame, matrix, or dist object.")
  if(!class(X) == "dist" & k > nrow(X) - 1) stop("k must be at most equal to nrow(X) - 1.")

  # Case 1: X is a point cloud data set
  if(is.data.frame(X) | is.matrix(X)){
    kNN <- get.knn(X, k)
    adj <- cbind(from=rep(1:nrow(X), each=k), to=as.integer(t(kNN[[1]])), weight=as.numeric(t(kNN[[2]])))
    I <- adj[,1] > adj[,2]
    adj[I, c(1, 2)] <- adj[I, c(2, 1)]
    adj <- adj[!duplicated(adj[,c(1, 2)]),]
    G <- graph.data.frame(adj, directed=FALSE)
    V(G)$name <- rownames(X)[as.integer(V(G)$name)]
  }

  # Case 2: X is a distance object
  else if(class(X) == "dist"){
    if(k > attr(X, "Size") - 1) stop("k must be at least equal to nrow(X) - 1")
    adj <- cbind(from=rep(1:attr(X, "Size"), each=k),
                  to=as.integer(apply(as.matrix(X), 1, function(r) order(r)[2:(k + 1)])))
    I <- adj[,1] > adj[,2]
    adj[I, c(1, 2)] <- adj[I, c(2, 1)]
    adj <- adj[!duplicated(adj[,c(1, 2)]),]
    G <- graph.data.frame(adj, directed=FALSE)
    V(G)$name <- attr(X, "Labels")[as.integer(V(G)$name)]
    E(G)$weight <- as.matrix(X)[adj]
  }

  # Return kNN graph
  return(G)
}
