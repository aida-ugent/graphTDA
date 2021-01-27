#' Compute a backbone in a tree graph G through Constrained Leaves Optimal subForest (CLOF)
#'
#' @param G          forest graph from which a backbone is to be retrieved
#' @param leaves     (Optional) Maximal number of leaves to be included in the backbone.
#'                   May also be a vector with a length equal to the number of components of G.
#'                   Note that this number may not be achieved, e.g., if there are not enough leaves in G.
#' @param CLOF_vcost (Optional) vertex valued function used to solve CLOF.
#'                   Standard is vertex betweenness.
#' @param CLOF_ecost (Optional) edge valued function used to solve CLOF.
#'                   Standard is NULL, and a vertex valued function is used.
#'                   If not NULL, CLOF_vcost is ignored.
#' @param preprune   (Optional) Whether to preprune the f-pine, i.e., discard all its leaves once if the cost for CLOF is constant on leaves.
#'                   This may significantly improve computation cost, and one may directly obtain the solution in the original pine by adding arbitrary leaves.
#'                   Standard is TRUE.
#'
#' @return A list containing the following items
#'                - B: the solution, i.e., backbone, obtained through CLOF
#'                - membership: the membership of the nodes in G according to the components of G
#'                - cost: a data.frame corresponding to the obtained cost according to the number of leaves
#'                - full_cost: the full cost of G
#'                - (if CLOF_vcost is used) includedV: a list marking the nodes in the sequential solutions for CLOF
#'                - (if CLOF_ecost is used) includedE: a list marking the nodes in the sequential solutions for CLOF

CLOF_tree <- function(G,
                      leaves=NA,
                      CLOF_vcost=betweenness,
                      CLOF_ecost=NULL,
                      preprune=TRUE)
  {
  # Preliminary checks
  if(!is.igraph(G)) stop("G must be an igraph object")
  if(length(V(G)) == 0) stop("G may not be the empty graph")
  if(length(E(G)) != length(V(G)) - 1 | any(is.na(bfs(G, 1, unreachable=FALSE)$order))) stop("G must be a tree graph, i.e., be connected and contain no cycles")
  if(!class(CLOF_vcost) == "function") stop("CLOF_vcost must be a valid vertex-valued function")
  if(!is.null(CLOF_ecost) & !class(CLOF_ecost) == "function") stop("CLOF_ecost must be a valid edge-valued function, or NULL to use CLOF_vcost")

  # Make the necessary variables
  if(is.na(leaves)) leaves <- Inf
  if(!is.null(CLOF_ecost)){
    vcost <- NULL
    ecost <- CLOF_ecost(G)
    full_cost <- sum(ecost)
  }
  else{
    vcost <- CLOF_vcost(G)
    ecost <- NULL
    full_cost <- sum(vcost)
  }
  pruned_once <- FALSE
  ends <- V(G)$name[degree(G) == 1]
  to <- NULL

  # Check for preprune
  if(preprune & length(V(G)) > 2 & ((is.null(ecost) & is.null(vcost)) | length(unique(vcost[ends])) == 1)){
    pruned_once <- TRUE
    value <- if(is.null(vcost)) 1 else vcost[ends[1]]
    to <- setdiff(V(G)$name, ends)
    H <- delete_vertices(G, ends)
    ends <- V(H)$name[degree(H) == 1]
    rm(H)
  }

  # We check for special cases for which we manually define the output
  if((length(V(G)) == 2 & !pruned_once) | length(ends) == 0 | leaves < 2){ # special cases
    l <- 0
    includedE <- NULL
    includedV <- NULL
    if(is.null(vcost) | full_cost == 0) cost <- 0
    else{
      i <- which.max(vcost[if(is.null(to)) V(G) else V(G)[to]])
      includedV <- list((if(is.null(to)) V(G)$name else V(G)[to]$name)[i])
      cost <- (vcost[if(is.null(to)) V(G) else V(G)[to]])[i] / full_cost
    }
    if(leaves > 1 & length(ends) > 0){
      l <- c(l, 2)
      cost <- c(cost, 1)
      if(is.null(vcost)) includedE <- list(integer(0), t(ends(G, E(G))))
      else includedV[[2]] <- V(G)$name
    }

    return(list(B=if(is.null(vcost)) subgraph.edges(G, get.edge.ids(G, t(ends(G, includedE[[length(includedE)]])))) else
      induced.subgraph(G, includedV[[length(includedV)]]), cost=data.frame("leaves"=l, "cost"=cost), full_cost=full_cost,
      includedV=includedV, includedE=includedE))
  }

  # Obtain the costs defined by all paths from nodes to leaves
  if(is.null(vcost)){
    store_nodes <- FALSE
    vcost <- integer(length(V(G))) # initialize by zeros to always obtain the added cost from D and vcost below
    D <- distances(G, v=ends, to=if(is.null(to)) V(G) else to, weight=ecost)
    rm(ecost)
  }
  else{
    store_nodes <- TRUE
    D <- (distances(G, v=ends, to=if(is.null(to)) V(G) else to, weight=sum_weights(G, vcost)) +
            sapply(vcost[if(is.null(to)) V(G) else to], function(c) c + vcost[ends])) / 2
  }
  names(vcost) <- names(V(G))

  # Determine the initial optimal path (solution for 2 leaves)
  rowMaxInd <- apply(D[ends, ends], 1, which.max)
  u <- which.max(D[ends, ends][cbind(seq(length(ends)), rowMaxInd)])
  v <- ends[rowMaxInd[u]]
  u <- ends[u]
  P <- shortest_paths(G, from=u, to=v, output="both")
  if(store_nodes){
    i <- which.max(vcost[if(is.null(to)) V(G) else V(G)[to]])
    cost <- c((vcost[if(is.null(to)) V(G) else V(G)[to]])[i] / full_cost, D[u, v] / full_cost)
    includedE <- NULL
    includedV <- list((if(is.null(to)) V(G)$name else V(G)[to]$name)[i], P$vpath[[1]]$name)
  } else{
    cost <- c(0, D[u, v] / full_cost)
    includedE <- list(integer(0), t(ends(G, P$epath[[1]])))
    includedV <- NULL
  }

  # Determine the next solutions for 2 < leaves if required
  l <- 2
  add_bif <- FALSE
  leaves_to_add <- setdiff(rownames(D)[!is.na(D[,u])], c(u, v))
  if(length(leaves_to_add) > 0 & l < leaves){
    add_bif <- TRUE
    if(length(leaves_to_add) == 1) closest_point <- names(P$vpath[[1]][which.min(D[leaves_to_add, names(P$vpath[[1]])])])
    else closest_point <- names(P$vpath[[1]][apply(D[leaves_to_add, names(P$vpath[[1]])], 1, which.min)])
    names(closest_point) <- leaves_to_add
    while(l < leaves & length(leaves_to_add) > 0){
      x <- leaves_to_add[which.max(D[cbind(leaves_to_add, closest_point[leaves_to_add])] - vcost[closest_point[leaves_to_add]])]
      y <- closest_point[x]
      P <- shortest_paths(G, from=x, to=y,  output="both")
      if(l == 2) cost <- c(cost[1], max(D[u, y], D[v, y]) / full_cost, cost[2]) # fictional cost for the unachievable case of one leaf
      l <- l + 1
      cost <- c(cost, cost[length(cost)] + (D[x, y] - vcost[y]) / full_cost)
      if(store_nodes) includedV[[l]] <- P$vpath[[1]]$name[-length(P$vpath[[1]])]
      else includedE[[l]] <- t(ends(G, P$epath[[1]]))
      leaves_to_add <- setdiff(leaves_to_add, x)
      for(x in leaves_to_add){
        y <- names(P$vpath[[1]][which.min(D[x, names(P$vpath[[1]])] - vcost[P$vpath[[1]]])])
        if(D[x, y] - vcost[y] < D[x, closest_point[x]] - vcost[closest_point[x]]) closest_point[x] <- y
      }
    }
  }

  # Construct final output object
  if(pruned_once){
    cost <- c(cost, cost[length(cost)] + value / full_cost)
    l <- l + 1
  }
  if(store_nodes) solution <- induced.subgraph(G, unlist(includedV)[-1])
  else solution <- subgraph.edges(G, get.edge.ids(G, unlist(includedE)))
  cost <- data.frame("leaves"=c(0, if(add_bif) 1 else NULL, 2:l), "cost"=cost)
  backbone <- list(B=solution, cost=cost, full_cost=full_cost, includedV=includedV, includedE=includedE)

  # Return solution
  return(backbone)
}
