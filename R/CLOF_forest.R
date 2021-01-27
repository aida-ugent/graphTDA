#' Compute a backbone in a forest graph G through Constrained Leaves Optimal subForest (CLOF)
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
#'                - cost: a data.frame corresponding to the obtained cost according to the number of leaves and the component
#'                - full_cost: a vector of the full cost of each component
#'                - (if CLOF_vcost is used) includedV: a list of lists, one for each component, marking the nodes in the sequential solutions for CLOF
#'                - (if CLOF_ecost is used) includedE: a list of lists, one for each component, marking the edges in the sequential solutions for CLOF

CLOF_forest <- function(G,
                        leaves=NA,
                        CLOF_vcost=betweenness,
                        CLOF_ecost=NULL,
                        preprune=TRUE)
{
  # Preliminary checks
  if(!is.igraph(G)) stop("G must be an igraph object")
  C <- components(G)
  if(!length(E(G)) == length(V(G)) - C$no) stop("G must be a forest graph, i.e., contain no cycles")
  if(!class(CLOF_vcost) == "function") stop("CLOF_vcost must be a valid vertex-valued function")
  if(!is.null(CLOF_ecost) & !class(CLOF_ecost) == "function") stop("CLOF_ecost must be a valid edge-valued function, or NULL to use CLOF_vcost")

  # Make the necessary variables
  if(length(leaves) == 1) leaves <- rep(leaves, C$no)
  else if(length(leaves) != C$no) stop("Length of 'leaves' must be 1 or the number of components of G")
  nodes_to_add <- character()
  edges_included <- integer()
  cost <- data.frame(integer(0), numeric(0), integer(0))
  full_cost <- numeric(C$no)
  if(is.null(CLOF_ecost)){
    includedV <- rep(list(NULL), C$no)
    includedE <- NULL
  }
  else{
    includedV <- NULL
    includedE <- rep(list(NULL), C$no)
  }

  # Perform CLOF on each tree component separately
  for(c in seq(C$no)){
    H <- induced.subgraph(G, which(C$membership == c))
    H <- CLOF_tree(H, leaves=leaves[c], CLOF_vcost=betweenness, CLOF_ecost=CLOF_ecost, preprune=preprune)
    if(!is.null(H$B)){
      nodes_to_add <- c(nodes_to_add, V(H$B)$name[which(degree(H$B) == 0)])
      edges_included <- c(edges_included, get.edge.ids(G, t(ends(H$B, E(H$B)))))
      if(!is.null(H$cost)) cost <- rbind(cost, cbind(H$cost, data.frame(component=c)))
      if(!is.null(H$full_cost)) full_cost[c] <- H$full_cost
      if(!is.null(H$includedV)) includedV[[c]] <- H$includedV
      if(!is.null(H$includedE)) includedE[[c]] <- H$includedE
    }
  }

  # Construct final output object
  cost[,"component"] <- factor(cost[,"component"])
  solution <- list(B=subgraph.edges(G, edges_included), cost=cost, full_cost=full_cost,
                   membership=C$membership, includedV=includedV, includedE=includedE)
  if(length(nodes_to_add) > 0){ # add isolated nodes to the backbone
    solution$B <- add_vertices(solution$B, length(nodes_to_add))
    V(solution$B)$name[seq(length(V(solution$B)) - length(nodes_to_add) + 1, length(V(solution$B)))] <- nodes_to_add
  }

  # Return solution
  return(solution)
}
