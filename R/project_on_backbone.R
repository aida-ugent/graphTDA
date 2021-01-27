#' Project a given graph G on a backbone subgraph G.
#'
#' @param G  igraph object.
#' @param B  Backbone subgraph of G.
#' @param dG (Optional) pairwise distance matrix of G.
#'           Computed if not provided.
#'
#' @return  The projection of G on B (list/igraph object), with two added entries:
#'               - addedV: the set difference between the nodes of G and B.
#'               - connectTo: for each node in addedV, the point in B on which it is projected.

project_on_backbone <- function(G, B, dG=NULL){
  # Preliminary checks
  if(!is_igraph(G) | !is_igraph(B)) stop("G and B must be igraph objects.")
  if(!all(V(B)$name %in% V(G)$name)) stop("B must be a subgraph of G.")

  # Make necessary variables
  VtoAdd <- setdiff(names(V(G)), names(V(B)))
  if(is.null(dG)) dG <- distances(G, VtoAdd, names(V(B)))

  # Determine edges to in the projection
  connectTo <- names(V(B))[apply(dG[VtoAdd, names(V(B))], 1, which.min)]
  EtoAdd <- unique(unlist(lapply(1:length(VtoAdd), function(n)
    shortest_paths(G, VtoAdd[n], connectTo[n], output="epath")$epath[[1]])))

  # Add nodes and edges to the backbone graph
  GprojB <- add.edges(add_vertices(B, length(VtoAdd), name=VtoAdd), t(ends(G, EtoAdd)))
  if(!is.null(E(B)$weight)) E(GprojB)$weight[(length(E(B)) + 1):(length(E(GprojB)))] <- E(G)$weight[EtoAdd]
  GprojB$addedV <- VtoAdd
  GprojB$connectTo <- connectTo

  # Return the projection of G on B
  return(GprojB)
}
