#' Decompose a bacbone according to branches, i.e., maximal paths with no multifurcations.
#'
#' @param backbone A backbone object.
#' @param G        (Optional) the original graph on top of the backbone.
#'                 Must be included if not contained in the backbone object.
#'                 Standard is NULL.
#' @param prob     (Optional & if assign) whether to assign nodes in G to paths according to a probability measure
#'                 if FALSE, then nodes are assigned to exactly one path through a one-hot-encoding.
#'                 Standard is TRUE.
#'
#' @return The backbone object (list) with the following (new) entries:
#'                - branch: a vector of integers clustering the edges of backbone into paths.
#'                - prob: a sparse n x k matrix, where n = |V(G)| and k is the number of branches in the backbone, clustering points of G according to branches of the backbone.
#'                - Palette: a vector of hex code colors, one for each branch.
#'                - Col: a vector of |V(G)| hex code colors, interpolating between the palette colors according to prob.

decompose <- function(backbone, G=NULL, prob=TRUE)
{
  # Preliminary checks
  if(is.null(backbone$G) & is.null(G)) stop("backbone does not include original graph G, please provide G.")

  # Check which nodes lie either on one or on multiple branches
  decompV <- integer(length(V(backbone$B)))
  degrees <- degree(backbone$B)
  arm_nodes <- which(degrees < 3)
  decompV[arm_nodes] <- as.integer(components(induced_subgraph(backbone$B, arm_nodes))$membership)

  # Compute all distances to the backbone
  distToB <- distances(if(is.null(G)) backbone$G else G, v=V(backbone$B)$name)

  # Current branch assignment of backbone nodes according to connected component
  # Only needs to be changed if there is at least one multifurcation point
  if(length(arm_nodes) != length(V(backbone$B))){
    decompV[-arm_nodes] <- sapply(which(decompV == 0), function(v)
      decompV[arm_nodes[which.min(distToB[v, V(backbone$B)$name[arm_nodes]])]])
  }

  # Assign edges according to current assignment of backbone nodes
  # This corresponds to the final branches
  decompE <- apply(get.edges(backbone$B, E(backbone$B)), 1, function(e){
    c <- decompV[e[degrees[e] %in% c(1, 2)]][1]
    if(is.na(c)) c <- decompV[e[1]]
    return(c)
  })

  # Assign all nodes in G exactly to one branch
  BA <- sapply(V(backbone$pine), function(v) decompV[which.min(distToB[,v])])

  # Modify the one-hot-encoding for a more probabilistic interpretation, if required
  if(prob){
    BA_prob <- matrix(numeric(length(unique(BA)) * ncol(distToB)), nrow = ncol(distToB))
    for(i in seq(ncol(distToB))){
      if(!is.null(backbone$G)) f <- prop.table(table(BA[neighbors(backbone$G, V(backbone$G)[i])]))
      else f <- prop.table(table(BA[neighbors(G, V(G)[i])]))
      BA_prob[i, as.integer(names(f))] <- f
    }
  }

  # If no probabilistic interpretation is required, modify the branch assignment to a one-hot-encoding
  else BA_prob <- as.matrix(diag(max(BA))[BA,])

  # Provide a color assignment for visualization
  set.seed(42)
  palette <- distinctColorPalette(ncol(BA_prob))
  col <- rgb(t(sapply(seq(nrow(BA_prob)), function(i)
    apply((BA_prob[i,] * t(sapply(palette, function(c) col2rgb(c))) / 255), 2, sum))))
  names(col) <- if(is.null(G)) names(V(BCB$G)) else names(V(G))

  # Add the assignment objects to the backbone
  backbone$branch <- decompE
  backbone$prob <- as.spam(BA_prob)
  backbone$palette <- palette
  backbone$col <- col

  # Return backbone with branch assignment
  return(backbone)
}
