#' Get a new backbone with a different amount of leaves
#'
#' @param backbone A backbone object to be converted to a new backbone object with a different amount of leaves
#' @param leaves   Vector of integers, denoting the number of leaves to be included for each component in the pine.
#'                 Must have length 1 or length equal to the number of connected components included in backbone.
#'                 If length 1, we optimize for the total number of leaves equal to the number specified.
#'                 If length equal to the number of connected components, we extract the specified number of leaves for each component.
#' @param G        (Optional) the original graph on top of the backbone.
#'                 Must be included if not contained in the backbone object, and either assign is TRUE or the backbone already had a branch assignment.
#'                 Standard is NULL.
#' @param assign   (Optional) Whether to produce a branch assignment of the resulting new backbone.
#'                 Will always be treated as TRUE if the backbone already had a branch assignment.
#'                 Standard is FALSE.
#' @param prob     (Optional) Whether to assign nodes in G to paths according to a 'probability'.
#'                 If FALSE, then nodes are assigned to exactly one path through a one-hot-encoding.
#'                 Standard is TRUE.
#' @param stdize   (Optional) Whether to standardize the cost of each component according to its full cost.
#'                 Standard is FALSE, which may add bias towards representing larger components in the backbone.
#'
#' @return a backbone object with the new number of leaves

get_new_leaves <- function(backbone,
                           leaves,
                           G=NULL,
                           assign=FALSE,
                           prob=TRUE,
                           stdize=FALSE)
{
  # Preliminary checks
  if(is.null(backbone$G) & is.null(G) & (assign | !is.null(backbone$branch))) stop("backbone does not include original graph G, please provide G.")
  if(!is.null(backbone$includedV) & !length(leaves) %in% c(1, length(backbone$includedV))) stop("length of leaves must 1 or the number of components in the original graph.")
  if(!is.null(backbone$includedE) & !length(leaves) %in% c(1, length(backbone$includedE))) stop("length of leaves must 1 or the number of components in the original graph.")

  # Check if the relative cost (which is already standardized) needs to be adjusted
  if(stdize) adjust <- rep(1, length(backbone$includedV))
  else adjust <- backbone$full_cost

  # Case 1: CLOF was optimized for a vertex-valued function
  if(!is.null(backbone$includedV)){

    # Case 1a: if length(leaves) < components, we optimize over all possible combinations of leaves per component which sum to the specified number
    if(length(leaves) < length(backbone$includedV)){
      comp_rows <- lapply(1:length(backbone$includedV), function(c) which(backbone$cost[,"component"]==c))
      leave_rows <- lapply(1:leaves, function(l) which(backbone$cost[,"leaves"]==(if(l==1) 0 else l)))
      names(leave_rows) <- 1:leaves
      leavecombs <- get_leave_combs(length(backbone$includedV), leaves) # all possible combinations
      cost <- apply(leavecombs, 1, function(r){
        sum(sapply(1:length(r), function(i){
          cr <- intersect(comp_rows[[i]], leave_rows[[as.character(r[i])]])
          if(length(cr) == 0) return(0)
          backbone$cost[cr, "cost"] * adjust[i]
        }))
      })
      leaves <- leavecombs[which.max(cost),]
      rm(comp_rows, leave_rows, leavecombs, cost)
    }

    # Case 1(a->)b: now length(leaves) == components, the result may be directly obtained
    nodes <- character()
    for(c in seq(length(leaves))){
      if(is.na(leaves[c])) nodes <- c(nodes, unlist(backbone$includedV[[c]][min(length(backbone$includedV[[c]]), 2)]))
      else if(leaves[c] == 1) nodes <- c(nodes, unlist(backbone$includedV[[c]][1]))
      else if(leaves[c] > 1) nodes <- c(nodes, unlist(backbone$includedV[[c]][2:leaves[c]]))
    }

    # Construct the new backbone object
    backbone$B <- induced.subgraph(backbone$pine, nodes)
  }

  # Case 2: CLOF was optimized for an edge-valued function
  else{

    # Case 2a: if length(leaves) < components, we optimize over all possible combinations of leaves per component which sum to the specified number
    if(length(leaves) < length(backbone$includedE)){
      comp_rows <- lapply(1:length(backbone$includedE), function(c) which(backbone$cost[,"component"]==c))
      if(leaves %in% c(2, 3)) poss_leaves <- c(0, leaves)
      else poss_leaves <- c(0, 2:(leaves - 2), leaves)
      leave_rows <- lapply(poss_leaves, function(l) which(backbone$cost[,"leaves"]==l))
      names(leave_rows) <- poss_leaves
      leavecombs <- get_leave_combs(length(backbone$includedE), leaves)
      leavecombs <- leavecombs[!apply(leavecombs, 1, function(r) any(r==1)),]
      cost <- apply(leavecombs, 1, function(r){
        sum(sapply(1:length(r), function(i){
          cr <- intersect(comp_rows[[i]], leave_rows[[as.character(r[i])]])
          if(length(cr) == 0) return(0)
          backbone$cost[cr, "cost"] * adjust[i]
        }))
      })
      leaves <- leavecombs[which.max(cost),]
      rm(comp_rows, leave_rows, leavecombs, cost)
    }

    # Case 2(a->)b: now length(leaves) == components, the result may be directly obtained
    edges <- character()
    for(c in seq(length(leaves))){
      if(is.na(leaves[c]) & length(backbone$includedE[[c]]) > 1) edges <- c(edges, unlist(backbone$includedE[[c]][[1]]))
      else if(leaves[c] > 1) edges <- c(edges, unlist(backbone$includedE[[c]][1:leaves[c]]))
    }

    # Construct the new backbone object
    backbone$B <- subgraph.edges(backbone$pine, get.edge.ids(backbone$pine, edges))
  }

  # Perform (new) branch assignment if necessary
  if(assign | !is.null(backbone$branch)) backbone <- decompose(backbone, G=G, prob=prob)

  # Return the new backbone
  return(backbone)
}
