#' Compute a backbone, i.e., a simplified graph-structured topology, of a given graph or point cloud data set X.
#'
#' @param X             Data for which a backbone is to be computed in igraph/data.frame/matrix/dist format. or a 'matrix'
#'                      In all but the former case, a proximity graph G is constructed.
#'                      Otherwise G = X is used.
#' @param type          (Optional) Used if X is not in igraph format.
#'                      Must be either 'knn' if kNN graph is to be constructed, or 'rips' if Vietoris-Rips graph is to be constructed.
#'                      Standard is 'knn'
#' @param k             (Optional) Number of neighbors to be used when a kNN graph G is constructed from X.
#'                      Standard is 10.
#' @param eps           (Optional) Distance parameter to be used when a (Vietoris-)Rips graph G is constructed from X.
#'                      Must be provided if X is not in igraph format and type is rips.
#' @param f             (Optional) function values to compute f-pine.
#'                      Boundary coefficients are computed and used if not provided.
#' @param CLOF_vcost    (Optional) vertex valued function used to solve CLOF.
#'                      Standard is vertex betweenness.
#' @param CLOF_ecost    (Optional) edge valued function used to solve CLOF.
#'                      Standard is NULL, and a vertex valued function is used.
#'                      If not NULL, CLOF_vcost is ignored.
#' @param max_leaves    (Optional) Upper bound on the number of leaves to be included in the backbone.
#'                      May also be a vector with a length equal to the number of components of G.
#'                      A number of leaves of at most max_leaves is estimated by an elbow estimator.
#'                      Standard is NA (which is equivalent to +infinity).
#' @param leaves        (Optional) Exact number of leaves to be included in the backbone.
#'                      May also be a vector with a length equal to the number of components of G.
#'                      Note that this number may not be achieved, e.g., if there are not enough leaves in the f-pine.
#' @param preprune      (Optional) Whether to preprune the f-pine, i.e., discard all its leaves once if the cost for CLOF is constant on leaves.
#'                      This may significantly improve computation cost, and one may directly obtain the solution in the original pine by adding arbitrary leaves.
#'                      Standard is TRUE.
#' @param stdize        (Optional) Whether to standardize the cost of each component according to its full cost.
#'                      Standard is FALSE, which may add bias towards representing larger components in the backbone.
#' @param assign        (Optional) Whether to construct a branch assignment, marking for each node x which branch(es) represent x.
#'                      Standard is FALSE.
#' @param prob          (Optional & if assign) whether to assign nodes in G to paths according to a probability measure
#'                      if FALSE, then nodes are assigned to exactly one path through a one-hot-encoding.
#'                      Standard is TRUE.
#' @param constrain_mem (Optional & f is not provided) if TRUE, unnecessary distances are never computed.
#'                      Standard is FALSE, which currently results in faster computation (with cpp running in the background).
#'
#' @return A list storing the following items:
#'                - B: the backbone in X encoded as an igraph object.
#'                - (if constructed) G: the proximity constructed from X (if not stored, then G = X).
#'                - (if computed) f: the boundary coefficients of the (constructed) graph.
#'                - pine: an f-pine in G.
#'                - membership: the membership of the nodes according to the connected components in in G.
#'                - cost: a data.frame corresponding to the obtained cost according to the number of leaves and components.
#'                - full_cost: a vector of the maximal obtainable cost of each component.
#'                - (if CLOF optimizes with vertex-valued function) includedV: a list of lists (one for each connected component) containing the added nodes of G at each iteration.
#'                - (if CLOF optimizes with edge-valued function) includedE: a list of lists (one for each connected component) containing the added edges of G at each iteration.
#'                - (if assign) branch: a vector clustering the edges of B into paths.
#'                - (if assign) prob: a n x k matrix, where n = |V(G)| and k is the number of resulting branches in B, corresponding to a probabilistic assignment of the nodes in G according to these paths.
#'                - (if assign) palette: a vector of hex code colors, one for each path.
#'                - (if assign) col: a vector of |V(G)| hex code colors, interpolating between the palette colors according to prob.

backbone <- function(X,
                     type="knn",
                     k=10,
                     eps=NULL,
                     f=NULL,
                     CLOF_vcost=betweenness,
                     CLOF_ecost=NULL,
                     max_leaves=NA,
                     leaves=NA,
                     preprune=TRUE,
                     stdize=FALSE,
                     assign=FALSE,
                     prob=TRUE,
                     constrain_mem=FALSE)
{
  # Preliminary checks
  if(!(is_igraph(X) | is.data.frame(X) | is.matrix(X) | class(X)=="dist")) stop("X must be either an igraph object, a data frame, a matrix, or distance object.")
  if(!type %in% c("knn", "rips")) stop("type must be 'knn' or 'rips', standard is 'knn'.")
  if(type == "rips" & is.null(eps)) stop("type is 'rips', please provide a value for 'eps'.")
  if(!class(CLOF_vcost) == "function") stop("CLOF_vcost must be a valid vertex-valued function.")
  if(!is.null(CLOF_ecost) & !class(CLOF_ecost) == "function") stop("CLOF_ecost must be a valid edge-valued function, or NULL to use CLOF_vcost.")

  # Track complete computation time
  Old <- Sys.time()

  # Check if proximity graph needs to be constructed
  G <-  NULL
  store_G <- FALSE
  if(!is_igraph(X)){
    store_G <- TRUE
    old <- Sys.time()
    if(type == "knn") G <- kNNgraph(X, k)
    else G <- VRgraph(X, eps)
    new <- Sys.time() - old
    print(paste("Proximity graph constructed in", round(new, 3), attr(new, "units")))
  }
  else if(is.null(V(X)$name)) V(X)$name <- as.character(seq(length(V(X))))

  # Check if function f for computing f-pine is given
  # Otherwise compute boundary coefficients
  store_f <- FALSE
  if(is.null(f)){
    store_f <- TRUE
    old <- Sys.time()
    f <- boundary_coefficient(if(is.null(G)) X else G, constrain_mem=constrain_mem)
    new <- Sys.time() - old
    print(paste("Boundary coefficients obtained in", round(new, 3), attr(new, "units")))
  }

  # Compute f-pine
  old <- Sys.time()
  pine <- mst(if(is.null(G)) X else G, weights = sum_weights(if(is.null(G)) X else G, f))
  new <- Sys.time() - old
  print(paste("f-pine obtained in", round(new, 3), attr(new, "units")))

  # Compute backbone through CLOF
  old <- Sys.time()
  CLOF_leaves <- if(any(is.na(leaves))) max_leaves + 1 else leaves # max_leaves + 1 makes max_leaves possible as elbow
  backbone <- CLOF_forest(pine, leaves=CLOF_leaves, CLOF_vcost=CLOF_vcost, CLOF_ecost=CLOF_ecost, preprune=preprune)
  backbone$pine <- pine
  rm(pine)
  if(any(is.na(leaves))){
    backbone <- get_new_leaves(backbone, leaves=get_elbow(backbone$cost)[,1], G=if(is.null(G)) X else G, assign=FALSE, stdize=stdize)
    new <- Sys.time() - old
    print(paste("Backbone obtained in", round(new, 3), attr(new, "units")))
  }
  else if(length(leaves) == 1 & length(levels(factor(backbone$membership))) > 1){
    backbone <- get_new_leaves(backbone, leaves=leaves, G=if(is.null(G)) X else G, assign=FALSE, stdize=stdize)
    new <- Sys.time() - old
    }
  print(paste("Backbone obtained in", round(new, 3), attr(new, "units")))

  # Check for branch assignment
  if(assign){
    old <- Sys.time()
    backbone <- decompose(backbone, G=if(is.null(G)) X else G, prob=prob)
    new <- Sys.time() - old
    print(paste("Branch assignment obtained in", round(new, 3), attr(new, "units")))
  }

  # Construct final output object
  if(store_f) backbone$f <- f
  if(store_G) backbone$G <- G

  # Print total computation time
  new <- Sys.time() - Old
  print("--------------------------------------------")
  print(paste("Backbone pipeline conducted in", round(new, 3), attr(new, "units")))

  # Return backbone
  return(backbone)
}
