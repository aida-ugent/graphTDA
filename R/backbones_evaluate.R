#' Evaluate one or multiple backbones according to various metrics.
#'
#' @param BL A backbone or list of backbone or igraph objects.
#' @param G  (Optional) The original graph G, assumed to be equal for all backbones.
#'           Must be provided of no object in BL contains G.
#'           Standard is NULL.
#'
#' @return   A data frame with one row per backbone, evaluating it according to the following metrics:
#'                - V = fraction of included nodes
#'                - E = fraction of included edges
#'                - R = goodness of fit
#'                - sigma = smoothness
#'                - cor_act = average commute time preservation
#'                - leaves = number of leaves in the backbone

backbones_evaluate <- function(BL, G=NULL)
{
  # Preliminary checks
  for(idx in 1:length(BL)){
    if(!is.null(BL[[idx]]$G) & is.null(G)) G <- BL[[idx]]$G
    if(is_igraph(BL[[idx]])) BL[[idx]] <- list(B=BL[[idx]])
    else if(is.null(BL[[idx]]$B)) stop("List must consist of either backbone or igraph objects.")
  }
  if(is.null(G)) stop("No backbone contains original graph G. Please provide G.")

  # Determine longest path in original graph
  print("Computing longest path in original graph...")
  dG <- distances(G)
  uv <- which(dG == max(dG[dG != Inf]), arr.ind=TRUE)[1,]

  # Compute center(s) and total variance in original graph
  print("Computing center(s) and total variance in original graph...")
  Comp <- components(G)
  memberNodes <- lapply(1:Comp$no, function(C) V(G)$name[which(Comp$membership == C)])
  rm(Comp)
  Center <- sapply(memberNodes, function(N) N[which.min(apply(dG[N, N], 1, max, na.rm=TRUE))])
  totalVar <- sum(apply(as.matrix(dG[,Center], ncol=length(Center)), 1, min, na.rm=TRUE))

  # Compute average commute times in orginal graph
  print("Computing average commute times in original graph...")
  OmegaG <- matrix(rep(NaN, length(V(G))^2), nrow=length(V(G)))
  rownames(OmegaG) <- names(V(G))
  colnames(OmegaG) <- names(V(G))
  for(N in memberNodes) OmegaG[N, N] <- proxfun(induced_subgraph(G, N), N, N, method="act")
  diag(OmegaG) <- 0

  # Evaluate the backbones
  metrics <- do.call("rbind", lapply(1:length(BL), function(i){

    # Compute projection of G on current backbone
    print(paste0("Computing projection on backbone ", i, "..."))
    GprojB <- project_on_backbone(G, BL[[i]]$B, dG)

    # Compute goodness of fit
    print(paste0("Computing R: goodness of fit of backbone ", i, "..."))
    R <- 1 - sum(dG[cbind(GprojB$addedV, GprojB$connectTo)]) / totalVar

    # Compute smoothness
    print(paste0("Computing sigma: smoothness of backbone ", i, "..."))
    sigma <- dG[uv[1], uv[2]] / distances(GprojB, names(V(G))[uv[1]], names(V(G))[uv[2]])[1, 1]

    # Compute correlation to average commute times of projection on backbone
    print(paste0("Computing cor_act: correlation to average commute times of projection on backbone ", i, "..."))
    OmegaGprojB <- matrix(rep(NaN, length(V(GprojB))^2), nrow=length(V(G)))
    rownames(OmegaGprojB) <- names(V(G))
    colnames(OmegaGprojB) <- names(V(G))
    for(N in memberNodes) OmegaGprojB[N, N] <- proxfun(induced_subgraph(GprojB, N), N, N, method="act")
    diag(OmegaGprojB) <- 0
    cor_act <- cor(as.numeric(OmegaG), as.numeric(OmegaGprojB), use="complete.obs")

    # Return the metrics for this backbone
    return(data.frame(frac_V=length(V(BL[[i]]$B)) / length(V(G)),
                      frac_E=length(E(BL[[i]]$B)) / length(E(G)),
                      R=R,
                      sigma=sigma,
                      cor_act=cor_act,
                      leaves=sum(degree(BL[[i]]$B) == 1)))
  }))

  # Return the metrics for all backbones
  return(metrics)
}
