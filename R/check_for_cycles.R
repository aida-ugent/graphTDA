#' Check for any gaps in the obtained backbone by means of persistent homology.
#' As R is currently not ideal for persistent homology computation, this method can quikly lead to memory issues.
#'
#' @param backbone A backbone object for which gaps need to be identified.
#' @param G        (Optional) the original graph on top of the backbone.
#'                 Must be included if not contained in the backbone object.
#'                 Standard is NULL.
#'
#' @return The backbone object (list) with the following (new) entries:
#'              - diagram: the 0-dimensional and 1-dimensional persistence diagrams in a matrix.
#'              - repCycles: list of edges that represent the 1-dimensional holes in diagram.

check_for_cycles <- function(backbone, G=NULL)
{
  # Preliminary checks
  if(is.null(backbone$G) & is.null(G)) stop("backbone does not include original graph G, please provide G.")

  # Start time tracking
  old <- Sys.time()

  # Compute persistent homology
  DG <- distances(if(is.null(G)) backbone$G else G, v=V(backbone$B)$name, to=V(backbone$B)$name)
  persistence <- TDA::ripsDiag(X=DG, maxdimension=1, maxscale=Inf, dist="arbitrary", location=TRUE, library="Dionysus")

  # Stop time tracking and print time
  new <- Sys.time() - old
  print(paste("Persistent homology computation performed in", round(new, 3), attr(new, "units")))

  # Add diagram and cycles to backbone object
  backbone$diagram <- persistence$diagram
  backbone$repCycles <- lapply(persistence$cycleLocation[persistence$diagram[,"dimension"]==1], function(CL){
    return(cbind(V(backbone$B)$name[CL[,1]], V(backbone$B)$name[CL[,2]]))
  })

  # Return backbone object
  return(backbone)
}
