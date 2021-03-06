% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CLOF_forest.R
\name{CLOF_forest}
\alias{CLOF_forest}
\title{Compute a backbone in a forest graph G through Constrained Leaves Optimal subForest (CLOF)}
\usage{
CLOF_forest(
  G,
  leaves = NA,
  CLOF_vcost = betweenness,
  CLOF_ecost = NULL,
  preprune = TRUE
)
}
\arguments{
\item{G}{forest graph from which a backbone is to be retrieved}

\item{leaves}{(Optional) Maximal number of leaves to be included in the backbone.
May also be a vector with a length equal to the number of components of G.
Note that this number may not be achieved, e.g., if there are not enough leaves in G.}

\item{CLOF_vcost}{(Optional) vertex valued function used to solve CLOF.
Standard is vertex betweenness.}

\item{CLOF_ecost}{(Optional) edge valued function used to solve CLOF.
Standard is NULL, and a vertex valued function is used.
If not NULL, CLOF_vcost is ignored.}

\item{preprune}{(Optional) Whether to preprune the f-pine, i.e., discard all its leaves once if the cost for CLOF is constant on leaves.
This may significantly improve computation cost, and one may directly obtain the solution in the original pine by adding arbitrary leaves.
Standard is TRUE.}
}
\value{
A list containing the following items
               - B: the solution, i.e., backbone, obtained through CLOF
               - membership: the membership of the nodes in G according to the components of G
               - cost: a data.frame corresponding to the obtained cost according to the number of leaves and the component
               - full_cost: a vector of the full cost of each component
               - (if CLOF_vcost is used) includedV: a list of lists, one for each component, marking the nodes in the sequential solutions for CLOF
               - (if CLOF_ecost is used) includedE: a list of lists, one for each component, marking the edges in the sequential solutions for CLOF
}
\description{
Compute a backbone in a forest graph G through Constrained Leaves Optimal subForest (CLOF)
}
