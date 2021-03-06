% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bounded_hop_pairG.R
\name{bounded_hop_pairG}
\alias{bounded_hop_pairG}
\title{Compute the hop-k-approximation of a pairwise distance matrix of a given graph G.
We suggest a future implementation directly in (and called from) C++.}
\usage{
bounded_hop_pairG(G, k = 2, constrain_mem = FALSE)
}
\arguments{
\item{G}{The graph from which the distance matrix is to be approximated.}

\item{k}{(Optional) Maximum allowed hops between pairs of nodes for which the weighted distances are to be computed.
Standard is k = 2.}

\item{constrain_mem}{(Optional) If TRUE and pairD is missing, the entire distance matrix on G as is computed first, from which the hop-k-approximation is extracted afterwards.
If FALSE, entries that are not needed in pairD are never computed, requiring less storage.
Standard is FALSE, which leads to faster results with our current implementation.}
}
\value{
A sparse matrix storing the weighted distance D_uv between nodes u and v if the unweighted distance between u and v is less than k.
}
\description{
Compute the hop-k-approximation of a pairwise distance matrix of a given graph G.
We suggest a future implementation directly in (and called from) C++.
}
