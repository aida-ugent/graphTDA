% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_edges2D.R
\name{get_edges2D}
\alias{get_edges2D}
\title{Identify the endpoints of a proximity graph G in the 2D Euclidean plane through the data from which it was constructed.}
\usage{
get_edges2D(df, G)
}
\arguments{
\item{df}{Two-dimensional data.frame where the rownames correspond to the vertex names in G (not necessarily in order).}

\item{G}{A proximity graph constructed from the points (represented) in df.}
}
\value{
A four-dimensional data.frame containing the coordinates of the endpoints of in the form (x1, y1, x2, y2).
}
\description{
Identify the endpoints of a proximity graph G in the 2D Euclidean plane through the data from which it was constructed.
}
