# The purpose of this tutorial is to introduce persistent homology.
# The code can also be easily adapted to explain the method through custom point cloud data, e.g., for in academic papers.
# Note that is a very brief introduction to how and what we can compute with persistent homology in R.
# It does not serve as a replacement of any formal introduction to persistent homology.
# We first import the necessary libraries.

library("ggplot2") # plotting
library("latex2exp") # LaTeX text in figures
library("gridExtra") # arrange multiple plots

# We sample and plot a point cloud data, of which the underlying topological model equals two disconnected cycles.

set.seed(42)
n <- 50
theta <- runif(2 * n, 0, 2 * pi)
ycenters <- c(0, 3)
df <- data.frame(x=cos(theta) + rnorm(2 * n, sd=0.1), y=rep(ycenters, each=n) + sin(theta) + rnorm(2 * n, sd=0.1))
ggplot(df, aes(x=x, y=y)) +
  geom_point(size=1) +
  coord_fixed() +
  theme_bw()

# Being a two-dimensional point cloud, there are two types of topological features we may obtain from the Euclidean distances between points.
# The first identify zero-dimensional holes (H0), or connected components, in the underlying model.
# The second identify one-dimensional holes (H1), or cycles, in the underlying model.
# We can now compute persistent homology as follows

filtration <- ripsFiltration(df, maxdimension=1, maxscale=Inf)
diag <- filtrationDiag(filtration, maxdimension=1)

# We can visualize the result of persistent homology in two ways.
# The first is by means of persistence barcodes, where a bar from b to d represent a topological feature that persists from time b to d.

op <- par(mar = c(3.25, 1, 1, 1))
plot.diagram(diag[["diagram"]], diagLim=c(0,  2.5), barcode=TRUE)
legend(2, 100, legend=c("H0", "H1"), col=c("black", "red"), lty=1, lwd=2, box.lty=0)
par(op)

# The second is by means of persistence diagrams, where a point (b, d) represent a topological feature that persists from time b to d.

op <- par(mar = c(3.25, 3.25, 1, 1))
plot.diagram(diag[["diagram"]], diagLim=c(0, 2.5))
legend(2, 1, legend=c("H0", "H1"), col=c("black", "red"), pch=c(19, 2), pt.lwd=2, box.lty=0)
par(op)

# Through both visualizations two 'persisting' components for both connected components (H0) and cycles (H1).
# These are topological features (holes) of which the death time d is much larger than the birth time b.
# These observations are consistent with the underlying model of the data.
# To explain what is meant by these times, one requires more knowledge of what persistent homology is actually computing.
# More specifically, it computes the appearance and disappearance of holes in a time varying simplicial complex
# Such simplicial complex can be seen as a generalization of a graph, that also allows for triangles, tetrahedra, or in general: higher-dimensional simplices.
# In our case, the simplicial complex at time t is specified by adding a simplex on each subset of points with diameter at most t.
# Some examples of such complexes at different times, can be visualized as follows.
# Note that for larger diameters, it can take some time for the 'geom_polygon' function to finish.

epsilons <- c(0.01, 0.5, 0.75, 1, 1.6, 2.5)
simpPlots <- list()
for(idx in 1:length(epsilons)){
  print(paste("Constructing plot for time parameter epsilon =", epsilons[idx]))
  if(idx == 1){
    m <- -Inf
    nodes <- data.frame(x=numeric(0), y=numeric(0))
    edges <- data.frame(x1=numeric(0), y1=numeric(0), x2=numeric(0), y2=numeric(0))
    triangles  <- data.frame(id=integer(0), x=numeric(0), y=numeric(0))
    tId <- 0
  }
  for(simp_idx in 1:length(filtration$cmplx)){
    if(filtration$values[simp_idx] <= epsilons[idx] & filtration$values[simp_idx] > m){
      simp <- filtration$cmplx[[simp_idx]]
      if(length(simp)==1) nodes[nrow(nodes) + 1,] <- df[simp,]
      else if(length(simp)==2) edges[nrow(edges) + 1,] <- as.numeric(c(df[simp[1],], df[simp[2],]))
      else if(length(simp==3)){
        tId <- tId + 1
        triangles[(nrow(triangles) + 1):(nrow(triangles) + 3),] <- cbind(tId, df[simp,])
      }
    }
  }
  m <- epsilons[idx]
  simpPlots[[length(simpPlots) + 1]] <- ggplot(nodes, aes(x=x, y=y)) +
    geom_segment(data=edges, aes(x=x1, y=y1, xend=x2, yend=y2), color="black", size=0.5, alpha=0.5) +
    geom_point(size=1.5, alpha=0.75) +
    geom_polygon(data=triangles, aes(group=id), fill="green", alpha=0.25) +
    geom_segment(data=edges, aes(x=x1, y=y1, xend=x2, yend=y2), color="black", size=0.5, alpha=0.05) +
    geom_point(size=1.5, alpha=0.25) +
    coord_fixed() +
    theme_bw() +
    ggtitle(TeX(sprintf("$\\epsilon = %g$", epsilons[idx]))) +
    scale_x_continuous(breaks = c(-1, 0, 1)) +
    theme(plot.title = element_text(hjust = 0.5, size=22), text = element_text(size=20))
}

grid.arrange(grobs=simpPlots, ncol=length(epsilons))
