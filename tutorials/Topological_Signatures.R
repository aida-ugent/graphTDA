# The purpose of this tutorial is to show how we can compute and compare topological signatures of graphs in R.
# This is mostly performed through the TDA library in R.
# However, we provide code to modify igraph objects to be compatible with this library.
# We first import the necessary libraries.

devtools::load_all() # load graphTDA
library("ggplot2") # plotting
library("ggpubr") # plotting in grid

# We define a list of shapes, that will make up the underlying models of the point cloud data sets we will compare.
# These shapes will be metric trees, defined by a collection of segments between start and end points.

Ishape <- list(list(start=c(-1, -1), end=c(0, 0)),
               list(start=c(0, 0), end=c(1,  1)))

Yshape <- list(list(start=c(0, 0), end=c(1, 1)),
               list(start=c(0, 0), end=c(-1, 1)),
               list(start=c(0, 0), end=c(0, -sqrt(2))))

Xshape <- list(list(start=c(0, 0), end=c(1, 1)),
               list(start=c(0, 0), end=c(-1, 1)),
               list(start=c(0, 0), end=c(1, -1)),
               list(start=c(0, 0), end=c(-1, -1)))

Hshape <- list(list(start=c(-1/2, 0), end=c(0, 0)),
               list(start=c(1/2, 0), end=c(0, 0)),
               list(start=c(-1/2, 0), end=c(-1, 1)),
               list(start=c(-1/2, 0), end=c(-1, -1)),
               list(start=c(1/2, 0), end=c(1, 1)),
               list(start=c(1/2, 0), end=c(1, -1)))

shapes <- list("I"=Ishape, "Y"=Yshape, "X"=Xshape, "H"=Hshape)

# Before yhe models, we first construct their (ground-truth) topological signatures.
# These will be derived as the persistence diagrams of the sublevel filtrations defined by the model's normalized centrality functions.

trueGraphs <- lapply(shapes, function(shape){ # true discrete graph/tree-structured topological representations
  G <- t(sapply(shape, function(shape) do.call("c", shape)))
  weight <- apply(G, 1, function(r) norm(r[1:2] - r[3:4], type="2"))
  G <- graph_from_data_frame(data.frame(from=(apply(matrix(G[,1:2], ncol=2), 1, function(r) paste(r[1], r[2]))),
                                        to=apply(matrix(G[,3:4], ncol=2), 1, function(r) paste(r[1], r[2]))),
                             directed=FALSE)
  E(G)$weight <- weight
  G$NC <- normalized_centrality(G)
  return(G)
})

trueDiagrams <- lapply(trueGraphs, function(G) graph_persistence(S=G, f=G$NC)) # true topological signatures
names(trueDiagrams) <- names(shapes)

# To illustrate the topologies, we will first uniformly sample from them without noise.

npoints <- 600 # points per model, exact number of points may slightly differ due to rounding

set.seed(42)
cleanDataSets <- lapply(names(shapes), function(name){
  shape <- shapes[[name]]
  cleanData <- data.frame(x=numeric(), y=numeric(), NC=numeric())
  pointsPerBranch <- sapply(shape, function(l) norm(l$start - l$end, type="2"))
  pointsPerBranch <- round(pointsPerBranch / sum(pointsPerBranch) * npoints)
  for(idx in 1:length(shape)){
    branch <- shape[[idx]]
    t <- runif(pointsPerBranch[[idx]])
    pointsOnBranchX <- t * branch$start[1] + (1 - t) * branch$end[1]
    pointsOnBranchY <- t * branch$start[2] + (1 - t) * branch$end[2]
    NC <- t * trueGraphs[[name]]$NC[paste(branch$start[1], branch$start[2])] +
      (1 - t) * trueGraphs[[name]]$NC[paste(branch$end[1], branch$end[2])]
    cleanData <- rbind(cleanData, data.frame(x=pointsOnBranchX, y=pointsOnBranchY, NC=NC))
  }
  return(cleanData)
})
names(cleanDataSets) <- names(shapes)

# We can now visualize the clean data, and hence, the models, as well as the filtration defining functions, as follows

plotsOfCleanData <- lapply(names(cleanDataSets), function(name){
  ggplot(cleanDataSets[[name]], aes(x=x, y=y, col=NC)) +
    geom_point(size=2.5, col="black") +
    geom_point(size=2) +
    scale_colour_gradientn(colours=topo.colors(7)) +
    ggtitle(paste(name, "ground truth")) +
    labs(col="normalized centrality") +
    theme_bw() +
    theme(plot.title=element_text(hjust = 0.5)) +
    coord_fixed()
})
ggpubr::ggarrange(plotlist=plotsOfCleanData, nrow=2, ncol=2)

# From each model, we sample three noisy data sets, by adding random noise and rotations, and centering to (0, 0)

nDistortions <- 3
sigma <- 0.005

noisyDataSets <- lapply(cleanDataSets, function(cleanData){
  lapply(seq_len(nDistortions), function(distortionIdx){
    noise <- MASS::mvrnorm(nrow(cleanData), mu=c(0,0), Sigma=diag(2)*sigma)
    rotation <- mixAK::rRotationMatrix(1, 2)
    noisyData <- data.frame(scale(t(apply(cleanData[,1:2], 1, function(row) rotation %*% row)) + noise, scale=FALSE))
    colnames(noisyData) <- colnames(cleanData)[1:2]
    return(noisyData)
  })
})
names(noisyDataSets) <- names(shapes)

# We can visualize the noisy samples as follows

noisyLim <- c(min(unlist(noisyDataSets)), max(unlist(noisyDataSets)))
plotsOfNoisyData <- lapply(1:length(noisyDataSets), function(idx1){
  lapply(1:length(noisyDataSets[[idx1]]), function(idx2){
    ggplot(noisyDataSets[[idx1]][[idx2]], aes(x=x, y=y)) +
      geom_point(size=0.25) +
      xlim(noisyLim) +
      ylim(noisyLim) +
      theme_bw() +
      coord_fixed()
  })
})
ggpubr::ggarrange(plotlist=do.call("c", plotsOfNoisyData), common.legend=TRUE, nrow=length(shapes), ncol=nDistortions)

# From each noisy data set, we now construct the minimum spanning tree (MST) as well as the resulting normalized centralities

mstOfDataSets <- lapply(noisyDataSets, function(noisyDatasForShape){
  lapply(noisyDatasForShape, function(noisyData){
    G <- emstreeR::ComputeMST(noisyData[,1:2], verbose=FALSE)[,c("from", "to", "distance")]
    G <- data.frame(G)
    colnames(G)[3] <- "weight"
    G <- graph_from_data_frame(G, directed=FALSE)
    G$NC <- normalized_centrality(G)
    return(G)
  })
})

# We can now plot all clean data as well as the MSTs with normalized centralities as follows.
# Note that the edges of the MSTs may be difficult to visualize for smaller figure sizes.

plotsOfMST <- lapply(1:length(noisyDataSets), function(idx1){
  lapply(1:length(noisyDataSets[[idx1]]), function(idx2){
    euclideanEdges <- get_edges2D(noisyDataSets[[idx1]][[idx2]], mstOfDataSets[[idx1]][[idx2]])
    ggplot(noisyDataSets[[idx1]][[idx2]][names(V(mstOfDataSets[[idx1]][[idx2]])),], aes(x=x, y=y)) +
      geom_segment(data=euclideanEdges, aes(x=x1, y=y1, xend=x2, yend=y2), color="red", size=0.5) +
      geom_point(size=1.5, col="black") +
      geom_point(size=1, aes(col=mstOfDataSets[[idx1]][[idx2]]$NC)) +
      ggtitle(paste(names(noisyDataSets)[[idx1]], "sample", idx2)) +
      scale_colour_gradientn(colours=topo.colors(7)) +
      labs(col="Normalized Centrality") +
      xlim(noisyLim) +
      ylim(noisyLim) +
      theme_bw() +
      theme(plot.title=element_text(hjust = 0.5)) +
      coord_fixed()
  })
})
ggpubr::ggarrange(plotlist=do.call("c", lapply(1:length(shapes), function(idx) c(list(plotsOfCleanData[[idx]]), plotsOfMST[[idx]]))),
                  common.legend=TRUE, nrow=length(shapes), ncol=nDistortions + 1)

# We now computed the topological signatures, i.e., the zero-dimensional persistence diagram, for each MST.
# This through the sublevel filtrations defined by the normalized centrality functions.

diagramsOfMST <- list()
progress <- 0
for(idx1 in 1:length(mstOfDataSets)){
  diagramsOfMST[[idx1]] <- list()
  for(idx2 in 1:length(mstOfDataSets[[idx1]])){
    print(paste0("progress: ", round(100 * progress / (length(mstOfDataSets) * length(mstOfDataSets[[1]])), 2), "%"))
    G <- mstOfDataSets[[idx1]][[idx2]]
    diagramsOfMST[[idx1]][[idx2]] <- graph_persistence(S=G, f=G$NC)
    progress <- progress + 1
  }
  if(idx1==length(mstOfDataSets)) print(paste0("progress: 100%"))
}

# We can now plot all true and and empirical diagrams as follows.

colors <- distinctColorPalette(length(noisyDataSets))
opt <- par(mfrow=c(length(noisyDataSets), length(noisyDataSets[[1]]) + 1), mar=c(3, 4, 2, 2))
for(idx1 in 1:length(noisyDataSets)){
  TDA::plot.diagram(trueDiagrams[[idx1]], main=paste(names(shapes)[idx1], "ground truth"), diagLim=c(0,1.1), col=colors[idx1])
  for(idx2 in 1:length(noisyDataSets[[idx1]])){
    TDA::plot.diagram(diagramsOfMST[[idx1]][[idx2]], main=paste(names(shapes)[idx1], "sample", idx2), diagLim=c(0,1.1), col=colors[idx1])
  }
}; opt

# To compare the different topological signatures, we compute the bottleneck distances between them.

pairwiseBottleNecksMST <- matrix(numeric((length(noisyDataSets)*(length(noisyDataSets[[1]]) + 1))^2),
                                 nrow=length(noisyDataSets)*(length(noisyDataSets[[1]]) + 1))
for(idx1 in 1:(nrow(pairwiseBottleNecksMST) -  1)){
  for(idx2 in idx1:(nrow(pairwiseBottleNecksMST))){
    shapeIndex1 <- ceiling(idx1 / (nDistortions + 1))
    dataIndex1 <- idx1 - ((shapeIndex1 - 1) * (nDistortions + 1))
    if(dataIndex1 == 1) dgm1 <- trueDiagrams[[shapeIndex1]]
    else dgm1 <- diagramsOfMST[[shapeIndex1]][[dataIndex1 - 1]]
    shapeIndex2 <- ceiling(idx2 / (nDistortions + 1))
    dataIndex2 <- idx2 - ((shapeIndex2 - 1) * (nDistortions + 1))
    if(dataIndex2 == 1) dgm2 <- trueDiagrams[[shapeIndex2]]
    else dgm2 <- diagramsOfMST[[shapeIndex2]][[dataIndex2 - 1]]
    pairwiseBottleNecksMST[idx1, idx2] <- TDA::bottleneck(dgm1, dgm2, dimension=0)
    pairwiseBottleNecksMST[idx2, idx1] <- pairwiseBottleNecksMST[idx1, idx2]
  }
}

# We visualize the bottleneck distances through a heatmap of the matrix, as follows.

ggplot(reshape2::melt(pairwiseBottleNecksMST), aes(Var1, Var2, fill=value)) +
  geom_raster() +
  scale_x_continuous(breaks=1 + seq(0, nrow(pairwiseBottleNecksMST) - length(shapes), length.out=length(shapes)), labels=names(shapes)) +
  scale_y_continuous(breaks=1 + seq(0, nrow(pairwiseBottleNecksMST) - length(shapes), length.out=length(shapes)),
                     trans="reverse", labels=names(shapes)) +
  xlab("") +
  ylab("") +
  labs(fill="BN Dist.") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5), text=element_text(size=15)) +
  coord_fixed()

# Finally, we can also chart the data sets through a MDS plot of this distance matrix.

fitBottleNecksMST <- data.frame(cmdscale(pairwiseBottleNecksMST, k=2))
fitBottleNecksMST[,3] <- rep(names(noisyDataSets), each=length(noisyDataSets[[1]]) + 1)
colnames(fitBottleNecksMST) <- c("x", "y", "shape")
ggplot(fitBottleNecksMST, aes(x=x, y=y, col=shape)) +
  geom_point(size=4.5, col="white") +
  geom_point(size=4) +
  geom_point(data=fitBottleNecksMST[1 + seq(0, nrow(fitBottleNecksMST) - length(shapes), length.out=length(shapes)),], size=5, col="black") +
  geom_point(data=fitBottleNecksMST[1 + seq(0, nrow(fitBottleNecksMST) - length(shapes), length.out=length(shapes)),], size=4) +
  scale_color_manual(values=colors[order(names(noisyDataSets))]) +
  theme_bw() +
  coord_fixed()
