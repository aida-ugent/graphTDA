# GraphTDA

Topological Data Analysis of graphs. 


# Info

To start:
- Install the dependencies below.
- Open 'graphTDA.Rproj' in Rstudio.
- Build --> Clean and Rebuild.
- Open a tutorial from Rstudio and you are good to go!

# Tutorials

We provide following tutorials to investigate the functionality of our code:


- Persistence.R:

  The purpose of this tutorial is to introduce persistent homology. 
  The code can also be easily adapted to explain the method through custom point cloud data, e.g., for in academic papers.
  Note that is a very brief introduction to how and what we can compute with persistent homology in R.
  It does not serve as a replacement of any formal introduction to persistent homology.
  
  
- Basic_Backbone.R:

  The purpose of this tutorial is to show the basic functionality of the backbone function.
  
  
- Disconnected_with_Cycles.R:

  The purpose of this tutorial is to show how to infer backbones in disconnected graphs, and how to identify 'gaps' in the backbone.
  We recommend to go through the 'Basic_Backbone' and 'Persistence' tutorials first.
  
  
- Cell_Trajectory_Inference.R:

  The purpose of this tutorial is to show how to infer backbones in real-world cell trajectory data.
  We recommend to go through the 'Basic_Backbone' and tutorial first.
  
  
 - Topological_Signatures.R:

  The purpose of this tutorial is to show how we can compute and compare topological signatures of graphs in R.
  This is mostly performed through the TDA library in R.
  However, we provide code to modify igraph objects to be compatible with this library.


# Dependencies

most of these can be installed directly through the following command in R: install.packages("<package_name>") 

Package dependencies:
- "igraph" (https://igraph.org/r/)
- "spam" (https://cran.r-project.org/web/packages/spam/)
- "FNN" (https://cran.r-project.org/web/packages/FNN/)
- "pdist" (https://cran.r-project.org/web/packages/pdist/)
- "dplyr" (https://cran.r-project.org/web/packages/dplyr/)
- "TDA" (https://cran.r-project.org/web/packages/TDA/)
- "linkprediction" (https://cran.r-project.org/web/packages/linkprediction/)
- "randomcoloR" (https://cran.r-project.org/web/packages/randomcoloR/)

Tutorial dependencies:
- "ggplot2"
- "latex2exp"
- "gridExtra"
- "ggpubr"
- "dyndimred"


# Contact

Robin.Vandaele@UGent.be