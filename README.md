This repository contains scripts and markdown files used in a graduate research project in collaboration with the Knight Campus Graduate Internship Program and OHSU (Dr. Maxson and Dr. Braun labs). The goal of this project is to use trajectory analysis and pseudotime to identify differences in monocyte differentiation between healthy patients and those with chronic myelomonocytic leukemia (CMML), as well as acute myeloid leukemia (AML). This is a work in progress. 

## KNN_search.R

This is a *k* nearest neighbor search script which improves speeds by mapping reference cell UMAP values to a matrix and using the query UMAP values to search the matrix for the nearest reference neighbors in a spiral search pattern. This approach is fast but memory intensive. Time depends on the density of reference cells. Low density will be slower as it will take longer to find the desired number of neighbors. A more dense reference set will be faster since neighbors are close together. 

Functions:
generate.map()
find.neighbors()


## cmml_cord_bio_reps.Rmd      

This is a markdown file summarizing the simulation of biological replicates from normal and poisson distributions based on the CMML and cord pseudotime values estimated from the KNN_search.R script. It also includes simulated biological replicates from the pseudotime values themselves. 
