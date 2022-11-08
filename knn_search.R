library(tidyverse)
library(Signac)
#library(comprehenr)

# adjust these for scaling
a <- 20 # point addition. 10 uses 2.7 GB, 20 uses 6.5 GB. 
m <- 1000 # point multiplier
neighbors <- 10 # number of neighbors to find per cell. 

# cell technique
# cell must be first column, followed by
# umap1 and umap 2. multiply all dimensions by 1000 to convert from decimal
# to integer.  

generate.map <- function(ref) {
  
  print("Generating reference map...")
  dim1 <- (max(ref[,1]) + a) * m + 10
  dim2 <- (max(ref[,2]) + a) * m + 10
  map <- matrix("X", nrow=dim1, ncol=dim2)
  #map <- matrix("X", nrow=dim1, ncol=dim2)
  r <- nrow(ref)
  
  print("Populating map with reference cells...")
  for (row in 1:r) {
    if (row%%1000 == 0) {
      print(paste("Added", row, "cells to map."))
    }
    x <- ref[row, 1] + a
    y <- ref[row, 2] + a
    x = x * m
    y = y * m
    map[x, y] <- rownames(ref)[row] #ref[row,1]
    
  }
  print("Done.")
  print(paste("Added", row, "cells to map."))
  return(map)
}


# find the k nearest neighbors of each query cell in the reference map. 
# currently using about 4GB of memory.  
find.neighbors <- function(ref.map, query, k) {
  
  print("Finding neighbors...")
  r <- nrow(query)
  query_neighbors <- data.frame(query_cell=c("query_cells"),
                                ref_neighbor=c("ref_neighbors"), 
                                euclidean_dist=c("e_dist"),
                                )
  
  for (row in 1:r) {
    print(paste("Working on cell",row))
    x <- (query[row, 1] + a) * m
    y <- (query[row, 2] + a) * m
    i = 1
    
    
    k_found = 0
    current_pos = ref.map[x, y]
    sign = -1
    is_x = FALSE
    
    while (k_found < k) {
      
      if (length(current_pos) == 0){
        # somehow ended up out of bounds from the ref_map, return.
        print("out of bounds of reference map.")
        print(c(x,y))
        break
      }
      
      if (current_pos != "X" & k_found == 0) {
        # this is here for if a query cell and reference cell happen to already
        # have the same coordinates. 
        # convert back to original umap values to calc. euclidean dist. 
        x_real = x/m - a
        y_real = y/m - a
        e_dist = sqrt((x_real-query[row,1])^2 + (y_real-query[row,2])^2)
        query_neighbors[nrow(query_neighbors)+1,] <- c(rownames(query)[row], ref.map[x,y], e_dist)
        k_found = k_found + 1
      }

      # begin searching around the starting point. 
      z = i*2
      for (c in 1:z) {
        if (k_found == k) {
          break
        }
        

        if (is_x == TRUE) {
          x = x + 1*sign
          current_pos = ref.map[x, y]
        }
        else {
          y = y + 1*sign
          current_pos = ref.map[x, y]
        }
        
        if (c == i) {
          # switch between x/y being updated
          is_x = TRUE
        }

        if (length(current_pos) == 0){
          # somehow ended up out of bounds from the ref_map, return.
          print("out of bounds of reference map.")
          print(c(x,y))
          break
        }
        
        if (current_pos != "X") {
          # convert back to original umap values to calc. euclidean dist. 
          x_real = x/m - a
          y_real = y/m - a
          e_dist = sqrt((x_real-query[row,1])^2 + (y_real-query[row,2])^2)
          query_neighbors[nrow(query_neighbors)+1,] <- c(rownames(query)[row], ref.map[x,y], e_dist)
          k_found = k_found + 1
        }
        
      }
      
      sign = sign * -1 #switch sign/direction 
      i = i + 1
      is_x = FALSE
      
    }
    print(paste("Found neighbors for", row, "cell(s)."))  
  }
  return(query_neighbors)
}



query_data <- readRDS("../raw_data/from_garth/cord-cmml-granja-mapped.rds")
ref_data <- readRDS("../raw_data/from_garth/hematopoiesis-5000bp.rds")
query_to_ref_umaps <- query_data$ref.umap@cell.embeddings
ref_umaps <- ref_data$umap@cell.embeddings


# BELOW, THESE ARE JUST FOR TESTING 
# test_ref <- ref_umaps[1:4,]
# test_query <- query_to_ref_umaps[1:4,]

# ref_simple <- data.frame(umap1 = c(1.000, 1.000, 1.002, 1.005, 1.005),
#                          umap2 = c(1.000, 1.001, 1.000, 1.005, 1.003))

# query_simple <- data.frame(umap1 = c(1.001, 1.005, 1.003),
#                            umap2 = c(1.001, 1.005, 1.002))

# row.names(query_simple) <- c("cell a", "cell b", "cell c")
# row.names(ref_simple) <- to_vec(for(i in 1:5) paste("cell",as.character(i)))


map <- generate.map(ref_umaps)
query_neighbors <- find.neighbors(map, query_to_ref_umaps, neighbors)
write.csv(query_neighbors,"./query_neighbors.csv", row.names = FALSE)
