library(tidyverse)
library(Signac)

# run info: takes about 8 minutes for 20 neighbors, uses about 10GB total. 

# PARAMETERS #
a <- 20 # Value to add to each UMAP value so they are all positive. Affects reference matrix size. 10 uses 2.7 GB, 20 uses 6.5 GB. 
m <- 1000 # point multiplier after addition. Use 1000 for 3 significant figures. Example: 3.123 becomes 3123. 
neighbors <- 20 # number of neighbors to find per cell. Will continue searching until all are found. 
write_to_CSV = TRUE # write results to CSV file. 
map_dims = c(12, 12) # based on the largest absolute value observed in UMAP. 


# DATA #
query_data <- readRDS("../../raw_data/from_garth/cord-cmml-granja-mapped.rds")
ref_data <- readRDS("../../raw_data/from_garth/hematopoiesis-monocle3-no-inf.rds")
query_to_ref_umaps <- query_data$ref.umap@cell.embeddings
ref_umaps <- ref_data$umap@cell.embeddings
ref_pt <- ref_data@meta.data$pseudotime 


generate.map <- function(map_dims, ref, ref_pt) {
  # Generates the reference matrix for fast look up of query cells. 
  # Generates a matrix with dimensions max(reference UMAP value) + a * m. 
  # Returns the reference map. 
  
  print("Generating reference map...")
  dim1 <- (map_dims[1] + a) * m
  dim2 <- (map_dims[2] + a) * m
  map <- matrix("X", nrow=dim1, ncol=dim2)
  r <- length(ref_pt)
  
  print("Populating map with reference cells...")
  for (row in 1:r) {
    if (row%%1000 == 0) {
      print(paste("Added", row, "cells to map."))
    }
    x <- ref[row, 1] + a
    y <- ref[row, 2] + a
    x = x * m
    y = y * m

    map[x, y] <- ref_pt[row] #[1] # rownames(ref)[row],# cell name, pseudotime
    
  }
  print("Done.")
  print(paste("Added", row, "cells to map."))
  return(map)
}


find.neighbors <- function(ref.map, query, k, map_dims) {
  # Search for k nearest neighbors. 
  # query UMAP values are transformed in the same manner as the reference values in generate.map().
  # Transformed query UMAP values are mapped to the reference map. A spiral matrix traversal algorithm 
  # searches outward from the initial query location until all k neighbors are found. 
  # This approach is fast if the reference map is dense and can slow down significantly if the
  # reference map is sparse.  
  print("Finding neighbors...")
  r <- nrow(query)

  dim1 <- (map_dims[1] + a) * m
  dim2 <- (map_dims[2] + a) * m

  query_neighbors <- data.frame(query_cell=c("query_cells"),
                                qUMAP_1=c("umap1"),
                                qUMAP_2=c("umap2"),  
                                ref_pseudotime=c("ref_pt"),
                                rUMAP_1=c("umap1"),
                                rUMAP_2=c("umap2"), 
                                euclidean_dist=c("e_dist")
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
    out_of_bounds = FALSE
    
    while (k_found < k) {

      if (out_of_bounds == TRUE) {
        break
      }
      
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
        query_neighbors[nrow(query_neighbors)+1,] <- c(rownames(query)[row], query[row,1], query[row,2], ref.map[x,y], x_real, y_real, e_dist)
        k_found = k_found + 1
      }

      # begin searching around the starting point. 
      z = i*2
      for (c in 1:z) {

        if (x > dim1 | y > dim2) {
          print(paste("Only", k_found, "cells found for query cell", row, "but search out of bounds. Moving to next query cell. "))
          out_of_bounds = TRUE
          break
        }
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
          query_neighbors[nrow(query_neighbors)+1,] <- c(rownames(query)[row], query[row,1], query[row,2], ref.map[x,y], x_real, y_real, e_dist)
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



# FUNCTION CALLS #
map <- generate.map(map_dims, ref_umaps, ref_pt)
query_neighbors <- find.neighbors(map, query_to_ref_umaps, neighbors, map_dims)

if (write_to_CSV == TRUE) {
  write.csv(query_neighbors,"./query_neighbors.csv", row.names = FALSE)
}

