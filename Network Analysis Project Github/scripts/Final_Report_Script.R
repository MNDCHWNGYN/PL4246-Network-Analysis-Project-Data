setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Amanda/School/NUS/1. Psychology/PL4246/Network Analysis Project")

library(igraph)

# Load the brain connectivity matrix and labels ----
matrix <- read.delim('AveragedMatrix.txt', header = FALSE, sep = ' ')
# if need to use short labels: short_labels <- read.delim('region_names_abbrev_file.txt', header = FALSE)
long_labels <- read.delim('region_names_full_file.txt', header = FALSE)

# Pre-process the matrix for network analysis ----
matrix <- abs(matrix) 
matrix <- as.matrix(matrix)
colnames(matrix) <- long_labels$V1

# Setting thresholds to range from 10% to 50%, with an interval of 5%
thresholds <- seq(0.1, 0.5, by = 0.05) 
all_vals <- c(as.matrix(matrix))

results <- data.frame()
hub_list <- list()

# Threshold analysis ----
for (t in thresholds) {
  
  th_value <- quantile(all_vals, probs = 1 - t)
  
  mat_temp <- matrix
  mat_temp[mat_temp < th_value] <- 0
  
  net_temp <- graph_from_adjacency_matrix(mat_temp, 
                                          weighted = TRUE, 
                                          mode = "undirected")
  
  ## Real Network Metrics ----
  ### Global Clustering Coefficient ----
  C_real <- transitivity(net_temp, type = "global")
  ### Path Length ----
  L_real <- mean_distance(net_temp, 
                          weights = 1/E(net_temp)$weight, 
                          unconnected = TRUE)
  
  ## Random Network for SWI ----
  er_graphs <- list()
  for(i in 1:1000){  
    er_graphs[[i]] <- sample_gnm(
      n = gorder(net_temp), 
      m = gsize(net_temp), 
      directed = FALSE
    )
  }
  
  C_rand <- mean(sapply(er_graphs, transitivity, type = "global"))
  L_rand <- mean(sapply(er_graphs, mean_distance, unconnected = TRUE))
  
  ### SWI ----
  SWI <- (C_real / C_rand) / (L_real / L_rand)
  
  ## Modularity__Community Structure ---- 
  ### Louvain Community Detection ----
  comm_louvain <- cluster_louvain(net_temp, 
                                  weights = E(net_temp)$weight)
  modularity_louvain <- modularity(comm_louvain)
  
  ## Hub Structure__Node Centrality 
  ### Betweenness Centrality ----
  btw_centrality <- betweenness(net_temp, 
                                weights = 1/E(net_temp)$weight, 
                                normalized = TRUE)
  hub_list[[as.character(t)]] <- btw_centrality
  ### Hub Summary Metrics ----
  hub_index <- max(btw_centrality) / mean(btw_centrality)
  top_nodes <- names(sort(btw_centrality, decreasing = TRUE)[1:5])
  
  ## Assortativity__Connectivity Patterns
  ### Assortative mixing by degree ----
  assort_deg <- assortativity_degree(net_temp, directed = F)
  
  ### Assortative mixing by node attributes ----
  #### by hemisphere ----
  V(net_temp)$Hemisphere <- ifelse(
    grepl("^Left", V(net_temp)$name), "L", 
    ifelse(grepl("^Right", V(net_temp)$name), "R", "Other"))
  
  V(net_temp)$hemisphere <- as.factor(V(net_temp)$Hemisphere)
  
  assort_hemi <- assortativity_nominal(net_temp, 
                                       types = factor(V(net_temp)$hemisphere), 
                                       directed = F)
  
  #### by different lobes ----
  V(net_temp)$lobe <- ifelse(
    grepl("Frontal", V(net_temp)$name), "Frontal", 
    ifelse(grepl("Parietal", V(net_temp)$name), "Parietal", 
           ifelse(grepl("Temporal", V(net_temp)$name), "Temporal", 
                  ifelse(grepl("Occipital", V(net_temp)$name), "Occipital", "Other"))))
  
  V(net_temp)$Lobe <- as.factor(V(net_temp)$lobe)
  
  assort_lobe <- assortativity_nominal(net_temp, 
                                       types = factor(V(net_temp)$Lobe), 
                                       directed = F)
  
  ## Global Structure Controls ---- 
  density <- edge_density(net_temp)
  degree <- mean(degree(net_temp))
  
  components <- components(net_temp)$no
  
  diameter <- diameter(net_temp,
                       directed = FALSE,
                       weights = 1/E(net_temp)$weight,
                       unconnected = TRUE)
  
  ## Store Results ----
  results <- rbind(results, data.frame(
    threshold = t,
    
    # Small-world
    SWI = SWI,
    clustering = C_real,
    path_length = L_real,
    
    # Modularity
    modularity_louvain = modularity_louvain,
    
    # Hubs
    hub_index = hub_index,
    top_hubs = paste(top_nodes, collapse = ", "),
    
    # Connectivity
    assort_degree = assort_deg,
    assort_hemisphere = assort_hemi,
    assort_lobe = assort_lobe,
    
    # Global structure
    density = density,
    mean_degree = degree,
    components = components,
    diameter = diameter
  ))
}

results

# Hub stability ----
hub_table <- table(unlist(strsplit(results$top_hubs, ", ")))
sort(hub_table, decreasing = TRUE)

# Identifying Detached Node when Threshold = 10% ----
## Setting threshold to 10% ----
matrix[matrix < 0.255] <- 0 

## Clean up ---- 
matrix <- as.matrix(matrix)
colnames(matrix) <- long_labels$V1 # add node labels 

## Generating the network ----
network <- graph_from_adjacency_matrix(matrix, weighted = T, mode = "undirected")
summary(network)

## Creating table to show brain areas and their corresponding component in the network ----
comp_df <- data.frame(
  node = names(components(network)$membership), 
  component = components(network)$membership
  )
comp_df

## Specify the name of the detached node ---
detached_node = paste(names(
                      components(network)$membership[components(network)$membership != which.max(components(network)$csize)]
                      ))
detached_node # The node is the Right Temporal Pole

# Plots ----
## Plots for Small-World across thresholds ----
plot(results$threshold, results$SWI, type = "b",
     xlab = "Threshold",
     ylab = "SWI",
     main = "SWI Across Thresholds")
plot(results$threshold, results$clustering, type = "b",
     xlab = "Threshold",
     ylab = "Clustering",
     main = "Global Clustering Across Thresholds")
plot(results$threshold, results$path_length, type = "b",
     xlab = "Threshold",
     ylab = "Path Length",
     main = "Path Length Across Thresholds")

## Plots for Modularity across thresholds ---- 
plot(results$threshold, results$modularity_louvain, type = "b", 
     xlab = "Threshold",
     ylab = "Modularity",
     main = "Modularity Across Thresholds")

## Plot for hub_index across thresholds ----
plot(results$threshold, results$hub_index, type = "b",
     xlab = "Threshold",
     ylab = "Hub Index",
     main = "Hub Dominance Across Thresholds")

## Plot for Connectivity across thresholds ----
plot(results$threshold, results$assort_hemisphere, type = "b",
     xlab = "Threshold",
     ylab = "Hemisphere Assortativity",
     main = "Hemisphere Assortativity Across Thresholds")

plot(results$threshold, results$assort_lobe, type = "b",
     xlab = "Threshold",
     ylab = "Lobe Assortativity",
     main = "Lobe Assortativity Across Thresholds")

plot(results$threshold, results$assort_degree, type = "b",
     xlab = "Threshold",
     ylab = "Degree Assortativity",
     main = "Degree Assortativity Across Thresholds")


write.csv(results, "Final Results.csv", row.names = TRUE)
save.image(file = "Final_Report.RData")
