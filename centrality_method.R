# continued from marker_genes_network.R

library(CINNA)
library(igraph)
  
corr_genes1 <- abs(corr_genes)    # take absolute value of correlation
corr_genes1[is.na(corr_genes1)] <- 0   # convert NA values to 0, if any
diag(corr_genes1) <- 0    # make diagonal values 0

corr_graph <- graph_from_adjacency_matrix(corr_genes1, weighted = TRUE, mode = "undirected")

corr_graph <- simplify(corr_graph)

comp_corr_graph <- components(corr_graph)  
mod1 <- modularity(corr_graph, membership(comp_corr_graph))  #0

total_centralities <- proper_centralities(corr_graph)  # 50 measures recommended

centrality_values <- calculate_centralities(corr_graph, include = total_centralities[c(3,4,7,9,15,22,25)]) 
View(centrality_values)

centrality_pca <- pca_centralities(centrality_values, scale.unit = TRUE)  # PCA to find most effective centrality
centrality_pca$data

centrality_tsne <- tsne_centralities(centrality_values, dims = 2, perplexity = 2) # to find the centrality that contributes most
centrality_tsne$data

#visualize_graph(corr_graph, centrality.type = "Lin Centrality")

top_central <- data.frame(sort(centrality_values$`Lin Centrality`, decreasing = TRUE)[1:10])  # taking nodes with top 10 centralities as community centers

top.markers.intersect1 <- data.frame(top.markers.intersect)
  
top.markers.intersect1$cen_values <- centrality_values$`Lin Centrality`

top.markers.intersect1$visited <- 0   # make new column "visited", initialise it to 0
  
# nodes that are community centers are made visited=1
for (gene in rownames(top_central)) {
    gene_row = rownames(top.markers.intersect1)[which(top.markers.intersect1[,1]==gene)]
    top.markers.intersect1[gene_row,"visited"] = 1
  }
  
# make new empty graph  
new_graph <- make_empty_graph(n=dim(corr_genes)[1], directed = FALSE)  
new_graph <- set_vertex_attr(new_graph,"name",index = V(new_graph),value = top.markers.intersect)
comp_new_graph <- components(new_graph)

# visit all other nodes that are not community centers and assign them to one of the community centers
for (other_nodes in 1:length(top.markers.intersect)) {

  get_nodes <- data.frame(matrix(ncol = 3)) # empty matrix to save nodes

  if(top.markers.intersect1[other_nodes,"visited"]==0){
    top.markers.intersect1[other_nodes,"visited"] = 1 

    for (k in 1:dim(top_central)[1]) {
      # get correlation between two vertices
      get_nodes[k,] = c(top.markers.intersect1[other_nodes,1], rownames(top_central)[k], corr_graph[top.markers.intersect1[other_nodes,1],rownames(top_central)[k]])

    }

    max_value <- max(get_nodes[,3])
    if(max_value!=0){
      max_node <- get_nodes[which(get_nodes[,3]==max(get_nodes[,3])),]  # get pair of nodes with max correlation

      node1 = unlist(max_node[1])
      node2 = unlist(max_node[2])

      new_graph[node1,node2] = corr_graph[node1,node2]  # SUB-GRAPH

    }


  }
}


# GET MODULARITY OF SUB-GRAPH
comp_new_graph <- components(new_graph) 
mod1_new_graph <- modularity(new_graph, membership(comp_new_graph))  # save modularity of the sub-graph

# MERGE SUB-COMMUNITIES
 
  diff_mod = mod1_new_graph - mod1  # change in modularity
  
  org_graph = new_graph  # save the graph with initial sub-communities
  comp_org_graph = comp_new_graph
  merge_community = data.frame(matrix(ncol = 3)) # to save which communities to be merged
  
  
  for (m in 1:((comp_new_graph$no)-1)) {
    print(paste0("community1 -- ",m))
    merging_nodes = data.frame(matrix(ncol = 3))  # matrix to take maximum merging factor
    
    inter_sub_correlation = data.frame(matrix(ncol = 3)) # matrix to save inter cluster correlation
    
    com_nodes1 <- communities(comp_new_graph)[m] #get nodes of a community
    sub_graph1 <- induced_subgraph(corr_graph, unlist(com_nodes1)) # get the first individual cluster 
    
    edge_sum1 = sum(E(sub_graph1)$weight)  #sum of correlation edges of own cluster
    avg_sum1 = edge_sum1/length(com_nodes1[[1]])  # average of correlation within community
    merge_cen1 = intersect(unlist(com_nodes1),top_central) # get the center node 1
    
    for (n in (m+1):(comp_new_graph$no)) {
      print(paste0("community2 -- ",n))
      sum_edge_cut = 0
      com_nodes2 <- communities(comp_new_graph)[n]  #get nodes of second community
      sub_graph2 <- induced_subgraph(corr_graph, unlist(com_nodes2)) # get the second individual cluster 
      
      edge_sum2 = sum(E(sub_graph2)$weight)  #sum of edges of own cluster - second
      avg_sum1 = edge_sum1/length(com_nodes2[[1]])  # average of correlation within community - second
      #merge_cen2 = intersect(unlist(com_nodes2),cluster_centers) # get the center node 2
      
      # sum of correlation edges shared by the two clusters
      for (p in 1: length(com_nodes1[[1]])) {
        for (q in 1:length(com_nodes2[[1]])) {
          sum_edge_cut = sum_edge_cut + corr_graph[com_nodes1[[1]][p],com_nodes2[[1]][q]]
          new_entry <- c(com_nodes1[[1]][p],com_nodes2[[1]][q],corr_graph[com_nodes1[[1]][p],com_nodes2[[1]][q]]) 
          inter_sub_correlation = rbind(inter_sub_correlation, new_entry) 
        }
      }
      
      # MERGING FACTOR *************************
      merging_factor = (sum_edge_cut/(gsize(sub_graph1)+gsize(sub_graph2)))/((edge_sum1/gsize(sub_graph1)) + (edge_sum2/gsize(sub_graph2)))
      merging_nodes = rbind(merging_nodes, c(m,n,merging_factor)) # cluster1, cluster2, value
      
    }
    
    max_merge = max(merging_nodes[,3],na.rm = TRUE)  # maximum merging factor among two communities
    merge_com1 = m
    merge_community = rbind(merge_community,merging_nodes[which(merging_nodes[,3]==max_merge),])  # get the two communities
    
  }
  
  sort_merge_community <- merge_community[order(merge_community[,3], decreasing = TRUE),]
  sort_merge_community <- cbind(sort_merge_community, visited=0)
  sort_merge_community <- na.omit(sort_merge_community)
  sort_merge_community <- subset(sort_merge_community, sort_merge_community[,3]!=Inf)
  
  
  org_graph = new_graph
  comp_org_graph = comp_new_graph
  
  
# Merge pairs of communities and check modularity
  for (r in 1:dim(sort_merge_community)[1]) { 
    temp_graph = new_graph
    comp_temp_graph = comp_new_graph
    diff_mod1 = diff_mod
    mod_temp = mod1_new_graph
    
    com_nodes1 <- communities(comp_org_graph)[sort_merge_community[r,1]] #get nodes of a community from original
    sub_graph1 <- induced_subgraph(corr_graph, unlist(com_nodes1)) # get the individual cluster
    
    com_nodes2 <- communities(comp_org_graph)[sort_merge_community[r,2]] #get nodes of a community from original
    sub_graph2 <- induced_subgraph(corr_graph, unlist(com_nodes2)) # get the individual cluster 
    
    for (e1 in 1:length(unlist(com_nodes1))) {
      merge_node1 = unlist(com_nodes1)[e1]
      for (e2 in 1:length(unlist(com_nodes2))) {
        merge_node2 = unlist(com_nodes2)[e2]
        
        new_graph[merge_node1, merge_node2] = corr_graph[merge_node1, merge_node2]  # CREATING FINAL GRAPH
      }
      
      sort_merge_community[r,"visited"]=1
    }
    
    comp_new_graph <- components(new_graph)  # membership, csize(cluster sizes), no(no. of clusters)
    mod1_new_graph <- modularity(new_graph, membership(comp_new_graph))
    diff_mod = mod1_new_graph - mod_temp  # check difference in  modularity of current and previous graph
    print(paste0("diff_mod---  ",diff_mod,"----mod1_new_graph----",mod1_new_graph))
    
    if(diff_mod < diff_mod1){
      print(paste0(diff_mod," < ",diff_mod1))
      #new_graph[merge_node1, merge_node2] = 0 # remove the edge
      new_graph = temp_graph
      comp_new_graph = comp_temp_graph
      sort_merge_community[r,"visited"]=0
      
      #comp_new_graph <- components(new_graph)  # membership, csize(cluster sizes), no(no. of clusters)
      mod1_new_graph <- modularity(new_graph, membership(comp_new_graph))
    
      print(paste0("diff_mod---  ",diff_mod,"----mod1_new_graph----",mod1_new_graph))
    }
  }
  
  
  comp_new_graph$csize
  comp_new_graph$no
  mod1_new_graph
