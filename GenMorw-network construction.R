## Construction of the GenMorw-network and identification of its strongly connected components and cliques
# To construct the GenMorw-network, it is necessary to set the "retain_mmPPI" parameter to "True" in main.R 
# program to ensure that GenMorw retains gene interaction information within the heterogeneous network. 
# Below is the parameter statement for the program:
# Cancer_name: TCGA Abbreviation of the cancer for which the GenMorw-network is to be constructed
# min_num_of_component: Minimum number of genes contained in the smallest strongly connected component
# min_size_clique: Minimum number of genes contained in the smallest clique.
# thre: A confidence threshold (from 0 to 1), where 0.5 means that interactions must be supported by at least 50% of the heterogeneous networks.


Cancer_name <- "LAML"
min_num_of_component <- 5
min_size_clique <- 3
thre <- 0.5

library(igraph)
all_mut_genes <- scan(paste0("Data/", Cancer_name, "/all_mut_genes"), what = "c", quiet = T)
global_matrix <- matrix(rep(0, length(all_mut_genes) * length(all_mut_genes)), ncol = length(all_mut_genes))
colnames(global_matrix) <- all_mut_genes
rownames(global_matrix) <- all_mut_genes
mpath <- paste0("Models/Data/", Cancer_name, "/mmPPI")
all_group_files <- list.files(mpath)
i <- 0
for (file in all_group_files)
{
  i <- i + 1
  print(paste0(Cancer_name, " ", i, " / ", length(all_group_files)))
  demo <- readRDS(paste0(mpath, "/", file))
  diag(demo) <- 0
  demo <- apply(demo, 2, function(x) {
    thre <- quantile(x, 0.95)
    ifelse(x > thre, 1, 0)
  })

  global_matrix[rownames(demo), colnames(demo)] <- demo + global_matrix[rownames(demo), colnames(demo)]
}
global_matrix <- global_matrix[, colSums(global_matrix) != 0]
global_matrix <- global_matrix[rowSums(global_matrix) != 0, ]
saveRDS(global_matrix, paste0("Results/Data/", Cancer_name, "/", Cancer_name, "_GenMorw-network.rds"))

global_matrix <- log2(1+global_matrix)
global_matrix[global_matrix < log2(1+length(all_group_files)*thre)] <- 0
global_matrix <- global_matrix[, colSums(global_matrix) != 0]
global_matrix <- global_matrix[rowSums(global_matrix) != 0, ]

all_genes <- c(colnames(global_matrix), rownames(global_matrix)) %>% unique()
gm <- matrix(rep(0, length(all_genes) * length(all_genes)), ncol = length(all_genes))
rownames(gm) <- all_genes
colnames(gm) <- all_genes
gm[rownames(global_matrix), colnames(global_matrix)] <- global_matrix
global_matrix <- gm
rm(gm)
diag(global_matrix) <- 0

my_graph_directed <- graph_from_adjacency_matrix(t(global_matrix), mode = "directed", weighted = TRUE, diag = FALSE)
global_matrix2 <- (global_matrix + t(global_matrix)) / 2
my_graph_undirected <- graph_from_adjacency_matrix(global_matrix2, mode = "undirected", weighted = TRUE, diag = FALSE)

# 寻找强连通分量
components_strong <- components(my_graph_directed, mode = "strong")
# 找到节点数大于10的所有连通分量
large_components <- which(components_strong$csize > min_num_of_component)
# 在无向图中找到团结构
my_cliques <- max_cliques(my_graph_undirected, min = min_size_clique)

components_record <- list()
i <- 0
for(com in large_components)
{
  i <- i + 1
  subgraphs <- induced_subgraph(my_graph_directed, which(components_strong$membership == com))
  tc <- ends(subgraphs, E(subgraphs))
  colnames(tc) <- c("from","to")
  components_record[[paste0("Component ",i)]] <- tc
}

i <- 0
cliques_record <- list()
for(clique in my_cliques)
{
  i <- i + 1
  clique_graph <- induced_subgraph(my_graph_undirected, clique)
  cliques_record[[paste0("Clique ",i)]] <- V(clique_graph)$name
}

saveRDS(components_record, paste0("Results/Data/", Cancer_name, "/", Cancer_name, "_GenMorw-network_strongly_connected_components.rds"))
saveRDS(cliques_record, paste0("Results/Data/", Cancer_name, "/", Cancer_name, "_GenMorw-network_cliques.rds"))
