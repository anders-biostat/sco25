library( Matrix )
library( sparseMatrixStats )


## Read count matrix

Seurat::ReadMtx( 
  "~/tmp/ifnagrko/ifnagrko_raw_counts.mtx.gz",
  "~/tmp/ifnagrko/ifnagrko_obs.csv.gz",
  "~/tmp/ifnagrko/ifnagrko_var.csv", 
  feature.column=2, skip.cell=1, skip.feature=1, 
  cell.sep=",", feature.sep=",",
  mtx.transpose=TRUE ) -> count_matrix

count_matrix[ 1:5, 1:5 ]


## Calculate fraction

fractions <- t( t(count_matrix) / colSums(count_matrix) )


## Get highly variable genes

hvg <- names( head( sort( rowVars(fractions) / rowMeans(fractions), decreasing=TRUE ), 2000 ) )

plot( rowMeans(fractions), rowVars(fractions) / rowMeans(fractions), 
      cex=.3, log="xy", 
      col = scales::alpha( ifelse( rownames(fractions) %in% hvg, "red", "black" ), .3 ) )
abline( h = mean( 1/colSums(count_matrix)), col="orange")


# Log transform

expr <- log1p( fractions * 1e4 )


# PCA

pca <- irlba::prcomp_irlba( t( expr[hvg,] ), n=20, center=TRUE, scale.=TRUE )
rownames(pca$rotation) <- hvg
rownames(pca$x) <- colnames(expr)


# Distance matrix

dm <- as.matrix( dist( pca$x ) )

t(sapply( 1:ncol(expr), function(cell) 
  order( dm[cell,] )[1:20] )) -> nn
head( nn )

# Distance matrix, faster (using kd-tree algorithm: https://en.wikipedia.org/wiki/K-d_tree)
FNN::get.knn( pca$x, 20 ) -> nn2


# Make UMAP

uwot::umap( pca$x )-> ump
colnames(ump) <- c( "U1", "U2" )

plot( ump, cex=.1, asp=1, col="#00000050" )

sleepwalk::sleepwalk( ump, pca$x )


# Feature plots
library( tidyverse )
as_tibble( ump ) %>%
add_column( expr=expr["Mki67",] ) %>%
ggplot + geom_point( aes( x=U1, y=U2, col=expr ), size=.3 ) + coord_equal()


# Make a k-NN graph

library( igraph )

# Construct edge list: each cell to it's k-th nearest neighbor
k <- 3
#  start of the edge list:
( rbind( 1:nrow(nn), nn[,k+1] ) )[ , 1:20 ]
# the same, flattened (this is the format that igraph wants)
as.vector( rbind( 1:nrow(nn), nn[,k+1] ) )[ 1:40 ]
# now for all values of k
do.call( c, lapply( 1:9, function(k) as.vector( rbind( 1:nrow(nn), nn[,k+1] ) ) ) ) -> edgelist
# make a graph
make_graph( edges=edgelist, n=nrow(nn), directed=FALSE ) -> nn_graph

# Check the vertex degrees
hist( degree(nn_graph) )


# Leiden clustering
cluster_leiden( nn_graph, objective_function="modularity", resolution_parameter=.2 ) -> cm

table( membership(cm) )

as_tibble( ump ) %>%
add_column( cluster = factor( membership(cm) ) ) %>%
ggplot + geom_point( aes( x=U1, y=U2, col=cluster ), size=.3 ) + coord_equal()


# Check modularity score
modularity( nn_graph, membership(cm) )

# Calculate by hand:

# How many of the edges are connecting two vertices in the same group?
get.data.frame( nn_graph, what = "edges" ) %>%
as_tibble() %>%
mutate( 
  from_cluster = membership(cm)[from],
  to_cluster = membership(cm)[to] ) %>%
summarise( mean( from_cluster == to_cluster ) )
# -> 0.923

# And how many would we expect under random permutation of connections?
get.data.frame( nn_graph, what = "edges" ) %>%
as_tibble() %>%
mutate( to = sample(to) ) %>%
mutate(
  from_cluster = membership(cm)[from],
  to_cluster = membership(cm)[to] ) %>%
summarise( mean( from_cluster == to_cluster ) )
# -> 0.077

# The modularity score is the difference between these two:
0.923 - 0.077

# Can we calculated the expectation after permutation without
# actually performing the permutations?

# We can see for each community how many edge stubs it contains
tibble( 
  vertex = 1:vcount(nn_graph),
  degree = degree(nn_graph),
  group = membership(cm) ) %>%
group_by( group ) %>%
summarise( n_stubs = sum(degree) ) -> stub_counts
stub_counts

# Let's write s_g for the number of stubs in cluster g. If we
# randomly pick to stubs to connect, what is the probability
# for them to be in the same group?

sum(stub_counts$n_stubs^2) / sum(stub_counts$n_stubs)^2

# Here, we have the 0.077 from above


# Calculation with matrix

get.data.frame( nn_graph, what = "edges" ) %>%
as_tibble() %>%
mutate(
  from_cluster = membership(cm)[from],
  to_cluster = membership(cm)[to] ) %>%
group_by( from_cluster, to_cluster ) %>%
summarise( n=n() ) %>%
mutate(
  from_cluster = factor( from_cluster, 1:max(from_cluster) ),
  to_cluster = factor( to_cluster, 1:max(to_cluster) ) ) %>%
pivot_wider( id_cols = "from_cluster", names_from = "to_cluster", values_from = "n", values_fill = 0, names_sort=TRUE ) %>%
column_to_rownames("from_cluster") %>% as.matrix -> ccm  # cluster connection matrix
ccm + t(ccm) -> ccm   # both directions


# Fraction of within-group edges
sum(diag(ccm)) / sum(ccm)

# Expected fraction after permuting edge conenctions
colSums(ccm) -> vertex_sums
sum( vertex_sums^2 ) / sum(vertex_sums)^2

# Modularity
sum(diag(ccm)) / sum(ccm) - sum( colSums(ccm)^2 ) / sum( colSums(ccm) )^2


sum(diag(ccm)) / sum(ccm) - sum( ( vertex_sums/sum(vertex_sums) )^2 )

# Compare to igraph's calculation
modularity( nn_graph, membership(cm) )
