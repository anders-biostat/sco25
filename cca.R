library( tidyverse )
library( Seurat )

# Load `ifnb` data set
ifnb <- SeuratData::LoadData("ifnb")
ifnb

# Split by `stim`:
ifnbs <- ifnb
ifnbs[["RNA"]] <- split( ifnbs[["RNA"]], f = ifnbs$stim )
ifnbs

# Run standard analysis workflow
ifnbs %>%
NormalizeData() %>%
FindVariableFeatures() %>%
ScaleData() %>%
RunPCA() %>%
FindNeighbors( dims=1:30, reduction="pca" ) %>%
FindClusters( resolution=2, cluster.name="unintegrated_clusters" ) %>%
RunUMAP( dims=1:30, reduction="pca", reduction.name="umap.unintegrated" ) -> ifnbs

# Plot UMAP, colour by `stim` and by Leiden clusters
DimPlot( ifnbs, reduction = "umap.unintegrated", group.by = "stim" ) + coord_equal()
DimPlot( ifnbs, reduction = "umap.unintegrated", group.by = "seurat_clusters" ) + coord_equal()

# Also show UMAP using colours for annotation from previous analysis
DimPlot( ifnbs, reduction = "umap.unintegrated", group.by = "seurat_annotations", label=TRUE ) + coord_equal()

# Now, we use Sleepwalk to compare the two cell subsets:
sleepwalk::sleepwalk(
  list( 
    Embeddings(ifnbs,"umap.unintegrated")[ ifnbs$stim=="CTRL", ],
    Embeddings(ifnbs,"umap.unintegrated")[ ifnbs$stim=="STIM", ] ),
  list( 
    Embeddings(ifnbs,"pca")[ ifnbs$stim=="CTRL", ],
    Embeddings(ifnbs,"pca")[ ifnbs$stim=="STIM", ] ),
  same="features" )


# Perform CCA integration (This takes a few minutes)
ifnbs <- IntegrateLayers( object = ifnbs, method = CCAIntegration, 
            orig.reduction = "pca", new.reduction = "integrated.cca" )

# Save result of this lengthy calculation
#save( ifnbs, file="~/tmp/ifnbs_integrated.rda" )
load( "~/tmp/ifnbs_integrated.rda" )

# Re-join layers after integration
ifnb <- ifnbs
ifnb[["RNA"]] <- JoinLayers( ifnb[["RNA"]] )

# Redo neighbor search, clustering, and UMAP, now based on CCA instead of PCA:
ifnb %>% 
FindNeighbors( reduction = "integrated.cca", dims = 1:30 ) %>%
FindClusters( resolution = 1 ) %>%
RunUMAP( dims = 1:30, reduction = "integrated.cca" ) -> ifnb

# Plot UMAP again, coloured as before by `stim` covariate, clusters (from new clusterin),
# and cell types (from original analysis):
DimPlot( ifnb, reduction = "umap", group.by = "stim" ) + coord_equal()
DimPlot( ifnb, reduction = "umap", group.by = "seurat_clusters" ) + coord_equal()
DimPlot( ifnb, reduction = "umap", group.by = "seurat_annotations", label=TRUE ) + coord_equal()

# The plot coloured by `stim` is misleading, as the blue points are plotted after 
# (and hence on top of) the red points.
# Let's redo this plot manually and use `sample_frac` to permute the table rows and
# thus shuffle the order in which the points are placed
Embeddings(ifnb,"umap") %>%
as_tibble() %>%
add_column( cond = ifnb$stim ) %>%
sample_frac() %>%
ggplot() +
  geom_point(aes(x=umap_1,y=umap_2,col=cond),size=.01) +
  coord_equal()

## Now we perform a simplified CCA manually, to compare

# First, here's the number of cells of each condition:
table( ifnb$stim )

# Let's calculate a cross correlation between the control and the
# the stimulated cells. The following matrix has one row for each
# control cell and one column for each stimulated cells. The matrix entries
# are the dot products of the corresponding cells' expression vectors,
# as formed from the log-normalized expressions for the 2000 most variable genes.

crosscor <- as.matrix(
  t(LayerData(ifnb)[,ifnb$stim=="CTRL"][VariableFeatures(ifnb),]) %*% 
    LayerData(ifnb)[,ifnb$stim=="STIM"][VariableFeatures(ifnb),] ) 

# Strictly speaking, this is not a cross-correlation yet, as we should have
# standardized the two expression matrices before multiplying them by subtracting
# for each gene its mean and dividing by its standard deviation. I'll add this when
# I go through this code once more. 

# Now, do an SVD on the cross-correlation matrix
# svd <- irlba::svdr( crosscor, 20 )

#save( svd, file="~/tmp/svd.rda" )
load( "~/tmp/svd.rda" )

# We use the left and right singular vectors as new coordinates, to replace the
# PCA coordinates. The diagonal matrix in the middle is split and each square root
# multiplied to one side. This gives us better distances, as can be seen with Sleepwalk:

sleepwalk::sleepwalk(
  list( 
    Embeddings(ifnb,"umap.unintegrated")[ ifnb$stim=="CTRL", ],
    Embeddings(ifnb,"umap.unintegrated")[ ifnb$stim=="STIM", ] ),
  list(
    t( t(svd$u) * sqrt(svd$d) ),
    t( t(svd$v) * sqrt(svd$d) ) ),
  same="features" )


# Therefore, we use these as new reduction coordinates

cca <- rbind( 
  t( t(svd$u) * sqrt(svd$d) ),
  t( t(svd$v) * sqrt(svd$d) ) )
rownames(cca) <- c( colnames(ifnb)[ifnb$stim=="CTRL"], colnames(ifnb)[ifnb$stim=="STIM"] )
all( rownames(cca) == colnames(ifnb) )

# Add these to the Seurat object
ifnb@reductions$svd <- CreateDimReducObject( cca, assay="RNA", key="svd_" )

# Let's make a UMAP from these:
RunUMAP( ifnb, dims=1:20, reduction="svd", reduction.name = "umap.svd" ) -> ifnb

# Here's a plot for the new UMAP:
perm <- sample.int(ncol(ifnb))
plot( Embeddings(ifnb,"umap.svd")[perm,], cex=.2, asp=1, col=1+as.integer(factor(ifnb$stim[perm])) )

# And the usual colourings
DimPlot( ifnb, reduction = "umap.svd", group.by = "seurat_clusters" ) + coord_equal()
DimPlot( ifnb, reduction = "umap.svd", group.by = "seurat_annotations", label=TRUE ) + coord_equal()

# Of course, the clustering does not fit that well, as it was done using another embedding.

# And here a comparison between Seurat's CCA integration and our simplified one:
sleepwalk::sleepwalk(
  list( Embeddings(ifnb,"umap"), Embeddings(ifnb,"umap.svd") ),
  list( Embeddings(ifnb,"integrated.cca"), Embeddings(ifnb,"svd"))
)


plot( density( Embeddings(ifnb,"integrated.cca")[ifnb$stim=="CTRL",1] ) )
plot( density( Embeddings(ifnb,"integrated.cca")[ifnb$stim=="STIM",1] ) )


### Integration via mutual nearest neighbours

# (FROM HERE ON: only draft)

# We start with the PCA that was computed at the very beginning, using all cells
# and ignoring the stim covariate
Embeddings( ifnb, "pca" ) -> pca

# Here ios the number of cells 
table( ifnb$stim )

# cross-NN search:
nn_cs <- FNN::get.knnx( pca[ifnb$stim=="CTRL",], pca[ifnb$stim=="STIM",], k=25 )
nn_sc <- FNN::get.knnx( pca[ifnb$stim=="STIM",], pca[ifnb$stim=="CTRL",], k=25 )

bind_rows(
  nn_cs$nn.index %>%
    `colnames<-`( 1:ncol(.) ) %>%
    as_tibble() %>%
    add_column( center_cell = names(which(ifnb$stim=="STIM")) ) %>%
    pivot_longer( -`center_cell`, names_to="ngb_index", values_to="ngb_cell" ) %>%
    mutate( ngb_index = as.integer(ngb_index) ) %>%
    mutate( ngb_cell = names(which(ifnb$stim=="CTRL"))[ ngb_cell ] ) %>%
    add_column( center_is = "STIM" ),
  nn_sc$nn.index %>%
    `colnames<-`( 1:ncol(.) ) %>%
    as_tibble() %>%
    add_column( center_cell = names(which(ifnb$stim=="CTRL")) ) %>%
    pivot_longer( -`center_cell`, names_to="ngb_index", values_to="ngb_cell" ) %>%
    mutate( ngb_index = as.integer(ngb_index) ) %>%
    mutate( ngb_cell = names(which(ifnb$stim=="STIM"))[ ngb_cell ] ) %>%
    add_column( center_is = "CTRL" ) ) -> knn_tbl

head( knn_tbl )

# How many of the neighbours are mutual?
knn_tbl %>%
mutate( stim_cell = ifelse( center_is=="STIM", center_cell, ngb_cell ) ) %>%
mutate( ctrl_cell = ifelse( center_is=="CTRL", center_cell, ngb_cell ) ) %>%
group_by( stim_cell, ctrl_cell ) %>%
summarise( n=n() ) %>% ungroup() %>%
count( n ) 

# Make table of MNNs
knn_tbl %>%
mutate( stim_cell = ifelse( center_is=="STIM", center_cell, ngb_cell ) ) %>%
mutate( ctrl_cell = ifelse( center_is=="CTRL", center_cell, ngb_cell ) ) %>%
group_by( stim_cell, ctrl_cell ) %>%
summarise( n=n(), .groups="drop" ) %>%
filter( n==2 ) %>% select(-n) -> mnn_tbl
head( mnn_tbl )

# How many MNNs?
mnn_tbl %>%
group_by( ctrl_cell ) %>%
summarise(n())



## Save image
#save.image( file="~/tmp/cca.rda" )
