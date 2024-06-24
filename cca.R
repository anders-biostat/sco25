library( tidyverse )
library( Seurat )

ifnb <- SeuratData::LoadData("ifnb")

ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)
ifnb
# run standard anlaysis workflow
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb)
ifnb <- ScaleData(ifnb)
ifnb <- RunPCA(ifnb)
ifnb <- FindNeighbors(ifnb, dims = 1:30, reduction = "pca")
ifnb <- FindClusters(ifnb, resolution = 2, cluster.name = "unintegrated_clusters")

ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(ifnb, reduction = "umap.unintegrated", group.by = c("stim", "seurat_clusters"))

ifnb <- IntegrateLayers(object = ifnb, method = CCAIntegration, 
            orig.reduction = "pca", new.reduction = "integrated.cca" )

# re-join layers after integration
ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])

ifnb <- FindNeighbors(ifnb, reduction = "integrated.cca", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 1)

ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(ifnb, reduction = "umap.unintegrated", group.by = c("stim", "seurat_clusters"))

ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.cca")

# Visualization
DimPlot(ifnb, reduction = "umap", group.by = c("stim", "seurat_annotations"))

save( ifnb, file="~/tmp/ifnb_integrated.rda" )


Embeddings(ifnb,"umap") %>%
as_tibble() %>%
add_column( cond = ifnb$stim ) %>%
sample_frac() %>%
ggplot() +
  geom_point(aes(x=umap_1,y=umap_2,col=cond),size=.01) +
  coord_equal()









### Manual

# PCA performend on all cells together, regardless of origin
Embeddings( ifnb, "pca" ) -> pca

# Origin of the cells:
table( ifnb$stim )

# NN search:
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

# How many of the neighbours are mutual?
  
knn_tbl %>%
mutate( stim_cell = ifelse( center_is=="STIM", center_cell, ngb_cell ) ) %>%
mutate( ctrl_cell = ifelse( center_is=="CTRL", center_cell, ngb_cell ) ) %>%
group_by( stim_cell, ctrl_cell ) %>%
summarise( n=n() ) %>% ungroup() %>%
count( n ) 

knn_tbl %>%
mutate( stim_cell = ifelse( center_is=="STIM", center_cell, ngb_cell ) ) %>%
mutate( ctrl_cell = ifelse( center_is=="CTRL", center_cell, ngb_cell ) ) %>%
group_by( stim_cell, ctrl_cell ) %>%
summarise( n=n(), .groups="drop" ) %>%
filter( n==2 ) %>% select(-n) -> mnn_tbl

Embeddings(ifnb,"umap.unintegrated") %>%
as_tibble( rownames="cell" ) %>%
add_column( cond = ifnb$stim ) %>%
mutate( is_among_mnn = cell %in% c( mnn_tbl$ctrl_cell, mnn_tbl$stim_cell ) ) %>%
sample_frac() %>%
ggplot( aes( x=umapunintegrated_1, y=umapunintegrated_2, col=cond ) ) +
  geom_point( size=.1, alpha=.2 ) +
  geom_point( size=.1, alpha=1, data=function(x) filter(x,is_among_mnn) ) +
  coord_equal()



mnn_tbl %>%
group_by( ctrl_cell ) %>%
summarise(n())

sleepwalk::sleepwalk(
  list( 
    Embeddings(ifnb,"umap.unintegrated")[ ifnb$stim=="CTRL", ],
    Embeddings(ifnb,"umap.unintegrated")[ ifnb$stim=="STIM", ] ),
  list( 
    Embeddings(ifnb,"pca")[ ifnb$stim=="CTRL", ],
    Embeddings(ifnb,"pca")[ ifnb$stim=="STIM", ] ),
  same="features" )



#### SVD

crosscor <- as.matrix(
   t(LayerData(ifnb)[,ifnb$stim=="CTRL"][VariableFeatures(ifnb),]) %*% 
      LayerData(ifnb)[,ifnb$stim=="STIM"][VariableFeatures(ifnb),] ) 

svd <- irlba::svdr( crosscor, 20 )

sleepwalk::sleepwalk(
  list( 
    Embeddings(ifnb,"umap.unintegrated")[ ifnb$stim=="CTRL", ],
    Embeddings(ifnb,"umap.unintegrated")[ ifnb$stim=="STIM", ] ),
  list(
    t( t(svd$u) * sqrt(svd$d) ),
    t( t(svd$v) * sqrt(svd$d) ) ),
  same="features" )

cca <- rbind( 
  t( t(svd$u) * sqrt(svd$d) ),
  t( t(svd$v) * sqrt(svd$d) ) )
rownames(cca) <- c( colnames(ifnb)[ifnb$stim=="CTRL"], colnames(ifnb)[ifnb$stim=="STIM"] )
ump_cca <- uwot::umap(cca)

all( rownames(ump_cca) == colnames(ifnb) )

perm <- sample.int(nrow(ump_cca))
plot( ump_cca[perm,], cex=.2, asp=1, col=1+as.integer(factor(ifnb$stim[perm])) )

sleepwalk::sleepwalk(
  list( Embeddings(ifnb,"umap"), ump_cca ),
  list( Embeddings(ifnb,"integrated.cca"), cca)
)


qqplot( svd$u[,1], svd$v[,1] )

u1_sorted <- sort( svd$u[,1] )
v1_sorted <- sort( svd$v[,1] )
plot(
  u1_sorted[ ceiling( seq(0,1,length.out=300)*length(u1_sorted) ) ],
  v1_sorted[ ceiling( seq(0,1,length.out=300)*length(v1_sorted) ) ] )

plot( 
    cor( t(as.matrix(LayerData(ifnb)[,ifnb$stim=="CTRL"][VariableFeatures(ifnb),])), svd$u[,1] ),
    cor( t(as.matrix(LayerData(ifnb)[,ifnb$stim=="STIM"][VariableFeatures(ifnb),])), svd$v[,1] ) )
