library( tidyverse )
library( Seurat )

readRDS( "~/w/sco24/cmss.rda" ) %>%
CreateSeuratObject() %>%
NormalizeData() %>%
FindVariableFeatures() %>%
ScaleData() %>%
RunPCA( npcs=20 ) %>%
FindNeighbors( dims=1:20 ) %>%
FindClusters( resolution=0.5 ) %>%
RunUMAP( dims=1:20 ) -> seu  

UMAPPlot( seu, label=TRUE )
FeaturePlot( seu, "Aqp4")
FeaturePlot( seu, "Mki67")
FeaturePlot( seu, "Dcx")

sleepwalk::sleepwalk( Embeddings(seu,"umap"), Embeddings(seu,"pca") )

lineage_clusters <- c( 7, 1, 3, 5, 4, 0 )

x <- Embeddings(seu,"pca")[ seu$seurat_clusters %in% lineage_clusters, ]

library( princurve )
pc <- principal_curve( x, smoother="smooth_spline", df=10 )

pc$lambda %>%
enframe("cell","pt") %>%
left_join( 
  Embeddings(seu,"umap") %>% 
  as_tibble(rownames="cell") ) %>%
ggplot +
  geom_point( aes( x=umap_1, y=umap_2, col=pt ) )


cbind(
  Embeddings(seu,"umap"),
  Embeddings(seu,"pca") ) %>%
as_tibble( rownames="cell" ) %>%
ggplot +
  geom_point( aes( x=umap_1, y=umap_2, col=PC_3 ) )


CellCycleScoring( seu, cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes ) -> seu

orth <- AnnotationDbi::select(org.Hs.eg.db, keys=gene_symbols, columns=c("SYMBOL", "MOUSE"), keytype="SYMBOL")

library(biomaRt)
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
getLDS(attributes = c("hgnc_symbol", "mmusculus_homolog_associated_gene_name
"), filters = "hgnc_symbol", values = cc.genes.updated.2019$s.genes, mart = ensembl, attributesL = c("hgnc_symbol", "mgi_symbol"), martL = useMart("ensembl", dataset = "mmusculus_gene_ensembl"))



read_tsv("~/Downloads/mart_export(4).txt") %>%
rename( human=`Gene name`, mouse=`Mouse gene name` ) -> ortholog_table

ortholog_table %>%
filter( human %in% cc.genes$s.genes ) %>% 
filter( !is.na(mouse) ) %>% pull(mouse) -> mouse_s_genes

ortholog_table %>%
filter( human %in% cc.genes$g2m.genes ) %>% 
filter( !is.na(mouse) ) %>% pull(mouse) -> mouse_g2m_genes

Embeddings(seu,"umap") %>%
as_tibble( rownames="cell" ) %>%
add_column( s = colSums( LayerData(seu)[ mouse_s_genes, ] ) ) %>%
ggplot +
  geom_point(aes(x=umap_1,y=umap_2,col=s))

Embeddings(seu,"umap") %>%
as_tibble( rownames="cell" ) %>%
add_column( g2m = colSums( LayerData(seu)[ mouse_g2m_genes, ] ) ) %>%
  ggplot +
  geom_point(aes(x=umap_1,y=umap_2,col=g2m))

ScaleData()

seu %>%
CellCycleScoring( mouse_s_genes, mouse_g2m_genes ) %>%
ScaleData( vars.to.regress = c( "S.Score", "G2M.Score" ) ) %>%
RunPCA( npcs=20 ) %>%
FindNeighbors( dims=1:20 ) %>%
FindClusters( resolution=0.5 ) %>%
RunUMAP( dims=1:20 ) -> seu2
UMAPPlot(seu2,label=TRUE)
FeaturePlot(seu2,"Dcx")


lineage_clusters <- c(6, 1, 4, 8, 2, 0 )

x <- Embeddings(seu2,"pca")[ seu2$seurat_clusters %in% lineage_clusters, ]

library( princurve )
pc <- principal_curve( x, smoother="smooth_spline", df=10 )

pc$lambda %>%
  enframe("cell","pt") %>%
  left_join( 
    Embeddings(seu2,"umap") %>% 
      as_tibble(rownames="cell") ) %>%
  ggplot +
  geom_point( aes( x=umap_1, y=umap_2, col=pt ) ) + scale_color_viridis_c()
