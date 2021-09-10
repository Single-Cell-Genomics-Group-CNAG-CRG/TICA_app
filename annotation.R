### Define functions for the annotation tab

# function to create jaccard index plot from matchSCore

plot_heatmap <- function(atlas, atlas_markers, marker_genes, annotation, organism){
  
  # load selected annotation level markers
  
  atlas_markers <- filter(atlas_markers, level == annotation)
  atlas_markers <- split(atlas_markers$gene, atlas_markers$cluster)
  
  # read query marker genes
  
  marker_genes <- read.delim(marker_genes,
                             header = TRUE,
                             sep = ",")
  
  if(organism == "mmus"){
    marker_genes$gene <- stringr::str_to_upper(marker_genes$gene)
  }
  
  marker_genes <- split(marker_genes$gene, marker_genes$cluster)
  
  ms <- matchSCore2::matchSCore2(gene_cl.ref = atlas_markers[sort(names(atlas_markers))], 
                                 gene_cl.obs = marker_genes[sort(names(marker_genes))], 
                                 ylab = "Atlas cell types", 
                                 xlab = "Query clusters")
  return(ms$ggplot)
  
}