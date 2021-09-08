DimPlot_fun <- function(atlas, group, x, y) {
  
  tmp_df <- data.frame(
    x = atlas@meta.data[, x],
    y = atlas@meta.data[, y],
    color = atlas@meta.data[, group])
  
  dim_plot <- ggplot2::ggplot(tmp_df,
                              ggplot2::aes(x = x, y = y, color = color)) +
    ggplot2::geom_point() +
    ggplot2::theme_classic() +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5))) +
    ggplot2::labs(
      title = "",
      x = "DIM-1",
      y = "DIM-2",
      color = group)
  
  return(dim_plot)
}

FeatPlot_fun(atlas = atlas, gene = "CD3D",  x = "umap_x", y = "umap_y")

FeatPlot_fun <- function(atlas, gene, x, y) {
  
  tmp_df <- data.frame(
    x = atlas@meta.data[, x],
    y = atlas@meta.data[, y],
    gene = atlas@assays$RNA@data[gene, ])
  
  feat_plot <- ggplot2::ggplot(tmp_df,
                              ggplot2::aes(x = x, y = y, color = gene)) +
    ggplot2::geom_point() +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = "",
      x = "DIM-1",
      y = "DIM-2",
      color = "Expression") +
    ggplot2::scale_color_gradient(
      low = "lightgrey",
      high = "blue",
      limits = c(min(tmp_df[, "gene"]), 
                 max(tmp_df[, "gene"])))
  
  return(feat_plot)
}


VlnPlot_fun <- function(atlas, group, gene) {
  
  tmp_df <- data.frame(
    gene = atlas@assays$RNA@data[gene, ],
    color = atlas@meta.data[, group])
  
  vln_plot <- ggplot2::ggplot(tmp_df,
                               ggplot2::aes(x = color, y = gene,
                                            fill = color, color = color)) +
    ggplot2::geom_violin() +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1)) +
    ggplot2::labs(
      title = "",
      x = "",
      y = "Expression",
      fill = group,
      color = group)
  
  return(vln_plot)
}
