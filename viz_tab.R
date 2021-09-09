complete_palette <- readRDS("www/complete_palette.rds")

# create function to plot UMAP

DimPlot_fun <- function(atlas, group, x, y) {
  
  tmp_df <- data.frame(
    x = atlas@meta.data[, x],
    y = atlas@meta.data[, y],
    color = atlas@meta.data[, group])
  
  dim_plot <- ggplot2::ggplot(tmp_df,
                              ggplot2::aes(x = x, y = y, color = color)) +
    ggplot2::geom_point(size = 1) +
    ggplot2::theme_classic() +
    ggplot2::scale_color_manual(values = complete_palette[atlas@meta.data[, group] %>% unique %>% sort]) +
    ggplot2::labs(
      title = "",
      x = "DIM-1",
      y = "DIM-2",
      color = "") +
    ggplot2::theme(text = element_text(size = 20),
                   #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                   legend.text = element_text(size = 21),
                   legend.margin = margin(0, 0, 0, 0)) +
    ggplot2::guides(color = guide_legend(ncol = 1, override.aes = list(size = 4.5)))
  
  return(dim_plot)
}

#FeatPlot_fun(atlas = atlas, gene = "CD3D",  x = "umap_x", y = "umap_y")

FeatPlot_fun <- function(atlas, gene, x, y) {
  
  tmp_df <- data.frame(
    x = atlas@meta.data[, x],
    y = atlas@meta.data[, y],
    gene = atlas@assays$RNA@data[gene, ])
  
  feat_plot <- ggplot2::ggplot(tmp_df,
                               ggplot2::aes(x = x, y = y, color = gene)) +
    ggplot2::geom_point(size = 1) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = "",
      x = "DIM-1",
      y = "DIM-2",
      color = "Expression") +
    ggplot2::scale_color_gradient(
      low = "#D3D3D34D",
      high = "red",
      limits = c(min(tmp_df[, "gene"]), 
                 max(tmp_df[, "gene"]))) +
    # ggplot2::theme(text = element_text(size = 20),
    #                legend.text = element_text(size = 16),
    #                axis.text = element_text(size = 20))
    ggplot2::theme(text = element_text(size = 20),
                   #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                   legend.text = element_text(size = 18))
  
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
    ggplot2::scale_color_manual(values = complete_palette[atlas@meta.data[, group] %>% unique %>% sort]) +
    ggplot2::scale_fill_manual(values = complete_palette[atlas@meta.data[, group] %>% unique %>% sort]) +
    ggplot2::labs(
      title = "",
      x = "",
      y = "Expression",
      fill = group,
      color = group) +
    ggplot2::theme(text = element_text(size = 25),
                   axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                   legend.text = element_text(size = 16))
  
  return(vln_plot)
}