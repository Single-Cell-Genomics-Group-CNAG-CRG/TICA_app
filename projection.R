### Define functions for the projection tab

# load_palette
complete_palette <- readRDS("www/complete_palette.rds")

# function to project queries on the atlas to predict cell-wise annotation

projection <- function(atlas, query, annotation, ndims) {
  
  query <- CreateSeuratObject(counts = read.delim(query, 
                                                  sep = ",",
                                                  row.names = 1))
  
  # project and make predictions
  #### Modifiy annotation here!
  predictions <- TransferData(
    anchorset = FindTransferAnchors(
      reference = atlas, 
      query = query, 
      dims = 1:ndims, 
      normalization.method = "LogNormalize",
      reference.assay = "RNA", 
      query.assay = "RNA",
      verbose = FALSE), 
    refdata = atlas@meta.data[, annotation])
  
  # get predictions
  
  predictions <- predictions[, c("predicted.id", "prediction.score.max")]
  colnames(predictions) <- c("predicted cell type", "confidence score")
  
  return(predictions)
}


# function to obtain a frequency barplot for the predictions

plot_projections <- function(results){
  bplot <- ggplot(results, aes(x = `predicted cell type`, 
                               fill = `predicted cell type`)) +  
    geom_bar() +
    theme_minimal() +
    scale_fill_manual(values = complete_palette[results$`predicted cell type` %>% unique %>% sort]) +
    theme(text = element_text(size = 25),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.text = element_text(size = 16)) +
    guides(fill = guide_legend(ncol = 1, override.aes = list(size = 3)))
  
  return(bplot)
}