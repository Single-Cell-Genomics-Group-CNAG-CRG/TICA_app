### Define functions for the deconvolution tab

# load_palette
complete_palette <- readRDS("www/complete_palette.rds")

# function to use SPOTlight for deconvolution

deconvolution <- function(atlas, atlas_markers, query, annotation, organism) {
  
  # load selected annotation level
  
  atlas_markers <- filter(atlas_markers, level == annotation)
  
  # read query to deconvolute
  query <- read.delim(query,
                      header = TRUE,
                      sep = ",",
                      row.names = 1)
  
  if(organism == "mmus"){
    rownames(query) <- stringr::str_to_upper(rownames(query))
  }
  
  # Run SPOTlight
  spotlight_results <- spotlight_deconvolution(
    se_sc = atlas,
    counts_spatial = query,
    clust_vr = paste0(annotation, "_annot"),
    cluster_markers = atlas_markers,
    cl_n = 50,
    hvg = 5000,
    ntop = NULL
  )
  
  spotlight_results <- spotlight_results[[2]]
  return(spotlight_results)
}

plot_spotlight <- function(spotlight_results){
  
  # extract deconvolution results
  
  spotlight_results <- spotlight_results[,which(colnames(spotlight_results) != "res_ss")]
  spotlight_results <- as.data.frame(spotlight_results)
  rownames(spotlight_results) <- paste0("sample_", 1:nrow(spotlight_results))
  spotlight_results$sample <- rownames(spotlight_results)
  spotlight_results <- pivot_longer(as.data.frame(spotlight_results), cols = colnames(spotlight_results)[colnames(spotlight_results) != "sample"], names_to = "cell_type")
  spotlight_results <- dplyr::filter(spotlight_results, value > 0)
  spotlight_results$cell_type <- str_replace_all(spotlight_results$cell_type, pattern = "\\.", " ")
  
  spotlight_plot <- ggplot(spotlight_results, 
                           aes(x = sample,
                               y = value,
                               fill = cell_type)) +  
    geom_bar(stat="identity") +
    theme_minimal() +
    scale_fill_manual(values = complete_palette[spotlight_results$cell_type %>% unique %>% sort]) +
    theme(text = element_text(size = 25),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.text = element_text(size = 16)) +
    guides(fill = guide_legend(ncol = 1, override.aes = list(size = 3))) +
    labs(title = "SPOTlight deconvolution",
         x = "samples",
         y = "composition",
         fill = "")
  
  return(spotlight_plot)
  
}