
library(BiocManager)
library(shiny)
library(shinythemes)
library(Seurat)
library(tidyverse)
library(matchSCore2)

# cell_type_palette <- readRDS("www/complete_cell_type_palette.rds")
# cancer_type_palette <- readRDS("www/cancer_type_palette.rds")

# load selected annotation level
atlas <- readRDS("www/TICAtlas_complete_subset.rds")

# Load markers
atlas_markers1 <- readRDS("www/TICAtlas_markers_lv1.rds")
atlas_markers2 <- readRDS("www/TICAtlas_markers_lv2.rds")

complete_palette <- readRDS("www/complete_palette.rds")

#######################################################

### Load Visualization functions
source("viz_tab.R")

### Define functions

# function to project queries on the atlas to predict cell-wise annotation

projection <- function(query, annotation, ndims) {

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
    # refdata = atlas$annotation,
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


# function to create jaccard index plot from matchSCore

plot_heatmap <- function(marker_genes, annotation, organism){
  
  # load selected annotation level
  if(annotation == "lv1"){
    atlas_markers <- atlas_markers1
  }
  if(annotation == "lv2"){
    atlas_markers <- atlas_markers2
  }
  
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

#############################################################

### Define UI 

ui <- 
  fluidPage(
    theme = shinytheme("united"),
    
    ### TITLE
    
    #titlePanel(title=div(img(src="www/cnag-crg-logo.jpg", width="300"), "Tumor Immune Cell Atlas")),
    
    titlePanel(
      fluidRow(
        column(5, img(height = 160, src = "tica.jpg")),
        column(5, "Tumor Immune Cell Atlas"), 
        column(2, img(width = 300, src = "cnag-crg-logo.jpg"))
      )
      
    ),
    
    ### TABS
    
    tabsetPanel(
      type = "pills",
      
      ### TAB HOME
      
      tabPanel("Home", 
               br(),
               br(), 
               fluidRow(
                 column(2, ""),
                 column(8, 
                        includeHTML("home.html")
                        )
                 ),
               br(),
               hr(),
               div(
                 class = "footer",
                 includeHTML("footer.html")
               )
               #
               #img(src="www/tica.jpg", height="90%", width="90%", align="center")
               
      ),
      
      ### TAB HOW TO USE
      
      tabPanel("How to use",
               br(),
               br(), 
               fluidRow(
                 column(2, ""),
                 column(8, 
                        includeHTML("how_to.html")
                 )
               ),
               br(),
               hr(),
               div(
                 class = "footer",
                 includeHTML("footer.html")
               )
      ),
      
      ### TAB Visualization
      
      tabPanel("Visualization",
               h2("Visualize genes on the TICA", align = "center"),
               br(),
               # h5("Using Seurat LabelTransfer we project TICA's cell-types onto new datasets"),
               # br(),
               ##
               sidebarLayout(
                 sidebarPanel(
                   # fileInput(inputId = "projectionFile", "Upload your file", multiple = FALSE, accept = ".csv", width = NULL, buttonLabel = "Browse...", placeholder = "No file selected"),
                   # actionButton(inputId = "projectionCheckFile", "Verify format"),
                   # br(), br(),
                   selectizeInput(inputId = "annotLevel", "Annotation level to use", choices = c("level 1" = "lv1_annot", "level 2" = "lv2_annot", "cancer type" = "subtype"), multiple = FALSE),
                   selectizeInput(
                     inputId = "geneMarker",
                     "Gene to visualize",
                     selected = "CD3D",
                     choices = NULL,
                     options = list(create = TRUE, maxOptions = 5),
                     size = 1,
                     multiple = FALSE),
                   # sliderInput(inputId = "projectionDims", 
                   #             label = "Select number of dimensions to use", 
                   #             min = 1, 
                   #             max = 50, 
                   #             value = 25, 
                   #             step = 1,
                   #             round = TRUE
                   # ),
                   # uiOutput(outputId = "projectionRunAppear"),
                   # br(), br(),
                   # uiOutput(outputId = "projectionDownloadAppear")
                 ),
                 mainPanel(
                   plotOutput(outputId = "VizPlts")
                 )
               ),
               br(),
               hr(),
               div(
                 class = "footer",
                 includeHTML("footer.html")
               )
      ),
      
      ### TAB PROJECTION
      
      tabPanel("Projection",
               h2("Projection of new datasets onto the TICA", align = "center"),
               br(),
               h5("Using Seurat LabelTransfer we project TICA's cell-types onto new datasets"),
               br(),
               ##
               sidebarLayout(
                 sidebarPanel(
                   fileInput(inputId = "projectionFile", "Upload your file", multiple = FALSE, accept = ".csv", width = NULL, buttonLabel = "Browse...", placeholder = "No file selected"),
                   actionButton(inputId = "projectionCheckFile", "Verify format"),
                   br(), br(),
                   selectInput(inputId = "projectionLevel", "Annotation level to use", choices = c("level 1" = "lv1_annot", "level 2" = "lv2_annot"), multiple = FALSE),
                   sliderInput(inputId = "projectionDims", 
                               label = "Select number of dimensions to use", 
                               min = 1, 
                               max = 50, 
                               value = 25, 
                               step = 1,
                               round = TRUE
                               ),
                   uiOutput(outputId = "projectionRunAppear"),
                   br(), br(),
                   uiOutput(outputId = "projectionDownloadAppear")
                 ),
                 mainPanel(plotOutput(outputId = "projectionBarplot")
                           )
               ),
               br(),
               hr(),
               div(
                 class = "footer",
                 includeHTML("footer.html")
               )
      ),
      
      ### TAB ANNOTATION
      
      tabPanel("Annotation",
               h2("Annotation of new datasets using the TICA", align = "center"),
               br(),
               h5("Using matchSCore we compare input cluster markers to the TICA's cell-type markers"),
               br(),
               sidebarLayout(
                 sidebarPanel(
                   fileInput("annotationFile", "Upload your markers", multiple = FALSE, accept = ".csv", width = NULL, buttonLabel = "Browse...", placeholder = "No file selected"),
                   actionButton(inputId = "annotationCheckFile", "Verify format"),
                   br(), br(),
                   radioButtons("annotationOrg", "Organism:",
                                c("Human" = "hsap",
                                  "Mouse" = "mmus"
                                  )
                                ),
                   br(),
                   radioButtons("annotationLevel", "Annotation level:",
                                c("Level 1" = "lv1",
                                  "Level 2" = "lv2"
                                )
                   ),
                   br(),
                   uiOutput(outputId = "annotationRunAppear")
                 ),
                 mainPanel(
                   plotOutput("annotationHeatmap")
                 )
               ),
               br(),
               hr(),
               div(
                 class = "footer",
                 includeHTML("footer.html")
               )
      ),
      
      ### TAB DECONVOLUTION
      
      tabPanel("Deconvolution",
               tabPanel("Annotation",
                        h2("Deconvolution using the TICA", align = "center"),
                        br(),
               #         h5("Using the TICA to deconvolute spot-mixtures in spatial transcriptomic datasets"),
                        br(),
                        h3("This functionality is under development and will be available soon.", align = "center"),
                        h3("Sorry for the inconvenience!", align = "center"),
                        br(),
                        br()
               #          sidebarLayout(
               #            sidebarPanel(
               #              fileInput("deconvFile", "Upload your spots", multiple = FALSE, accept = NULL, width = NULL, buttonLabel = "Browse...", placeholder = "No file selected"),
               #              #sliderInput("degInput", "Differentially Expressed Genes", 0, 10000, 3000, step = 500),
               #              #radioButtons("normInput", "Normalization method", choices = c("SCT", "LogNormalize", "None"), selected = "None"),
               #              selectInput("orgInput", "Organism", choices = c("Human", "Mouse")),
               #              actionButton("button", "Run")
               #            ),
               #            mainPanel(img(src="www/spatial_stratification.PNG", height="100%", width="100%", align="right"))
                        ),
               br(),
               hr(),
               div(
                 class = "footer",
                 includeHTML("footer.html")
               )
      ),
      
      
      ### TAB CONTACT 
      
      tabPanel("Contact", 
               br(),
               br(), 
               fluidRow(
                 column(2, ""),
                 column(8, 
                        includeHTML("contact.html")
                 )
               ),
               br(),
               hr(),
               div(
                 class = "footer",
                 includeHTML("footer.html")
               )
      )
    )
  )


#############################################################

### SERVER

server <- function(input, output, session) {
  # Setting maximum file size to 8GB
  options(shiny.maxRequestSize = 8000 * 1024 ^ 2)
  
  # Visualization
  
  # Use this trick from the website below to be able to load all the rownames
  # https://shiny.rstudio.com/articles/selectize.html
  # You may use choices = NULL to create an empty selectize instance, so that it will load quickly initially, then use updateSelectize(server = TRUE) to pass the choices
  updateSelectizeInput(session, 'geneMarker', choices = rownames(atlas), server = TRUE)
  
  output$VizPlts <- renderPlot({
    
    # Build the pices of the final plot
    pt1 <- DimPlot_fun(atlas = atlas, group = input$annotLevel, x = "UMAP_1", y = "UMAP_2")
    pt2 <- FeatPlot_fun(atlas = atlas, gene = input$geneMarker, x = "UMAP_1", y = "UMAP_2")
    pt3 <- VlnPlot_fun(atlas = atlas, group = input$annotLevel, gene = input$geneMarker) +
      ggplot2::theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))
    legend <- ggpubr::get_legend(pt1)
    pt1 <- pt1 + Seurat::NoLegend()
    
    # Put them together with patchwork
    pt_grid <- (( pt1 | pt2 ) / pt3) | ggpubr::as_ggplot(legend)
    
    # Return plot
    pt_grid
    }, height = 750)
  
  #### Projection tab
  
  # validate input file and create run button
  
  observeEvent(input$projectionCheckFile,{
    
    if(is.null(input$projectionFile$datapath)){
      showModal(modalDialog(
        title = "Warning",
        "Please upload a csv file",
        easyClose = TRUE
      ))
    } else if(tools::file_ext(input$projectionFile$datapath) != "csv"){
      showModal(modalDialog(
        title = "Warning",
        "Please upload a csv file",
        easyClose = TRUE
      ))
    } else if(testit::has_error(data.matrix(read.delim(input$projectionFile$datapath, 
                                                       sep = ",")))){
      showModal(modalDialog(
        title = "Warning",
        "Please upload a valid csv file",
        easyClose = TRUE
      ))
      
    } else{
      output$projectionRunAppear <- renderUI(
        actionButton("projectionRun", "Run",
                     style = "color: #fff; background-color: #e95420; border-color: #c34113;
                                border-radius: 10px; 
                                border-width: 2px"),
      )
    }
  })
  
  
  # run projections and plot
  
  observeEvent(input$projectionRun, {

    file <- input$projectionFile
    path <- file$datapath

    level <- input$projectionLevel
    dimensions <- input$projectionDims

    projected_data <- projection(query = path,
                                 annotation = level,
                                 ndims = dimensions
                                 )
    
    bplot <- plot_projections(projected_data)
    
    output$projectionBarplot <- renderPlot({bplot}, height = 1000)
    
    output$projectionRunAppear <- renderUI(
      downloadButton("projectionDownload","Download"),
    )
    
  })
  
  
  

  # Annotation tab
  
  # validate input file and create run button
  
  observeEvent(input$annotationCheckFile,{
    
    if(is.null(input$annotationFile$datapath)){
      showModal(modalDialog(
        title = "Warning",
        "Please upload a csv file",
        easyClose = TRUE
      ))
    } else if(tools::file_ext(input$annotationFile$datapath) != "csv"){
      showModal(modalDialog(
        title = "Warning",
        "Please upload a csv file",
        easyClose = TRUE
      ))
    } else if(!c("cluster", "gene") %in% colnames(read.delim(input$annotationFile$datapath,
                                                             header = TRUE,
                                                             sep = ","))){
      showModal(modalDialog(
        title = "Warning",
        "Please upload a valid csv file",
        easyClose = TRUE
      ))
      
    } else{
      output$annotationRunAppear <- renderUI(
        actionButton("annotationRun","Run",
                     style="color: #fff; background-color: #e95420; border-color: #c34113;
                                border-radius: 10px; 
                                border-width: 2px"),
      )
    }
  })
  
  
  # run MatchSCore2 and plot
  
  observeEvent(input$annotationRun, {
    
    file <- input$annotationFile
    path <- file$datapath
    
    level <- input$annotationLevel
    org <- input$annotationOrg
    
    hmplot <- plot_heatmap(path, level, org)
    
    output$annotationHeatmap <- renderPlot({hmplot}, height = 800)
    
  })
  
}

# shinyApp
shinyApp(ui, server)
