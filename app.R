
library(BiocManager)
library(shiny)
library(shinythemes)
library(Seurat)
library(tidyverse)
library(matchSCore2)

cell_type_palette <- readRDS("www/complete_cell_type_palette.rds")

#######################################################

### Define functions

# function to project querys on the atlas to predict cell-wise annotation

projection <- function(query, annotation, ndims) {

  query <- CreateSeuratObject(counts = read.delim(query, 
                                                  sep = ","))
  
  # load selected annotation level
  
  if(annotation == "level 1"){
    atlas <- readRDS("www/TICAtlas_lv1_subset.rds")
  }
  if(annotation == "level 2"){
    atlas <- readRDS("www/TICAtlas_lv2_subset.rds")
  }
  
  # project and make predictions
  
  predictions <- TransferData(anchorset = FindTransferAnchors(reference = atlas, 
                                                              query = query, 
                                                              dims = 1:ndims, 
                                                              normalization.method = "LogNormalize",
                                                              reference.assay = "RNA", 
                                                              query.assay = "RNA",
                                                              verbose = FALSE), 
                              refdata = atlas$annotation)
  
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
    scale_fill_manual(values = cell_type_palette) +
    theme(text = element_text(size = 25),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.text=element_text(size=16)) +
    guides(fill=guide_legend(ncol=1, override.aes = list(size=3)))

    
  return(bplot)
}


# function to create jaccard index plot from matchSCore

plot_heatmap <- function(marker_genes, annotation, organism){
  
  # load selected annotation level
  
  if(annotation == "lv1"){
    atlas_markers <- readRDS("www/TICAtlas_markers_lv1.rds")
  }
  if(annotation == "lv2"){
    atlas_markers <- readRDS("www/TICAtlas_markers_lv2.rds")
  }
  
  # read query marker genes
  
  marker_genes <- read.delim(marker_genes,
                             header = TRUE,
                             sep = ",")
  
  if(organism == "mmus"){
    marker_genes$gene <- str_to_upper(marker_genes$gene)
  }
  
  marker_genes <- split(marker_genes$gene, marker_genes$cluster)
  
  ms <- matchSCore2(gene_cl.ref = atlas_markers[sort(names(atlas_markers))], 
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
                   selectInput(inputId = "projectionLevel", "Annotation level to use", choices = c("level 1", "level 2"), multiple = FALSE),
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
                        h3("This functionality is under development and will be available soon. Sorry for the inconvenience!", align = "center"),
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
               HTML("<p>You can find more information about the Tumor Immune Cell Atlas in <a href='https://github.com/Single-Cell-Genomics-Group-CNAG-CRG/Tumor-Immune-Cell-Atlas'>our github repository</a>! where you can also directly contact us or ask your questions directly. We appreciate every kind of suggestion, comment or discussion and we thank you for using our tool. </p>"),
               br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
               HTML("<p>This web App is currently actively under development. Thanks for understanding!</p>"),
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

server <- function(input, output) {
  
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
        actionButton("projectionRun","Run",
                     style="color: #fff; background-color: #e95420; border-color: #c34113;
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
