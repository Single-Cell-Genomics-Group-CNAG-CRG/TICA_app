
library(BiocManager)
library(shiny)
library(shinythemes)
library(Seurat)
library(tidyverse)
library(matchSCore2)
library(SPOTlight)

# load atlas
atlas <- readRDS("www/TICAtlas_complete_subset.rds")
# load markers 
atlas_markers <- readRDS("www/TICAtlas_markers.rds")

#######################################################

### Load functions
source("viz_tab.R")
source("projection.R")
source("annotation.R")
source("deconvolution.R")

#############################################################

### Define UI 

ui <- 
  fluidPage(
    theme = shinytheme("united"),
    
    ### TITLE
    
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
               sidebarLayout(
                 sidebarPanel(
                   selectizeInput(inputId = "annotLevel", "Variable to use for coloring", choices = c("level 1" = "lv1_annot", "level 2" = "lv2_annot", "cancer type" = "subtype"), multiple = FALSE),
                   selectizeInput(
                     inputId = "geneMarker",
                     "Gene to visualize",
                     selected = "CD3D",
                     choices = NULL,
                     options = list(create = TRUE, maxOptions = 5),
                     size = 1,
                     multiple = FALSE),
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
               sidebarLayout(
                 sidebarPanel(
                   fileInput(inputId = "projectionFile", "Upload your file", multiple = FALSE, accept = ".csv", width = NULL, buttonLabel = "Browse...", placeholder = "No file selected"),
                   actionButton(inputId = "projectionCheckFile", "Verify format"),
                   br(), 
                   br(),
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
                   br(), 
                   br(),
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
               h2("Deconvoltion using SPOTlight and the TICA", align = "center"),
               br(),
               h5("Using SPOTlight's NNMF-based deconvolution and the atlas as reference, we deconvolute mixture-like RNA-seq data into its proportions"),
               br(),
               ##
               sidebarLayout(
                 sidebarPanel(
                   fileInput(inputId = "deconvolutionFile", "Upload your file", multiple = FALSE, accept = ".csv", width = NULL, buttonLabel = "Browse...", placeholder = "No file selected"),
                   actionButton(inputId = "deconvolutionCheckFile", "Verify format"),
                   br(), br(),
                   selectInput(inputId = "deconvolutionLevel", "Annotation level to use", choices = c("level 1" = "lv1", "level 2" = "lv2_annot"), multiple = FALSE),
                   radioButtons("deconvolutionOrg", "Organism:",
                                c("Human" = "hsap",
                                  "Mouse" = "mmus"
                                )
                   ),
                   uiOutput(outputId = "deconvolutionRunAppear"),
                   br(), br(),
                   uiOutput(outputId = "deconvolutionDownloadAppear")
                 ),
                 mainPanel(plotOutput(outputId = "deconvolutionBarplot")
                 )
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
  
  # Visualization tab
  
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
    # Add & theme(strip.placement = NULL) bc patchwork has issues with certain themes
    # https://github.com/thomasp85/patchwork/issues/132
    pt <- (( pt1 | pt2 ) / pt3) | ggpubr::as_ggplot(legend) & ggplot2::theme(strip.placement = NULL)
    pt_grid <- pt + patchwork::plot_layout(widths = c(3, 1))
    # pt_grid <- (( pt1 | pt2 ) / pt3) | ggpubr::as_ggplot(legend) + patchwork::plot_layout(widths = c(2, 1))
    
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
    
    projected_data <- projection(atlas,
                                 query = path,
                                 annotation = level,
                                 ndims = dimensions
    )
    
    bplot <- plot_projections(projected_data)
    
    output$projectionBarplot <- renderPlot({bplot}, height = 1000)
    
    output$projectionRunAppear <- renderUI(
      downloadButton("projectionDownload","Download"),
    )
    
    output$projectionDownload <- downloadHandler(
        filename = function() {
          paste('TICA_projection_', Sys.Date(), '.csv', sep='')
        },
        content = function(con) {
          write.csv(projected_data, 
                    con,
                    quote = FALSE,
                    row.names = TRUE,
                    col.names = TRUE
                    )
        }
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
    
    hmplot <- plot_heatmap(atlas, atlas_markers, path, level, org)
    
    output$annotationHeatmap <- renderPlot({hmplot}, height = 800)
    
  })
  
  
  #### Deconvolution tab
  
  # validate input file and create run button
  
  observeEvent(input$deconvolutionCheckFile,{
    
    if(is.null(input$deconvolutionFile$datapath)){
      showModal(modalDialog(
        title = "Warning",
        "Please upload a csv file",
        easyClose = TRUE
      ))
    } else if(tools::file_ext(input$deconvolutionFile$datapath) != "csv"){
      showModal(modalDialog(
        title = "Warning",
        "Please upload a csv file",
        easyClose = TRUE
      ))
    } else if(testit::has_error(data.matrix(read.delim(input$deconvolutionFile$datapath, 
                                                       sep = ",")))){
      showModal(modalDialog(
        title = "Warning",
        "Please upload a valid csv file",
        easyClose = TRUE
      ))
      
    } else{
      output$deconvolutionRunAppear <- renderUI(
        actionButton("deconvolutionRun", "Run",
                     style = "color: #fff; background-color: #e95420; border-color: #c34113;
                                border-radius: 10px; 
                                border-width: 2px"),
      )
    }
  })
  
  
  # run deconvolution and plot
  
  observeEvent(input$deconvolutionRun, {
    
    file <- input$deconvolutionFile
    path <- file$datapath
    
    level <- input$deconvolutionLevel
    org <- input$deconvolutionOrg
    
    deconvoluted_data <- deconvolution(atlas,
                                       atlas_markers,
                                       query = path,
                                       annotation = level,
                                       organism = org
    )
    
    deconv_plot <- plot_spotlight(deconvoluted_data)
    
    output$deconvolutionBarplot <- renderPlot({deconv_plot}, height = 1000)
    
    output$deconvolutionRunAppear <- renderUI(
      downloadButton("deconvolutionDownload","Download"),
    )
    
    output$deconvolutionDownload <- downloadHandler(
      filename = function() {
        paste('TICA_deconvolution_', Sys.Date(), '.csv', sep='')
      },
      content = function(con) {
        write.csv(deconvoluted_data, 
                  con,
                  quote = FALSE,
                  row.names = TRUE,
                  col.names = TRUE
        )
      }
    )
    
    
  })
  
  
  
  
  
  
  
}

# shinyApp
shinyApp(ui, server)
