library(shiny)
library(shinyjs)
library(Seurat)
library(ggplot2)
library(dplyr)
#library(patchwork)
library(bslib)
library(fontawesome)

# Define the directory containing Seurat objects
data_dir <- "data/"
available_files <- list.files(data_dir, pattern = "\\.rds$", full.names = FALSE)

# Set the correct password
password_env <- Sys.getenv("SHINY_APP_PASSWORD")

# Define UI with password authentication
ui <- fluidPage(
  theme = bs_theme(bootswatch = "flatly"),  # Use Flatly Bootstrap theme
  useShinyjs(),  # Enable JavaScript for UI actions
  
  uiOutput("auth_ui"),  # This will dynamically show either login modal or the app
)

server <- function(input, output, session) {
  
  # Reactive value to track login status
  user_auth <- reactiveVal(FALSE)
  
  # Show password modal on startup
  observe({
    showModal(
      modalDialog(
        title = "Enter Password",
        passwordInput("password", "Password:"),
        footer = tagList(
          modalButton("Cancel"),
          actionButton("submit_password", "Submit")
        ),
        easyClose = FALSE  # Prevent closing without entering password
      )
    )
  })
  
  # Check password when user submits
  observeEvent(input$submit_password, {
    if (input$password == password_env) {
      user_auth(TRUE)
      removeModal()  # Remove the modal if the password is correct
    } else {
      showNotification("Incorrect password! Try again.", type = "error")
    }
  })
  
  # Dynamically display the app UI if logged in
  output$auth_ui <- renderUI({
    if (user_auth()) {
      # Main App UI
      div(class = "container-custom",
          div(class = "card",
              fluidRow(
                column(6, 
                       selectInput("seurat_file", "Select Seurat Object:", choices = available_files)),
                column(6, 
                       uiOutput("metadata_selector"))
              ),
              
              fluidRow(
                column(8, 
                       selectizeInput("gene_select", "Select Gene(s):", choices = NULL, multiple = TRUE,
                                      options = list(placeholder = "Type to search...", maxOptions = 100))),
                column(4, 
                       actionButton("update", "Update Visualization", class = "btn btn-primary btn-custom"))
              ),
              
              # Status Message
              div(class = "loading-status",
                  uiOutput("status_message"))
          ),
          
          # Visualization Panel
          div(class = "card",
              tabsetPanel(
                tabPanel("UMAP Plot", plotOutput("umapPlot", height = "500px")),
                tabPanel("Feature Plot", uiOutput("featurePlots")),
                tabPanel("Dot Plot", plotOutput("dotPlot", height = "500px"))
              )
          ),
          
          footer = div(class = "text-center mt-5 p-3", 
                       tags$p("© 2025 Seurat Data Visualization (XiaofanZhao, XiaLab OHSU). Powered by Shiny."))
      )
    } else {
      div()  # Empty UI until password is entered
    }
  })
  
  # Load the selected Seurat object
  seurat_data <- reactive({
    req(input$seurat_file)
    output$status_message <- renderUI({ icon("spinner", class = "fa-spin") })  # Show loading icon
    file_path <- file.path(data_dir, input$seurat_file)
    data <- readRDS(file_path)
    output$status_message <- renderUI({ icon("check-circle", style = "color: green;") })  # Indicate success
    return(data)
  })
  
  # Update UI elements based on loaded data
  observeEvent(input$seurat_file, {
    seurat_obj <- seurat_data()
    
    # Update gene selection dropdown with default selection "GNLY"
    default_gene <- if ("GNLY" %in% rownames(seurat_obj)) "GNLY" else NULL
    
    updateSelectizeInput(session, "gene_select", 
                         choices = rownames(seurat_obj), 
                         selected = default_gene,  
                         server = TRUE)
    
    # Filter metadata to include only categorical variables (factors or characters)
    metadata_columns <- colnames(seurat_obj@meta.data)
    categorical_columns <- metadata_columns[sapply(seurat_obj@meta.data, function(x) is.factor(x) || is.character(x))]
    
    if (length(categorical_columns) == 0) {
      showNotification("No categorical metadata found.", type = "warning")
      categorical_columns <- NULL
    }
    
    output$metadata_selector <- renderUI({
      selectInput("metadata", "Select Grouping Metadata:", 
                  choices = categorical_columns, 
                  selected = ifelse(length(categorical_columns) > 0, categorical_columns[1], NULL))
    })
  })
  
  # Render UMAP plot
  output$umapPlot <- renderPlot({
    req(seurat_data(), input$metadata)
    
    seurat_obj <- seurat_data()
    
    # Check if 'umap' exists, otherwise use 'umap.unintegrated'
    reduction_to_use <- if ("umap" %in% names(seurat_obj@reductions)) {
      "umap"
    } else if ("umap.unintegrated" %in% names(seurat_obj@reductions)) {
      "umap.unintegrated"
    } else {
      showNotification("No valid UMAP reduction found.", type = "error")
      return(NULL)
    }
    
    DimPlot(seurat_obj, reduction = reduction_to_use, group.by = input$metadata) + 
      ggtitle(paste("UMAP Visualization using", reduction_to_use))
  })
  
  # Render feature plots for all selected genes
  output$featurePlots <- renderUI({
    req(seurat_data(), input$gene_select)
    
    gene_plots <- lapply(input$gene_select, function(gene) {
      if (gene %in% rownames(seurat_data())) {
        plotOutput(outputId = paste0("featurePlot_", gene), height = "500px")
      } else {
        showNotification(paste("Gene", gene, "not found in dataset."), type = "error")
      }
    })
    
    do.call(tagList, gene_plots)
  })
  
  # Generate feature plots dynamically
  observe({
    req(seurat_data(), input$gene_select)
    
    for (gene in input$gene_select) {
      local({
        gene_name <- gene  # Capture the current gene
        output[[paste0("featurePlot_", gene_name)]] <- renderPlot({
          FeaturePlot(seurat_data(), features = gene_name) + 
            ggtitle(paste("Feature Plot for", gene_name))
        })
      })
    }
  })
  
  # Render dot plot for selected genes
  output$dotPlot <- renderPlot({
    req(seurat_data(), input$gene_select)
    DotPlot(seurat_data(), features = input$gene_select, group.by = input$metadata) + 
      scale_color_gradient2(low = "#9a031e", mid = "#EDEDED", high = "#90a955", midpoint = 0.5) +
      ggtitle("Dot Plot for Selected Genes")
  })
}

# Run the application
shinyApp(ui = ui, server = server)
