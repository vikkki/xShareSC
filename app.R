# Load required libraries
library(shiny)
library(ggplot2)
library(DT)

# Define UI
ui <- fluidPage(
  titlePanel("Basic Shiny App"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("num", 
                  "Select Number of Points:", 
                  min = 10, max = 1000, value = 500),
      selectInput("color", 
                  "Select Color:", 
                  choices = c("red", "blue", "green"), 
                  selected = "blue"),
      actionButton("refresh", "Refresh Plot")
    ),
    
    mainPanel(
      plotOutput("scatterPlot"),
      dataTableOutput("dataTable")
    )
  )
)

# Define Server
server <- function(input, output, session) {
  
  # Reactive data based on input
  dataset <- reactive({
    input$refresh  # Trigger when refresh button is clicked
    isolate({
      data.frame(
        x = rnorm(input$num),
        y = rnorm(input$num)
      )
    })
  })
  
  # Generate plot
  output$scatterPlot <- renderPlot({
    ggplot(dataset(), aes(x = x, y = y)) +
      geom_point(color = input$color) +
      theme_minimal() +
      labs(title = "Scatter Plot", x = "X-axis", y = "Y-axis")
  })
  
  # Generate table
  output$dataTable <- renderDataTable({
    datatable(dataset(), options = list(pageLength = 10))
  })
}

# Run the app
shinyApp(ui = ui, server = server)
