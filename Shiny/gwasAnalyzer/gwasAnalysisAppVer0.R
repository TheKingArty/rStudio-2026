# With help from tutorial by Data Professor

# Load R packages and datasets
library(shiny) # shiny
library(shinythemes) # themes
library(DT) # datatables
library(manhattanly)
library(plotly)

#library(bslib) #xtra
#library(palmerpenguins) #xtra
#library(dplyr) #xtra

data("airquality")

bmixdrinking <- read.csv("C:/Users/arthu/Documents/Research/localGWASdata/GCST90681941.h.tsv/GCST90681941.h.tsv", sep = "\t")

bchr1 <- subset(bmixdrinking, chromosome %in% 1)

bchr1ovrp5 <- subset(bchr1, p_value < 1e-5)


# Define UI
ui <- fluidPage(theme = shinytheme("cerulean"),
                navbarPage(
                  #theme = "cerulean"
                  "My first app",
                  tabPanel("Navbar 1",
                           sidebarPanel(
                             tags$h2("Input:"),
                             textInput("txt1", "First Name:", ""),
                             textInput("txt2", "Last Name:", ""),
                             
                           ), # sidebarPanel
                           mainPanel(
                             h1("Name"),
                             
                             h4("Full Name"),
                             verbatimTextOutput("txtout"),
                             
                           ) # mainPanel
                           
                  ), # Navbar 1, tabPanel
                  tabPanel("Navbar 2",
                           titlePanel("Interactive GWAS Manhattan Plot"),
                           
                           mainPanel(
                             # Use plotlyOutput to handle the htmlwidget
                             plotlyOutput("manhattanPlot", height = "600px")
                           )
                  ),
                  tabPanel("Navbar 3", 
                           dataTableOutput("table"),
                           "Hello")
                  
                  
                ) # navbarPage
) # fluidPage


# Define server function  
server <- function(input, output) {
  
  output$txtout <- renderText({
    paste( input$txt1, input$txt2, sep = " " )
  })
  output$manhattanPlot <- renderPlotly({
    # Use renderPlotly to wrap the manhattanly function call
    manhattanly(
      bchr1ovrp5,          # Swap with your actual dataset variable
      chr = "chromosome",
      bp = "base_pair_location",
      p = "p_value",
      snp = "rsid"
    )
  })
  output$table <- renderDataTable({datatable(bchr1ovrp5)})
} # server


# Create Shiny object
shinyApp(ui = ui, server = server)