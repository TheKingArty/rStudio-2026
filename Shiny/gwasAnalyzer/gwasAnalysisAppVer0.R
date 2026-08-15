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

# bmixdrinking <- read.csv("C:/Users/arthu/Documents/Research/localGWASdata/GCST90681941.h.tsv/GCST90681941.h.tsv", sep = "\t")

# fread

# bchr1 <- subset(bmixdrinking, chromosome %in% 1)

# bchr1ovrp5 <- subset(bchr1, p_value < 1e-5)

# sys1 <- read.csv("C:/Users/arthu/Documents/Research/localGWASdata/gwasTestSystems/gwas_sys1.B1.glm.logistic", sep = "\t")

# sys1chr1 <- subset(sys1, CHROM %in% 1)

# sys1chr1ovrp5 <- subset(sys1chr1, P < 1e-5)

#sys2 <- read.csv("C:/Users/arthu/Documents/Research/localGWASdata/gwasTestSystems/gwas_sys2.B1.glm.logistic", sep = "\t")

#sys2chr1 <- subset(sys2, CHROM %in% 1)

#sys2chr1ovrp5 <- subset(sys2chr1, P < 1e-5)

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
                             plotlyOutput("manhattanPlot1", height = "600px"),
                             plotlyOutput("manhattanPlot2", height = "600px")
                           )
                  ),
                  tabPanel("Navbar 3", 
                           dataTableOutput("table1"),
                           dataTableOutput("table2"),
                           "Hello")
                  
                  
                ) # navbarPage
) # fluidPage


# Define server function  
server <- function(input, output) {
  
  output$txtout <- renderText({
    paste( input$txt1, input$txt2, sep = " " )
  })
  output$manhattanPlot1 <- renderPlotly({
    # Use renderPlotly to wrap the manhattanly function call
    manhattanly(
      sys1chr1ovrp5,          # Swap with your actual dataset variable
      chr = "CHROM",
      bp = "POS",
      p = "P",
      snp = "ID"
    )
  })
  output$manhattanPlot2 <- renderPlotly({
    # Use renderPlotly to wrap the manhattanly function call
    manhattanly(
      sys2chr1ovrp5,          # Swap with your actual dataset variable
      chr = "CHROM",
      bp = "POS",
      p = "P",
      snp = "ID"
    )
  })
  output$table1 <- renderDataTable({datatable(sys1chr1ovrp5)})
  output$table2 <- renderDataTable({datatable(sys2chr1ovrp5)})
} # server


# Create Shiny object
shinyApp(ui = ui, server = server)