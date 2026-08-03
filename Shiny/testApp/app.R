# With help from tutorial by Data Professor

# Load R packages and datasets
library(shiny)
library(shinythemes)
data("airquality")


# Define UI
ui <- fluidPage(theme = shinytheme("darkly"),
                navbarPage(
                  #theme = "cerulean",  # <--- To use a theme, uncomment this
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
                           sidebarLayout(
                             
                             # Sidebar panel for inputs
                             sidebarPanel(
                               
                               # Input: Slider for the number of bins
                               sliderInput(inputId = "bins",
                                           label = "Number of bins:",
                                           min = 1,
                                           max = 50,
                                           value = 30)
                               
                             ),
                             
                             # Main panel for displaying outputs
                             mainPanel(
                               
                               # Output: Histogram ----
                               plotOutput(outputId = "distPlot")
                               
                             )
                           )
                           ),
                  tabPanel("Navbar 3", "Hello")
                  
                ) # navbarPage
) # fluidPage


# Define server function  
server <- function(input, output) {
  
  output$txtout <- renderText({
    paste( input$txt1, input$txt2, sep = " " )
  })
  output$distPlot <- renderPlot({
    
    x    <- airquality$Ozone
    x    <- na.omit(x)
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    hist(x, breaks = bins, col = "#75AADB", border = "black",
         xlab = "Ozone level",
         main = "Histogram of Ozone level")
  })
} # server


# Create Shiny object
shinyApp(ui = ui, server = server)