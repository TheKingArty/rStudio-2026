# With help from tutorial by Data Professor

# Load R packages and datasets
library(shiny) # shiny
library(shinythemes) # themes
library(DT) # datatables
library(manhattanly)
library(data.table)
library(plotly)

#set max upload size 1GB

options(shiny.maxRequestSize = 1000 * 1024^2)

#library(bslib) #xtra
#library(palmerpenguins) #xtra
#library(dplyr) #xtra

#bmixdrinking <- read.csv("C:/Users/arthu/Documents/Research/localGWASdata/GCST90681941.h.tsv/GCST90681941.h.tsv", sep = "\t")

# fread

#bchr1 <- subset(bmixdrinking, chromosome %in% 1)

#bchr1ovrp5 <- subset(bchr1, p_value < 1e-5)

#sys1 <- fread("C:/Users/arthu/Documents/Research/localGWASdata/gwasTestSystems/gwas_sys1.B1.glm.logistic", sep = "\t")

#sys1chr1 <- subset(sys1, CHROM %in% 1)

#sys1chr1ovrp5 <- subset(sys1chr1, P < 1e-5)

#sys2 <- fread("C:/Users/arthu/Documents/Research/localGWASdata/gwasTestSystems/gwas_sys2.B1.glm.logistic", sep = "\t")

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
                           fluidRow(
                             column(
                               width = 6,
                               fileInput("sys1_file", "Upload Sys 1 Dataset", accept = c(".tsv", ".logistic", ".txt", ".csv"))
                             ),
                             column(
                               width = 6,
                               fileInput("sys2_file", "Upload Sys 2 Dataset", accept = c(".tsv", ".logistic", ".txt", ".csv"))
                             )
                           ),
                           
                           hr(), #seperator
                           
                           # Row 2: Dynamic Column Mapping & Filters
                           fluidRow(
                             # System 1 Controls
                             column(6,
                                    wellPanel(
                                      h4("System 1 Settings"),
                                      fluidRow(
                                        column(6, textInput("sys1_chr_col", "Chromosome Column Name", value = "CHROM")),
                                        column(6, textInput("sys1_bp_col", "Position (BP) Column Name", value = "POS"))
                                      ),
                                      fluidRow(
                                        column(6, textInput("sys1_p_col", "P-Value Column Name", value = "P")),
                                        column(6, textInput("sys1_snp_col", "SNP ID Column Name", value = "ID"))
                                      ),
                                      fluidRow(
                                        column(6, textInput("sys1_chr_filter", "Chromosome Filter (e.g. 1, or leave blank for All)", value = "1")),
                                        column(6, numericInput("sys1_p_thresh", "Max P-Value Threshold", value = 1e-5, step = 1e-6))
                                      )
                                    )
                             ),
                             
                             # System 2 Controls
                             column(6,
                                    wellPanel(
                                      h4("System 2 Settings"),
                                      fluidRow(
                                        column(6, textInput("sys2_chr_col", "Chromosome Column Name", value = "CHROM")),
                                        column(6, textInput("sys2_bp_col", "Position (BP) Column Name", value = "POS"))
                                      ),
                                      fluidRow(
                                        column(6, textInput("sys2_p_col", "P-Value Column Name", value = "P")),
                                        column(6, textInput("sys2_snp_col", "SNP ID Column Name", value = "ID"))
                                      ),
                                      fluidRow(
                                        column(6, textInput("sys2_chr_filter", "Chromosome Filter (e.g. 1, or leave blank for All)", value = "1")),
                                        column(6, numericInput("sys2_p_thresh", "Max P-Value Threshold", value = 1e-5, step = 1e-6))
                                      )
                                    )
                             )
                           ),
                           
                           hr(),
                           
                           fluidRow(
                             column(
                               width = 6,
                               plotlyOutput("manhattanPlot1",
                                            height = "600px")
                             ),
                             column(
                               width = 6,
                               plotlyOutput("manhattanPlot2",
                                            height = "600px")
                             )
                           )
                           
                  ),
                  tabPanel("Navbar 3", 
                           dataTableOutput("table1"),
                           dataTableOutput("table2"),
                           "Hello")
                  
                  
                ) # navbarPage
) # fluidPage



# Define server function  
server <- function(input, output, session) {
  
  output$txtout <- renderText({
    paste( input$txt1, input$txt2, sep = " " )
  })
  
  # --- SYSTEM 1 LOGIC ---
  
  # 1. Read raw file 1
  sys1_raw <- reactive({
    req(input$sys1_file)
    fread(input$sys1_file$datapath)
  })
  
  # 2. Automatically update column choices when File 1 is uploaded
  observeEvent(sys1_raw(), {
    cols <- names(sys1_raw())
    
    # Try to pick sensible default selections if column names exist
    chr_default <- grep("chr|chrom", cols, ignore.case = TRUE, value = TRUE)[1]
    bp_default  <- grep("pos|bp", cols, ignore.case = TRUE, value = TRUE)[1]
    p_default   <- grep("^p$|p_val|p.val|pval", cols, ignore.case = TRUE, value = TRUE)[1]
    snp_default <- grep("snp|id|rs", cols, ignore.case = TRUE, value = TRUE)[1]
    
    updateSelectInput(session, "sys1_chr_col", choices = cols, selected = ifelse(is.na(chr_default), cols[1], chr_default))
    updateSelectInput(session, "sys1_bp_col", choices = cols, selected = ifelse(is.na(bp_default), cols[1], bp_default))
    updateSelectInput(session, "sys1_p_col", choices = cols, selected = ifelse(is.na(p_default), cols[1], p_default))
    updateSelectInput(session, "sys1_snp_col", choices = cols, selected = ifelse(is.na(snp_default), cols[1], snp_default))
  })
  
  # 3. Process and filter Dataset 1 reactively
  inp1 <- reactive({
    req(sys1_raw(), input$sys1_chr_col, input$sys1_p_col)
    
    df <- sys1_raw()
    chr_col <- input$sys1_chr_col
    p_col <- input$sys1_p_col
    
    # Filter Chromosome if user typed something in
    if (nzchar(input$sys1_chr_filter)) {
      df <- df[get(chr_col) %in% input$sys1_chr_filter]
    }
    
    # Filter P-Value
    if (!is.null(input$sys1_p_thresh) && !is.na(input$sys1_p_thresh)) {
      df <- df[get(p_col) < as.numeric(input$sys1_p_thresh)]
    }
    
    df
  })
  
  # --- SYSTEM 2 LOGIC ---
  
  # 1. Read raw file 2
  sys2_raw <- reactive({
    req(input$sys2_file)
    fread(input$sys2_file$datapath)
  })
  
  # 2. Automatically update column choices when File 2 is uploaded
  observeEvent(sys2_raw(), {
    cols <- names(sys2_raw())
    
    chr_default <- grep("chr|chrom", cols, ignore.case = TRUE, value = TRUE)[1]
    bp_default  <- grep("pos|bp", cols, ignore.case = TRUE, value = TRUE)[1]
    p_default   <- grep("^p$|p_val|p.val|pval", cols, ignore.case = TRUE, value = TRUE)[1]
    snp_default <- grep("snp|id|rs", cols, ignore.case = TRUE, value = TRUE)[1]
    
    updateSelectInput(session, "sys2_chr_col", choices = cols, selected = ifelse(is.na(chr_default), cols[1], chr_default))
    updateSelectInput(session, "sys2_bp_col", choices = cols, selected = ifelse(is.na(bp_default), cols[1], bp_default))
    updateSelectInput(session, "sys2_p_col", choices = cols, selected = ifelse(is.na(p_default), cols[1], p_default))
    updateSelectInput(session, "sys2_snp_col", choices = cols, selected = ifelse(is.na(snp_default), cols[1], snp_default))
  })
  
  # 3. Process and filter Dataset 2 reactively
  inp2 <- reactive({
    req(sys2_raw(), input$sys2_chr_col, input$sys2_p_col)
    
    df <- sys2_raw()
    chr_col <- input$sys2_chr_col
    p_col <- input$sys2_p_col
    
    # Filter Chromosome if user typed something in
    if (nzchar(input$sys2_chr_filter)) {
      df <- df[get(chr_col) %in% input$sys2_chr_filter]
    }
    
    # Filter P-Value
    if (!is.null(input$sys2_p_thresh) && !is.na(input$sys2_p_thresh)) {
      df <- df[get(p_col) < as.numeric(input$sys2_p_thresh)]
    }
    
    df
  })
  
  # --- PLOTS & TABLES ---
  
  output$manhattanPlot1 <- renderPlotly({
    req(inp1(), input$sys1_chr_col, input$sys1_bp_col, input$sys1_p_col, input$sys1_snp_col)
    p1 <- manhattanly(
      inp1(),
      chr = input$sys1_chr_col,
      bp = input$sys1_bp_col,
      p = input$sys1_p_col,
      snp = input$sys1_snp_col
    )
    # Assign a source ID to trace zoom/pan events
    p1$x$source <- "manhattanPlot1"
    p1
  })
  
  output$manhattanPlot2 <- renderPlotly({
    req(inp2(), input$sys2_chr_col, input$sys2_bp_col, input$sys2_p_col, input$sys2_snp_col)
    p2 <- manhattanly(
      inp2(),
      chr = input$sys2_chr_col,
      bp = input$sys2_bp_col,
      p = input$sys2_p_col,
      snp = input$sys2_snp_col
    )
    # Assign a source ID to trace zoom/pan events
    p2$x$source <- "manhattanPlot2"
    p2
  })
  
  # --- SYNC ZOOM/PAN BETWEEN PLOTS ---
  
  # 1. Listen for Zoom/Pan on Plot 1 and update Plot 2
  observeEvent(event_data("plotly_relayout", source = "manhattanPlot1"), {
    relayout_data <- event_data("plotly_relayout", source = "manhattanPlot1")
    
    # Check if x-axis or y-axis bounds were modified
    if (!is.null(relayout_data)) {
      plotlyProxy("manhattanPlot2", session) %>%
        plotlyProxyInvoke("relayout", relayout_data)
    }
  }, ignoreInit = TRUE)
  
  # 2. Listen for Zoom/Pan on Plot 2 and update Plot 1
  observeEvent(event_data("plotly_relayout", source = "manhattanPlot2"), {
    relayout_data <- event_data("plotly_relayout", source = "manhattanPlot2")
    
    if (!is.null(relayout_data)) {
      plotlyProxy("manhattanPlot1", session) %>%
        plotlyProxyInvoke("relayout", relayout_data)
    }
  }, ignoreInit = TRUE)
  
  output$table1 <- renderDataTable({
    req(inp1())
    datatable(inp1())
  })
  
  output$table2 <- renderDataTable({
    req(inp2())
    datatable(inp2())
  })
} # server


# Create Shiny object
shinyApp(ui = ui, server = server)