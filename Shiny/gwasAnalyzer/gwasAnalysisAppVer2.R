# With help from tutorial by Data Professor

# Load R packages and datasets
library(shiny) # shiny
library(shinythemes) # themes
library(DT) # datatables
library(manhattanly)
library(data.table)
library(plotly)
library(colourpicker) # color input

# Set max upload size to 1GB
options(shiny.maxRequestSize = 1000 * 1024^2)

# Helper function to parse inputs like "1, 3, 5-7"
parse_chr_input <- function(input_str) {
  if (!nzchar(trimws(input_str))) return(NULL)
  
  parts <- unlist(strsplit(input_str, ","))
  parts <- trimws(parts)
  parts <- parts[parts != ""]
  
  parsed_chrs <- c()
  
  for (part in parts) {
    if (grepl("-", part)) {
      range_bounds <- unlist(strsplit(part, "-"))
      range_bounds <- trimws(range_bounds)
      num_start <- as.numeric(range_bounds[1])
      num_end <- as.numeric(range_bounds[2])
      
      if (!is.na(num_start) && !is.na(num_end)) {
        parsed_chrs <- c(parsed_chrs, as.character(seq(num_start, num_end)))
      } else {
        parsed_chrs <- c(parsed_chrs, part)
      }
    } else {
      parsed_chrs <- c(parsed_chrs, part)
    }
  }
  
  return(unique(parsed_chrs))
}

# Define UI
ui <- fluidPage(
  theme = shinytheme("cerulean"),
  navbarPage(
    "My first app",
    tabPanel("Navbar 1",
             sidebarPanel(
               tags$h2("Input:"),
               textInput("txt1", "First Name:", ""),
               textInput("txt2", "Last Name:", "")
             ),
             mainPanel(
               h1("Name"),
               h4("Full Name"),
               verbatimTextOutput("txtout")
             )
    ),
    tabPanel("Navbar 2",
             titlePanel("Interactive GWAS Manhattan Plot"),
             fluidRow(
               column(width = 6, fileInput("sys1_file", "Upload Sys 1 Dataset", accept = c(".tsv", ".logistic", ".txt", ".csv"))),
               column(width = 6, fileInput("sys2_file", "Upload Sys 2 Dataset", accept = c(".tsv", ".logistic", ".txt", ".csv")))
             ),
             
             hr(),
             
             # Row 2: Dynamic Column Mapping & Filters
             fluidRow(
               # System 1 Controls
               column(6,
                      wellPanel(
                        h4("System 1 Settings"),
                        fluidRow(
                          column(6, selectInput("sys1_chr_col", "Chromosome Column", choices = NULL)),
                          column(6, selectInput("sys1_bp_col", "Position (BP) Column", choices = NULL))
                        ),
                        fluidRow(
                          column(6, selectInput("sys1_p_col", "P-Value Column", choices = NULL)),
                          column(6, selectInput("sys1_snp_col", "SNP ID Column", choices = NULL))
                        ),
                        fluidRow(
                          column(6, textInput("sys1_chr_filter", "Chromosome Filter (e.g. 1, or blank)", value = "1")),
                          column(6, numericInput("sys1_p_thresh", "Max P-Value Threshold", value = 1e-5, step = 1e-6))
                        ),
                        fluidRow(
                          column(12, colourInput("sys1_col", "System 1 Point Color", value = "#1F77B4"))
                        )
                      )
               ),
               
               # System 2 Controls
               column(6,
                      wellPanel(
                        h4("System 2 Settings"),
                        fluidRow(
                          column(6, selectInput("sys2_chr_col", "Chromosome Column", choices = NULL)),
                          column(6, selectInput("sys2_bp_col", "Position (BP) Column", choices = NULL))
                        ),
                        fluidRow(
                          column(6, selectInput("sys2_p_col", "P-Value Column", choices = NULL)),
                          column(6, selectInput("sys2_snp_col", "SNP ID Column", choices = NULL))
                        ),
                        fluidRow(
                          column(6, textInput("sys2_chr_filter", "Chromosome Filter (e.g. 1, or blank)", value = "1")),
                          column(6, numericInput("sys2_p_thresh", "Max P-Value Threshold", value = 1e-5, step = 1e-6))
                        ),
                        fluidRow(
                          column(12, colourInput("sys2_col", "System 2 Point Color", value = "#FF7F0E"))
                        )
                      )
               )
             ),
             
             hr(),
             
             fluidRow(
               column(12,
                      wellPanel(
                        radioButtons("plot_mode", "Plot Display Mode:",
                                     choices = c("Side-by-Side" = "side", "Overlay Both Datasets" = "overlay"),
                                     selected = "side", inline = TRUE)
                      )
               )
             ),
             
             # Dynamic container handles both side-by-side and overlay modes
             uiOutput("plot_container")
    ),
    
    tabPanel("Navbar 3", 
             dataTableOutput("table1"),
             dataTableOutput("table2")
    )
  )
)

# Define server function   
server <- function(input, output, session) {
  
  output$txtout <- renderText({
    paste(input$txt1, input$txt2, sep = " ")
  })
  
  # --- SYSTEM 1 LOGIC ---
  
  sys1_raw <- reactive({
    req(input$sys1_file)
    fread(input$sys1_file$datapath)
  })
  
  observeEvent(sys1_raw(), {
    cols <- names(sys1_raw())
    
    chr_default <- grep("chr|chrom", cols, ignore.case = TRUE, value = TRUE)[1]
    bp_default  <- grep("pos|bp", cols, ignore.case = TRUE, value = TRUE)[1]
    p_default   <- grep("^p$|p_val|p.val|pval", cols, ignore.case = TRUE, value = TRUE)[1]
    snp_default <- grep("snp|id|rs", cols, ignore.case = TRUE, value = TRUE)[1]
    
    updateSelectInput(session, "sys1_chr_col", choices = cols, selected = ifelse(is.na(chr_default), cols[1], chr_default))
    updateSelectInput(session, "sys1_bp_col", choices = cols, selected = ifelse(is.na(bp_default), cols[1], bp_default))
    updateSelectInput(session, "sys1_p_col", choices = cols, selected = ifelse(is.na(p_default), cols[1], p_default))
    updateSelectInput(session, "sys1_snp_col", choices = cols, selected = ifelse(is.na(snp_default), cols[1], snp_default))
  })
  
  inp1 <- reactive({
    req(sys1_raw(), input$sys1_chr_col, input$sys1_p_col)
    
    df <- sys1_raw()
    chr_col <- input$sys1_chr_col
    p_col <- input$sys1_p_col
    
    target_chrs <- parse_chr_input(input$sys1_chr_filter)
    if (!is.null(target_chrs)) {
      df <- df[as.character(get(chr_col)) %in% target_chrs]
    }
    
    if (!is.null(input$sys1_p_thresh) && !is.na(input$sys1_p_thresh)) {
      df <- df[get(p_col) < as.numeric(input$sys1_p_thresh)]
    }
    
    df
  })
  
  # --- SYSTEM 2 LOGIC ---
  
  sys2_raw <- reactive({
    req(input$sys2_file)
    fread(input$sys2_file$datapath)
  })
  
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
  
  inp2 <- reactive({
    req(sys2_raw(), input$sys2_chr_col, input$sys2_p_col)
    
    df <- sys2_raw()
    chr_col <- input$sys2_chr_col
    p_col <- input$sys2_p_col
    
    target_chrs <- parse_chr_input(input$sys2_chr_filter)
    if (!is.null(target_chrs)) {
      df <- df[as.character(get(chr_col)) %in% target_chrs]
    }
    
    if (!is.null(input$sys2_p_thresh) && !is.na(input$sys2_p_thresh)) {
      df <- df[get(p_col) < as.numeric(input$sys2_p_thresh)]
    }
    
    df
  })
  
  # Dynamic UI Layout for Plot Mode
  output$plot_container <- renderUI({
    if (input$plot_mode == "side") {
      fluidRow(
        column(6, plotlyOutput("manhattanPlot1", height = "600px")),
        column(6, plotlyOutput("manhattanPlot2", height = "600px"))
      )
    } else {
      fluidRow(
        column(12, plotlyOutput("overlayPlot", height = "650px"))
      )
    }
  })
  
  # --- SYSTEM 1 INDIVIDUAL PLOT ---
  output$manhattanPlot1 <- renderPlotly({
    req(inp1(), input$sys1_chr_col, input$sys1_bp_col, input$sys1_p_col, input$sys1_snp_col)
    
    p1 <- manhattanly(
      inp1(),
      chr = input$sys1_chr_col,
      bp = input$sys1_bp_col,
      p = input$sys1_p_col,
      snp = input$sys1_snp_col,
      col = c(input$sys1_col, input$sys1_col)
    )
    p1$x$source <- "manhattanPlot1"
    p1
  })
  
  # --- SYSTEM 2 INDIVIDUAL PLOT ---
  output$manhattanPlot2 <- renderPlotly({
    req(inp2(), input$sys2_chr_col, input$sys2_bp_col, input$sys2_p_col, input$sys2_snp_col)
    
    p2 <- manhattanly(
      inp2(),
      chr = input$sys2_chr_col,
      bp = input$sys2_bp_col,
      p = input$sys2_p_col,
      snp = input$sys2_snp_col,
      col = c(input$sys2_col, input$sys2_col)
    )
    p2$x$source <- "manhattanPlot2"
    p2
  })
  
  # --- OVERLAY PLOT LOGIC ---
  output$overlayPlot <- renderPlotly({
    req(inp1(), inp2(), input$sys1_chr_col, input$sys1_bp_col, input$sys1_p_col, input$sys1_snp_col)
    req(input$sys2_chr_col, input$sys2_bp_col, input$sys2_p_col, input$sys2_snp_col)
    
    # Copy datasets and create uniform column names for plotting
    d1 <- copy(inp1())
    d1_chr <- input$sys1_chr_col
    d1_bp  <- input$sys1_bp_col
    d1_p   <- input$sys1_p_col
    d1_snp <- input$sys1_snp_col
    
    d1[, `:=`(
      CHR_plot = as.character(get(d1_chr)),
      BP_plot  = as.numeric(get(d1_bp)),
      P_plot   = -log10(as.numeric(get(d1_p))),
      SNP_plot = as.character(get(d1_snp))
    )]
    
    d2 <- copy(inp2())
    d2_chr <- input$sys2_chr_col
    d2_bp  <- input$sys2_bp_col
    d2_p   <- input$sys2_p_col
    d2_snp <- input$sys2_snp_col
    
    d2[, `:=`(
      CHR_plot = as.character(get(d2_chr)),
      BP_plot  = as.numeric(get(d2_bp)),
      P_plot   = -log10(as.numeric(get(d2_p))),
      SNP_plot = as.character(get(d2_snp))
    )]
    
    # Construct overlay plot with separate traces for each system
    plot_ly() %>%
      add_trace(
        data = d1,
        x = ~BP_plot,
        y = ~P_plot,
        type = 'scatter',
        mode = 'markers',
        name = 'System 1',
        marker = list(color = input$sys1_col, size = 6, opacity = 0.7),
        text = ~paste("SNP:", SNP_plot, "<br>CHR:", CHR_plot, "<br>BP:", BP_plot, "<br>-log10(P):", round(P_plot, 3)),
        hoverinfo = "text"
      ) %>%
      add_trace(
        data = d2,
        x = ~BP_plot,
        y = ~P_plot,
        type = 'scatter',
        mode = 'markers',
        name = 'System 2',
        marker = list(color = input$sys2_col, size = 6, opacity = 0.7),
        text = ~paste("SNP:", SNP_plot, "<br>CHR:", CHR_plot, "<br>BP:", BP_plot, "<br>-log10(P):", round(P_plot, 3)),
        hoverinfo = "text"
      ) %>%
      layout(
        title = "Overlay Manhattan Plot",
        xaxis = list(title = "Base Pair Position (BP)"),
        yaxis = list(title = "-log10(p-value)"),
        legend = list(title = list(text = '<b>Dataset</b>'))
      )
  })
  
  # --- SYNC ZOOM/PAN BETWEEN PLOTS ---
  observeEvent(event_data("plotly_relayout", source = "manhattanPlot1"), {
    relayout_data <- event_data("plotly_relayout", source = "manhattanPlot1")
    if (!is.null(relayout_data)) {
      plotlyProxy("manhattanPlot2", session) %>%
        plotlyProxyInvoke("relayout", relayout_data)
    }
  }, ignoreInit = TRUE)
  
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
}

# Create Shiny object
shinyApp(ui = ui, server = server)