# Load R packages
library(shiny)
library(shinythemes)
library(DT)
library(manhattanly)
library(data.table)
library(plotly)
library(colourpicker)
library(httr2)
library(jsonlite)

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

# gnomad url helper
build_gnomad_url <- function(query, dataset) {
  q <- trimws(query)
  if (!nzchar(q)) return(NULL)
  
  # Strip "chr" prefix if present
  q_clean <- gsub("^chr", "", q, ignore.case = TRUE)
  
  # 1. Ensembl Gene ID (e.g., ENSG00000012048)
  if (grepl("^ENSG[0-9]+$", q, ignore.case = TRUE)) {
    return(paste0("https://gnomad.broadinstitute.org/gene/", toupper(q), "?dataset=", dataset))
  }
  
  # 2. Variant ID (e.g., 1-55516888-G-GA, 1:55516888:G:GA)
  var_match <- regmatches(q_clean, regexec("^([0-9]{1,2}|X|Y|MT|M)[-:_]([0-9]+)[-:_]([A-Za-z]+)[-:_]([A-Za-z]+)$", q_clean, ignore.case = TRUE))[[1]]
  if (length(var_match) == 5) {
    formatted_var <- paste(toupper(var_match[2]), var_match[3], toupper(var_match[4]), toupper(var_match[5]), sep = "-")
    return(paste0("https://gnomad.broadinstitute.org/variant/", formatted_var, "?dataset=", dataset))
  }
  
  # 3. Genomic Region / Locus (e.g., 1:55516888-55520000)
  reg_match <- regmatches(q_clean, regexec("^([0-9]{1,2}|X|Y|MT|M)[-:_]([0-9]+)[-:_]([0-9]+)$", q_clean, ignore.case = TRUE))[[1]]
  if (length(reg_match) == 4) {
    formatted_reg <- paste(toupper(reg_match[2]), reg_match[3], reg_match[4], sep = "-")
    return(paste0("https://gnomad.broadinstitute.org/region/", formatted_reg, "?dataset=", dataset))
  }
  
  # 4. Universal Search Endpoint
  return(paste0("https://gnomad.broadinstitute.org/search?dataset=", dataset, "&q=", URLencode(q, reserved = TRUE)))
}

# Helper to extract raw Variant ID from hover text/customdata
clean_variant_id <- function(raw_input) {
  if (is.null(raw_input) || length(raw_input) == 0) return(NULL)
  raw_str <- as.character(raw_input[1])
  
  if (grepl("SNP:", raw_str)) {
    match <- regmatches(raw_str, regexec("SNP:\\s*([^<\\s]+)", raw_str))[[1]]
    if (length(match) >= 2) return(trimws(match[2]))
  }
  
  var_match <- regmatches(raw_str, regexec("([0-9]{1,2}|X|Y|MT|M)[-:_]([0-9]+)[-:_]([A-Za-z]+)[-:_]([A-Za-z]+)", raw_str, ignore.case = TRUE))[[1]]
  if (length(var_match) >= 1) return(var_match[1])
  
  return(trimws(raw_str))
}

# Helper to fetch Live Allele Frequency from gnoMAD GraphQL API
fetch_gnomad_af <- function(variant_id, dataset = "gnomad_r4") {
  cleaned_id <- clean_variant_id(variant_id)
  if (is.null(cleaned_id) || !nzchar(cleaned_id)) return("N/A")
  
  no_chr_id <- gsub("^chr", "", cleaned_id, ignore.case = TRUE)
  formatted_id <- gsub("[:_]", "-", no_chr_id)
  
  is_canonical_var <- grepl("^([0-9]{1,2}|X|Y|MT|M)-[0-9]+-[A-Za-z]+-[A-Za-z]+$", formatted_id, ignore.case = TRUE)
  is_rsid <- grepl("^rs[0-9]+$", formatted_id, ignore.case = TRUE)
  
  # Construct GraphQL query based on variant format type
  if (is_rsid) {
    query_string <- sprintf('{
      rsid(rsid: "%s") {
        variants(dataset: %s) {
          variant_id
          genome { af }
          exome { af }
          joint { af }
        }
      }
    }', formatted_id, dataset)
  } else if (is_canonical_var) {
    query_string <- sprintf('{
      variant(variant_id: "%s", dataset: %s) {
        genome { af }
        exome { af }
        joint { af }
      }
    }', formatted_id, dataset)
  } else {
    message("gnoMAD Fetch Skipped: ID '", formatted_id, "' is not in CHROM-POS-REF-ALT or rsID format.")
    return("No REF/ALT Alleles")
  }
  
  tryCatch({
    resp <- request("https://gnomad.broadinstitute.org/api") %>%
      req_headers("Content-Type" = "application/json") %>%
      req_body_json(list(query = query_string)) %>%
      req_timeout(10) %>%
      req_error(is_error = function(resp) FALSE) %>%
      req_perform()
    
    if (resp_status(resp) >= 400) {
      message("gnoMAD API HTTP Error Code: ", resp_status(resp))
      return("API Error")
    }
    
    res_data <- resp_body_json(resp)
    
    if (!is.null(res_data$errors)) {
      message("gnoMAD GraphQL returned error: ", jsonlite::toJSON(res_data$errors))
      return("Not Found")
    }
    
    # Extract variant object
    if (is_rsid) {
      var_list <- res_data$data$rsid$variants
      if (is.null(var_list) || length(var_list) == 0) return("Not Found")
      var_res <- var_list[[1]]
    } else {
      var_res <- res_data$data$variant
    }
    
    if (is.null(var_res)) return("Not Found")
    
    af_val <- NULL
    if (!is.null(var_res$joint) && !is.null(var_res$joint$af)) {
      af_val <- var_res$joint$af
    } else if (!is.null(var_res$genome) && !is.null(var_res$genome$af)) {
      af_val <- var_res$genome$af
    } else if (!is.null(var_res$exome) && !is.null(var_res$exome$af)) {
      af_val <- var_res$exome$af
    }
    
    if (is.null(af_val)) return("Not Found")
    return(formatC(as.numeric(af_val), format = "e", digits = 3))
    
  }, error = function(e) {
    message("gnoMAD Fetch Exception: ", e$message)
    return("Error Fetching")
  })
}

# Store clicked annotations for each plot
clicked_point_sys1 <- reactiveVal(NULL)
clicked_point_sys2 <- reactiveVal(NULL)
clicked_point_overlay <- reactiveVal(NULL)

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
                                     choices = c("Individual Plots" = "individual", "Overlay Both Datasets" = "overlay"),
                                     selected = "individual", inline = TRUE),
                        
                        conditionalPanel(
                          condition = "input.plot_mode == 'individual'",
                          radioButtons("individual_layout", "Individual Plot Layout:",
                                       choices = c("Side-by-Side" = "side", "Top-to-Bottom" = "stacked"),
                                       selected = "side", inline = TRUE)
                        )
                      )
               )
             ),
             
             uiOutput("plot_container"),
             
             hr(),
             
             fluidRow(
               column(12,
                      wellPanel(
                        h4("gnoMAD Gene & Variant Direct Link"),
                        fluidRow(
                          column(4, textInput("gnomad_query", "Enter Gene Symbol, RSID, or Locus (e.g., PCSK9 or rs1234)", value = "")),
                          column(3, selectInput("gnomad_build", "Genome Build", choices = c("GRCh38" = "gnomad_r4", "GRCh37" = "gnomad_r2_1"))),
                          column(5, 
                                 br(),
                                 uiOutput("gnomad_link_button")
                          )
                        )
                      )
               )
             )
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
  
  output$plot_container <- renderUI({
    if (input$plot_mode == "individual") {
      if (input$individual_layout == "side") {
        fluidRow(
          column(6, plotlyOutput("manhattanPlot1", height = "600px")),
          column(6, plotlyOutput("manhattanPlot2", height = "600px"))
        )
      } else {
        fluidRow(
          column(12, plotlyOutput("manhattanPlot1", height = "500px")),
          column(12, br()),
          column(12, plotlyOutput("manhattanPlot2", height = "500px"))
        )
      }
    } else {
      fluidRow(
        column(12, plotlyOutput("overlayPlot", height = "650px"))
      )
    }
  })
  
  output$manhattanPlot1 <- renderPlotly({
    req(inp1(), input$sys1_chr_col, input$sys1_bp_col, input$sys1_p_col, input$sys1_snp_col)
    
    p1 <- manhattanly(
      inp1(),
      chr = input$sys1_chr_col,
      bp = input$sys1_bp_col,
      p = input$sys1_p_col,
      snp = input$sys1_snp_col,
      annotation1 = input$sys1_snp_col,
      annotation2 = input$sys1_bp_col,
      col = c(input$sys1_col, input$sys1_col)
    )
    p1$x$source <- "manhattanPlot1"
    layout(p1, uirevision = "manhattanPlot1_state")
  })
  
  output$manhattanPlot2 <- renderPlotly({
    req(inp2(), input$sys2_chr_col, input$sys2_bp_col, input$sys2_p_col, input$sys2_snp_col)
    
    p2 <- manhattanly(
      inp2(),
      chr = input$sys2_chr_col,
      bp = input$sys2_bp_col,
      p = input$sys2_p_col,
      snp = input$sys2_snp_col,
      annotation1 = input$sys2_snp_col,
      annotation2 = input$sys2_bp_col,
      col = c(input$sys2_col, input$sys2_col)
    )
    p2$x$source <- "manhattanPlot2"
    layout(p2, uirevision = "manhattanPlot2_state")
  })
  
  output$overlayPlot <- renderPlotly({
    req(inp1(), inp2(), input$sys1_chr_col, input$sys1_bp_col, input$sys1_p_col, input$sys1_snp_col)
    req(input$sys2_chr_col, input$sys2_bp_col, input$sys2_p_col, input$sys2_snp_col)
    
    d1 <- copy(inp1())
    d1[, `:=`(
      CHR_plot = as.character(get(input$sys1_chr_col)),
      BP_plot  = as.numeric(get(input$sys1_bp_col)),
      P_plot   = -log10(as.numeric(get(input$sys1_p_col))),
      SNP_plot = as.character(get(input$sys1_snp_col))
    )]
    
    d2 <- copy(inp2())
    d2[, `:=`(
      CHR_plot = as.character(get(input$sys2_chr_col)),
      BP_plot  = as.numeric(get(input$sys2_bp_col)),
      P_plot   = -log10(as.numeric(get(input$sys2_p_col))),
      SNP_plot = as.character(get(input$sys2_snp_col))
    )]
    
    p_overlay <- plot_ly() %>%
      add_trace(
        data = d1,
        x = ~BP_plot,
        y = ~P_plot,
        type = 'scatter',
        mode = 'markers',
        name = 'System 1',
        marker = list(color = input$sys1_col, size = 6, opacity = 0.7),
        customdata = ~SNP_plot,
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
        customdata = ~SNP_plot,
        text = ~paste("SNP:", SNP_plot, "<br>CHR:", CHR_plot, "<br>BP:", BP_plot, "<br>-log10(P):", round(P_plot, 3)),
        hoverinfo = "text"
      ) %>%
      layout(
        title = "Overlay Manhattan Plot",
        xaxis = list(title = "Base Pair Position (BP)"),
        yaxis = list(title = "-log10(p-value)"),
        legend = list(title = list(text = '<b>Dataset</b>')),
        uirevision = "overlayPlot_state"
      )
    
    p_overlay$x$source <- "overlayPlot"
    p_overlay
  })
  
  output$gnomad_link_button <- renderUI({
    req(input$gnomad_query)
    query_trimmed <- trimws(input$gnomad_query)
    if (!nzchar(query_trimmed)) return(NULL)
    
    target_url <- build_gnomad_url(query_trimmed, input$gnomad_build)
    if (is.null(target_url)) return(NULL)
    
    tags$a(
      href = target_url,
      target = "_blank",
      class = "btn btn-primary",
      icon("external-link-alt"),
      paste("Open", query_trimmed, "in gnoMAD")
    )
  })
  
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
  
  # Click Observers
  observeEvent(event_data("plotly_click", source = "manhattanPlot1"), {
    click_data <- event_data("plotly_click", source = "manhattanPlot1")
    if (!is.null(click_data)) {
      current <- clicked_point_sys1()
      
      if (!is.null(current) && identical(current$x, click_data$x) && identical(current$y, click_data$y)) {
        clicked_point_sys1(NULL)
        plotlyProxy("manhattanPlot1", session) %>%
          plotlyProxyInvoke("relayout", list(annotations = list()))
      } else {
        raw_val <- if (!is.null(click_data$customdata)) click_data$customdata else click_data$text
        snp_id <- clean_variant_id(raw_val)
        if (is.null(snp_id) || !nzchar(snp_id)) snp_id <- paste0(click_data$x)
        
        af_val <- fetch_gnomad_af(snp_id, input$gnomad_build)
        
        hover_txt <- paste0(
          "<b>SNP:</b> ", snp_id, "<br>",
          "<b>BP:</b> ", click_data$x, "<br>",
          "<b>-log10(P):</b> ", round(click_data$y, 3), "<br>",
          "<b>gnoMAD AF:</b> ", af_val
        )
        
        annotation <- list(
          x = click_data$x,
          y = click_data$y,
          text = hover_txt,
          showarrow = FALSE,
          xanchor = "left",
          yanchor = "bottom",
          bgcolor = "rgba(255, 255, 255, 0.95)",
          bordercolor = "#444444",
          borderwidth = 1,
          borderpad = 6,
          font = list(size = 12, color = "#000000")
        )
        clicked_point_sys1(annotation)
        plotlyProxy("manhattanPlot1", session) %>%
          plotlyProxyInvoke("relayout", list(annotations = list(annotation)))
      }
    }
  })
  
  observeEvent(event_data("plotly_click", source = "manhattanPlot2"), {
    click_data <- event_data("plotly_click", source = "manhattanPlot2")
    if (!is.null(click_data)) {
      current <- clicked_point_sys2()
      
      if (!is.null(current) && identical(current$x, click_data$x) && identical(current$y, click_data$y)) {
        clicked_point_sys2(NULL)
        plotlyProxy("manhattanPlot2", session) %>%
          plotlyProxyInvoke("relayout", list(annotations = list()))
      } else {
        raw_val <- if (!is.null(click_data$customdata)) click_data$customdata else click_data$text
        snp_id <- clean_variant_id(raw_val)
        if (is.null(snp_id) || !nzchar(snp_id)) snp_id <- paste0(click_data$x)
        
        af_val <- fetch_gnomad_af(snp_id, input$gnomad_build)
        
        hover_txt <- paste0(
          "<b>SNP:</b> ", snp_id, "<br>",
          "<b>BP:</b> ", click_data$x, "<br>",
          "<b>-log10(P):</b> ", round(click_data$y, 3), "<br>",
          "<b>gnoMAD AF:</b> ", af_val
        )
        
        annotation <- list(
          x = click_data$x,
          y = click_data$y,
          text = hover_txt,
          showarrow = FALSE,
          xanchor = "left",
          yanchor = "bottom",
          bgcolor = "rgba(255, 255, 255, 0.95)",
          bordercolor = "#444444",
          borderwidth = 1,
          borderpad = 6,
          font = list(size = 12, color = "#000000")
        )
        clicked_point_sys2(annotation)
        plotlyProxy("manhattanPlot2", session) %>%
          plotlyProxyInvoke("relayout", list(annotations = list(annotation)))
      }
    }
  })
  
  observeEvent(event_data("plotly_click", source = "overlayPlot"), {
    click_data <- event_data("plotly_click", source = "overlayPlot")
    if (!is.null(click_data)) {
      current <- clicked_point_overlay()
      
      if (!is.null(current) && identical(current$x, click_data$x) && identical(current$y, click_data$y)) {
        clicked_point_overlay(NULL)
        plotlyProxy("overlayPlot", session) %>%
          plotlyProxyInvoke("relayout", list(annotations = list()))
      } else {
        raw_val <- if (!is.null(click_data$customdata)) click_data$customdata else click_data$text
        snp_id <- clean_variant_id(raw_val)
        if (is.null(snp_id) || !nzchar(snp_id)) snp_id <- paste0(click_data$x)
        
        af_val <- fetch_gnomad_af(snp_id, input$gnomad_build)
        
        hover_txt <- paste0(
          "<b>SNP:</b> ", snp_id, "<br>",
          "<b>BP:</b> ", click_data$x, "<br>",
          "<b>-log10(P):</b> ", round(click_data$y, 3), "<br>",
          "<b>gnoMAD AF:</b> ", af_val
        )
        
        annotation <- list(
          x = click_data$x,
          y = click_data$y,
          text = hover_txt,
          showarrow = FALSE,
          xanchor = "left",
          yanchor = "bottom",
          bgcolor = "rgba(255, 255, 255, 0.95)",
          bordercolor = "#444444",
          borderwidth = 1,
          borderpad = 6,
          font = list(size = 12, color = "#000000")
        )
        clicked_point_overlay(annotation)
        plotlyProxy("overlayPlot", session) %>%
          plotlyProxyInvoke("relayout", list(annotations = list(annotation)))
      }
    }
  })
  
  output$table1 <- renderDataTable({
    req(inp1())
    datatable(inp1())
  })
  
  output$table2 <- renderDataTable({
    req(inp2())
    datatable(inp2())
  })
}

shinyApp(ui = ui, server = server)