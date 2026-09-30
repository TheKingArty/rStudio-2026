# Load required R packages
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

# Helper function to parse chromosome range inputs like "1, 3, 5-7"
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

# gnoMAD URL builder helper
build_gnomad_url <- function(query, dataset) {
  q <- trimws(query)
  if (!nzchar(q)) return(NULL)
  
  q_clean <- gsub("^chr", "", q, ignore.case = TRUE)
  
  # 1. Ensembl Gene ID
  if (grepl("^ENSG[0-9]+$", q, ignore.case = TRUE)) {
    return(paste0("https://gnomad.broadinstitute.org/gene/", toupper(q), "?dataset=", dataset))
  }
  
  # 2. Variant ID
  var_match <- regmatches(q_clean, regexec("^([0-9]{1,2}|X|Y|MT|M)[-:_]([0-9]+)[-:_]([A-Za-z0-9*]+)[-:_]([A-Za-z0-9*]+)$", q_clean, ignore.case = TRUE))[[1]]
  if (length(var_match) == 5) {
    formatted_var <- paste(toupper(var_match[2]), var_match[3], toupper(var_match[4]), toupper(var_match[5]), sep = "-")
    return(paste0("https://gnomad.broadinstitute.org/variant/", formatted_var, "?dataset=", dataset))
  }
  
  # 3. Genomic Region / Locus
  reg_match <- regmatches(q_clean, regexec("^([0-9]{1,2}|X|Y|MT|M)[-:_]([0-9]+)[-:_]([0-9]+)$", q_clean, ignore.case = TRUE))[[1]]
  if (length(reg_match) == 4) {
    formatted_reg <- paste(toupper(reg_match[2]), reg_match[3], reg_match[4], sep = "-")
    return(paste0("https://gnomad.broadinstitute.org/region/", formatted_reg, "?dataset=", dataset))
  }
  
  # 4. Universal Search Endpoint
  return(paste0("https://gnomad.broadinstitute.org/search?dataset=", dataset, "&q=", URLencode(q, reserved = TRUE)))
}

clean_variant_id <- function(raw_input, chr_fallback = NULL, bp_fallback = NULL) {
  if (is.null(raw_input) || length(raw_input) == 0) {
    if (!is.null(chr_fallback) && !is.null(bp_fallback)) {
      c_val <- toupper(gsub("^chr", "", as.character(chr_fallback), ignore.case = TRUE))
      b_val <- as.character(bp_fallback)
      if (nzchar(c_val) && nzchar(b_val)) return(paste(c_val, b_val, sep = "-"))
    }
    return(NULL)
  }
  
  raw_str <- unlist(raw_input)[1]
  if (is.na(raw_str) || !nzchar(raw_str)) {
    if (!is.null(chr_fallback) && !is.null(bp_fallback)) {
      c_val <- toupper(gsub("^chr", "", as.character(chr_fallback), ignore.case = TRUE))
      b_val <- as.character(bp_fallback)
      if (nzchar(c_val) && nzchar(b_val)) return(paste(c_val, b_val, sep = "-"))
    }
    return(NULL)
  }
  
  raw_str <- as.character(raw_str)
  
  # Strip HTML tags & newlines
  clean_str <- gsub("<[^>]+>", " ", raw_str)
  clean_str <- gsub("[\r\n]", " ", clean_str)
  
  # 1. Try matching full 4-part CHROM-POS-REF-ALT
  m_var <- regmatches(clean_str, regexec("(?:chr)?([0-9]{1,2}|X|Y|MT|M)[-:_]([0-9]+)[-:_]([A-Za-z0-9*<>-]+)[-:_/]([A-Za-z0-9*<>-]+)", clean_str, ignore.case = TRUE))[[1]]
  if (length(m_var) == 5 && nzchar(m_var[1])) {
    return(toupper(paste(m_var[2], m_var[3], m_var[4], m_var[5], sep = "-")))
  }
  
  # 2. Try matching rsID
  m_rs <- regmatches(clean_str, regexpr("rs[0-9]+", clean_str, ignore.case = TRUE))
  if (length(m_rs) > 0 && nzchar(m_rs[1])) {
    return(tolower(m_rs[1]))
  }
  
  # 3. Try matching 2-part Position-Only (CHROM-POS)
  m_pos <- regmatches(clean_str, regexec("(?:chr)?([0-9]{1,2}|X|Y|MT|M)[-:_]([0-9]+)", clean_str, ignore.case = TRUE))[[1]]
  if (length(m_pos) == 3 && nzchar(m_pos[1])) {
    return(toupper(paste(m_pos[2], m_pos[3], sep = "-")))
  }
  
  # 4. Fallback to passed chromosome & BP if available
  if (!is.null(chr_fallback) && !is.null(bp_fallback)) {
    c_val <- toupper(gsub("^chr", "", as.character(chr_fallback), ignore.case = TRUE))
    b_val <- as.character(bp_fallback)
    if (nzchar(c_val) && nzchar(b_val) && !is.na(c_val) && !is.na(b_val)) {
      return(paste(c_val, b_val, sep = "-"))
    }
  }
  
  trimmed <- trimws(clean_str)
  if (!grepl("\\s", trimmed) && nzchar(trimmed)) {
    return(toupper(trimmed))
  }
  
  return(NULL)
}

# Robust gnoMAD AF Fetcher with Header & Fallback Fixes
fetch_gnomad_af <- function(variant_id, dataset = "gnomad_r4") {
  cleaned_id <- clean_variant_id(variant_id)
  if (is.null(cleaned_id) || !nzchar(cleaned_id)) return("N/A")
  
  no_chr_id <- gsub("^chr", "", cleaned_id, ignore.case = TRUE)
  formatted_id <- toupper(gsub("[:_]", "-", no_chr_id))
  if (startsWith(formatted_id, "RS")) formatted_id <- tolower(formatted_id)
  
  is_rsid          <- grepl("^rs[0-9]+$", formatted_id, ignore.case = TRUE)
  is_canonical_var <- grepl("^([0-9]{1,2}|X|Y|MT|M)-[0-9]+-[^-]+-[^-]+$", formatted_id, ignore.case = TRUE)
  is_pos_only      <- grepl("^([0-9]{1,2}|X|Y|MT|M)-[0-9]+$", formatted_id, ignore.case = TRUE)
  
  # Conditionally request joint AF only for v4 dataset
  af_fields <- if (identical(dataset, "gnomad_r4")) {
    "genome { af } exome { af } joint { af }"
  } else {
    "genome { af } exome { af }"
  }
  
  if (is_rsid) {
    query_string <- sprintf('{
      rsid(rsid: "%s") {
        variants(dataset: %s) {
          variant_id
          %s
        }
      }
    }', formatted_id, dataset, af_fields)
    
  } else if (is_canonical_var) {
    query_string <- sprintf('{
      variant(variant_id: "%s", dataset: %s) {
        %s
      }
    }', formatted_id, dataset, af_fields)
    
  } else if (is_pos_only) {
    parts <- unlist(strsplit(formatted_id, "-"))
    chr_val <- gsub("^chr", "", parts[1], ignore.case = TRUE)
    pos_val <- as.numeric(parts[2])
    
    ref_genome <- if (grepl("r2", dataset)) "GRCh37" else "GRCh38"
    query_string <- sprintf('{
      region(chrom: "%s", start: %d, stop: %d, reference_genome: %s) {
        variants(dataset: %s) {
          variant_id
          %s
        }
      }
    }', chr_val, pos_val, pos_val, ref_genome, dataset, af_fields)
    
  } else {
    message("gnoMAD Fetch Skipped: Could not parse ID '", formatted_id, "'")
    return("Invalid Format")
  }
  
  tryCatch({
    # Primary Request
    resp <- request("https://gnomad.broadinstitute.org/api") %>%
      req_user_agent("Mozilla/5.0 (Windows NT 10.0; Win64; x64) R-Shiny-App") %>%
      req_headers("Content-Type" = "application/json") %>%
      req_body_json(list(query = query_string)) %>%
      req_timeout(10) %>%
      req_error(is_error = function(resp) FALSE) %>%
      req_perform()
    
    # If joint query failed or returned HTTP status >= 400, retry without joint field
    if (resp_status(resp) >= 400) {
      err_msg <- resp_body_string(resp)
      message("gnoMAD API HTTP Error ", resp_status(resp), ": ", err_msg)
      
      if (grepl("joint", query_string, ignore.case = TRUE)) {
        fallback_query <- gsub(" joint \\{ af \\}", "", query_string)
        resp_fb <- request("https://gnomad.broadinstitute.org/api") %>%
          req_user_agent("Mozilla/5.0 (Windows NT 10.0; Win64; x64) R-Shiny-App") %>%
          req_headers("Content-Type" = "application/json") %>%
          req_body_json(list(query = fallback_query)) %>%
          req_timeout(10) %>%
          req_error(is_error = function(resp) FALSE) %>%
          req_perform()
        
        if (resp_status(resp_fb) < 400) {
          resp <- resp_fb
        } else {
          return("API Error")
        }
      } else {
        return("API Error")
      }
    }
    
    res_data <- resp_body_json(resp)
    
    if (!is.null(res_data$errors)) {
      message("gnoMAD GraphQL error: ", jsonlite::toJSON(res_data$errors))
      return("Not Found")
    }
    
    if (is.null(res_data$data)) return("Not Found")
    
    var_res <- NULL
    if (is_rsid) {
      var_list <- res_data$data$rsid$variants
      if (!is.null(var_list) && length(var_list) > 0) var_res <- var_list[[1]]
    } else if (is_canonical_var) {
      var_res <- res_data$data$variant
    } else if (is_pos_only) {
      var_list <- res_data$data$region$variants
      if (!is.null(var_list) && length(var_list) > 0) var_res <- var_list[[1]]
    }
    
    if (is.null(var_res)) return("Not Found")
    
    af_val <- NULL
    if (!is.null(var_res$joint$af)) {
      af_val <- var_res$joint$af
    } else if (!is.null(var_res$genome$af)) {
      af_val <- var_res$genome$af
    } else if (!is.null(var_res$exome$af)) {
      af_val <- var_res$exome$af
    }
    
    if (is.null(af_val) || is.na(af_val)) return("Not Found")
    return(formatC(as.numeric(af_val), format = "e", digits = 3))
    
  }, error = function(e) {
    message("gnoMAD Fetch Exception: ", e$message)
    return("Error Fetching")
  })
}

# Store clicked annotations
clicked_point_sys1 <- reactiveVal(NULL)
clicked_point_sys2 <- reactiveVal(NULL)
clicked_point_overlay <- reactiveVal(NULL)

# UI Definition
ui <- fluidPage(
  theme = shinytheme("flatly"),
  navbarPage(
    "GAnalyzer",
    id = "navbar",
    tabPanel("Name",
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
    tabPanel("Plot",
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
                        
                        checkboxInput("sys1_has_ref_alt", "Dataset includes REF and ALT columns", value = FALSE),
                        conditionalPanel(
                          condition = "input.sys1_has_ref_alt == true",
                          fluidRow(
                            column(6, selectInput("sys1_ref_col", "REF Allele Column", choices = NULL)),
                            column(6, selectInput("sys1_alt_col", "ALT Allele Column", choices = NULL))
                          )
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
                        
                        checkboxInput("sys2_has_ref_alt", "Dataset includes REF and ALT columns", value = FALSE),
                        conditionalPanel(
                          condition = "input.sys2_has_ref_alt == true",
                          fluidRow(
                            column(6, selectInput("sys2_ref_col", "REF Allele Column", choices = NULL)),
                            column(6, selectInput("sys2_alt_col", "ALT Allele Column", choices = NULL))
                          )
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
                        fluidRow(
                          column(8,
                                 radioButtons("plot_mode", "Plot Display Mode:",
                                              choices = c("Individual Plots" = "individual", "Overlay Both Datasets" = "overlay"),
                                              selected = "individual", inline = TRUE),
                                 
                                 conditionalPanel(
                                   condition = "input.plot_mode == 'individual'",
                                   radioButtons("individual_layout", "Individual Plot Layout:",
                                                choices = c("Side-by-Side" = "side", "Top-to-Bottom" = "stacked"),
                                                selected = "side", inline = TRUE)
                                 )
                          ),
                          column(4, class = "text-right",
                                 actionButton("reset_plots", "Reset Plot View & Axes", 
                                              icon = icon("rotate-left"), 
                                              class = "btn-warning", 
                                              style = "margin-top: 15px;")
                          )
                        )
                      )
               )
             ),
             
             fluidRow(
               column(12,
                      actionButton("go_to_nav3", "Download data", class = "btn-primary", style = "margin-bottom: 15px;")
               )
             ),
             
             uiOutput("plot_container"),
             
             hr(),
             
             # Region Data Preview Tables (< 50 points threshold)
             fluidRow(
               column(6,
                      wellPanel(
                        h4("System 1 Region Preview"),
                        uiOutput("sys1_preview_ui")
                      )
               ),
               column(6,
                      wellPanel(
                        h4("System 2 Region Preview"),
                        uiOutput("sys2_preview_ui")
                      )
               )
             ),
             
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
    
    tabPanel("Table", 
             h3("Full Datasets"),
             dataTableOutput("table1"),
             dataTableOutput("table2"),
             
             hr(),
             
             h3("Region-Filtered Datasets"),
             downloadButton("download_tables", "Download tables", class = "btn-success"),
             br(), br(),
             dataTableOutput("filtered_table1"),
             dataTableOutput("filtered_table2")
    )
  )
)

# Server Logic
server <- function(input, output, session) {
  # Reactive counter to force plot axis resets
  reset_counter <- reactiveVal(0)
  
  # Reset Button Handler
  observeEvent(input$reset_plots, {
    # 1. Clear clicked annotations
    clicked_point_sys1(NULL)
    clicked_point_sys2(NULL)
    clicked_point_overlay(NULL)
    
    # 2. Increment counter (forces Plotly to auto-scale axes)
    reset_counter(reset_counter() + 1)
    
    # 3. Clear annotations via proxy
    plotlyProxy("manhattanPlot1", session) %>% plotlyProxyInvoke("relayout", list(annotations = list()))
    plotlyProxy("manhattanPlot2", session) %>% plotlyProxyInvoke("relayout", list(annotations = list()))
    plotlyProxy("overlayPlot", session) %>% plotlyProxyInvoke("relayout", list(annotations = list()))
  })
  # Store the last relayout state to prevent infinite zooming loops
  last_relayout <- reactiveVal(NULL)
  
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
    ref_default <- grep("^ref$|reference|a1|allele1", cols, ignore.case = TRUE, value = TRUE)[1]
    alt_default <- grep("^alt$|alternate|a2|allele2", cols, ignore.case = TRUE, value = TRUE)[1]
    
    updateSelectInput(session, "sys1_chr_col", choices = cols, selected = ifelse(is.na(chr_default), cols[1], chr_default))
    updateSelectInput(session, "sys1_bp_col", choices = cols, selected = ifelse(is.na(bp_default), cols[1], bp_default))
    updateSelectInput(session, "sys1_p_col", choices = cols, selected = ifelse(is.na(p_default), cols[1], p_default))
    updateSelectInput(session, "sys1_snp_col", choices = cols, selected = ifelse(is.na(snp_default), cols[1], snp_default))
    updateSelectInput(session, "sys1_ref_col", choices = cols, selected = ifelse(is.na(ref_default), cols[1], ref_default))
    updateSelectInput(session, "sys1_alt_col", choices = cols, selected = ifelse(is.na(alt_default), cols[1], alt_default))
  })
  
  inp1 <- reactive({
    req(sys1_raw(), input$sys1_chr_col, input$sys1_p_col, input$sys1_bp_col)
    
    df <- copy(sys1_raw())
    chr_col <- input$sys1_chr_col
    bp_col  <- input$sys1_bp_col
    p_col   <- input$sys1_p_col
    
    target_chrs <- parse_chr_input(input$sys1_chr_filter)
    if (!is.null(target_chrs)) {
      df <- df[as.character(get(chr_col)) %in% target_chrs]
    }
    
    if (!is.null(input$sys1_p_thresh) && !is.na(input$sys1_p_thresh)) {
      df <- df[get(p_col) < as.numeric(input$sys1_p_thresh)]
    }
    
    if (isTRUE(input$sys1_has_ref_alt) && !is.null(input$sys1_ref_col) && !is.null(input$sys1_alt_col)) {
      ref_c <- input$sys1_ref_col
      alt_c <- input$sys1_alt_col
      
      df[, gnomad_var_id := paste(
        toupper(gsub("^chr", "", get(chr_col), ignore.case = TRUE)),
        get(bp_col),
        toupper(get(ref_c)),
        toupper(get(alt_c)),
        sep = "-"
      )]
    } else {
      snp_c <- input$sys1_snp_col
      if (!is.null(snp_c) && snp_c %in% names(df)) {
        df[, gnomad_var_id := as.character(get(snp_c))]
      } else {
        df[, gnomad_var_id := paste(
          toupper(gsub("^chr", "", get(chr_col), ignore.case = TRUE)),
          get(bp_col), 
          sep = "-"
        )]
      }
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
    ref_default <- grep("^ref$|reference|a1|allele1", cols, ignore.case = TRUE, value = TRUE)[1]
    alt_default <- grep("^alt$|alternate|a2|allele2", cols, ignore.case = TRUE, value = TRUE)[1]
    
    updateSelectInput(session, "sys2_chr_col", choices = cols, selected = ifelse(is.na(chr_default), cols[1], chr_default))
    updateSelectInput(session, "sys2_bp_col", choices = cols, selected = ifelse(is.na(bp_default), cols[1], bp_default))
    updateSelectInput(session, "sys2_p_col", choices = cols, selected = ifelse(is.na(p_default), cols[1], p_default))
    updateSelectInput(session, "sys2_snp_col", choices = cols, selected = ifelse(is.na(snp_default), cols[1], snp_default))
    updateSelectInput(session, "sys2_ref_col", choices = cols, selected = ifelse(is.na(ref_default), cols[1], ref_default))
    updateSelectInput(session, "sys2_alt_col", choices = cols, selected = ifelse(is.na(alt_default), cols[1], alt_default))
  })
  
  inp2 <- reactive({
    req(sys2_raw(), input$sys2_chr_col, input$sys2_p_col, input$sys2_bp_col)
    
    df <- copy(sys2_raw())
    chr_col <- input$sys2_chr_col
    bp_col  <- input$sys2_bp_col
    p_col   <- input$sys2_p_col
    
    target_chrs <- parse_chr_input(input$sys2_chr_filter)
    if (!is.null(target_chrs)) {
      df <- df[as.character(get(chr_col)) %in% target_chrs]
    }
    
    if (!is.null(input$sys2_p_thresh) && !is.na(input$sys2_p_thresh)) {
      df <- df[get(p_col) < as.numeric(input$sys2_p_thresh)]
    }
    
    if (isTRUE(input$sys2_has_ref_alt) && !is.null(input$sys2_ref_col) && !is.null(input$sys2_alt_col)) {
      ref_c <- input$sys2_ref_col
      alt_c <- input$sys2_alt_col
      
      df[, gnomad_var_id := paste(
        toupper(gsub("^chr", "", get(chr_col), ignore.case = TRUE)),
        get(bp_col),
        toupper(get(ref_c)),
        toupper(get(alt_c)),
        sep = "-"
      )]
    } else {
      snp_c <- input$sys2_snp_col
      if (!is.null(snp_c) && snp_c %in% names(df)) {
        df[, gnomad_var_id := as.character(get(snp_c))]
      } else {
        df[, gnomad_var_id := paste(
          toupper(gsub("^chr", "", get(chr_col), ignore.case = TRUE)),
          get(bp_col), 
          sep = "-"
        )]
      }
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
    req(inp1(), input$sys1_chr_col, input$sys1_bp_col, input$sys1_p_col)
    
    p1 <- manhattanly(
      inp1(),
      chr = input$sys1_chr_col,
      bp = input$sys1_bp_col,
      p = input$sys1_p_col,
      snp = "gnomad_var_id",
      annotation1 = "gnomad_var_id",
      annotation2 = input$sys1_bp_col,
      col = c(input$sys1_col, input$sys1_col)
    )
    p1$x$source <- "manhattanPlot1"
    layout(p1, uirevision = reset_counter())
  })
  
  output$manhattanPlot2 <- renderPlotly({
    req(inp2(), input$sys2_chr_col, input$sys2_bp_col, input$sys2_p_col)
    
    p2 <- manhattanly(
      inp2(),
      chr = input$sys2_chr_col,
      bp = input$sys2_bp_col,
      p = input$sys2_p_col,
      snp = "gnomad_var_id",
      annotation1 = "gnomad_var_id",
      annotation2 = input$sys2_bp_col,
      col = c(input$sys2_col, input$sys2_col)
    )
    p2$x$source <- "manhattanPlot2"
    layout(p2, uirevision = reset_counter())
  })
  
  output$overlayPlot <- renderPlotly({     req(inp1(), inp2(), input$sys1_chr_col, input$sys1_bp_col, input$sys1_p_col)
    req(input$sys2_chr_col, input$sys2_bp_col, input$sys2_p_col)
    
    d1 <- copy(inp1())
    d1[, `:=`(
      CHR_plot = as.character(get(input$sys1_chr_col)),
      BP_plot  = as.numeric(get(input$sys1_bp_col)),
      P_plot   = -log10(as.numeric(get(input$sys1_p_col))),
      SNP_plot = as.character(gnomad_var_id)
    )]
    
    d2 <- copy(inp2())
    d2[, `:=`(
      CHR_plot = as.character(get(input$sys2_chr_col)),
      BP_plot  = as.numeric(get(input$sys2_bp_col)),
      P_plot   = -log10(as.numeric(get(input$sys2_p_col))),
      SNP_plot = as.character(gnomad_var_id)
    )]
    
    p_overlay <- plot_ly() %>%
      add_trace(
        data = d1,
        x = ~BP_plot,
        y = ~P_plot,
        type = 'scatter',
        mode = 'markers',
        name = 'System 1',
        key = ~SNP_plot,         # <--- ADD THIS
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
        key = ~SNP_plot,         # <--- ADD THIS
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
        uirevision = reset_counter()
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
  
  # Helper to extract numeric X-axis range from relayout events
  extract_x_range <- function(relayout_data) {
    if (is.null(relayout_data)) return(NULL)
    x0 <- relayout_data[["xaxis.range[0]"]]
    x1 <- relayout_data[["xaxis.range[1]"]]
    if (!is.null(x0) && !is.null(x1)) return(c(as.numeric(x0), as.numeric(x1)))
    if (!is.null(relayout_data[["xaxis.range"]])) return(as.numeric(relayout_data[["xaxis.range"]]))
    return(NULL)
  }
  
  last_synced_range <- reactiveVal(NULL)
  
  # Synchronize Zoom/Pan with Echo Suppression
  observeEvent(event_data("plotly_relayout", source = "manhattanPlot1"), {
    relayout_data <- event_data("plotly_relayout", source = "manhattanPlot1")
    req(relayout_data)
    
    rng <- extract_x_range(relayout_data)
    last_rng <- last_synced_range()
    
    # Ignore echo if the numerical range is already identical
    if (!is.null(rng) && !is.null(last_rng)) {
      if (isTRUE(all.equal(rng, last_rng, tolerance = 1e-4))) return()
    }
    
    if (!is.null(rng)) last_synced_range(rng)
    
    plotlyProxy("manhattanPlot2", session) %>%
      plotlyProxyInvoke("relayout", relayout_data)
  }, ignoreInit = TRUE)
  
  observeEvent(event_data("plotly_relayout", source = "manhattanPlot2"), {
    relayout_data <- event_data("plotly_relayout", source = "manhattanPlot2")
    req(relayout_data)
    
    rng <- extract_x_range(relayout_data)
    last_rng <- last_synced_range()
    
    # Ignore echo if the numerical range is already identical
    if (!is.null(rng) && !is.null(last_rng)) {
      if (isTRUE(all.equal(rng, last_rng, tolerance = 1e-4))) return()
    }
    
    if (!is.null(rng)) last_synced_range(rng)
    
    plotlyProxy("manhattanPlot1", session) %>%
      plotlyProxyInvoke("relayout", relayout_data)
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
        raw_val <- click_data$key
        if (is.null(raw_val)) raw_val <- click_data$customdata
        if (is.null(raw_val)) raw_val <- click_data$text
        chr_fb  <- parse_chr_input(input$sys1_chr_filter)[1]
        bp_fb   <- click_data$x
        
        snp_id <- clean_variant_id(raw_val, chr_fallback = chr_fb, bp_fallback = bp_fb)
        if (is.null(snp_id) || !nzchar(snp_id)) {
          snp_id <- paste(ifelse(is.null(chr_fb), "1", chr_fb), bp_fb, sep = "-")
        }
        
        af_val <- fetch_gnomad_af(snp_id, input$gnomad_build)
        
        hover_txt <- paste0(
          "<b>SNP/Variant:</b> ", snp_id, "<br>",
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
        raw_val <- click_data$key
        if (is.null(raw_val)) raw_val <- click_data$customdata
        if (is.null(raw_val)) raw_val <- click_data$text
        chr_fb  <- parse_chr_input(input$sys2_chr_filter)[1]
        bp_fb   <- click_data$x
        
        snp_id <- clean_variant_id(raw_val, chr_fallback = chr_fb, bp_fallback = bp_fb)
        if (is.null(snp_id) || !nzchar(snp_id)) {
          snp_id <- paste(ifelse(is.null(chr_fb), "1", chr_fb), bp_fb, sep = "-")
        }
        
        af_val <- fetch_gnomad_af(snp_id, input$gnomad_build)
        
        hover_txt <- paste0(
          "<b>SNP/Variant:</b> ", snp_id, "<br>",
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
        raw_val <- click_data$key
        if (is.null(raw_val)) raw_val <- click_data$customdata
        
        # Safely extract CHR from hover text if present
        txt_content <- click_data$text
        chr_extracted <- NULL
        if (!is.null(txt_content)) {
          m_chr <- regmatches(txt_content, regexec("CHR:\\s*([0-9A-Za-z]+)", txt_content))[[1]]
          if (length(m_chr) >= 2) chr_extracted <- m_chr[2]
        }
        
        if (is.null(chr_extracted)) {
          chr_extracted <- parse_chr_input(input$sys1_chr_filter)[1]
        }
        
        bp_fb <- click_data$x
        
        snp_id <- clean_variant_id(raw_val, chr_fallback = chr_extracted, bp_fallback = bp_fb)
        if (is.null(snp_id) || !nzchar(snp_id)) {
          snp_id <- paste(ifelse(is.null(chr_extracted), "1", chr_extracted), bp_fb, sep = "-")
        }
        
        af_val <- fetch_gnomad_af(snp_id, input$gnomad_build)
        
        hover_txt <- paste0(
          "<b>SNP/Variant:</b> ", snp_id, "<br>",
          "<b>CHR:</b> ", ifelse(is.null(chr_extracted), "N/A", chr_extracted), "<br>",
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
  # --- REGION FILTERING REACTIVES ---
  sys1_range <- reactive({
    relayout <- event_data("plotly_relayout", source = if (identical(input$plot_mode, "overlay")) "overlayPlot" else "manhattanPlot1")
    if (is.null(relayout)) return(NULL)
    
    x0 <- relayout[["xaxis.range[0]"]]
    x1 <- relayout[["xaxis.range[1]"]]
    if (!is.null(x0) && !is.null(x1)) return(c(as.numeric(x0), as.numeric(x1)))
    if (!is.null(relayout[["xaxis.range"]])) return(as.numeric(relayout[["xaxis.range"]]))
    return(NULL)
  })
  
  inp1_region <- reactive({
    df <- inp1()
    req(df, input$sys1_bp_col)
    rng <- sys1_range()
    if (!is.null(rng)) {
      bp_col <- input$sys1_bp_col
      df <- df[get(bp_col) >= rng[1] & get(bp_col) <= rng[2]]
    }
    df
  })
  
  sys2_range <- reactive({
    relayout <- event_data("plotly_relayout", source = if (identical(input$plot_mode, "overlay")) "overlayPlot" else "manhattanPlot2")
    if (is.null(relayout)) return(NULL)
    
    x0 <- relayout[["xaxis.range[0]"]]
    x1 <- relayout[["xaxis.range[1]"]]
    if (!is.null(x0) && !is.null(x1)) return(c(as.numeric(x0), as.numeric(x1)))
    if (!is.null(relayout[["xaxis.range"]])) return(as.numeric(relayout[["xaxis.range"]]))
    return(NULL)
  })
  
  inp2_region <- reactive({
    df <- inp2()
    req(df, input$sys2_bp_col)
    rng <- sys2_range()
    if (!is.null(rng)) {
      bp_col <- input$sys2_bp_col
      df <- df[get(bp_col) >= rng[1] & get(bp_col) <= rng[2]]
    }
    df
  })
  
  # --- NAVIGATION BUTTON ---
  observeEvent(input$go_to_nav3, {
    updateNavbarPage(session, "navbar", selected = "Navbar 3")
  })
  
  # --- NAVBAR 2 PREVIEW TABLES (< 50 POINTS CHECK) ---
  output$sys1_preview_ui <- renderUI({
    df <- inp1_region()
    if (is.null(df)) return(p("No dataset loaded."))
    if (nrow(df) >= 50) {
      p(sprintf("Region contains %d points. Zoom in until there are fewer than 50 points to show preview.", nrow(df)))
    } else {
      DTOutput("sys1_preview_table")
    }
  })
  
  output$sys1_preview_table <- renderDT({
    req(inp1_region())
    datatable(inp1_region(), options = list(pageLength = 5, scrollX = TRUE))
  })
  
  output$sys2_preview_ui <- renderUI({
    df <- inp2_region()
    if (is.null(df)) return(p("No dataset loaded."))
    if (nrow(df) >= 50) {
      p(sprintf("Region contains %d points. Zoom in until there are fewer than 50 points to show preview.", nrow(df)))
    } else {
      DTOutput("sys2_preview_table")
    }
  })
  
  output$sys2_preview_table <- renderDT({
    req(inp2_region())
    datatable(inp2_region(), options = list(pageLength = 5, scrollX = TRUE))
  })
  
  # --- NAVBAR 3 FILTERED TABLES & DOWNLOAD ---
  output$filtered_table1 <- renderDataTable({
    req(inp1_region())
    datatable(inp1_region(), options = list(scrollX = TRUE))
  })
  
  output$filtered_table2 <- renderDataTable({
    req(inp2_region())
    datatable(inp2_region(), options = list(scrollX = TRUE))
  })
  
  output$download_tables <- downloadHandler(
    filename = function() {
      paste0("region_filtered_tables_", Sys.Date(), ".zip")
    },
    content = function(file) {
      tmpdir <- tempdir()
      f1 <- file.path(tmpdir, "sys1_region_filtered.csv")
      f2 <- file.path(tmpdir, "sys2_region_filtered.csv")
      
      write.csv(inp1_region(), f1, row.names = FALSE)
      write.csv(inp2_region(), f2, row.names = FALSE)
      
      zip(file, files = c(f1, f2), flags = "-j")
    }
  )
}

shinyApp(ui = ui, server = server)