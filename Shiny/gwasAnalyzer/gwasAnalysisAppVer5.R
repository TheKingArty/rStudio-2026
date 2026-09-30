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

# Fallback operator for safe column selection
`%||%` <- function(x, y) {
  if (length(x) > 0 && !is.na(x[1]) && nzchar(as.character(x[1]))) x else y
}

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
    return("Invalid Format")
  }
  
  tryCatch({
    resp <- request("https://gnomad.broadinstitute.org/api") %>%
      req_user_agent("Mozilla/5.0 (Windows NT 10.0; Win64; x64) R-Shiny-App") %>%
      req_headers("Content-Type" = "application/json") %>%
      req_body_json(list(query = query_string)) %>%
      req_timeout(10) %>%
      req_error(is_error = function(resp) FALSE) %>%
      req_perform()
    
    if (resp_status(resp) >= 400) {
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
    if (!is.null(res_data$errors) || is.null(res_data$data)) return("Not Found")
    
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
    return("Error Fetching")
  })
}

# Store clicked annotations for up to 4 datasets + overlay
clicked_point_sys1 <- reactiveVal(NULL)
clicked_point_sys2 <- reactiveVal(NULL)
clicked_point_sys3 <- reactiveVal(NULL)
clicked_point_sys4 <- reactiveVal(NULL)
clicked_point_overlay <- reactiveVal(NULL)

# UI Definition
ui <- fluidPage(
  theme = shinytheme("cerulean"),
  navbarPage(
    "My first app",
    id = "navbar",
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
             
             # Dataset count selector slider
             fluidRow(
               column(6,
                      wellPanel(
                        sliderInput("num_datasets", "Number of Datasets to Analyze:", 
                                    min = 1, max = 4, value = 2, step = 1)
                      )
               )
             ),
             
             # File Uploads (Dynamic based on num_datasets)
             fluidRow(
               column(3, conditionalPanel(condition = "input.num_datasets >= 1", fileInput("sys1_file", "Upload Sys 1 Dataset", accept = c(".tsv", ".logistic", ".txt", ".csv")))),
               column(3, conditionalPanel(condition = "input.num_datasets >= 2", fileInput("sys2_file", "Upload Sys 2 Dataset", accept = c(".tsv", ".logistic", ".txt", ".csv")))),
               column(3, conditionalPanel(condition = "input.num_datasets >= 3", fileInput("sys3_file", "Upload Sys 3 Dataset", accept = c(".tsv", ".logistic", ".txt", ".csv")))),
               column(3, conditionalPanel(condition = "input.num_datasets >= 4", fileInput("sys4_file", "Upload Sys 4 Dataset", accept = c(".tsv", ".logistic", ".txt", ".csv"))))
             ),
             
             hr(),
             
             # Settings Panels Row 1 (System 1 & System 2)
             fluidRow(
               column(6,
                      conditionalPanel(
                        condition = "input.num_datasets >= 1",
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
                      )
               ),
               
               column(6,
                      conditionalPanel(
                        condition = "input.num_datasets >= 2",
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
               )
             ),
             
             # Settings Panels Row 2 (System 3 & System 4)
             fluidRow(
               column(6,
                      conditionalPanel(
                        condition = "input.num_datasets >= 3",
                        wellPanel(
                          h4("System 3 Settings"),
                          fluidRow(
                            column(6, selectInput("sys3_chr_col", "Chromosome Column", choices = NULL)),
                            column(6, selectInput("sys3_bp_col", "Position (BP) Column", choices = NULL))
                          ),
                          fluidRow(
                            column(6, selectInput("sys3_p_col", "P-Value Column", choices = NULL)),
                            column(6, selectInput("sys3_snp_col", "SNP ID Column", choices = NULL))
                          ),
                          checkboxInput("sys3_has_ref_alt", "Dataset includes REF and ALT columns", value = FALSE),
                          conditionalPanel(
                            condition = "input.sys3_has_ref_alt == true",
                            fluidRow(
                              column(6, selectInput("sys3_ref_col", "REF Allele Column", choices = NULL)),
                              column(6, selectInput("sys3_alt_col", "ALT Allele Column", choices = NULL))
                            )
                          ),
                          fluidRow(
                            column(6, textInput("sys3_chr_filter", "Chromosome Filter (e.g. 1, or blank)", value = "1")),
                            column(6, numericInput("sys3_p_thresh", "Max P-Value Threshold", value = 1e-5, step = 1e-6))
                          ),
                          fluidRow(
                            column(12, colourInput("sys3_col", "System 3 Point Color", value = "#2CA02C"))
                          )
                        )
                      )
               ),
               
               column(6,
                      conditionalPanel(
                        condition = "input.num_datasets >= 4",
                        wellPanel(
                          h4("System 4 Settings"),
                          fluidRow(
                            column(6, selectInput("sys4_chr_col", "Chromosome Column", choices = NULL)),
                            column(6, selectInput("sys4_bp_col", "Position (BP) Column", choices = NULL))
                          ),
                          fluidRow(
                            column(6, selectInput("sys4_p_col", "P-Value Column", choices = NULL)),
                            column(6, selectInput("sys4_snp_col", "SNP ID Column", choices = NULL))
                          ),
                          checkboxInput("sys4_has_ref_alt", "Dataset includes REF and ALT columns", value = FALSE),
                          conditionalPanel(
                            condition = "input.sys4_has_ref_alt == true",
                            fluidRow(
                              column(6, selectInput("sys4_ref_col", "REF Allele Column", choices = NULL)),
                              column(6, selectInput("sys4_alt_col", "ALT Allele Column", choices = NULL))
                            )
                          ),
                          fluidRow(
                            column(6, textInput("sys4_chr_filter", "Chromosome Filter (e.g. 1, or blank)", value = "1")),
                            column(6, numericInput("sys4_p_thresh", "Max P-Value Threshold", value = 1e-5, step = 1e-6))
                          ),
                          fluidRow(
                            column(12, colourInput("sys4_col", "System 4 Point Color", value = "#D62728"))
                          )
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
                                              choices = c("Individual Plots" = "individual", "Overlay All Datasets" = "overlay"),
                                              selected = "individual", inline = TRUE),
                                 
                                 conditionalPanel(
                                   condition = "input.plot_mode == 'individual'",
                                   radioButtons("individual_layout", "Individual Plot Layout:",
                                                choices = c("Side-by-Side/Grid" = "side", "Top-to-Bottom Stacked" = "stacked"),
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
               column(3, conditionalPanel(condition = "input.num_datasets >= 1", wellPanel(h4("Sys 1 Preview"), uiOutput("sys1_preview_ui")))),
               column(3, conditionalPanel(condition = "input.num_datasets >= 2", wellPanel(h4("Sys 2 Preview"), uiOutput("sys2_preview_ui")))),
               column(3, conditionalPanel(condition = "input.num_datasets >= 3", wellPanel(h4("Sys 3 Preview"), uiOutput("sys3_preview_ui")))),
               column(3, conditionalPanel(condition = "input.num_datasets >= 4", wellPanel(h4("Sys 4 Preview"), uiOutput("sys4_preview_ui"))))
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
    
    tabPanel("Navbar 3", 
             h3("Full Datasets"),
             conditionalPanel(condition = "input.num_datasets >= 1", dataTableOutput("table1")),
             conditionalPanel(condition = "input.num_datasets >= 2", dataTableOutput("table2")),
             conditionalPanel(condition = "input.num_datasets >= 3", dataTableOutput("table3")),
             conditionalPanel(condition = "input.num_datasets >= 4", dataTableOutput("table4")),
             
             hr(),
             
             h3("Region-Filtered Datasets"),
             downloadButton("download_tables", "Download tables", class = "btn-success"),
             br(), br(),
             conditionalPanel(condition = "input.num_datasets >= 1", dataTableOutput("filtered_table1")),
             conditionalPanel(condition = "input.num_datasets >= 2", dataTableOutput("filtered_table2")),
             conditionalPanel(condition = "input.num_datasets >= 3", dataTableOutput("filtered_table3")),
             conditionalPanel(condition = "input.num_datasets >= 4", dataTableOutput("filtered_table4"))
    )
  )
)

# Server Logic
server <- function(input, output, session) {
  reset_counter <- reactiveVal(0)
  
  # Reset Button Handler
  observeEvent(input$reset_plots, {
    clicked_point_sys1(NULL)
    clicked_point_sys2(NULL)
    clicked_point_sys3(NULL)
    clicked_point_sys4(NULL)
    clicked_point_overlay(NULL)
    
    reset_counter(reset_counter() + 1)
    
    plotlyProxy("manhattanPlot1", session) %>% plotlyProxyInvoke("relayout", list(annotations = list()))
    plotlyProxy("manhattanPlot2", session) %>% plotlyProxyInvoke("relayout", list(annotations = list()))
    plotlyProxy("manhattanPlot3", session) %>% plotlyProxyInvoke("relayout", list(annotations = list()))
    plotlyProxy("manhattanPlot4", session) %>% plotlyProxyInvoke("relayout", list(annotations = list()))
    plotlyProxy("overlayPlot", session) %>% plotlyProxyInvoke("relayout", list(annotations = list()))
  })
  
  last_synced_range <- reactiveVal(NULL)
  
  output$txtout <- renderText({
    paste(input$txt1, input$txt2, sep = " ")
  })
  
  # --- DATASET 1 LOGIC ---
  sys1_raw <- reactive({ req(input$sys1_file); fread(input$sys1_file$datapath) })
  observeEvent(sys1_raw(), {
    cols <- names(sys1_raw())
    updateSelectInput(session, "sys1_chr_col", choices = cols, selected = grep("chr|chrom", cols, ignore.case = TRUE, value = TRUE)[1] %||% cols[1])
    updateSelectInput(session, "sys1_bp_col", choices = cols, selected = grep("pos|bp", cols, ignore.case = TRUE, value = TRUE)[1] %||% cols[1])
    updateSelectInput(session, "sys1_p_col", choices = cols, selected = grep("^p$|p_val|p.val|pval", cols, ignore.case = TRUE, value = TRUE)[1] %||% cols[1])
    updateSelectInput(session, "sys1_snp_col", choices = cols, selected = grep("snp|id|rs", cols, ignore.case = TRUE, value = TRUE)[1] %||% cols[1])
    updateSelectInput(session, "sys1_ref_col", choices = cols, selected = grep("^ref$|reference|a1|allele1", cols, ignore.case = TRUE, value = TRUE)[1] %||% cols[1])
    updateSelectInput(session, "sys1_alt_col", choices = cols, selected = grep("^alt$|alternate|a2|allele2", cols, ignore.case = TRUE, value = TRUE)[1] %||% cols[1])
  })
  
  inp1 <- reactive({
    req(sys1_raw(), input$sys1_chr_col, input$sys1_p_col, input$sys1_bp_col)
    df <- copy(sys1_raw())
    target_chrs <- parse_chr_input(input$sys1_chr_filter)
    if (!is.null(target_chrs)) df <- df[as.character(get(input$sys1_chr_col)) %in% target_chrs]
    if (!is.null(input$sys1_p_thresh) && !is.na(input$sys1_p_thresh)) df <- df[get(input$sys1_p_col) < as.numeric(input$sys1_p_thresh)]
    
    if (isTRUE(input$sys1_has_ref_alt) && !is.null(input$sys1_ref_col) && !is.null(input$sys1_alt_col)) {
      df[, gnomad_var_id := paste(toupper(gsub("^chr", "", get(input$sys1_chr_col), ignore.case = TRUE)), get(input$sys1_bp_col), toupper(get(input$sys1_ref_col)), toupper(get(input$sys1_alt_col)), sep = "-")]
    } else {
      snp_c <- input$sys1_snp_col
      if (!is.null(snp_c) && snp_c %in% names(df)) df[, gnomad_var_id := as.character(get(snp_c))]
      else df[, gnomad_var_id := paste(toupper(gsub("^chr", "", get(input$sys1_chr_col), ignore.case = TRUE)), get(input$sys1_bp_col), sep = "-")]
    }
    df
  })
  
  # --- DATASET 2 LOGIC ---
  sys2_raw <- reactive({ req(input$sys2_file); fread(input$sys2_file$datapath) })
  observeEvent(sys2_raw(), {
    cols <- names(sys2_raw())
    updateSelectInput(session, "sys2_chr_col", choices = cols, selected = grep("chr|chrom", cols, ignore.case = TRUE, value = TRUE)[1] %||% cols[1])
    updateSelectInput(session, "sys2_bp_col", choices = cols, selected = grep("pos|bp", cols, ignore.case = TRUE, value = TRUE)[1] %||% cols[1])
    updateSelectInput(session, "sys2_p_col", choices = cols, selected = grep("^p$|p_val|p.val|pval", cols, ignore.case = TRUE, value = TRUE)[1] %||% cols[1])
    updateSelectInput(session, "sys2_snp_col", choices = cols, selected = grep("snp|id|rs", cols, ignore.case = TRUE, value = TRUE)[1] %||% cols[1])
    updateSelectInput(session, "sys2_ref_col", choices = cols, selected = grep("^ref$|reference|a1|allele1", cols, ignore.case = TRUE, value = TRUE)[1] %||% cols[1])
    updateSelectInput(session, "sys2_alt_col", choices = cols, selected = grep("^alt$|alternate|a2|allele2", cols, ignore.case = TRUE, value = TRUE)[1] %||% cols[1])
  })
  
  inp2 <- reactive({
    req(sys2_raw(), input$sys2_chr_col, input$sys2_p_col, input$sys2_bp_col)
    df <- copy(sys2_raw())
    target_chrs <- parse_chr_input(input$sys2_chr_filter)
    if (!is.null(target_chrs)) df <- df[as.character(get(input$sys2_chr_col)) %in% target_chrs]
    if (!is.null(input$sys2_p_thresh) && !is.na(input$sys2_p_thresh)) df <- df[get(input$sys2_p_col) < as.numeric(input$sys2_p_thresh)]
    
    if (isTRUE(input$sys2_has_ref_alt) && !is.null(input$sys2_ref_col) && !is.null(input$sys2_alt_col)) {
      df[, gnomad_var_id := paste(toupper(gsub("^chr", "", get(input$sys2_chr_col), ignore.case = TRUE)), get(input$sys2_bp_col), toupper(get(input$sys2_ref_col)), toupper(get(input$sys2_alt_col)), sep = "-")]
    } else {
      snp_c <- input$sys2_snp_col
      if (!is.null(snp_c) && snp_c %in% names(df)) df[, gnomad_var_id := as.character(get(snp_c))]
      else df[, gnomad_var_id := paste(toupper(gsub("^chr", "", get(input$sys2_chr_col), ignore.case = TRUE)), get(input$sys2_bp_col), sep = "-")]
    }
    df
  })
  
  # --- DATASET 3 LOGIC ---
  sys3_raw <- reactive({ req(input$sys3_file); fread(input$sys3_file$datapath) })
  observeEvent(sys3_raw(), {
    cols <- names(sys3_raw())
    updateSelectInput(session, "sys3_chr_col", choices = cols, selected = grep("chr|chrom", cols, ignore.case = TRUE, value = TRUE)[1] %||% cols[1])
    updateSelectInput(session, "sys3_bp_col", choices = cols, selected = grep("pos|bp", cols, ignore.case = TRUE, value = TRUE)[1] %||% cols[1])
    updateSelectInput(session, "sys3_p_col", choices = cols, selected = grep("^p$|p_val|p.val|pval", cols, ignore.case = TRUE, value = TRUE)[1] %||% cols[1])
    updateSelectInput(session, "sys3_snp_col", choices = cols, selected = grep("snp|id|rs", cols, ignore.case = TRUE, value = TRUE)[1] %||% cols[1])
    updateSelectInput(session, "sys3_ref_col", choices = cols, selected = grep("^ref$|reference|a1|allele1", cols, ignore.case = TRUE, value = TRUE)[1] %||% cols[1])
    updateSelectInput(session, "sys3_alt_col", choices = cols, selected = grep("^alt$|alternate|a2|allele2", cols, ignore.case = TRUE, value = TRUE)[1] %||% cols[1])
  })
  
  inp3 <- reactive({
    req(sys3_raw(), input$sys3_chr_col, input$sys3_p_col, input$sys3_bp_col)
    df <- copy(sys3_raw())
    target_chrs <- parse_chr_input(input$sys3_chr_filter)
    if (!is.null(target_chrs)) df <- df[as.character(get(input$sys3_chr_col)) %in% target_chrs]
    if (!is.null(input$sys3_p_thresh) && !is.na(input$sys3_p_thresh)) df <- df[get(input$sys3_p_col) < as.numeric(input$sys3_p_thresh)]
    
    if (isTRUE(input$sys3_has_ref_alt) && !is.null(input$sys3_ref_col) && !is.null(input$sys3_alt_col)) {
      df[, gnomad_var_id := paste(toupper(gsub("^chr", "", get(input$sys3_chr_col), ignore.case = TRUE)), get(input$sys3_bp_col), toupper(get(input$sys3_ref_col)), toupper(get(input$sys3_alt_col)), sep = "-")]
    } else {
      snp_c <- input$sys3_snp_col
      if (!is.null(snp_c) && snp_c %in% names(df)) df[, gnomad_var_id := as.character(get(snp_c))]
      else df[, gnomad_var_id := paste(toupper(gsub("^chr", "", get(input$sys3_chr_col), ignore.case = TRUE)), get(input$sys3_bp_col), sep = "-")]
    }
    df
  })
  
  # --- DATASET 4 LOGIC ---
  sys4_raw <- reactive({ req(input$sys4_file); fread(input$sys4_file$datapath) })
  observeEvent(sys4_raw(), {
    cols <- names(sys4_raw())
    updateSelectInput(session, "sys4_chr_col", choices = cols, selected = grep("chr|chrom", cols, ignore.case = TRUE, value = TRUE)[1] %||% cols[1])
    updateSelectInput(session, "sys4_bp_col", choices = cols, selected = grep("pos|bp", cols, ignore.case = TRUE, value = TRUE)[1] %||% cols[1])
    updateSelectInput(session, "sys4_p_col", choices = cols, selected = grep("^p$|p_val|p.val|pval", cols, ignore.case = TRUE, value = TRUE)[1] %||% cols[1])
    updateSelectInput(session, "sys4_snp_col", choices = cols, selected = grep("snp|id|rs", cols, ignore.case = TRUE, value = TRUE)[1] %||% cols[1])
    updateSelectInput(session, "sys4_ref_col", choices = cols, selected = grep("^ref$|reference|a1|allele1", cols, ignore.case = TRUE, value = TRUE)[1] %||% cols[1])
    updateSelectInput(session, "sys4_alt_col", choices = cols, selected = grep("^alt$|alternate|a2|allele2", cols, ignore.case = TRUE, value = TRUE)[1] %||% cols[1])
  })
  
  inp4 <- reactive({
    req(sys4_raw(), input$sys4_chr_col, input$sys4_p_col, input$sys4_bp_col)
    df <- copy(sys4_raw())
    target_chrs <- parse_chr_input(input$sys4_chr_filter)
    if (!is.null(target_chrs)) df <- df[as.character(get(input$sys4_chr_col)) %in% target_chrs]
    if (!is.null(input$sys4_p_thresh) && !is.na(input$sys4_p_thresh)) df <- df[get(input$sys4_p_col) < as.numeric(input$sys4_p_thresh)]
    
    if (isTRUE(input$sys4_has_ref_alt) && !is.null(input$sys4_ref_col) && !is.null(input$sys4_alt_col)) {
      df[, gnomad_var_id := paste(toupper(gsub("^chr", "", get(input$sys4_chr_col), ignore.case = TRUE)), get(input$sys4_bp_col), toupper(get(input$sys4_ref_col)), toupper(get(input$sys4_alt_col)), sep = "-")]
    } else {
      snp_c <- input$sys4_snp_col
      if (!is.null(snp_c) && snp_c %in% names(df)) df[, gnomad_var_id := as.character(get(snp_c))]
      else df[, gnomad_var_id := paste(toupper(gsub("^chr", "", get(input$sys4_chr_col), ignore.case = TRUE)), get(input$sys4_bp_col), sep = "-")]
    }
    df
  })
  
  # --- DYNAMIC PLOT CONTAINER ---
  output$plot_container <- renderUI({
    n <- input$num_datasets
    if (input$plot_mode == "individual") {
      if (input$individual_layout == "side") {
        if (n == 1) {
          fluidRow(column(12, plotlyOutput("manhattanPlot1", height = "600px")))
        } else if (n == 2) {
          fluidRow(column(6, plotlyOutput("manhattanPlot1", height = "600px")), column(6, plotlyOutput("manhattanPlot2", height = "600px")))
        } else if (n == 3) {
          fluidRow(column(4, plotlyOutput("manhattanPlot1", height = "600px")), column(4, plotlyOutput("manhattanPlot2", height = "600px")), column(4, plotlyOutput("manhattanPlot3", height = "600px")))
        } else {
          fluidRow(
            column(6, plotlyOutput("manhattanPlot1", height = "500px")), column(6, plotlyOutput("manhattanPlot2", height = "500px")),
            column(6, plotlyOutput("manhattanPlot3", height = "500px")), column(6, plotlyOutput("manhattanPlot4", height = "500px"))
          )
        }
      } else {
        # Stacked layout
        res <- list()
        for (i in seq_len(n)) {
          res[[length(res) + 1]] <- column(12, plotlyOutput(paste0("manhattanPlot", i), height = "500px"))
          if (i < n) res[[length(res) + 1]] <- column(12, br())
        }
        do.call(fluidRow, res)
      }
    } else {
      fluidRow(column(12, plotlyOutput("overlayPlot", height = "650px")))
    }
  })
  
  # Render individual manhattan plots
  output$manhattanPlot1 <- renderPlotly({
    req(inp1(), input$sys1_chr_col, input$sys1_bp_col, input$sys1_p_col)
    p1 <- manhattanly(inp1(), chr = input$sys1_chr_col, bp = input$sys1_bp_col, p = input$sys1_p_col, snp = "gnomad_var_id", annotation1 = "gnomad_var_id", annotation2 = input$sys1_bp_col, col = c(input$sys1_col, input$sys1_col))
    p1$x$source <- "manhattanPlot1"
    layout(p1, uirevision = reset_counter())
  })
  
  output$manhattanPlot2 <- renderPlotly({
    req(input$num_datasets >= 2, inp2(), input$sys2_chr_col, input$sys2_bp_col, input$sys2_p_col)
    p2 <- manhattanly(inp2(), chr = input$sys2_chr_col, bp = input$sys2_bp_col, p = input$sys2_p_col, snp = "gnomad_var_id", annotation1 = "gnomad_var_id", annotation2 = input$sys2_bp_col, col = c(input$sys2_col, input$sys2_col))
    p2$x$source <- "manhattanPlot2"
    layout(p2, uirevision = reset_counter())
  })
  
  output$manhattanPlot3 <- renderPlotly({
    req(input$num_datasets >= 3, inp3(), input$sys3_chr_col, input$sys3_bp_col, input$sys3_p_col)
    p3 <- manhattanly(inp3(), chr = input$sys3_chr_col, bp = input$sys3_bp_col, p = input$sys3_p_col, snp = "gnomad_var_id", annotation1 = "gnomad_var_id", annotation2 = input$sys3_bp_col, col = c(input$sys3_col, input$sys3_col))
    p3$x$source <- "manhattanPlot3"
    layout(p3, uirevision = reset_counter())
  })
  
  output$manhattanPlot4 <- renderPlotly({
    req(input$num_datasets >= 4, inp4(), input$sys4_chr_col, input$sys4_bp_col, input$sys4_p_col)
    p4 <- manhattanly(inp4(), chr = input$sys4_chr_col, bp = input$sys4_bp_col, p = input$sys4_p_col, snp = "gnomad_var_id", annotation1 = "gnomad_var_id", annotation2 = input$sys4_bp_col, col = c(input$sys4_col, input$sys4_col))
    p4$x$source <- "manhattanPlot4"
    layout(p4, uirevision = reset_counter())
  })
  
  # Overlay Plot with up to 4 datasets
  output$overlayPlot <- renderPlotly({
    n <- input$num_datasets
    req(inp1(), input$sys1_chr_col, input$sys1_bp_col, input$sys1_p_col)
    if (n >= 2) req(inp2(), input$sys2_chr_col, input$sys2_bp_col, input$sys2_p_col)
    if (n >= 3) req(inp3(), input$sys3_chr_col, input$sys3_bp_col, input$sys3_p_col)
    if (n >= 4) req(inp4(), input$sys4_chr_col, input$sys4_bp_col, input$sys4_p_col)
    
    p_overlay <- plot_ly()
    
    # Dataset 1 trace
    d1 <- copy(inp1())
    d1[, `:=`(CHR_plot = as.character(get(input$sys1_chr_col)), BP_plot = as.numeric(get(input$sys1_bp_col)), P_plot = -log10(as.numeric(get(input$sys1_p_col))), SNP_plot = as.character(gnomad_var_id))]
    p_overlay <- p_overlay %>% add_trace(data = d1, x = ~BP_plot, y = ~P_plot, type = 'scatter', mode = 'markers', name = 'System 1', key = ~SNP_plot, marker = list(color = input$sys1_col, size = 6, opacity = 0.7), customdata = ~SNP_plot, text = ~paste("SNP:", SNP_plot, "<br>CHR:", CHR_plot, "<br>BP:", BP_plot, "<br>-log10(P):", round(P_plot, 3)), hoverinfo = "text")
    
    if (n >= 2) {
      d2 <- copy(inp2())
      d2[, `:=`(CHR_plot = as.character(get(input$sys2_chr_col)), BP_plot = as.numeric(get(input$sys2_bp_col)), P_plot = -log10(as.numeric(get(input$sys2_p_col))), SNP_plot = as.character(gnomad_var_id))]
      p_overlay <- p_overlay %>% add_trace(data = d2, x = ~BP_plot, y = ~P_plot, type = 'scatter', mode = 'markers', name = 'System 2', key = ~SNP_plot, marker = list(color = input$sys2_col, size = 6, opacity = 0.7), customdata = ~SNP_plot, text = ~paste("SNP:", SNP_plot, "<br>CHR:", CHR_plot, "<br>BP:", BP_plot, "<br>-log10(P):", round(P_plot, 3)), hoverinfo = "text")
    }
    
    if (n >= 3) {
      d3 <- copy(inp3())
      d3[, `:=`(CHR_plot = as.character(get(input$sys3_chr_col)), BP_plot = as.numeric(get(input$sys3_bp_col)), P_plot = -log10(as.numeric(get(input$sys3_p_col))), SNP_plot = as.character(gnomad_var_id))]
      p_overlay <- p_overlay %>% add_trace(data = d3, x = ~BP_plot, y = ~P_plot, type = 'scatter', mode = 'markers', name = 'System 3', key = ~SNP_plot, marker = list(color = input$sys3_col, size = 6, opacity = 0.7), customdata = ~SNP_plot, text = ~paste("SNP:", SNP_plot, "<br>CHR:", CHR_plot, "<br>BP:", BP_plot, "<br>-log10(P):", round(P_plot, 3)), hoverinfo = "text")
    }
    
    if (n >= 4) {
      d4 <- copy(inp4())
      d4[, `:=`(CHR_plot = as.character(get(input$sys4_chr_col)), BP_plot = as.numeric(get(input$sys4_bp_col)), P_plot = -log10(as.numeric(get(input$sys4_p_col))), SNP_plot = as.character(gnomad_var_id))]
      p_overlay <- p_overlay %>% add_trace(data = d4, x = ~BP_plot, y = ~P_plot, type = 'scatter', mode = 'markers', name = 'System 4', key = ~SNP_plot, marker = list(color = input$sys4_col, size = 6, opacity = 0.7), customdata = ~SNP_plot, text = ~paste("SNP:", SNP_plot, "<br>CHR:", CHR_plot, "<br>BP:", BP_plot, "<br>-log10(P):", round(P_plot, 3)), hoverinfo = "text")
    }
    
    p_overlay <- p_overlay %>% layout(
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
    
    tags$a(href = target_url, target = "_blank", class = "btn btn-primary", icon("external-link-alt"), paste("Open", query_trimmed, "in gnoMAD"))
  })
  
  extract_x_range <- function(relayout_data) {
    if (is.null(relayout_data)) return(NULL)
    x0 <- relayout_data[["xaxis.range[0]"]]
    x1 <- relayout_data[["xaxis.range[1]"]]
    if (!is.null(x0) && !is.null(x1)) return(c(as.numeric(x0), as.numeric(x1)))
    if (!is.null(relayout_data[["xaxis.range"]])) return(as.numeric(relayout_data[["xaxis.range"]]))
    return(NULL)
  }
  
  # Zoom synchronization across individual plots
  lapply(1:4, function(i) {
    observeEvent(event_data("plotly_relayout", source = paste0("manhattanPlot", i)), {
      relayout_data <- event_data("plotly_relayout", source = paste0("manhattanPlot", i))
      req(relayout_data)
      rng <- extract_x_range(relayout_data)
      last_rng <- last_synced_range()
      if (!is.null(rng) && !is.null(last_rng) && isTRUE(all.equal(rng, last_rng, tolerance = 1e-4))) return()
      if (!is.null(rng)) last_synced_range(rng)
      
      for (j in 1:input$num_datasets) {
        if (i != j) {
          plotlyProxy(paste0("manhattanPlot", j), session) %>% plotlyProxyInvoke("relayout", relayout_data)
        }
      }
    }, ignoreInit = TRUE)
  })
  
  # Click observers generator for systems 1-4
  local({
    gen_click_observer <- function(sys_id) {
      observeEvent(event_data("plotly_click", source = paste0("manhattanPlot", sys_id)), {
        click_data <- event_data("plotly_click", source = paste0("manhattanPlot", sys_id))
        if (!is.null(click_data)) {
          reactive_val_getter <- get(paste0("clicked_point_sys", sys_id))
          reactive_val_setter <- get(paste0("clicked_point_sys", sys_id))
          current <- reactive_val_getter()
          
          if (!is.null(current) && identical(current$x, click_data$x) && identical(current$y, click_data$y)) {
            reactive_val_setter(NULL)
            plotlyProxy(paste0("manhattanPlot", sys_id), session) %>% plotlyProxyInvoke("relayout", list(annotations = list()))
          } else {
            raw_val <- click_data$key
            if (is.null(raw_val)) raw_val <- click_data$customdata
            if (is.null(raw_val)) raw_val <- click_data$text
            chr_fb  <- parse_chr_input(input[[paste0("sys", sys_id, "_chr_filter")]])[1]
            bp_fb   <- click_data$x
            
            snp_id <- clean_variant_id(raw_val, chr_fallback = chr_fb, bp_fallback = bp_fb)
            if (is.null(snp_id) || !nzchar(snp_id)) snp_id <- paste(ifelse(is.null(chr_fb), "1", chr_fb), bp_fb, sep = "-")
            
            af_val <- fetch_gnomad_af(snp_id, input$gnomad_build)
            hover_txt <- paste0("<b>SNP/Variant:</b> ", snp_id, "<br><b>BP:</b> ", click_data$x, "<br><b>-log10(P):</b> ", round(click_data$y, 3), "<br><b>gnoMAD AF:</b> ", af_val)
            
            annotation <- list(x = click_data$x, y = click_data$y, text = hover_txt, showarrow = FALSE, xanchor = "left", yanchor = "bottom", bgcolor = "rgba(255, 255, 255, 0.95)", bordercolor = "#444444", borderpad = 6, font = list(size = 12, color = "#000000"))
            reactive_val_setter(annotation)
            plotlyProxy(paste0("manhattanPlot", sys_id), session) %>% plotlyProxyInvoke("relayout", list(annotations = list(annotation)))
          }
        }
      })
    }
    for (i in 1:4) gen_click_observer(i)
  })
  
  observeEvent(event_data("plotly_click", source = "overlayPlot"), {
    click_data <- event_data("plotly_click", source = "overlayPlot")
    if (!is.null(click_data)) {
      current <- clicked_point_overlay()
      if (!is.null(current) && identical(current$x, click_data$x) && identical(current$y, click_data$y)) {
        clicked_point_overlay(NULL)
        plotlyProxy("overlayPlot", session) %>% plotlyProxyInvoke("relayout", list(annotations = list()))
      } else {
        raw_val <- click_data$key
        if (is.null(raw_val)) raw_val <- click_data$customdata
        txt_content <- click_data$text
        chr_extracted <- NULL
        if (!is.null(txt_content)) {
          m_chr <- regmatches(txt_content, regexec("CHR:\\s*([0-9A-Za-z]+)", txt_content))[[1]]
          if (length(m_chr) >= 2) chr_extracted <- m_chr[2]
        }
        if (is.null(chr_extracted)) chr_extracted <- parse_chr_input(input$sys1_chr_filter)[1]
        bp_fb <- click_data$x
        
        snp_id <- clean_variant_id(raw_val, chr_fallback = chr_extracted, bp_fallback = bp_fb)
        if (is.null(snp_id) || !nzchar(snp_id)) snp_id <- paste(ifelse(is.null(chr_extracted), "1", chr_extracted), bp_fb, sep = "-")
        
        af_val <- fetch_gnomad_af(snp_id, input$gnomad_build)
        hover_txt <- paste0("<b>SNP/Variant:</b> ", snp_id, "<br><b>CHR:</b> ", ifelse(is.null(chr_extracted), "N/A", chr_extracted), "<br><b>BP:</b> ", click_data$x, "<br><b>-log10(P):</b> ", round(click_data$y, 3), "<br><b>gnoMAD AF:</b> ", af_val)
        
        annotation <- list(x = click_data$x, y = click_data$y, text = hover_txt, showarrow = FALSE, xanchor = "left", yanchor = "bottom", bgcolor = "rgba(255, 255, 255, 0.95)", bordercolor = "#444444", borderpad = 6, font = list(size = 12, color = "#000000"))
        clicked_point_overlay(annotation)
        plotlyProxy("overlayPlot", session) %>% plotlyProxyInvoke("relayout", list(annotations = list(annotation)))
      }
    }
  })
  
  # Tables full output
  output$table1 <- renderDataTable({ req(inp1()); datatable(inp1()) })
  output$table2 <- renderDataTable({ req(input$num_datasets >= 2, inp2()); datatable(inp2()) })
  output$table3 <- renderDataTable({ req(input$num_datasets >= 3, inp3()); datatable(inp3()) })
  output$table4 <- renderDataTable({ req(input$num_datasets >= 4, inp4()); datatable(inp4()) })
  
  # --- REGION FILTERING REACTIVES ---
  gen_region_reactives <- function(sys_id) {
    r_source <- reactive({
      relayout <- event_data("plotly_relayout", source = if (identical(input$plot_mode, "overlay")) "overlayPlot" else paste0("manhattanPlot", sys_id))
      if (is.null(relayout)) return(NULL)
      x0 <- relayout[["xaxis.range[0]"]]
      x1 <- relayout[["xaxis.range[1]"]]
      if (!is.null(x0) && !is.null(x1)) return(c(as.numeric(x0), as.numeric(x1)))
      if (!is.null(relayout[["xaxis.range"]])) return(as.numeric(relayout[["xaxis.range"]]))
      return(NULL)
    })
    
    r_filtered <- reactive({
      df <- get(paste0("inp", sys_id))()
      req(df, input[[paste0("sys", sys_id, "_bp_col")]])
      rng <- r_source()
      if (!is.null(rng)) {
        bp_col <- input[[paste0("sys", sys_id, "_bp_col")]]
        df <- df[get(bp_col) >= rng[1] & get(bp_col) <= rng[2]]
      }
      df
    })
    return(r_filtered)
  }
  
  inp1_region <- gen_region_reactives(1)
  inp2_region <- gen_region_reactives(2)
  inp3_region <- gen_region_reactives(3)
  inp4_region <- gen_region_reactives(4)
  
  observeEvent(input$go_to_nav3, {
    updateNavbarPage(session, "navbar", selected = "Navbar 3")
  })
  
  # Preview tables UI & rendering
  lapply(1:4, function(i) {
    output[[paste0("sys", i, "_preview_ui")]] <- renderUI({
      if (input$num_datasets < i) return(NULL)
      df <- get(paste0("inp", i, "_region"))()
      if (is.null(df)) return(p("No dataset loaded."))
      if (nrow(df) >= 50) {
        p(sprintf("Region contains %d points. Zoom in (< 50 pts) for preview.", nrow(df)))
      } else {
        DTOutput(paste0("sys", i, "_preview_table"))
      }
    })
    
    output[[paste0("sys", i, "_preview_table")]] <- renderDT({
      req(input$num_datasets >= i, get(paste0("inp", i, "_region"))())
      datatable(get(paste0("inp", i, "_region"))(), options = list(pageLength = 5, scrollX = TRUE))
    })
    
    output[[paste0("filtered_table", i)]] <- renderDataTable({
      req(input$num_datasets >= i, get(paste0("inp", i, "_region"))())
      datatable(get(paste0("inp", i, "_region"))(), options = list(scrollX = TRUE))
    })
  })
  
  output$download_tables <- downloadHandler(
    filename = function() { paste0("region_filtered_tables_", Sys.Date(), ".zip") },
    content = function(file) {
      tmpdir <- tempdir()
      file_paths <- c()
      n <- input$num_datasets
      for (i in seq_len(n)) {
        df_reg <- get(paste0("inp", i, "_region"))()
        if (!is.null(df_reg)) {
          f_path <- file.path(tmpdir, sprintf("sys%d_region_filtered.csv", i))
          write.csv(df_reg, f_path, row.names = FALSE)
          file_paths <- c(file_paths, f_path)
        }
      }
      if (length(file_paths) > 0) {
        zip(file, files = file_paths, flags = "-j")
      }
    }
  )
}

shinyApp(ui = ui, server = server)