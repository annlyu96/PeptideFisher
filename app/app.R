library(shiny)
library(dplyr)
library(readr)
library(readxl)
library(DT)
library(shinydashboard)
library(scales)

options(shiny.maxRequestSize = 100*1024^2)

# ==========================================================
# Peptide statistics
# ==========================================================
compute_peptide_stats <- function(data, control_cols, drug_cols) {
  
  control_mat <- as.matrix(data[, control_cols])
  drug_mat    <- as.matrix(data[, drug_cols])
  
  mean_control <- rowMeans(control_mat, na.rm = TRUE)
  mean_drug    <- rowMeans(drug_mat, na.rm = TRUE)
  
  pvals <- apply(
    cbind(drug_mat, control_mat),
    1,
    function(x) {
      d <- x[1:ncol(drug_mat)]
      c <- x[(ncol(drug_mat)+1):length(x)]
      tryCatch(t.test(d, c)$p.value, error = function(e) NA_real_)
    }
  )
  
  data %>%
    mutate(
      FoldChange = mean_drug / mean_control,
      log2FoldChange = log2(FoldChange),
      `p-value` = pvals,
      `-log10p-value` = -log10(`p-value`)
    )
}

# ==========================================================
# Fisher function (dual behaviour)
# ==========================================================
compute_fisher_by_topN_p_noimpute <- function(data, N, mode) {
  
  data %>%
    group_by(Accession) %>%
    arrange(desc(`-log10p-value`), .by_group = TRUE) %>%
    slice_head(n = N) %>%
    summarise(
      `Gene name` = first(`Gene name`),
      `Protein description` = first(`Protein description`),
      
      fisher_X2 = -2 * sum(log(`p-value`), na.rm = TRUE),
      df_used = 2 * n(),
      fisher_p = pchisq(fisher_X2, df = df_used, lower.tail = FALSE),
      `-log10fisher_p` = -log10(fisher_p),
      
      log2_mean_FC = if (mode == "partial") {
        pos_count <- sum(log2FoldChange > 0, na.rm = TRUE)
        neg_count <- sum(log2FoldChange < 0, na.rm = TRUE)
        direction <- ifelse(pos_count >= neg_count, 1, -1)
        direction * mean(abs(log2FoldChange), na.rm = TRUE)
      } else {
        NA_real_
      },
      
      .groups = "drop"
    )
}

# ==========================================================
# UI 
# ==========================================================
ui <- dashboardPage(
  skin = "blue",
  
  dashboardHeader(
    title = tags$img(
      src = "RZlab_fisher.jpg",
      height = "55px",
      style = "padding-top:0px;"
    ),
    titleWidth = 250,
    tags$li(
      class = "dropdown",
      tags$span(
        "Fisher’s method for peptide-to-protein data aggregation",
        style = "color: white; font-size: 25px; font-weight: bold; 
               line-height: 50px; margin-left: 15px;"
      )
    )
  ),

  
  dashboardSidebar(
    width = 250,
    
    radioButtons(
      "data_mode",
      "Data Type",
      choices = c(
        "Expression / Solubility" = "expression",
        "Partial proteolysis (AFDIP / HOLSER)" = "partial"
      )
    ),
    
    fileInput("peptide_file", "Upload Peptide File"),
    uiOutput("peptide_selector"),
    hr(),
    
    conditionalPanel(
      condition = "input.data_mode == 'expression'",
      fileInput("protein_file", "Upload Protein File"),
      uiOutput("protein_selector"),
      hr()
    ),
    
    numericInput(
      "N_value",
      "Top N peptides",
      value = 4,
      min = 1,
      step = 1
    ),
    
    actionButton("run", "Run Analysis"),
    br(), br(),
    
    
    div(
      style = "padding-left: 15px; padding-right: 15px;",
      downloadButton(
        "download",
        "Download Result",
        width = "100%",
        class = "btn btn-primary"
      )
    )
  ),
  
  dashboardBody(
    tags$style(HTML('
      
      /* Main background */
      body, .content-wrapper, .main-sidebar {
        background-color: #EDF4F4 !important;
        color: #4F0433 !important;
      }
      
      /* Header background */
      .skin-blue .main-header .logo, .skin-blue .main-header .navbar {
      background-color: #4F0433 !important;
      }
      .main-header .navbar {
      display: flex !important;
      align-items: center !important;
      }
      
      /* sidebar font */
      .skin-blue .main-sidebar,
      .skin-blue .main-sidebar label,
      .skin-blue .main-sidebar span,
      .skin-blue .main-sidebar a {
        color: #4F0433 !important;
      }



      /* button */
      .btn-primary {
        background-color: #FF876F !important;
        color: #4F0433 !important;
        border-color: #870052 !important;
      }

      .btn-primary:hover {
        background-color: #FF7050 !important;
      }

      /* result box */
      .custom-box {
        border: 1px solid #4F0433;
      }

      .custom-box .box-header {
        background-color: #4F0433 !important;
        color: white !important;
      }

    ')),
    
    fluidRow(
      box(
        title = "Fisher Result",
        width = 12,
        solidHeader = TRUE,
        DTOutput("result_table"),
        class = "custom-box"
      )
    )
  )
)

# ==========================================================
# SERVER
# ==========================================================
server <- function(input, output, session) {
  
  read_user_file <- function(file) {
    ext <- tools::file_ext(file$datapath)
    if (ext %in% c("tsv","txt")) read_tsv(file$datapath)
    else if (ext == "csv") read_csv(file$datapath)
    else if (ext %in% c("xlsx","xls")) read_excel(file$datapath)
  }
  
  peptide_data <- reactive({
    req(input$peptide_file)
    read_user_file(input$peptide_file)
  })
  
  protein_data <- reactive({
    req(input$protein_file)
    read_user_file(input$protein_file)
  })
  
  output$peptide_selector <- renderUI({
    req(peptide_data())
    cols <- colnames(peptide_data())
    tagList(
      selectizeInput("peptide_control","Peptide Control",
                     choices = cols, multiple = TRUE),
      selectizeInput("peptide_drug","Peptide Drug",
                     choices = cols, multiple = TRUE)
    )
  })
  
  output$protein_selector <- renderUI({
    req(protein_data())
    cols <- colnames(protein_data())
    tagList(
      selectizeInput("protein_control","Protein Control",
                     choices = cols, multiple = TRUE),
      selectizeInput("protein_drug","Protein Drug",
                     choices = cols, multiple = TRUE)
    )
  })
  
  
  final_result <- eventReactive(input$run, {
    
    # 1. Initialize progress bar
    withProgress(message = 'Running analysis...', value = 0, {
      
      # ----------------------------------------------------
      # stage 1: data validation (Progress: 10%)
      incProgress(0.1, detail = "Validating data...")
      
      validate(
        need(input$peptide_file, "Please upload peptide file"),
        need(length(input$peptide_control) > 0, "Select peptide control columns"),
        need(length(input$peptide_drug) > 0, "Select peptide drug columns"),
        need("Accession" %in% colnames(peptide_data()),
             "Peptide file must contain column: Accession")
      )
      
      if (input$data_mode == "partial") {
        required_cols <- c("Accession","Gene name","Protein description")
        validate(
          need(all(required_cols %in% colnames(peptide_data())),
               "Partial mode requires: Accession, Gene name, Protein description")
        )
      }
      
      
      withProgress(message = 'Running analysis...', value = 0, {
        
        # ----------------------------------------------------
        # Stage 1: Data Validation (Progress: 10%)
        incProgress(0.1, detail = "Validating data...")
        
        # 1. validation on number of peptide replicates (both modes)
        validate(
          need(input$peptide_file, "Please upload peptide file"),
          need(length(input$peptide_control) >= 2, 
               "Insufficient Peptide Control replicates: Please select at least 2 columns."),
          need(length(input$peptide_drug) >= 2, 
               "Insufficient Peptide Drug replicates: Please select at least 2 columns."),
          need("Accession" %in% colnames(peptide_data()),
               "Peptide file must contain column: Accession")
        )
        
        # 2. validations for different modes
        if (input$data_mode == "partial") {
          # Partial mode: check column names
          required_cols <- c("Accession", "Gene name", "Protein description")
          validate(
            need(all(required_cols %in% colnames(peptide_data())),
                 "Partial mode requires: Accession, Gene name, Protein description")
          )
        } else if (input$data_mode == "expression") {
          # Expression mode: check protein column names and number of replicates
          validate(
            need(input$protein_file, "Please upload protein file for Expression mode"),
            need(length(input$protein_control) >= 2, "Insufficient Protein Control: select at least 2 columns"),
            need(length(input$protein_drug) >= 2, "Insufficient Protein Drug: select at least 2 columns"),
            need(all(required_info_cols %in% colnames(protein_data())),
                 "Protein file must contain: Accession, Gene name, Protein description")
          )
        }
      
      
      
      # ----------------------------------------------------
      # stage 2: peptide calculation (Progress: 30%)
      incProgress(0.2, detail = "Computing peptide statistics...")
      
      pep_stats <- compute_peptide_stats(
        peptide_data(),
        input$peptide_control,
        input$peptide_drug
      )
      
      # ----------------------------------------------------
      # stage 3: peptide-to-protein aggregation with Fisher's method (Progress: 60%)
      incProgress(0.3, detail = "Performing Fisher's aggregation...")
      
      fisher_res <- compute_fisher_by_topN_p_noimpute(
        pep_stats,
        input$N_value,
        input$data_mode
      )
      
      # ----------------------------------------------------
      # stage 4: result integration and ranks (Progress: 90%)
      incProgress(0.3, detail = "Finalizing results...")
      
      if (input$data_mode == "expression") {
        # Extra calculation of fold change on the protein level for the Expression/Solubility mode
        validate(need(input$protein_file, "Please upload protein file"))
        
        prot_stats <- compute_peptide_stats(
          protein_data(),
          input$protein_control,
          input$protein_drug
        )
        
        final <- left_join(
          fisher_res,
          prot_stats %>% select(Accession, log2FoldChange),
          by = "Accession"
        ) %>%
          mutate(
            score = abs(log2FoldChange) * `-log10fisher_p`,
            rank = rank(-score)
          )
        
      } else {
        # For the Partial proteolysis mode
        final <- fisher_res %>%
          mutate(
            score = abs(log2_mean_FC) * `-log10fisher_p`,
            rank = rank(-score)
          )
      }
      
      res <- final %>% arrange(rank)
      
      # Finishing
      incProgress(0.1, detail = "Done!")
      return(res)
    })
  })
  

  output$result_table <- renderDT({
    datatable(final_result(), rownames = FALSE, options = list(pageLength = 20))
  })
  
  output$download <- downloadHandler(
    filename = function() "Fisher_Result.csv",
    content = function(file) write_csv(final_result(), file)
  )
}

shinyApp(ui, server)
