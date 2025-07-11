setwd("C:/Users/Bmark/Desktop/PCR")

# Required packages
required_packages <- c(
  "shiny", "readxl", "dplyr", "tidyr",
  "writexl", "tibble", "ggplot2", "knitr"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
lapply(required_packages, library, character.only = TRUE)

ui <- fluidPage(
tags$head(
  tags$link(href = "https://fonts.googleapis.com/css2?family=Inter&display=swap", rel = "stylesheet"),
  tags$style(HTML("
    body {
      font-family: 'Inter', sans-serif;
      background-color: #e1ecf7; /* soft elegant blue */
      transition: background-color 0.3s ease;
    }

    h2 {
      font-weight: 600;
      color: #1f2937;
    }

    .sidebarPanel {
      background-color: white;
      padding: 20px;
      border-radius: 15px;
      box-shadow: 0px 4px 12px rgba(0,0,0,0.05);
      transition: box-shadow 0.3s ease;
    }

    .sidebarPanel:hover {
      box-shadow: 0px 6px 16px rgba(0,0,0,0.1);
    }

    .mainPanel {
      background-color: white;
      padding: 20px;
      border-radius: 15px;
      box-shadow: 0px 4px 12px rgba(0,0,0,0.05);
      transition: box-shadow 0.3s ease;
    }

    .mainPanel:hover {
      box-shadow: 0px 6px 16px rgba(0,0,0,0.1);
    }

    .tab-pane {
      padding-top: 20px;
    }

    .shiny-input-container {
      margin-bottom: 15px;
      transition: all 0.2s ease-in-out;
    }

    .btn {
      border-radius: 8px;
      padding: 8px 14px;
      transition: background-color 0.2s ease, transform 0.2s ease;
    }

    .btn:hover {
      background-color: #d1d5db;
      transform: scale(1.03);
    }

    .form-control {
      border-radius: 8px;
      transition: border 0.2s ease;
    }

    .form-control:focus {
      border: 2px solid #6366f1;
      outline: none;
    }

    .well {
      border-radius: 15px;
    }

    /* --- Tab Styling Glow-Up --- */
    .nav-tabs {
      border-bottom: none;
      background-color: #ffffff;
      border-radius: 12px;
      overflow: hidden;
      box-shadow: 0px 2px 6px rgba(0,0,0,0.05);
      margin-bottom: 20px;
    }

    .nav-tabs > li {
      margin-bottom: -1px;
    }

    .nav-tabs > li > a {
  background-color: transparent;
  border: none;
  color: #1f2937;
  font-weight: 500;
  border-radius: 0;
  padding: 10px 20px;
  transition: all 0.3s ease;
}

    .nav-tabs > li > a:hover {
      background-color: #dceeff;
      color: #111827;
    }

    .nav-tabs > li.active > a,
    .nav-tabs > li.active > a:focus,
    .nav-tabs > li.active > a:hover {
      background-color: #e1ecf7;
      border: none;
      color: #111827;
      font-weight: 600;
      box-shadow: inset 0 -3px 0 #6366f1;
    }
  "))
),
  tags$div(style = "display: flex; align-items: center; margin-bottom: 20px;",
           tags$img(src = "him.jpg", height = "60px", style = "margin-right: 15px; border-radius: 50%; box-shadow: 0px 2px 8px rgba(0,0,0,0.1); transition: transform 0.3s ease;"),
           tags$h2("qRTea Time", style = "margin: 0;")
  ),

  sidebarLayout(
    sidebarPanel(
  fileInput("file", "Upload raw PCR Excel file", accept = ".xlsx"),
  uiOutput("ref_selector"),
  textInput("control_sample", "Control Sample Name", value = ""),
  actionButton("analyze", "Crunch the Numbers"),

  br(), br(),

  selectInput(
    "download_file",
    "Choose a file to download:",
    choices = c(
      "Summary Table" = "summary",
      "M-values Table" = "m_values",
      "Per-Gene M-values" = "per_gene_m",
      "Delta Fold Change" = "delta_fold",
      "HTML Report" = "html_report"
    ),
    selected = "summary"
  ),
  downloadButton("download_selected", "Download Selected File"),

  hr(),


  fileInput("crisp_file", "Upload CRISPR Validation Data (CSV or XLSX)", accept = c(".xlsx", ".csv")),

  hr(),


  fileInput("mir_target_file", "Upload miR Target Prediction Data (CSV or XLSX)", accept = c(".csv", ".xls", ".xlsx")),

  br(), br(),
      
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Summary Table", tableOutput("summary_table")),
        tabPanel("Per-sample M-values", tableOutput("m_value_table")),
        tabPanel("Per-gene M-values", tableOutput("per_gene_m_table")),
        tabPanel("Delta Ct & Fold Change", tableOutput("delta_fold_table")),
        tabPanel("Fold Change Plot", plotOutput("fold_plot")),
        tabPanel("CRISPR Validation", tableOutput("crisp_table"), plotOutput("crisp_fold_plot")),
        tabPanel("miR Target Prediction", tableOutput("mir_target_table"))
      )
    )
  )
)

server <- function(input, output, session) {
  all_possible_refs <- c(
    "U6", "RNU44", "RNU48", "SNORD44", "SNORD48",
    "WARS2", "TBP", "P4HA2", "GAPDH", "miR_21", "NCEH1", "MTMR4", "LRP12", "C8ORF33", "BLMH"
  )
  
crispr_data <- reactive({
  req(input$crisp_file)
  ext <- tools::file_ext(input$crisp_file$name)
  if (ext == "csv") {
    read.csv(input$crisp_file$datapath)
  } else {
    readxl::read_excel(input$crisp_file$datapath)
  }
})

crispr_merged <- reactive({
  req(crispr_data(), delta_fold_table())
  pcr_data <- delta_fold_table()
  crispr <- crispr_data()
  
  merged <- dplyr::left_join(pcr_data, crispr, by = "TargetGene")
  merged
})

output$crisp_table <- renderTable({
  crispr_merged()
})

output$crisp_fold_plot <- renderPlot({
  df <- crispr_merged()
  req(nrow(df) > 0)
  
  ggplot(df, aes(x = FoldChange, y = CRISPR_FoldChange, label = TargetGene)) +
    geom_point(color = "#4a90e2", size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "darkgray") +
    geom_text(nudge_y = 0.2, size = 3.5) +
    theme_minimal() +
    labs(
      x = "PCR Fold Change (2^-ΔΔCt)",
      y = "CRISPR Fold Change",
      title = "CRISPR vs PCR Validation"
    )
})

mir_targets <- reactive({
  req(input$mir_target_file)
  ext <- tools::file_ext(input$mir_target_file$name)
  if (ext == "csv") {
    read.csv(input$mir_target_file$datapath)
  } else {
    readxl::read_excel(input$mir_target_file$datapath)
  }
})

mir_prediction_overlap <- reactive({
  req(mir_targets(), delta_fold_table())
  pcr_data <- delta_fold_table()
  targets <- mir_targets()

  overlap <- dplyr::inner_join(pcr_data, targets, by = "TargetGene")
  overlap
})

output$mir_target_table <- renderTable({
  mir_prediction_overlap()
})

  # PCR data reactive
  data <- reactive({
    req(input$file)
    df <- readxl::read_excel(input$file$datapath)
    df$Cq <- as.numeric(df$Cq)
    df
  })
  
  final_summary <- eventReactive(input$analyze, {
    raw_data <- data()
    
    mean_ct <- raw_data %>%
      dplyr::filter(!is.na(Cq)) %>%
      dplyr::group_by(Sample, Target) %>%
      dplyr::summarise(Mean_Ct = mean(Cq), SD_Ct = sd(Cq), .groups = "drop")
    
    cv_df <- mean_ct %>%
      dplyr::mutate(CV_percent = (SD_Ct / Mean_Ct) * 100)
    
    mean_ct_wide <- mean_ct %>%
      dplyr::select(Sample, Target, Mean_Ct) %>%
      tidyr::pivot_wider(names_from = Target, values_from = Mean_Ct)
    
    stats_wide <- cv_df %>%
      dplyr::select(Sample, Target, SD_Ct, CV_percent) %>%
      tidyr::pivot_wider(names_from = Target, values_from = c(SD_Ct, CV_percent), names_glue = "{.value}_{Target}")
    
    combined <- dplyr::left_join(mean_ct_wide, stats_wide, by = "Sample") %>%
      dplyr::rename_with(~ gsub("-", "_", .x), dplyr::everything())
    
    combined
  })
  
  output$ref_selector <- renderUI({
  req(input$file)
  df <- readxl::read_excel(input$file$datapath)
  all_possible_refs <- c(
    "U6", "RNU44", "RNU48", "SNORD44", "SNORD48",
    "WARS2", "TBP", "P4HA2", "GAPDH", "miR_21", "NCEH1", "MTMR4", "LRP12", "C8ORF33", "BLMH"
  )


  df$Cq <- as.numeric(df$Cq)
  mean_ct <- df %>%
    filter(!is.na(Cq)) %>%
    group_by(Sample, Target) %>%
    summarise(Mean_Ct = mean(Cq), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Target, values_from = Mean_Ct)

  choices <- intersect(all_possible_refs, colnames(mean_ct))

  selectInput("selected_refs", "Choose Reference Gene(s)", choices = choices,
              selected = choices[1:min(2, length(choices))], multiple = TRUE)
})
  
  m_value_table <- reactive({
    df <- final_summary()
    req(input$selected_refs)
    ref_cols <- input$selected_refs
    if (length(ref_cols) < 2) return(NULL)
    
    m_list <- list()
    for (i in 1:(length(ref_cols) - 1)) {
      for (j in (i + 1):length(ref_cols)) {
        col1 <- ref_cols[i]
        col2 <- ref_cols[j]
        m_name <- paste("M", col1, col2, sep = "_")
        m_vals <- abs(df[[col1]] - df[[col2]])
        m_list[[m_name]] <- m_vals
      }
    }
    dplyr::bind_cols(Sample = df$Sample, tibble::as_tibble(m_list))
  })
  
  per_gene_m_values <- reactive({
    df <- final_summary()
    req(df, input$selected_refs)
    ref_cols <- input$selected_refs
    if(length(ref_cols) < 2) return(NULL)
    
    Q <- df %>%
      dplyr::select(Sample, dplyr::all_of(ref_cols)) %>%
      tibble::column_to_rownames("Sample") %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ 2^(-.)))
    
    pair_sd_list <- list()
    for(i in 1:(length(ref_cols)-1)) {
      for(j in (i+1):length(ref_cols)) {
        gene1 <- ref_cols[i]
        gene2 <- ref_cols[j]
        log_ratios <- log2(Q[[gene1]] / Q[[gene2]])
        sd_val <- stats::sd(log_ratios, na.rm = TRUE)
        pair_sd_list[[paste(gene1, gene2, sep = "_vs_")]] <- sd_val
      }
    }
    
    m_values <- sapply(ref_cols, function(gene) {
      related_sds <- sapply(names(pair_sd_list), function(pair_name) {
        parts <- unlist(strsplit(pair_name, "_vs_"))
        if(gene %in% parts) pair_sd_list[[pair_name]] else NA
      })
      mean(related_sds, na.rm = TRUE)
    })
    
    tibble::tibble(Reference_Gene = ref_cols, M_value = m_values) %>%
      dplyr::arrange(M_value)
  })
  
  delta_fold_table <- reactive({
    df <- final_summary()
    control <- input$control_sample
    req(input$selected_refs)
    present_refs <- input$selected_refs
    
    gene_cols <- colnames(df)
    mean_ct_cols <- gene_cols[!grepl("^(SD_Ct|CV_percent)_", gene_cols) & gene_cols != "Sample"]
    target_genes <- setdiff(mean_ct_cols, present_refs)
    
    if (length(present_refs) == 0 || length(target_genes) == 0) return(NULL)
    if (!(control %in% df$Sample)) {
      showNotification(paste("Control sample", control, "not found!"), type = "error")
      return(NULL)
    }
    
    delta_list <- list()
    for (ref in present_refs) {
      for (target in target_genes) {
        delta_df <- df %>%
          dplyr::select(Sample, !!rlang::sym(ref), !!rlang::sym(target)) %>%
          dplyr::mutate(
            RefGene = ref,
            TargetGene = target,
            DeltaCt = !!rlang::sym(target) - !!rlang::sym(ref)
          ) %>%
          dplyr::select(Sample, RefGene, TargetGene, DeltaCt)
        
        control_delta <- delta_df %>% dplyr::filter(Sample == control) %>% dplyr::pull(DeltaCt)
        
        if (length(control_delta) == 0 || is.na(control_delta)) {
          next
        }
        
        delta_df <- delta_df %>%
          dplyr::mutate(
            DeltaDeltaCt = DeltaCt - control_delta,
            FoldChange = 2^(-DeltaDeltaCt)
          )
        
        delta_list[[paste(ref, target, sep = "_")]] <- delta_df
      }
    }
    
    dplyr::bind_rows(delta_list)
  })
  
 
  output$summary_table <- renderTable({ final_summary() })
  output$m_value_table <- renderTable({ m_value_table() })
  output$per_gene_m_table <- renderTable({ per_gene_m_values() })
  output$delta_fold_table <- renderTable({ delta_fold_table() })
  
  output$fold_plot <- renderPlot({
    df <- delta_fold_table()
    req(df)
    target_to_plot <- df$TargetGene[1]
    plot_data <- df %>% dplyr::filter(TargetGene == target_to_plot)
    
    ggplot2::ggplot(plot_data, ggplot2::aes(x = Sample, y = FoldChange, fill = RefGene)) +
      ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste("Fold Change (2^-ΔΔCt) for", target_to_plot),
        y = "Fold Change",
        x = "Sample"
      ) +
      ggplot2::theme(legend.position = "bottom")
  })
  

  output$download_selected <- downloadHandler(
    filename = function() {
      paste0(input$download_file, ifelse(input$download_file == "html_report", ".html", ".xlsx"))
    },
    content = function(file) {
      switch(input$download_file,
        summary = writexl::write_xlsx(final_summary(), file),
        m_values = writexl::write_xlsx(m_value_table(), file),
        per_gene_m = writexl::write_xlsx(per_gene_m_values(), file),
        delta_fold = writexl::write_xlsx(delta_fold_table(), file),
        html_report = {
          summary <- final_summary()
          mvals <- m_value_table()
          per_gene <- per_gene_m_values()
          foldchange <- delta_fold_table()
          
          html <- paste0(
            "<html><head><title>qPCR Report</title>",
            "<style>",
            "body { font-family: 'Segoe UI', sans-serif; margin: 20px; }",
            "h2 { color: #4a4a4a; }",
            "table { border-collapse: collapse; width: 100%; margin-bottom: 40px; }",
            "th, td { border: 1px solid #ddd; padding: 8px; text-align: center; }",
            "th { background-color: #f2f2f2; }",
            "</style></head><body>",
            "<h2>Summary Table</h2>", knitr::kable(summary, format = "html"),
            "<h2>M-values Table</h2>", knitr::kable(mvals, format = "html"),
            "<h2>Per-Gene M-values</h2>", knitr::kable(per_gene, format = "html"),
            "<h2>Delta Ct & Fold Change</h2>", knitr::kable(foldchange, format = "html"),
            "<p style='font-size:12px; color:gray;'>Generated on: ", Sys.Date(), "</p>",
            "</body></html>"
          )
          
          writeLines(html, file)
        }
      )
    }
  )
}

shiny::runApp()
