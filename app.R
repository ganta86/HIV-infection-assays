# =========================================================
# app.R — HIV Assay Analysis
# =========================================================

library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(writexl)
library(ggrepel)
library(stringr)
library(cowplot)

# =========================================================
# IC50 Analysis Functions
# =========================================================

dose_response_hill <- function(x, bottom, top, ic50, hill_slope) {
  bottom + (top - bottom) / (1 + (ic50 / x)^hill_slope)
}

calculate_ic50_interpolation <- function(conc, resp) {
  df <- data.frame(conc = conc, resp = resp) %>% arrange(conc)
  lower <- tail(df[df$resp < 50, ], 1)
  upper <- head(df[df$resp > 50, ], 1)
  if (nrow(lower) == 0 || nrow(upper) == 0) return(NA)
  lower$conc + (50 - lower$resp) * (upper$conc - lower$conc) / (upper$resp - lower$resp)
}

bootstrap_ic50 <- function(conc, resp, n_bootstraps = 1000) {
  df <- data.frame(conc = conc, resp = resp)
  estimates <- numeric()
  for (i in 1:n_bootstraps) {
    sample_df <- df[sample(1:nrow(df), replace = TRUE), ]
    ic50 <- calculate_ic50_interpolation(sample_df$conc, sample_df$resp)
    if (!is.na(ic50)) estimates <- c(estimates, ic50)
  }
  if (!length(estimates)) return(NULL)
  list(
    n_valid = length(estimates),
    mean = mean(estimates),
    std = sd(estimates),
    se = sd(estimates) / sqrt(length(estimates)),
    cv = 100 * sd(estimates) / mean(estimates),
    ci_lower = quantile(estimates, 0.025),
    ci_upper = quantile(estimates, 0.975),
    estimates = estimates
  )
}

fit_hill_curve <- function(conc, resp) {
  hill_equation <- function(x, p) {
    bottom <- p[1]; top <- p[2]; ic50 <- p[3]; hill_slope <- p[4]
    bottom + (top - bottom) / (1 + (ic50 / x)^hill_slope)
  }
  objective_function <- function(par, x, y) sum((y - hill_equation(x, par))^2)
  tryCatch({
    init <- c(min(resp), max(resp), median(conc), 1)
    res <- optim(
      par = init, fn = objective_function,
      x = conc, y = resp, method = "L-BFGS-B",
      lower = c(0, 50, 0.01, 0.1),
      upper = c(20, 120, 10, 10),
      control = list(maxit = 10000)
    )
    preds <- hill_equation(conc, res$par)
    r2 <- 1 - sum((resp - preds)^2) / sum((resp - mean(resp))^2)
    rmse <- sqrt(mean((resp - preds)^2))
    list(parameters = setNames(res$par, c("bottom","top","ic50","hill_slope")),
         r2 = r2, rmse = rmse, predictions = preds)
  }, error = function(e) NULL)
}

plot_dose_response <- function(df, results, virus_name) {
  df <- df %>% filter(conc > 0) %>% mutate(log_conc = log10(conc))
  main_plot <- ggplot(df, aes(x = log_conc, y = inhibition)) +
    geom_point(color = "red", size = 2.5) +
    labs(x = "Log Concentration (μg/mL)", y = "% Inhibition",
         title = paste("Soluble CD4 (sCD4) Inhibition of", virus_name)) +
    geom_hline(yintercept = 50, color = "red", linetype = "dashed", linewidth = 1) +
    theme_minimal()
  
  if (!is.null(results$hill_fit)) {
    curve_data <- data.frame(log_conc = seq(min(df$log_conc), max(df$log_conc), length.out = 100))
    curve_data$conc <- 10^curve_data$log_conc
    params <- results$hill_fit$parameters
    curve_data$inhibition <- dose_response_hill(curve_data$conc, params['bottom'], params['top'], params['ic50'], params['hill_slope'])
    main_plot <- main_plot +
      geom_line(data = curve_data, aes(x = log_conc, y = inhibition), color = "blue", linewidth = 1.2) +
      annotate("text", x = min(df$log_conc), y = 15, hjust = 0,
               label = sprintf("R² = %.4f\nRMSE = %.4f", results$hill_fit$r2, results$hill_fit$rmse))
  }
  
  if (!is.na(results$ic50_interpolation)) {
    main_plot <- main_plot +
      annotate("text", x = min(df$log_conc), y = 5, hjust = 0,
               label = sprintf("IC50 = %.4f μg/mL", results$ic50_interpolation))
  }
  
  if (!is.null(results$bootstrap)) {
    ci <- results$bootstrap
    main_plot <- main_plot +
      annotate("text", x = min(df$log_conc), y = 95, hjust = 0,
               label = sprintf("Bootstrap mean = %.4f\n95%% CI = [%.4f, %.4f]\nCV = %.2f%%",
                               ci$mean, ci$ci_lower, ci$ci_upper, ci$cv), size = 3)
    
    hist_data <- data.frame(estimate = ci$estimates)
    hist_plot <- ggplot(hist_data, aes(x = estimate)) +
      geom_histogram(color = "black", fill = "skyblue", bins = 30) +
      geom_vline(xintercept = results$ic50_interpolation, color = "darkgreen", linetype = "dashed") +
      geom_vline(xintercept = ci$ci_lower, color = "blue", linetype = "dotted") +
      geom_vline(xintercept = ci$ci_upper, color = "blue", linetype = "dotted") +
      theme_minimal(base_size = 8) +
      labs(title = "Bootstrap IC50 Distribution", x = "IC50 (μg/mL)", y = "Count")
    
    return(ggdraw() +
             draw_plot(main_plot) +
             draw_plot(hist_plot, x = 0.55, y = 0.15, width = 0.4, height = 0.25))
  } else {
    return(main_plot)
  }
}

analyze_virus <- function(conc, r1, r2, virus_name) {
  rlu_avg <- (r1 + r2) / 2
  valid <- !(is.na(conc) | is.na(rlu_avg))
  df <- data.frame(conc = conc[valid], rlu = rlu_avg[valid]) %>% arrange(conc)
  blank <- min(df$rlu, na.rm = TRUE)
  df$rlu_bg <- df$rlu - blank
  control <- max(df$rlu_bg)
  df$poc <- df$rlu_bg / control * 100
  df$inhibition <- pmin(pmax(100 - df$poc, 0), 100)
  results <- list(); notes <- ""
  if (any(df$inhibition < 50) && any(df$inhibition > 50)) {
    results$ic50_interpolation <- calculate_ic50_interpolation(df$conc, df$inhibition)
    results$bootstrap <- bootstrap_ic50(df$conc, df$inhibition)
    results$hill_fit <- fit_hill_curve(df$conc, df$inhibition)
  } else {
    results$ic50_interpolation <- NA
    notes <- if (max(df$inhibition) < 50) "Max inhibition < 50%" else "Does not cross 50%"
  }
  list(summary_row = data.frame(
    Virus = virus_name,
    `IC50 (μg/mL)` = ifelse(is.na(results$ic50_interpolation) && notes != "",
                            ">50", results$ic50_interpolation),
    Mean_Bootstrap = if (!is.null(results$bootstrap)) results$bootstrap$mean else NA,
    CI_Lower = if (!is.null(results$bootstrap)) results$bootstrap$ci_lower else NA,
    CI_Upper = if (!is.null(results$bootstrap)) results$bootstrap$ci_upper else NA,
    CV = if (!is.null(results$bootstrap)) results$bootstrap$cv else NA,
    R2 = if (!is.null(results$hill_fit)) results$hill_fit$r2 else NA,
    RMSE = if (!is.null(results$hill_fit)) results$hill_fit$rmse else NA,
    Notes = notes
  ),
  plot = plot_dose_response(df, results, virus_name))
}
# =========================================================
# UI
# =========================================================
ui <- dashboardPage(
  dashboardHeader(title = "HIV Assay Analysis"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Upload",          tabName = "upload",         icon = icon("upload")),
      menuItem("Infection Assay",      tabName = "infection",      icon = icon("virus")),
      menuItem("Neutralization Assay", tabName = "neutralization", icon = icon("vial")),
      menuItem("Download Results",     tabName = "download",       icon = icon("download"))
    )
  ),
  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper, .right-side { background-color: #f4f4f4; }
    "))),
    tabItems(
      # -----------------------------
      # Upload tab
      # -----------------------------
      tabItem(tabName = "upload",
              fluidRow(
                box(title = "Upload Data Files", status = "primary", solidHeader = TRUE, width = 12,
                    h3("Infection Assay Data"),
                    fileInput("infection_file", "Choose Excel file for Infection Assay:", accept = c(".xlsx", ".xls")),
                    hr(),
                    h3("Neutralization Assay Data"),
                    fileInput("neutralization_file", "Choose Excel file for Neutralization Assay:", accept = c(".xlsx", ".xls")),
                    br(),
                    verbatimTextOutput("upload_status")
                )
              )
      ),
      
      # -----------------------------
      # Infection assay tab
      # -----------------------------
      tabItem(tabName = "infection",
              fluidRow(
                box(title = "Process Infection Data", status = "primary", solidHeader = TRUE, width = 12,
                    actionButton("process_infection", "Process Infection Data", class = "btn btn-success btn-lg"),
                    br(), br(),
                    verbatimTextOutput("infection_status")
                )
              ),
              fluidRow(
                box(title = "Data Preview", status = "info", solidHeader = TRUE, width = 12,
                    h4("Processed Data Structure:"),
                    DT::dataTableOutput("infection_preview"))
              ),
              fluidRow(
                box(title = "FFU/mL in TZM-bl", status = "info", solidHeader = TRUE, width = 12,
                    plotOutput("tzmbl_plot", height = "500px"))
              ),
              fluidRow(
                box(title = "FFU/mL by Sample and Cell Line", status = "warning", solidHeader = TRUE, width = 12,
                    plotOutput("ffu_plot", height = "600px"))
              ),
              fluidRow(
                box(title = "% Infectivity in JC53 relative to TZM-bl", status = "success", solidHeader = TRUE, width = 12,
                    plotOutput("jc53_infectivity_plot", height = "500px"))
              ),
              fluidRow(
                box(title = "% Infectivity in JC10 relative to TZM-bl", status = "success", solidHeader = TRUE, width = 12,
                    plotOutput("jc10_infectivity_plot", height = "500px"))
              ),
              fluidRow(
                box(title = "% Infectivity in RC49 relative to TZM-bl", status = "success", solidHeader = TRUE, width = 12,
                    plotOutput("rc49_infectivity_plot", height = "500px"))
              ),
              fluidRow(
                box(title = "% Infectivity Relative to TZM-bl (All Cell Lines)", status = "info", solidHeader = TRUE, width = 12,
                    plotOutput("infectivity_plot", height = "600px"))
              )
      ),
      
      # -----------------------------
      # Neutralization assay tab
      # -----------------------------
      tabItem(tabName = "neutralization",
              fluidRow(
                box(title = "Process Neutralization Data", status = "primary", solidHeader = TRUE, width = 12,
                    actionButton("process_neutralization", "Process Neutralization Data", class = "btn btn-success btn-lg"),
                    br(), br(),
                    verbatimTextOutput("neutralization_status"))
              ),
              fluidRow(
                box(title = "Raw Data Preview", status = "info", solidHeader = TRUE, width = 12,
                    DT::dataTableOutput("neutralization_raw_preview"))
              ),
              fluidRow(
                box(title = "IC50 Summary Table", status = "warning", solidHeader = TRUE, width = 12,
                    DT::dataTableOutput("ic50_table"))
              ),
              fluidRow(
                box(title = "Dose–Response Curves", status = "info", solidHeader = TRUE, width = 12,
                    uiOutput("dose_response_plots_ui"))
              ),
              fluidRow(
                box(title = "IC50 Distribution by Virus Category", status = "success", solidHeader = TRUE, width = 12,
                    plotOutput("ic50_category_plot", height = "600px"))
              ),
              fluidRow(
                box(title = "Category Summary Statistics", status = "warning", solidHeader = TRUE, width = 12,
                    DT::dataTableOutput("category_summary_table"))
              )
      ),
      
      # -----------------------------
      # Downloads tab
      # -----------------------------
      tabItem(tabName = "download",
              fluidRow(
                box(title = "Download Results", status = "primary", solidHeader = TRUE, width = 12,
                    h3("Infection Assay Results"),
                    downloadButton("download_infection", "Download Infection Results (.xlsx)", class = "btn btn-success btn-lg"),
                    br(), br(), hr(),
                    h3("Neutralization Assay Results"),
                    downloadButton("download_neutralization", "Download IC50 Summary (.xlsx)", class = "btn btn-success btn-lg"),
                    br(), br(),
                    downloadButton("download_category_summary", "Download Category Summary (.xlsx)", class = "btn btn-info btn-lg")
                )
              )
      )
    )
  )
)
# =========================================================
# SERVER
# =========================================================
 server <- function(input, output, session) {
  values <- reactiveValues(
    long_data = NULL,
    summary_data = NULL,
    percent_infectivity = NULL,
    neutralization_raw_data = NULL,
    ic50_summary = NULL,
    dose_response_data = NULL
  )
  
  # -----------------------------
  # Upload status
  # -----------------------------
  output$upload_status <- renderText({
    status <- c()
    if (!is.null(input$infection_file))     
      status <- c(status, paste("✓ Infection file uploaded:", input$infection_file$name))
    if (!is.null(input$neutralization_file)) 
      status <- c(status, paste("✓ Neutralization file uploaded:", input$neutralization_file$name))
    if (!length(status)) return("No files uploaded yet.")
    paste(status, collapse = "\n")
  })
  
  # -----------------------------
  # Infection assay processing
  # -----------------------------
  observeEvent(input$process_infection, {
    req(input$infection_file)
    tryCatch({
      df_raw <- read_excel(input$infection_file$datapath, col_names = FALSE)
      
      # Header rows
      cell_lines_raw <- unlist(df_raw[1, -1], use.names = FALSE)
      virus_ids_raw  <- unlist(df_raw[2, -1], use.names = FALSE)
      volumes_raw    <- as.numeric(unlist(df_raw[3, -1], use.names = FALSE))
      
      # Fill blanks in merged headers
      for (i in 2:length(cell_lines_raw)) if (is.na(cell_lines_raw[i])) cell_lines_raw[i] <- cell_lines_raw[i - 1]
      for (i in 2:length(virus_ids_raw))  if (is.na(virus_ids_raw[i]))  virus_ids_raw[i]  <- virus_ids_raw[i - 1]
      
      ncols    <- length(volumes_raw)
      cell_lines <- cell_lines_raw[seq_len(ncols)]
      virus_ids  <- virus_ids_raw[seq_len(ncols)]
      volumes    <- volumes_raw
      
      # Data body
      ffu_values <- df_raw[-c(1:3), ]
      colnames(ffu_values) <- c("Sample", paste0(cell_lines, "_", virus_ids, "_", volumes))
      
      # Long format
      long_data <- ffu_values %>%
        pivot_longer(cols = -Sample, names_to = "Cell_Virus_Volume", values_to = "FFU") %>%
        separate(Cell_Virus_Volume, into = c("CellLine", "Virus", "Volume"), sep = "_") %>%
        mutate(
          FFU    = as.numeric(FFU),
          Volume = as.numeric(Volume),
          FFU_per_mL = FFU * (1000 / Volume)
        )
      
      summary_data <- long_data %>%
        group_by(Sample, CellLine, Virus) %>%
        summarise(
          Mean_FFU_per_mL = mean(FFU_per_mL, na.rm = TRUE),
          SD_FFU_per_mL   = sd(FFU_per_mL,   na.rm = TRUE),
          .groups = "drop"
        )
      
      ref_tbl <- summary_data %>%
        filter(CellLine == "TZMBL") %>%
        select(Sample, Virus, Ref_FFU = Mean_FFU_per_mL)
      
      percent_infectivity <- summary_data %>%
        left_join(ref_tbl, by = c("Sample","Virus")) %>%
        mutate(Percent_Infectivity = 100 * Mean_FFU_per_mL / Ref_FFU)
      
      values$long_data          <- long_data
      values$summary_data       <- summary_data
      values$percent_infectivity<- percent_infectivity
      
      showNotification("Infection data processed successfully!", type = "success")
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  output$infection_status <- renderText({
    if (is.null(values$summary_data)) return("Click 'Process Infection Data' to analyze the uploaded file.")
    paste(
      "✓ Data processed successfully",
      paste0("\n• Samples: ", length(unique(values$summary_data$Sample))),
      paste0("\n• Cell lines: ", paste(unique(values$summary_data$CellLine), collapse = ", ")),
      paste0("\n• Volumes detected: ", paste(sort(unique(values$long_data$Volume)), collapse = ", ")),
      paste0("\n• Total data points: ", nrow(values$long_data))
    )
  })
  
  output$infection_preview <- DT::renderDataTable({
    req(values$long_data)
    preview_data <- values$long_data %>%
      select(Sample, CellLine, Virus, Volume, FFU, FFU_per_mL) %>%
      mutate(FFU_per_mL = round(FFU_per_mL, 2)) %>%
      head(100)
    datatable(preview_data, options = list(pageLength = 10, scrollX = TRUE),
              caption = "First 100 rows of processed data.")
  })
  output$tzmbl_plot <- renderPlot({
    req(values$long_data)
    tzmbl <- values$long_data %>% filter(CellLine == "TZMBL")
    ggplot(tzmbl, aes(x = Sample, y = FFU_per_mL)) +
      geom_bar(stat = "summary", fun = mean, fill = "forestgreen", color = "black") +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
      geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
      facet_wrap(~ Virus, scales = "free_y") +
      labs(title = "FFU/mL in TZM-bl", y = "FFU/mL", x = "Sample") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$ffu_plot <- renderPlot({
    req(values$summary_data, values$long_data)
    ggplot(values$summary_data,
           aes(x = Sample, y = Mean_FFU_per_mL, fill = CellLine)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
      geom_errorbar(aes(ymin = Mean_FFU_per_mL - SD_FFU_per_mL,
                        ymax = Mean_FFU_per_mL + SD_FFU_per_mL),
                    width = 0.2,
                    position = position_dodge(width = 0.9)) +
      geom_point(data = values$long_data,
                 aes(x = Sample, y = FFU_per_mL, color = CellLine, group = CellLine),
                 position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9),
                 size = 2, alpha = 0.7, inherit.aes = FALSE) +
      facet_wrap(~ CellLine, scales = "free_y") +
      labs(title = "FFU/mL by Sample and Cell Line", y = "FFU/mL", x = "Sample") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      guides(fill = "none", color = "none")
  })
  
  compute_rel_infectivity <- function(target_cell) {
    jd <- values$long_data %>% filter(CellLine == target_cell)
    ref <- values$long_data %>%
      filter(CellLine == "TZMBL") %>%
      group_by(Sample, Virus) %>%
      summarise(Ref_FFU = mean(FFU_per_mL, na.rm = TRUE), .groups = "drop")
    jd %>%
      left_join(ref, by = c("Sample","Virus")) %>%
      mutate(Percent_Infectivity = 100 * FFU_per_mL / Ref_FFU)
  }
  
  output$jc53_infectivity_plot <- renderPlot({
    req(values$long_data); if (!"JC53" %in% values$long_data$CellLine) return(NULL)
    dat <- compute_rel_infectivity("JC53")
    ggplot(dat, aes(x = Sample, y = Percent_Infectivity)) +
      geom_bar(stat = "summary", fun = mean, fill = "darkorange", color = "black") +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
      geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
      geom_hline(yintercept = 10, linetype = "dotted", color = "blue", size = 1) +
      geom_hline(yintercept = 20, linetype = "dotted", color = "gray40", size = 1) +
      geom_hline(yintercept = 30, linetype = "dotted", color = "red", size = 1) +
      labs(title = "% Infectivity in JC53 relative to TZM-bl", y = "% Infectivity", x = "Sample") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$jc10_infectivity_plot <- renderPlot({
    req(values$long_data); if (!"JC10" %in% values$long_data$CellLine) return(NULL)
    dat <- compute_rel_infectivity("JC10")
    ggplot(dat, aes(x = Sample, y = Percent_Infectivity)) +
      geom_bar(stat = "summary", fun = mean, fill = "skyblue", color = "black") +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
      geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
      geom_hline(yintercept = 10, linetype = "dotted", color = "blue", size = 1) +
      geom_hline(yintercept = 20, linetype = "dotted", color = "gray40", size = 1) +
      geom_hline(yintercept = 30, linetype = "dotted", color = "red", size = 1) +
      labs(title = "% Infectivity in JC10 relative to TZM-bl", y = "% Infectivity", x = "Sample") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$rc49_infectivity_plot <- renderPlot({
    req(values$long_data); if (!"RC49" %in% values$long_data$CellLine) return(NULL)
    dat <- compute_rel_infectivity("RC49")
    ggplot(dat, aes(x = Sample, y = Percent_Infectivity)) +
      geom_bar(stat = "summary", fun = mean, fill = "lightgreen", color = "black") +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
      geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
      geom_hline(yintercept = 10, linetype = "dotted", color = "blue", size = 1) +
      geom_hline(yintercept = 20, linetype = "dotted", color = "gray40", size = 1) +
      geom_hline(yintercept = 30, linetype = "dotted", color = "red", size = 1) +
      labs(title = "% Infectivity in RC49 relative to TZM-bl", y = "% Infectivity", x = "Sample") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$infectivity_plot <- renderPlot({
    req(values$percent_infectivity, values$long_data)
    raw_infectivity <- values$long_data %>%
      filter(CellLine != "TZMBL") %>%
      left_join(
        values$long_data %>%
          filter(CellLine == "TZMBL") %>%
          group_by(Sample, Virus) %>%
          summarise(Ref_FFU = mean(FFU_per_mL, na.rm = TRUE), .groups = "drop"),
        by = c("Sample","Virus")
      ) %>%
      mutate(Percent_Infectivity = 100 * FFU_per_mL / Ref_FFU)
    
    ggplot(raw_infectivity, aes(x = Sample, y = Percent_Infectivity)) +
      geom_bar(data = values$percent_infectivity %>% filter(CellLine != "TZMBL"),
               aes(y = Percent_Infectivity, fill = CellLine),
               stat = "identity", position = position_dodge(width = 0.9)) +
      geom_point(aes(color = CellLine, group = CellLine),
                 position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9),
                 size = 2, alpha = 0.7) +
      geom_hline(yintercept = 10, linetype = "dotted", color = "blue",  size = 1) +
      geom_hline(yintercept = 20, linetype = "dotted", color = "gray40", size = 1) +
      geom_hline(yintercept = 30, linetype = "dotted", color = "red",   size = 1) +
      facet_wrap(~ CellLine, scales = "free_y") +
      labs(title = "% Infectivity Relative to TZMBL", y = "% Infectivity", x = "Sample") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      guides(fill = "none", color = "none")
  })
  # -----------------------------
  # Neutralization assay processing
  # -----------------------------
  observeEvent(input$process_neutralization, {
    req(input$neutralization_file)
    tryCatch({
      sheet <- read_excel(input$neutralization_file$datapath)
      values$neutralization_raw_data <- sheet
      
      concentration <- as.numeric(sheet[-1, 1][[1]])
      virus_cols <- colnames(sheet)[-1]
      
      summary_df <- data.frame()
      plot_list  <- list()
      
      # Each virus: two adjacent columns for replicates
      for (i in seq(1, length(virus_cols), by = 2)) {
        col1 <- virus_cols[i]
        col2 <- if (i + 1 <= length(virus_cols)) virus_cols[i + 1] else NA
        if (is.na(col2)) next
        
        virus_name <- sub("\\.1$", "", col1)
        r1 <- as.numeric(sheet[-1, col1][[1]])
        r2 <- as.numeric(sheet[-1, col2][[1]])
        
        res <- analyze_virus(concentration, r1, r2, virus_name)
        summary_df <- bind_rows(summary_df, res$summary_row)
        plot_list[[virus_name]] <- res$plot
      }
      
      values$ic50_summary <- summary_df
      values$dose_response_data <- plot_list
      
      showNotification(paste("Neutralization data processed successfully! Analyzed",
                             nrow(summary_df), "viruses."), type = "success")
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  output$neutralization_status <- renderText({
    if (is.null(values$ic50_summary)) return("Click 'Process Neutralization Data' to analyze the uploaded file.")
    sum_success <- sum(!grepl("^>", values$ic50_summary$`IC50 (μg/mL)`))
    sum_failed  <- sum(grepl("^>",  values$ic50_summary$`IC50 (μg/mL)`))
    paste(
      "✓ Data processed successfully",
      paste0("\n• Viruses analyzed: ", nrow(values$ic50_summary)),
      paste0("\n• Successful IC50 calculations: ", sum_success),
      paste0("\n• Failed fits: ", sum_failed)
    )
  })
  
  output$neutralization_raw_preview <- DT::renderDataTable({
    req(values$neutralization_raw_data)
    datatable(values$neutralization_raw_data, options = list(pageLength = 10, scrollX = TRUE),
              caption = "Raw neutralization data from Excel file")
  })
  
  # Render list of per-virus plots
  output$dose_response_plots_ui <- renderUI({
    req(values$dose_response_data)
    lapply(names(values$dose_response_data), function(vn) {
      plotOutput(outputId = paste0("plot_", vn), height = "600px")
    })
  })
  
  observe({
    req(values$dose_response_data)
    for (vn in names(values$dose_response_data)) {
      local({
        virus_name <- vn
        output[[paste0("plot_", virus_name)]] <- renderPlot({
          values$dose_response_data[[virus_name]]
        })
      })
    }
  })
  
  # -----------------------------
  # IC50 summary table
  # -----------------------------
  output$ic50_table <- DT::renderDataTable({
    req(values$ic50_summary)
    datatable(values$ic50_summary, options = list(pageLength = 20, scrollX = TRUE))
  })
  
  # -----------------------------
  # Dynamic IC50 plot by virus
  # -----------------------------
  output$ic50_category_plot <- renderPlot({
    req(values$ic50_summary)
    df <- values$ic50_summary
    
    # Detect IC50 column dynamically
    ic50_col <- grep("IC50", names(df), value = TRUE, ignore.case = TRUE)[1]
    if (is.na(ic50_col) || ic50_col == "") {
      stop("No IC50 column found in data.")
    }
    
    # Convert IC50 to numeric, replace '>50' with 50
    df$IC50_numeric <- as.numeric(gsub("^>", "", df[[ic50_col]]))
    df$IC50_numeric[grepl("^>", df[[ic50_col]])] <- 50
    
    # Use Virus name directly as category
    df$Category <- as.factor(df$Virus)
    
    ggplot(df, aes(x = Category, y = IC50_numeric, color = Category)) +
      geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.9) +
      geom_text_repel(aes(label = Virus), size = 3.5, box.padding = 0.5, max.overlaps = 100) +
      geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray") +
      geom_hline(yintercept = 5.0, linetype = "dashed", color = "red", linewidth = 1) +
      annotate("text", x = n_distinct(df$Category), y = 1.5, label = "IC50 = 1 μg/mL", hjust = 1, size = 3.5) +
      annotate("text", x = n_distinct(df$Category), y = 5.5, label = "IC50 = 5 μg/mL", hjust = 1, size = 4, color = "red") +
      labs(title = "IC50 Distribution by Virus",
           x = "Virus", y = "IC50 (μg/mL)") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # -----------------------------
  # Dynamic category summary table
  # -----------------------------
  output$category_summary_table <- DT::renderDataTable({
    req(values$ic50_summary)
    df <- values$ic50_summary
    
    # Detect IC50 column dynamically
    ic50_col <- grep("IC50", names(df), value = TRUE, ignore.case = TRUE)[1]
    if (is.na(ic50_col) || ic50_col == "") {
      stop("No IC50 column found in data.")
    }
    
    # Convert IC50 to numeric, replace '>50' with 50
    df$IC50_numeric <- as.numeric(gsub("^>", "", df[[ic50_col]]))
    df$IC50_numeric[grepl("^>", df[[ic50_col]])] <- 50
    
    summary_tbl <- df %>%
      group_by(Virus) %>%
      summarise(
        Count = n(),
        Mean_IC50 = mean(IC50_numeric, na.rm = TRUE),
        Median_IC50 = median(IC50_numeric, na.rm = TRUE),
        SD = sd(IC50_numeric, na.rm = TRUE),
        Min_IC50 = min(IC50_numeric, na.rm = TRUE),
        Max_IC50 = max(IC50_numeric, na.rm = TRUE),
        .groups = "drop"
      )
    
    datatable(summary_tbl, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # -----------------------------
  # Download handlers
  # -----------------------------
  output$download_infection <- downloadHandler(
    filename = function() paste0("infection_results_", Sys.Date(), ".xlsx"),
    content  = function(file) {
      sheets <- list(
        "Raw_Data"            = values$long_data,
        "Summary"             = values$summary_data,
        "Percent_Infectivity" = values$percent_infectivity
      )
      write_xlsx(sheets, file)
    }
  )
  
  output$download_neutralization <- downloadHandler(
    filename = function() paste0("ic50_summary_", Sys.Date(), ".xlsx"),
    content  = function(file) {
      if (!is.null(values$ic50_summary)) {
        write_xlsx(values$ic50_summary, file)
      }
    }
  )
  
  output$download_category_summary <- downloadHandler(
    filename = function() paste0("IC50_category_summary_", Sys.Date(), ".xlsx"),
    content  = function(file) {
      if (!is.null(values$ic50_summary)) {
        df <- values$ic50_summary
        
        # Detect IC50 column dynamically
        ic50_col <- grep("IC50", names(df), value = TRUE, ignore.case = TRUE)[1]
        if (is.na(ic50_col) || ic50_col == "") {
          stop("IC50 column not found in summary data.")
        }
        
        # Convert to numeric, replace '>50' with 50
        df$IC50_numeric <- as.numeric(gsub("^>", "", df[[ic50_col]]))
        df$IC50_numeric[grepl("^>", df[[ic50_col]])] <- 50
        
        category_summary <- df %>%
          group_by(Virus) %>%
          summarise(
            Count       = n(),
            Mean_IC50   = mean(IC50_numeric, na.rm = TRUE),
            Median_IC50 = median(IC50_numeric, na.rm = TRUE),
            SD          = sd(IC50_numeric, na.rm = TRUE),
            Min_IC50    = min(IC50_numeric, na.rm = TRUE),
            Max_IC50    = max(IC50_numeric, na.rm = TRUE),
            .groups     = "drop"
          )
        
        write_xlsx(category_summary, file)
      }
    }
  )
 }
  # =============================
  # Run the Shiny app
  # =============================
  shinyApp(ui = ui, server = server)
  
