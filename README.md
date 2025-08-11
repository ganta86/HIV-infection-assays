# HIV-infection-assays
An interactive Shiny app for analyzing HIV infection and neutralization assay data. Features include data upload, automated processing of FFU/mL and IC50 calculations, bootstrapping, curve fitting, summary tables, and downloadable results.
# HIV Assay Analysis Shiny App

This repository contains a Shiny application for interactive analysis of HIV infection and neutralization assay data, designed for researchers working with data from FFU/mL and IC50 determinations.

# HIV Infection Assays - KL LAB - Dr GANTA Shiny App

This Shiny app allows users to perform analysis and visualization for HIV infection assays and neutralization assays (IC50 calculations). It supports file upload, dynamic virus and cell line selection, FFU conversion, IC50 curve fitting, and summary statistics/plots—all without requiring coding experience.

---

## 📁 File Structure
You can place this project in any folder, e.g., your Desktop. The folder should contain:

```
HIV_Infection_Assays_App/
  ├── app.R                 # Unified Shiny app (UI + server logic)
  ├── report_template.Rmd   # For downloadable PDF reports
  ├── www/                  # Optional: include CSS/logo/images
```

---

## ✅ Features
- Upload infection data or IC50 assay Excel files
- Select experiment type (Infection or Neutralization)
- Choose cell lines (e.g., TZMBL, Macrophage)
- Choose viruses (e.g., B33, B59, JRFL, LN8, LN10, JRCSF)
- View FFU/mL and % infectivity calculations
- Perform IC50 analysis with Hill fit and bootstrapped CI
- Generate summary IC50 plots with category-based coloring
- Interactive UI elements: tabs, sliders for dilution factors, virus/cell line filters
- Download PDF report with plots and results

---

## ▶️ Running the App
```R
library(shiny)
runApp("HIV_Infection_Assays_App")
```

You only need to run `app.R`. All logic is self-contained.

---

## 📄 report_template.Rmd (included in this folder)
```rmd
---
title: "HIV Assay Summary Report"
output: pdf_document
params:
  assay_type: ""
  virus: ""
  cell_line: ""
  summary_table: NULL
  summary_plot: NULL
---

```{r, echo=FALSE, message=FALSE}
cat("### Assay Type: ", params$assay_type, "\n")
cat("### Virus: ", params$virus, "\n")
cat("### Cell Line: ", params$cell_line, "\n")
```

```{r, echo=FALSE, results='asis'}
knitr::kable(params$summary_table)
```

```{r, echo=FALSE}
print(params$summary_plot)
```
```

---

## 📌 Instructions to Add Analysis Code

1. The **Infection Assay Module** is integrated with FFU calculation, dilution choice, and infectivity %.
2. The **Neutralization Assay Module** includes:
   - `fit_hill_curve`, `calculate_ic50_interpolation`, and `bootstrap_ic50`
   - `analyze_virus()` and `plot_dose_response()`
   - IC50 summary code for plotting and exporting results
3. All code is already embedded in `app.R` (no file name is hardcoded; user upload handles this).
4. The Hill fitting code is fully implemented as in the original R script.
5. The output of FFU input is made dynamic; you upload the FFU input Excel file via the UI.
6. You can download IC50 summary and reports via button clicks in the app.

---

## 🛠 Dependencies
Make sure these R packages are installed:
```R
install.packages(c("shiny", "readxl", "ggplot2", "dplyr", "tidyr", "writexl", "ggrepel", "cowplot", "openxlsx", "rmarkdown"))
```

Let me know when you’re ready to extend with viral load plots, plate layouts, or batch processing!
**For questions or suggestions, open an issue or contact [ganta86](https://github.com/ganta86).**
