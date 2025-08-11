# HIV-infection-assays
An interactive Shiny app for analyzing HIV infection and neutralization assay data. Features include data upload, automated processing of FFU/mL and IC50 calculations, bootstrapping, curve fitting, summary tables, and downloadable results.
# HIV Assay Analysis Shiny App

This repository contains a Shiny application for interactive analysis of HIV infection and neutralization assay data, designed for researchers working with data from FFU/mL and IC50 determinations.

## Features

- Upload and process infection assay data and neutralization assay data from Excel files.
- Automated calculation and plotting of FFU/mL for different cell lines.
- IC50 calculation with bootstrapping and nonlinear curve fitting.
- Downloadable results and summary tables.
- Interactive visualizations for quality control and exploration.

## Requirements

- **R** (version ≥ 4.0.0 recommended)
- **RStudio** (optional, but recommended for ease of use)

### Required R packages

This app uses the following R packages. Install them before running the app:

```r
install.packages(c(
  "shiny", "shinydashboard", "DT", "ggplot2", "dplyr",
  "tidyr", "readxl", "writexl", "ggrepel", "stringr", "cowplot"
))
```

## Usage

1. **Clone this repository:**

   ```sh
   git clone https://github.com/ganta86/my-shiny-app.git
   cd my-shiny-app
   ```

2. **Open R or RStudio in this directory.**

3. **Install the required packages** (see above).

4. **Run the app:**

   ```r
   shiny::runApp()
   ```

5. **Follow the app instructions to upload your Excel files and process your data.**

## Notes

- The app expects Excel files with a specific format for infection and neutralization assays. See the app’s upload page for details.
- Downloaded results will be in `.xlsx` format.

## License

[MIT](LICENSE) (or specify your preferred license)

---

**For questions or suggestions, open an issue or contact [ganta86](https://github.com/ganta86).**
