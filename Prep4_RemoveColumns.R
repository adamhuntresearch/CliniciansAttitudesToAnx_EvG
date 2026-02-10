# THIS IS THE THIRD PREPARATION R Script for the Clinician reactions to psychoeducation of Anxiety Study 
# by Hunt, Carpenter et al.  Pre-registration available at https://doi.org/10.17605/OSF.IO/7X98Y
# This is a simple R Script to Remove Final Unnecessary Columns of completion, duration, and consent.
# After this, the data is ready to analyse!

# ---
# 1. SETUP: LOAD LIBRARIES AND DEFINE FILE PATHS
# ---

# This section loads the R packages needed for data manipulation.
# If you don't have them installed, run: install.packages(c("tidyverse", "here"))
library(tidyverse)
library(here)

# Define the path to your 'cleaned_data' folder
cleaned_data_path <- here("cleaned_data")

# Define the input file and the final, analysis-ready output file name
input_file <- file.path(cleaned_data_path, "final_analysis_dataset.csv")
output_file <- file.path(cleaned_data_path, "FINALLY_analysis_ready_data.csv")


# ---
# 2. LOAD DATA AND REMOVE FINAL COLUMNS
# ---

# Check if the input file exists before trying to read it
if (file.exists(input_file)) {
  
  # Read the dataset
  final_data <- read_csv(input_file, col_types = cols(.default = "c"))
  
  # Remove the three specified columns
  # We use starts_with() for the long consent column name for simplicity
  analysis_ready_data <- final_data %>%
    select(
      -Progress_Progress,
      -Duration_in_seconds_Duration_in_seconds,
      -starts_with("Q31_I_consent")
    )
  
  # ---
  # 3. SAVE THE FINAL, ANALYSIS-READY DATA
  # ---
  
  # Save the analysis-ready data to a new file
  write_csv(analysis_ready_data, output_file)
  
  # Print a success message to the console
  cat("SUCCESS! Final columns removed. Your data is ready for analysis! 🎉\n")
  cat("The final dataset is saved at:", output_file, "\n")
  
} else {
  cat("ERROR: The input file was not found at:", input_file, "\n")
  cat("Please make sure the 'final_analysis_dataset.csv' file is in your 'cleaned_data' folder.\n")
}
