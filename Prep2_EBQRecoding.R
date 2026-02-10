# THIS IS THE SECOND PREPARATION R Script for the Clinician reactions to psychoeducation of Anxiety Study 
# by Hunt, Carpenter et al.  Pre-registration available at https://doi.org/10.17605/OSF.IO/7X98Y

#It is necessary to Recode EBQ Questionnaire Values from raw Prolific Data before running the analysis

# ---
# 1. SETUP: LOAD LIBRARIES AND DEFINE FILE PATHS
# ---

# This section loads the R packages needed for data manipulation.
library(tidyverse)
library(here)

# Define the path to the 'cleaned_data' folder
cleaned_data_path <- here("cleaned_data")

# Define the input file (the one we just created in stage one of cleaning) and the output file name
input_file <- file.path(cleaned_data_path, "cleaned_psych_data_FULL.csv")
output_file <- file.path(cleaned_data_path, "cleaned_psych_data_EBQRECODED.csv")

# ---
# 2. LOAD AND RECODE THE DATA
# ---

# Read the fully cleaned dataset

full_data <- read_csv(input_file, col_types = cols(.default = "c"))

# Recode the EBQ columns, because prolific outputted them with 1,2,3, 15,16,17 and 18, but they should be 1-7.
# 1→1, 2→2, 3→3, 15→4, 16→5, 17→6, 18→7

# We use mutate(across(...)) to apply the same logic to all columns starting with "EBQ"
# The case_when() function checks each value and replaces it according to your rules.
recoded_data <- full_data %>%
  mutate(across(starts_with("EBQ"), ~case_when(
    . == "15" ~ "4",
    . == "16" ~ "5",
    . == "17" ~ "6",
    . == "18" ~ "7",
    TRUE ~ . # This keeps all other values of 1, 2, 3 as they are.
  )))

# ---
# 3. SAVE THE FINAL, RECODED DATA
# ---

# Save the recoded data to a new file
write_csv(recoded_data, output_file)

# Print a success message to the console
cat("SUCCESS! EBQ values have been recoded. 🚀\n")
cat("The final dataset is saved at:", output_file, "\n")