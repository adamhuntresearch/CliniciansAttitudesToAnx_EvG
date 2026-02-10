# This is the first cleaning script for the Clinician reactions to psychoeducation of Anxiety Study 
# by Hunt, Carpenter et al. Pre-registration available at https://doi.org/10.17605/OSF.IO/7X98Y

# This cleaning script involves three functions which clean the many Prolific raw files and combine them 
# into a single datasheet. 
# The first function renames all the columns.
# The second function checks no names are duplicated and then deletes unnecessary columns to anonymise.
# The third function adds additional variables such as whether Evo/Gen, sessionID, participant ID, live/virtual
# The fourth stage combines all the functions and runs them on the raw_data folder.
# This cleaning script will be followed by a second 'EBQ recoding script', because 
# Qualtrics output incorrect numbers for the EBQ.
#
# After these first and second cleaning and prep scripts, analysis can start on the single dataframe.
# Note that non-consenters and non-completers have not yet been removed from the data, and will need to be dealt with.
#
# ---
# 1. SETUP: LOAD LIBRARIES AND DEFINE FILE PATHS
# ---

# This section loads the R packages needed for data manipulation.
# If you don't have them installed, run: install.packages(c("tidyverse", "here"))
library(tidyverse)
library(here)

# Define the paths for your data directories.
# The 'here()' function makes the script's paths relative to your R project root,
# which helps with reproducibility.
# Please create a 'raw_data' folder in your project directory and place all the CSV files there.
raw_data_path <- here("raw_data")
# This is where the final, clean dataset will be saved.
cleaned_data_path <- here("cleaned_data")

# Create the output directory if it doesn't exist yet.
if (!dir.exists(cleaned_data_path)) {
  dir.create(cleaned_data_path)
}

## FUNCTION 0: Check that no subject data lines are duplicated within 
## the renamed/manually edited files - run this first, checking ResponseID


check_duplicate_responses <- function() {
  
  # The path is specified directly inside the function.
  raw_data_path <- "./raw_data/"
  
  # Get a list of all the CSV files from the specified folder.
  file_paths <- list.files(
    path = raw_data_path,
    pattern = "^(Evo|Gen)_\\d+\\.csv$",
    full.names = TRUE
  )
  
  # Read the "ResponseId" column from each file and combine them.
  all_ids <- do.call(c, lapply(file_paths, function(file) {
    header <- read.csv(file, nrows = 1, header = FALSE, as.is = TRUE)
    response_id_col_index <- which(header == "ResponseId")
    
    if (length(response_id_col_index) > 0) {
      read.csv(file, skip = 3, header = FALSE, stringsAsFactors = FALSE)[[response_id_col_index]]
    } else {
      NULL
    }
  }))
  
  # Identify any ResponseIds that appear more than once.
  duplicated_ids <- unique(all_ids[duplicated(all_ids)])
  
  # If duplicates are found, stop. Otherwise, print a success message.
  if (length(duplicated_ids) > 0) {
    stop(paste("Execution halted: Duplicate ResponseIds detected ->", 
               paste(duplicated_ids, collapse = ", ")))
  } else {
    message("No duplicates found! Cleaning can continue!")
  }
}

## FUNCTION 1: RENAME COLUMNS
# This function reads a raw Qualtrics CSV file, combines the first two header rows
# into a single, clean column name for each column, and returns the data frame.

rename_qualtrics_columns <- function(file_path) {
  
  # Read the first two rows for header info
  header_row1 <- read_csv(file_path, n_max = 1, col_names = FALSE, show_col_types = FALSE)
  header_row2 <- read_csv(file_path, skip = 1, n_max = 1, col_names = FALSE, show_col_types = FALSE)
  
  # Read the main data, skipping the 3 header rows.
  # IMPORTANT: We force all columns to be read as 'character' type to prevent
  # conflicts when binding rows later.
  main_data <- read_csv(
    file_path, 
    skip = 3, 
    col_names = FALSE, 
    show_col_types = FALSE,
    col_types = cols(.default = "c") 
  )
  
  # Create and clean the new column names
  new_colnames <- paste(header_row1[1, ], header_row2[1, ], sep = "_")
  new_colnames <- gsub("\\s+", "_", new_colnames)
  new_colnames <- gsub("[^A-Za-z0-9_]", "", new_colnames)
  new_colnames[is.na(new_colnames)] <- "NA_Column"
  
  # Assign the new names
  colnames(main_data) <- new_colnames
  
  return(main_data)
}


## FUNCTION 2: DELETE UNNECESSARY COLUMNS
# This function takes the data frame and removes columns based on their starting text.
delete_unnecessary_columns <- function(data) {
  
  # Use a single starts_with() call to remove all specified columns.
  cleaned_data <- data %>%
    select(
      -starts_with(c(
        "StartDate", "EndDate", "Status", "RecordedDate", "ResponseId",
        "DistributionChannel", "UserLanguage", "Consent_detail",
        "Q33", "Q24", "Q32", "Q2", "Q3_", "Q28", "Finished"
      ))
    )
  
  return(cleaned_data)
}


## FUNCTION 3: ADD SESSION LABELS
# This function adds EduType, a unique SessionID, Delivery mode (live = 1, virtual =0), and ParticipantID.
add_session_labels <- function(data, file_path) {
  
  # Extract the base file name (e.g., "Evo_1") to use as the unique SessionID
  session_id_unique <- tools::file_path_sans_ext(basename(file_path))
  
  # Define which sessions were "Live"
  live_sessions <- c("Evo_11", "Evo_10", "Evo_7", "Evo_4", "Gen_3", "Gen_7" )
  
  labeled_data <- data %>%
    mutate(
      EduType = str_extract(session_id_unique, "^[A-Za-z]+"),
      SessionID = session_id_unique,
      Delivery = if_else(session_id_unique %in% live_sessions, 1, 0),
      .before = 1 # Adds these new columns at the beginning of the data frame
    ) %>%
    # Add a unique participant ID for each row (e.g., "Evo_1-1", "Evo_1-2")
    mutate(
      ParticipantID = paste(SessionID, row_number(), sep = "-"),
      .after = Delivery
    )
  
  return(labeled_data)
}


# ---
# 4. PROCESS ALL FILES AND COMBINE
# ---

# Get a list of all the CSV files from your raw data folder.
all_files <- list.files(
  path = raw_data_path,
  pattern = "^(Evo|Gen)_\\d+\\.csv$",
  full.names = TRUE
)

# Check if any files were found before trying to process them
if (length(all_files) > 0) {
  
  # Use 'lapply' to run the cleaning pipeline on every file.
  # This creates a list where each item is a fully cleaned data frame.
  list_of_cleaned_data <- lapply(all_files, function(file) {
    
    # The pipeline for each file:
    # 1. Rename columns
    # 2. Delete unnecessary columns
    # 3. Add session labels
    cleaned_df <- rename_qualtrics_columns(file) %>%
      delete_unnecessary_columns() %>%
      add_session_labels(file)
    
    return(cleaned_df)
  })
  
  # Combine the list of data frames into one single, master data frame.
  final_combined_data <- bind_rows(list_of_cleaned_data)
  
  # Save the final dataset to a new CSV file.
  output_filename <- file.path(cleaned_data_path, "cleaned_psych_data_FULL.csv")
  write_csv(final_combined_data, output_filename)
  
  cat("SUCCESS! All files have been processed and combined. 🚀\n")
  cat("The final dataset is saved at:", output_filename, "\n")
  cat("Total participants processed:", nrow(final_combined_data), "\n")
  
} else {
  cat("No data files found in the 'raw_data' directory. Please check the path and file names.\n")
}
