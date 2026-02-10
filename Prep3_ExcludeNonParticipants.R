# THIS IS THE THIRD PREPARATION R Script for the Clinician reactions to psychoeducation of Anxiety Study 
# by Hunt, Carpenter et al.  Pre-registration available at https://doi.org/10.17605/OSF.IO/7X98Y
# Its purpose is to exclude Participants and Analyze Dropouts. It prints how many dropouts and at what
# stage after opening the questionnaire they dropped out. It prints tables and creates graphs 
# to show dropout rates per session and per education type.
# it also excludes and lists the amount of people who didn't consent 
# or who failed the attention check.

# ---
# 1. SETUP: LOAD LIBRARIES AND DEFINE FILE PATHS
# ---

# This section loads the R packages needed for data manipulation.
# If you don't have them installed, run: install.packages(c("tidyverse", "here"))
library(tidyverse)
library(here)
library(BayesFactor)

# Define the path to your 'cleaned_data' folder
cleaned_data_path <- here("cleaned_data")

# Define the input file (the one you just created) and the final output file name
input_file <- file.path(cleaned_data_path, "cleaned_psych_data_EBQRECODED.csv")
output_file <- file.path(cleaned_data_path, "final_analysis_dataset.csv")


# ---
# 2. LOAD AND PRE-PROCESS THE DATA
# ---

# Read the fully cleaned and recoded dataset
# We read all columns as character to avoid type issues initially
full_data <- read_csv(input_file, col_types = cols(.default = "c"))
  
# Convert 'Progress' column to numeric for filtering
full_data <- full_data %>%
  mutate(Progress_Progress = as.numeric(Progress_Progress))
  
# ---
# 3. SEPARATE COMPLETERS AND ANALYZE NON-COMPLETERS
# ---
  
# Dataframe of participants who did NOT complete 100%
non_completers <- full_data %>%
  filter(Progress_Progress < 100)
  
# Dataframe of participants who completed 100% of the survey
completers <- full_data %>%
  filter(Progress_Progress == 100)
  
cat(sprintf("Initial data split into %d completers and %d non-completers.\n\n", nrow(completers), nrow(non_completers)))
  


###### Analyze the non-completers by time they left they dropped out
######    Create two points in the questionnaire for dropouts to be bracketed before/after
    # Find the exact names of key columns for dropout analysis
    continue_col_name <- names(non_completers)[str_starts(names(non_completers), "password_continue")]
    post_session_col_name <- names(non_completers)[str_starts(names(non_completers), "ExAnSt21")]
    
    # Categorize non-completers based on where they dropped out
    dropout_analysis <- non_completers %>%
      mutate(DropoutStage = case_when(
        is.na(Q4_Age) ~ "1. Early Dropout (Pre-Demographics)",
        !is.na(Q4_Age) & is.na(.data[[continue_col_name]]) ~ "2. Pre-Education Dropout",
        !is.na(.data[[continue_col_name]]) & is.na(.data[[post_session_col_name]]) ~ "3. Post-Education Dropout",
        TRUE ~ "4. Late Dropout (Post-Session)"
      ))
    
    # Create and print a summary table of the dropout analysis
    dropout_summary <- dropout_analysis %>%
      count(DropoutStage, name = "Number Of Participants") %>%
      arrange(DropoutStage)
    
    cat("--- Non-Completer Dropout Analysis ---\n")
    print(dropout_summary, n = Inf)
    cat("\n")
    
    # Create and print a summary table of dropout stage by education arm
    dropout_stage_by_arm <- dropout_analysis %>%
      count(EduType, DropoutStage, name = "number_of_participants") %>%
      arrange(EduType, DropoutStage)
    
    cat("--- Non-Completer Dropout Analysis by Stage and Education Arm ---\n")
    print(dropout_stage_by_arm, n = Inf)
    cat("\n")
    
##### MAKING A HISTOGRAM OF DROPOUT STAGE BY EDUCATION ARM    
    # Create an interaction variable to map both arm and stage to a unique color
    dropout_analysis <- dropout_analysis %>%
      mutate(fill_group = interaction(EduType, DropoutStage))
    
    # Define the order of levels from bottom-to-top for the stacking.
    # This ensures "1. Early Dropout" is at the bottom of the stack.
    bottom_to_top_levels <- c(
      "Evo.1. Early Dropout (Pre-Demographics)",
      "Evo.2. Pre-Education Dropout",
      "Evo.3. Post-Education Dropout",
      "Evo.4. Late Dropout (Post-Session)",
      "Gen.1. Early Dropout (Pre-Demographics)",
      "Gen.2. Pre-Education Dropout",
      "Gen.3. Post-Education Dropout",
      "Gen.4. Late Dropout (Post-Session)"
    )
    
    # Define two distinct color palettes (4 shades for each arm)
    stage_color_palette <- c(
      # Evo Arm Palette (shades of blue, dark to light)
      `Evo.1. Early Dropout (Pre-Demographics)` = "#c6dbef",
      `Evo.2. Pre-Education Dropout`             = "#6baed6",
      `Evo.3. Post-Education Dropout`            = "#2171b5",
      `Evo.4. Late Dropout (Post-Session)`       = "#084594",
      # Gen Arm Palette (shades of orange, dark to light)
      `Gen.1. Early Dropout (Pre-Demographics)` = "#fdd0a2",
      `Gen.2. Pre-Education Dropout`             = "#fd8d3c",
      `Gen.3. Post-Education Dropout`            = "#e6550d",
      `Gen.4. Late Dropout (Post-Session)`       = "#a63603"
    )
    
    # Create the ggplot object
    # Use position_stack(reverse = TRUE) so the bottom-to-top order matches the vector above
    # and keep the legend in the same order (Early → Late).
    dropout_stage_plot <- ggplot(
      dropout_analysis,
      aes(
        x = EduType,
        fill = factor(fill_group, levels = bottom_to_top_levels)
      )
    ) +
      geom_bar(position = position_stack(reverse = TRUE)) +
      scale_fill_manual(
        name = "Arm & Dropout Stage",
        values = stage_color_palette,
        labels = function(level) gsub("^Evo\\.|^Gen\\.", "", level),
        guide = guide_legend(reverse = FALSE)
      ) +
      labs(
        title = "Dropout Stage by Education Arm",
        x = "Education Arm",
        y = "Number of Non-Completers"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right"
      )
    
    # Print the plot
    cat("--- Displaying Dropout Stage Visualization ---\n")
    print(dropout_stage_plot)
    cat("Plot has been generated. Check the 'Plots' pane.\n\n")
    
##########    
# Analyse non-completers by whether they were in a particular education arm/session
    
    dropout_by_group <- full_data %>%
      
      # Step 1: Group data by the education arm and the session ID.
      # This ensures calculations are performed on each subgroup independently.
      group_by(EduType, SessionID) %>%
      
      # Step 2: Summarise each group to get key metrics.
      # We create new columns for our calculations.
      summarise(
        total_participants = n(),
        non_completers = sum(Progress_Progress < 100, na.rm = TRUE),
        .groups = 'drop' # Best practice to drop grouping after summarising.
      ) %>%
      
      # Step 3: Calculate the dropout rate as a percentage.
      # This is done in a separate step for clarity.
      mutate(
        dropout_rate_percent = (non_completers / total_participants) * 100
      ) %>%
      
      # Step 4: Arrange the final table for easy reading.
      arrange(EduType, SessionID)
    
    # Print the detailed dropout analysis table to the console.
    cat("--- Non-Completer Dropout Analysis by Group ---\n")
    print(dropout_by_group, n = Inf)
    cat("\n")
  
##### --- Analyse non-completers by education arm without separating by session ---
    
    # This code calculates the overall dropout rate for each education arm.
    dropout_by_arm <- full_data %>%
      # Group data only by the education arm.
      group_by(EduType) %>%
      # Summarise each group to get the total counts.
      summarise(
        total_participants = n(),
        non_completers = sum(Progress_Progress < 100, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      # Calculate the overall dropout rate as a percentage.
      mutate(
        dropout_rate_percent = (non_completers / total_participants) * 100
      )
    
    # Print the overall dropout analysis table to the console.
    cat("--- Non-Completer Dropout Analysis by Education Arm (Overall) ---\n")
    print(dropout_by_arm, n = Inf)
    cat("\n")
    
    
    # --- Visualize non-completers by education arm and session (Chart) ---
    
    # First, determine the correct logical order for the SessionIDs on the x-axis
    sorted_session_levels <- full_data %>%
      distinct(SessionID) %>%
      mutate(
        alpha_part = str_extract(SessionID, "^[A-Za-z]+"),
        numeric_part = as.integer(str_extract(SessionID, "\\d+"))
      ) %>%
      arrange(alpha_part, numeric_part) %>%
      pull(SessionID)
    
    # Now, create the dataframe to be used for the plot
    plot_data <- full_data %>%
      mutate(
        status = if_else(Progress_Progress == 100, "Completer", "Non-completer"),
        SessionID = factor(SessionID, levels = sorted_session_levels),
        grouping_variable = interaction(status, EduType)
      )
    
    # Define the color palette and legend labels for clarity
    color_palette <- c(
      "Completer.Evo" = "#a9cce3",      # Light Blue
      "Non-completer.Evo" = "#2980b9",  # Dark Blue
      "Completer.Gen" = "#f5cba7",      # Light Orange
      "Non-completer.Gen" = "#e67e22"   # Dark Orange
    )
    
    legend_labels <- c(
      "Completer.Evo" = "Evo - Completer",
      "Non-completer.Evo" = "Evo - Non-completer",
      "Completer.Gen" = "Gen - Completer",
      "Non-completer.Gen" = "Gen - Non-completer"
    )
    
    # Create the ggplot object
    participant_status_plot <- ggplot(plot_data, aes(x = SessionID, fill = grouping_variable)) +
      geom_bar(position = "stack") +
      scale_fill_manual(
        name = "Education Arm & Status",
        values = color_palette,
        labels = legend_labels
      ) +
      labs(
        title = "Participant Completion Status by Session and Education Arm",
        x = "Session ID",
        y = "Number of Participants"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "top"
      )
    
    # Print the plot to the RStudio plots pane
    cat("--- Displaying Non-Completer Visualization ---\n")
    print(participant_status_plot)
    cat("Plot has been generated. Check the 'Plots' pane.\n\n")
    
    # --- END: VISUALIZATION CODE ---
    
    
  
    
    
    # ---
    # 4. BAYESIAN ANALYSIS OF DROPOUT RATES
    # ---
    # This section performs a Bayesian analysis to determine the evidence for or against
    # a difference in dropout rates between the two education arms.
    
    
    cat("--- Bayesian Analysis of Dropout Rates by Arm ---\n")
    
    # To perform the analysis, we need a 2x2 contingency table:
    # Rows: Evo, Gen
    # Columns: Non-Completers, Completers
    
    # Step 1: Create a summary table that includes completer counts.
    # We reuse the 'dropout_by_arm' table and add a 'completers' column.
    contingency_summary <- dropout_by_arm %>%
      mutate(completers = total_participants - non_completers) %>%
      # Arrange by EduType to ensure a consistent row order ("Evo", then "Gen").
      arrange(EduType)
    
    # Step 2: Convert the summary into a 2x2 matrix format required by the function.
    # We select only the count columns.
    contingency_matrix <- contingency_summary %>%
      select(non_completers, completers) %>%
      as.matrix()
    
    # Step 3: Assign the education arm names as the row names for clarity in the output.
    rownames(contingency_matrix) <- contingency_summary$EduType
    
    cat("Contingency table for analysis:\n")
    print(contingency_matrix)
    cat("\n")
    
    # Step 4: Run the Bayes Factor analysis on the contingency table.
    # 'sampleType = "indepMulti"' specifies that we have two independent groups (the arms).
    # 'fixedMargin = "rows"' tells the function that the number of participants per arm (our rows) was fixed.
    bf_result <- contingencyTableBF(contingency_matrix, sampleType = "indepMulti", fixedMargin = "rows")
    
    # Step 5: Print and interpret the results.
    cat("Bayes Factor Result:\n")
    print(bf_result)
    cat("\n")
    
    # The output gives the Bayes Factor for the alternative hypothesis (BF10).
    # A value < 1 provides evidence FOR the null hypothesis (no difference).
    # To make it more intuitive, we can calculate the inverse (BF01), which is the
    # evidence in favor of the null.
    # We must first extract the numeric value from the 'S4' object using as.vector().
    bf_null <- 1 / as.vector(bf_result)
    
    cat("--- Interpretation ---\n")
    cat(sprintf(
      "The Bayes Factor in favor of the null hypothesis (BF₀₁) is approximately %.2f.\n",
      bf_null
    ))
    cat(sprintf(
      "This means the observed data are about %.2f times more likely under the null hypothesis (that there is no difference in dropout rates between arms) than under the alternative hypothesis (that there is a difference).\n",
      bf_null
    ))
    cat("This result provides strong evidence to support the statement that dropout rates did not differ significantly between the arms.\n\n")
    
    
    
  # ---
  # 5. FILTER THE COMPLETERS BASED ON QUALITY CHECKS AND NON-CONSENT
  # ---
  
  # Find the exact name of the consent column
  consent_col_name <- names(completers)[str_starts(names(completers), "Q31_I_consent")]
  
  # Find all attention check column names
  attention_check_cols <- names(completers)[str_starts(names(completers), "Attention")]
  
  # Apply filters for consent and attention checks
  final_data <- completers %>%
    # Filter 1: Keep only those who consented (value is "1")
    filter(.data[[consent_col_name]] == "1") %>%
    # Filter 2: Keep only those who passed all attention checks (value is "2")
    # Using if_all() is the modern, recommended syntax
    filter(if_all(all_of(attention_check_cols), ~ . == "2"))
  
  # Report on the exclusions
  non_consent_removed <- nrow(completers) - nrow(filter(completers, .data[[consent_col_name]] == "1"))
  attention_fail_removed <- nrow(filter(completers, .data[[consent_col_name]] == "1")) - nrow(final_data)
  
  cat("--- Completer Quality Check ---\n")
  cat(sprintf("Removed %d completers who did not consent.\n", non_consent_removed))
  cat(sprintf("Removed %d completers who failed the attention check.\n\n", attention_fail_removed))
  
  
  # ---
  # 6. SAVE THE FINAL DATASET
  # ---
  
  # Save the clean dataset (consenting, attentive completers only) to a new file
  write_csv(final_data, output_file)
  
  cat(sprintf("SUCCESS! Final dataset with %d participants saved to:\n%s\n", nrow(final_data), output_file))
  
  # ---
  # 6. FINAL PARTICIPANT SUMMARY
  # ---
  
  # This section calculates and prints the number of participants per arm
  # and the total number of participants in the final dataset.
  
  cat("\n--- Final Participant Summary ---\n")
  
  # Use dplyr's count() to get the number of participants in each arm
  participant_summary <- final_data %>%
    count(EduType, name = "NumberOfParticipants")
  
  # Print the summary table
  print(participant_summary)
  
  # Print the total number of participants
  cat(sprintf("\nTotal participants in final dataset: %d\n", nrow(final_data)))