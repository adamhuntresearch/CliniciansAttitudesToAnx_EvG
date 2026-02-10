# THIS IS AN EXPLORATORY ANALYSIS R Script for the Clinician reactions to psychoeducation of Anxiety Study 
# by Hunt, Carpenter et al.  Pre-registration available at https://doi.org/10.17605/OSF.IO/7X98Y
#
## It firstly looks at demographic data, including participant age, specialisation, and stage of training/position.

#-----------------------------------------------------------------------------
# STEP 1: SETUP - LOAD LIBRARIES AND DATA
#-----------------------------------------------------------------------------

library(tidyverse)

# Load the dataset.
file_path <- "FINALLY_analysis_ready_data.csv"
psy_data <- read_csv(file_path)

#-----------------------------------------------------------------------------
# STEP 2: CLEAN AND PREPARE THE AGE DATA
#-----------------------------------------------------------------------------
# The 'Q4_Age' column contains non-numeric characters and needs to be cleaned
# before it can be used for calculations.

age_data <- psy_data %>%
  mutate(
    # Use str_extract() to pull out only the numbers from the 'Q4_Age' column.
    # This correctly handles entries like "33years" by extracting "33".
    Age_Clean = as.numeric(str_extract(Q4_Age, "\\d+"))
  ) %>%
  # Remove rows where the cleaned age is NA (missing).
  filter(!is.na(Age_Clean))

#-----------------------------------------------------------------------------
# STEP 3: CALCULATE DESCRIPTIVE AGE STATISTICS
#-----------------------------------------------------------------------------
# Now that the data is clean, we can calculate summary statistics.

# Calculate the mean age
mean_age <- mean(age_data$Age_Clean)

# Calculate the age range
min_age <- min(age_data$Age_Clean)
max_age <- max(age_data$Age_Clean)

# Print the results to the console
print(paste("Mean Age:", round(mean_age, 1)))
print(paste("Age Range:", min_age, "-", max_age))

#-----------------------------------------------------------------------------
# STEP 4: CREATE AND ANALYZE AGE BRACKETS
#-----------------------------------------------------------------------------
# Grouping ages into brackets helps to visualize the distribution of the sample.
# Based on the typical distribution of clinicians, 10-year brackets are sensible.

age_brackets <- age_data %>%
  mutate(
    # Use the cut() function to create age brackets.
    # 'right = TRUE' means the interval is closed on the right (e.g., (20,30] includes 30).
    Age_Bracket = cut(Age_Clean, 
                      breaks = c(20, 29, 39, 49, 59, 69, 79),
                      labels = c("20-29", "30-39", "40-49", "50-59", "60-69", "69+"),
                      right = TRUE)
  )

# Count the number of participants in each bracket
bracket_counts <- age_brackets %>%
  count(Age_Bracket)

# Print the distribution table
print("Distribution of Participants by Age Bracket:")
print(bracket_counts)

# Create a simple bar plot to visualize this
ggplot(bracket_counts, aes(x = Age_Bracket, y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Participant Distribution by Age",
       x = "Age Bracket",
       y = "Number of Participants")


#-----------------------------------------------------------------------------
# STEP 5: ANALYZE COMPREHENSION RATINGS
#-----------------------------------------------------------------------------
# This section focuses on the comprehension question, which asks participants 
# to rate their understanding of the psychoeducational presentation they viewed.
# We first calculate descriptive statistics for each group (Evo vs. Gen)
# and then perform a Bayesian analysis to compare them.

# Select and prepare the data, using the full column name for clarity and
# accuracy. We rename the column to 'comprehension_rating' to make the
# subsequent code cleaner and more readable.
comprehension_data <- psy_data %>%
  select(
    EduType, 
    comprehension_rating = `Q38_How_would_you_rate_your_comprehension_of_the_presentation`
  ) %>%
  # Remove any rows where the rating is missing
  filter(!is.na(comprehension_rating)) %>%
  # Ensure the column is treated as a number for calculations
  mutate(comprehension_rating = as.numeric(comprehension_rating))

#-----------------------------------------------------------------------------
# STEP 5A: DESCRIPTIVE STATISTICS FOR COMPREHENSION
#-----------------------------------------------------------------------------

# Calculate summary statistics (count, mean, sd, median) grouped by education type
comprehension_summary <- comprehension_data %>%
  group_by(EduType) %>%
  summarise(
    count = n(),
    mean_rating = mean(comprehension_rating, na.rm = TRUE),
    sd_rating = sd(comprehension_rating, na.rm = TRUE),
    median_rating = median(comprehension_rating, na.rm = TRUE)
  )

# Print the resulting summary table
print("Summary of Comprehension Ratings by Group:")
print(comprehension_summary)

#-----------------------------------------------------------------------------
# STEP 5B: BAYESIAN ANALYSIS OF COMPREHENSION RATINGS
#-----------------------------------------------------------------------------
# To determine if there's a substantial difference in comprehension ratings
# between the 'Evo' and 'Gen' groups, a Bayesian approach is useful. It allows
# us to estimate the magnitude of the difference and our uncertainty in that
# estimate, which is often more informative than a simple p-value from a
# frequentist t-test.

# We will use a Bayesian regression model. This model will estimate the 
# comprehension rating based on which group the participant was in.
# The 'brms' package provides a user-friendly interface for this type of modeling.

# Define and run the Bayesian model.
# The formula now uses our cleaned column name, 'comprehension_rating'.
# brms uses sensible, weakly informative priors by default, which is a robust
# approach for many analyses. We set a seed for reproducibility.
comprehension_model <- brm(
  formula = comprehension_rating ~ EduType,
  data = comprehension_data,
  family = gaussian(), # Assumes the rating data is approximately normally distributed
  seed = 1234 # Ensures the analysis is reproducible
)

# Print the summary of the model results.
# The key section to examine is 'Population-Level Effects'.
# - The 'Intercept' represents the estimated mean rating for the baseline group (Gen).
# - The 'EduTypeEvo' estimate is the main parameter of interest. It represents 
#   the estimated difference in mean rating for the Evo group compared to the Gen group.
# - The '95% CrI' (Credible Interval) provides a range of plausible values for 
#   this difference. If this interval does not overlap substantially with zero, 
#   it suggests a credible difference between the groups.
print("Summary of the Bayesian Model:")
summary(comprehension_model)

# --- How to Interpret the Bayesian Model Output ---
# In the summary output, look for the 'EduTypeEvo' row.
# 1.  Estimate: This is the model's best guess for the average difference 
#    in rating points (Evo - Gen). A positive value means the Evo group's
#    average rating was higher.
# 2.  l-95% CrI & u-95% CrI: These are the lower and upper bounds of the 95%
#    credible interval. For example, if the interval is [0.2, 1.0], we can be 
#    95% certain that the true difference in means is between 0.2 and 1.0.
#    Since this interval is entirely positive, it provides strong evidence
#    that the Evo group rated their comprehension higher.
# 3.  If the interval contained zero (e.g., [-0.4, 0.5]), it would suggest
#    that there is no strong evidence for a meaningful difference between the groups.

