# THIS IS THE SECOND PRIMARY ANALYSIS R Script for the Clinician reactions to psychoeducation of Anxiety Study 
# by Hunt, Carpenter et al.  Pre-registration available at https://doi.org/10.17605/OSF.IO/7X98Y
#
# This script runs Bayesian cumulative-link ordinal regression models
# using the 'brms' package to test hypotheses H4 and H5.
#
# There are three types of models; logit and probit (which are essentially identical 
# but show odds or effect size, respectively)
# and lastly, monotonic, which is an altered version of incorporating the pre-score, used as a robustness check
#
# The script then proceeds to draw figures of the posterior probability and raw scores.

#-----------------------------------------------------------------------------
# STEP 1: SETUP - LOAD LIBRARIES AND DATA
#-----------------------------------------------------------------------------

# Load necessary R packages.
# 'brms' is for Bayesian regression modeling.
# 'tidyverse' is a collection of packages for data manipulation and visualization,
# particularly 'dplyr' for data wrangling and 'readr' for reading the CSV.

library(dplyr)
library(brms)
library(tidyverse)

# Set a seed for reproducibility. This ensures that the MCMC sampling process
# will produce the exact same results every time the script is run.
set.seed(1234)

# Load the dataset from the CSV file.
# Make sure the file "FINALLY_analysis_ready_data.csv" is in your R working directory,
# or provide the full file path.
file_path <- "FINALLY_analysis_ready_data.csv"
psy_data <- read_csv(file_path)

#-----------------------------------------------------------------------------
# STEP 2: DATA PREPARATION
#-----------------------------------------------------------------------------

# The dependent variables are responses to questions Q88 and Q90.
# These are on a Likert scale and should be treated as ordered factors for
# cumulative (ordinal) regression. We also need to ensure our predictors
# are correctly formatted as factors.

# Create a new, cleaned-up data frame for the analysis.
analysis_data <- psy_data %>%
  mutate(
    # Use backticks (`) around the names because they contain spaces or special characters.
    Q88_ord = factor(`Q88_Do_you_think_this_information_is_useful_for_patients`, ordered = TRUE, levels = 1:5),
    Q90_ord = factor(`Q90_Do_you_think_this_information_is_useful_for_clinicians`, ordered = TRUE, levels = 1:5),
    
    # Convert the predictors to factors.
    # 'EduType' is the main experimental condition.
    # 'Delivery' is the context (live vs. virtual).
    # 'SessionID' is the random effect grouping variable.
    EduType = factor(EduType, levels = c("Gen", "Evo")),,
    Delivery = as.factor(Delivery),
    SessionID = as.factor(SessionID)
  )


#-----------------------------------------------------------------------------
# STEP 3: HYPOTHESES 5 & 6 - PRIMARY ANALYSIS (LOGIT AND PROBIT MODELS)
#-----------------------------------------------------------------------------
# Here, we run the primary models for H5 (Q88) and H6 (Q90).
# We use a cumulative family with a 'logit' link. The logit link models the
# log-odds of moving to a higher category on the Likert scale.
# The coefficients (estimates) will be in log-odds ratios. Taking exp() of a
# coefficient gives the odds ratio (OR).
# OR > 1: Indicates that an increase in the predictor increases the odds of a
#         higher rating.
# OR < 1: Indicates that an increase in the predictor decreases the odds of a
#         higher rating.


# Define the set of regularizing priors
WR_priors <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 1.5), class = "Intercept")
)

# --- Model for H5 (Question 88) ---
# The formula is defined directly inside the brm call for clarity and to avoid evaluation errors.
# Formula: Q88_ord ~ 1 + EduType + Delivery + (1 | SessionID)
# 1: Intercept term.
# EduType: Main predictor of interest.
# Delivery: Covariate.
# (1 | SessionID): Random intercept for SessionID.

h5_logit_model <- brm(
  formula = bf(Q88_ord ~ 1 + EduType + Delivery + (1 | SessionID)),
  data = analysis_data,
  family = cumulative(link = "logit"),
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

saveRDS(h5_logit_model, file = "h5_logit_model.rds")

h5_logit_model_WRprior <- brm(
  formula = bf(Q88_ord ~ 1 + EduType + Delivery + (1 | SessionID)),
  data = analysis_data,
  family = cumulative(link = "logit"),
  prior = WR_priors,
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

saveRDS(h5_logit_model_WRprior, file = "h5_logit_model_WRprior.rds")


# --- Probit Model for H5 (Question 88) ---
h5_probit_model <- brm(
  formula = bf(Q88_ord ~ 1 + EduType + Delivery + (1 | SessionID)),
  data = analysis_data,
  family = cumulative(link = "probit"),
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

saveRDS(h5_probit_model, file = "h5_probit_model.rds")

h5_probit_model_WRprior <- brm(
  formula = bf(Q88_ord ~ 1 + EduType + Delivery + (1 | SessionID)),
  data = analysis_data,
  family = cumulative(link = "probit"),
  prior = WR_priors,
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

saveRDS(h5_probit_model_WRprior, file = "h5_probit_model_WRprior.rds")


# --- Model for H6 (Question 90) ---
# We use the same model structure for the second dependent variable, 'Q90_ord'.
h6_logit_model <- brm(
  formula = bf(Q90_ord ~ 1 + EduType + Delivery + (1 | SessionID)),
  data = analysis_data,
  family = cumulative(link = "logit"),
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

saveRDS(h6_logit_model, file = "h6_logit_model.rds")

h6_logit_model_WRprior <- brm(
  formula = bf(Q90_ord ~ 1 + EduType + Delivery + (1 | SessionID)),
  data = analysis_data,
  family = cumulative(link = "logit"),
  prior = WR_priors,
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

saveRDS(h6_logit_model_WRprior, file = "h6_logit_model_WRprior.rds")


# --- Probit Model for H6 (Question 90) ---
h6_probit_model <- brm(
  formula = bf(Q90_ord ~ 1 + EduType + Delivery + (1 | SessionID)),
  data = analysis_data,
  family = cumulative(link = "probit"),
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

saveRDS(h6_probit_model, file = "h6_probit_model.rds")

h6_probit_model_WRprior <- brm(
  formula = bf(Q90_ord ~ 1 + EduType + Delivery + (1 | SessionID)),
  data = analysis_data,
  family = cumulative(link = "probit"),
  prior = WR_priors,
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)


saveRDS(h6_probit_model_WRprior, file = "h6_probit_model_WRprior.rds")


#-----------------------------------------------------------------------------
# STEP 4: LOAD MODELS AND VIEW RESULTS
#-----------------------------------------------------------------------------

# Load the saved model objects from disk. This avoids re-running the models.
h5_logit_model <- readRDS("h5_logit_model.rds")
h5_logit_model_WRprior <- readRDS("h5_logit_model_WRprior.rds")
h5_probit_model <- readRDS("h5_probit_model.rds")
h5_probit_model_WRprior <- readRDS("h5_probit_model_WRprior.rds")

h6_logit_model <- readRDS("h6_logit_model.rds")
h6_logit_model_WRprior <- readRDS("h6_logit_model_WRprior.rds")
h6_probit_model <- readRDS("h6_probit_model.rds")
h6_probit_model_WRprior <- readRDS("h6_probit_model_WRprior.rds")


# --- View Primary Model Results ---
# The summary provides the posterior distribution for each parameter, including
# the mean estimate, standard error, and 95% credible intervals (CIs). It also 
# calculates the percentage of the posterior above 0

# --- H5 (Q88) Logit Model ---
print("--- H5 (Q88) Logit Model Summary ---")
summary(h5_logit_model)

# Extract draws into a data frame
draws_h5_logit <- as_draws_df(h5_logit_model)

# Use the data frame created by as_draws_df() for calculations
# Calculate and print the percentage of posterior draws > 0
cat("Percentage of posterior distribution for b_EduTypeEvo > 0:",
    round(mean(draws_h5_logit$b_EduTypeEvo > 0) * 100, 2), "%\n\n")


# --- H5 (Q88) Logit Model with Priors ---
print("--- H5 (Q88) Logit Model with Priors Summary ---")
summary(h5_logit_model_WRprior)

# Extract draws and print the percentage > 0
draws_h5_logit_wr <- as_draws_df(h5_logit_model_WRprior)
cat("Percentage of posterior distribution for b_EduTypeEvo > 0:",
    round(mean(draws_h5_logit_wr$b_EduTypeEvo > 0) * 100, 2), "%\n\n")


# --- H6 (Q90) Logit Model ---
print("--- H6 (Q90) Logit Model Summary ---")
summary(h6_logit_model)

# Extract draws and print the percentage > 0
draws_h6_logit <- as_draws_df(h6_logit_model)
cat("Percentage of posterior distribution for b_EduTypeEvo > 0:",
    round(mean(draws_h6_logit$b_EduTypeEvo > 0) * 100, 2), "%\n\n")


# --- H6 (Q90) Logit Model with Priors ---
print("--- H6 (Q90) Logit Model with Priors Summary ---")
summary(h6_logit_model_WRprior)

# Extract draws and print the percentage > 0
draws_h6_logit_wr <- as_draws_df(h6_logit_model_WRprior)
cat("Percentage of posterior distribution for b_EduTypeEvo > 0:",
    round(mean(draws_h6_logit_wr$b_EduTypeEvo > 0) * 100, 2), "%\n\n")


# --- View Probit Model Results ---
# --- H5 (Q88) Probit Model ---
print("--- H5 (Q88) Probit Model Summary (Effect Size) ---")
summary(h5_probit_model)

# Extract draws and calculate percentage > 0 for probit model
draws_h5_probit <- as_draws_df(h5_probit_model)
cat("Percentage of posterior distribution for b_EduTypeEvo > 0:",
    round(mean(draws_h5_probit$b_EduTypeEvo > 0) * 100, 2), "%\n\n")


# --- H5 (Q88) Probit Model with Priors ---
print("--- H5 (Q88) Probit Model with Priors Summary ---")
summary(h5_probit_model_WRprior)

# Extract draws and calculate percentage > 0 for probit model with priors
draws_h5_probit_wr <- as_draws_df(h5_probit_model_WRprior)
cat("Percentage of posterior distribution for b_EduTypeEvo > 0:",
    round(mean(draws_h5_probit_wr$b_EduTypeEvo > 0) * 100, 2), "%\n\n")


# --- H6 (Q90) Probit Model ---
print("--- H6 (Q90) Probit Model Summary (Effect Size) ---")
summary(h6_probit_model)

# Extract draws and calculate percentage > 0 for probit model
draws_h6_probit <- as_draws_df(h6_probit_model)
cat("Percentage of posterior distribution for b_EduTypeEvo > 0:",
    round(mean(draws_h6_probit$b_EduTypeEvo > 0) * 100, 2), "%\n\n")


# --- H6 (Q90) Probit Model with Priors ---
print("--- H6 (Q90) Probit Model with Priors Summary ---")
summary(h6_probit_model_WRprior)

# Extract draws and calculate percentage > 0 for probit model with priors
draws_h6_probit_wr <- as_draws_df(h6_probit_model_WRprior)
cat("Percentage of posterior distribution for b_EduTypeEvo > 0:",
    round(mean(draws_h6_probit_wr$b_EduTypeEvo > 0) * 100, 2), "%\n\n")

#-----------------------------------------------------------------------------
# STEP 5: DESCRIPTIVE STATISTICS
#-----------------------------------------------------------------------------

#SIMPLE MEAN LIKERT SCORES BETWEEN EVO AND GENGROUPS
mean_scores <- analysis_data %>%
  group_by(EduType) %>%
  summarise(
    mean_useful_for_patients = mean(as.numeric(Q88_ord), na.rm = TRUE),
    mean_useful_for_clinicians = mean(as.numeric(Q90_ord), na.rm = TRUE)
  )

# Print the resulting table
print(mean_scores)




################################################################################
# STEP 5 (Analysis of H5 and H6): POSTERIOR DISTRIBUTION PLOTS (ODDS RATIO SCALE)
################################################################################
#
# This script creates posterior distribution plots for the main logit models
# for H5 and H6. It converts the log-odds coefficients into odds ratios
# (by using exp()) to make the x-axis more interpretable.
#
# Finally, it uses the 'patchwork' package to combine the two plots into a
# single figure for a comprehensive overview.
#

# Ensure necessary libraries are loaded
library(brms)
library(tidyverse)
library(ggplot2)
library(patchwork)

# --- H5: Posterior Distribution Plot for Odds Ratio ---

# Step 1: Extract and convert posterior draws for H5
# We use the logit model with weakly regularizing priors.
log_odds_draws_h5 <- as_draws_df(h5_logit_model_WRprior)$b_EduTypeEvo
odds_ratio_draws_h5 <- exp(log_odds_draws_h5)

# Step 2: Calculate density data for H5
density_data_h5 <- density(odds_ratio_draws_h5)
density_df_h5 <- data.frame(x = density_data_h5$x, y = density_data_h5$y)

# Step 3: Calculate probability for H5
prob_positive_h5 <- mean(odds_ratio_draws_h5 > 1) * 100

# Step 4: Build and SAVE the plot for H5 to a variable
plot_h5 <- ggplot(density_df_h5, aes(x = x, y = y)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = subset(density_df_h5, x > 1), 
              aes(ymin = 0, ymax = y), 
              fill = "#D55E00", alpha = 0.8) +
  geom_ribbon(data = subset(density_df_h5, x < 1), 
              aes(ymin = 0, ymax = y), 
              fill = "#D55E00", alpha = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +
  labs(
    title = "Usefulness for Patients (H5)",
    subtitle = paste0("Probability of positive effect: ", round(prob_positive_h5, 2), "%"),
    x = "Odds Ratio",
    y = "Density"
  ) +
  # --- Constrain the x-axis to zoom in on the main density area ---
  coord_cartesian(xlim = c(0, 15)) +
  theme_minimal(base_size = 15)


# --- H6: Posterior Distribution Plot for Odds Ratio ---

# Step 1: Extract and convert posterior draws for H6
log_odds_draws_h6 <- as_draws_df(h6_logit_model_WRprior)$b_EduTypeEvo
odds_ratio_draws_h6 <- exp(log_odds_draws_h6)

# Step 2: Calculate density data for H6
density_data_h6 <- density(odds_ratio_draws_h6)
density_df_h6 <- data.frame(x = density_data_h6$x, y = density_data_h6$y)

# Step 3: Calculate probability for H6
prob_positive_h6 <- mean(odds_ratio_draws_h6 > 1) * 100

# Step 4: Build and SAVE the plot for H6 to a variable
plot_h6 <- ggplot(density_df_h6, aes(x = x, y = y)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = subset(density_df_h6, x > 1), 
              aes(ymin = 0, ymax = y), 
              fill = "#D55E00", alpha = 0.8) +
  geom_ribbon(data = subset(density_df_h6, x < 1), 
              aes(ymin = 0, ymax = y), 
              fill = "#D55E00", alpha = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +
  labs(
    title = "Usefulness for Clinicians (H6)",
    subtitle = paste0("Probability of positive effect: ", round(prob_positive_h6, 2), "%"),
    x = "Odds Ratio",
    y = "Density"
  ) +
  theme_minimal(base_size = 15)


# --- COMBINE THE TWO PLOTS ---

# Arrange the two saved plots side-by-side and add an overall title.
(plot_h5 + plot_h6) + 
  plot_annotation(
    title = 'Posterior Distributions for Perceived Usefulness of Evo. vs. Gen. Education',
    theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))
  )



########                                     #########
########    JITTER PLOTS FOR H5 & H6         #########
########                                     #########

# --- STEP 1: PREPARE DATA FOR BOTH H5 AND H6 ---

# We will prepare the data for both plots at once to be more efficient.
# We select the relevant columns, rename them, and pivot them into a single
# long-format dataframe with a new column to identify the hypothesis.
plot_data_h5h6_long <- psy_data %>%
  # Select the group variable and the two outcome variables
  select(
    EduType,
    H5_Patient_Usefulness = `Q88_Do_you_think_this_information_is_useful_for_patients`,
    H6_Clinician_Usefulness = `Q90_Do_you_think_this_information_is_useful_for_clinicians`
  ) %>%
  # Reshape the data from wide to long format
  tidyr::pivot_longer(
    cols = c(H5_Patient_Usefulness, H6_Clinician_Usefulness),
    names_to = "hypothesis",
    values_to = "score"
  ) %>%
  # Clean up the hypothesis names for plot titles
  mutate(
    hypothesis = gsub("_", " ", hypothesis),
    # Ensure EduType is a factor with the correct labels
    EduType = factor(EduType, levels = c("Gen", "Evo"), labels = c("Genetic", "Evolutionary"))
  )

# --- STEP 2: CALCULATE SUMMARY STATISTICS (MEAN AND 95% CI) ---

summary_data_h5h6 <- plot_data_h5h6_long %>%
  # Group by both hypothesis and education type
  group_by(hypothesis, EduType) %>%
  summarise(
    n = n(),
    mean_score = mean(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    se_score = sd_score / sqrt(n),
    ci_lower = mean_score - 1.96 * se_score,
    ci_upper = mean_score + 1.96 * se_score,
    .groups = 'drop' # Ungroup the data
  )

# --- STEP 3: BUILD THE PLOTS ---

# -- Plot for H5 --
plot_h5 <- ggplot(
  # Filter the data for H5
  data = filter(plot_data_h5h6_long, grepl("H5", hypothesis)), 
  aes(x = EduType, y = score, color = EduType)
) +
  # Add jittered raw data points
  geom_jitter(width = 0.2, alpha = 0.3, size = 2) +
  # Add points and error bars for the means, using the summary data
  geom_point(data = filter(summary_data_h5h6, grepl("H5", hypothesis)), aes(y = mean_score), size = 6) +
  geom_errorbar(
    data = filter(summary_data_h5h6, grepl("H5", hypothesis)),
    aes(ymin = ci_lower, ymax = ci_upper, y = mean_score),
    width = 0.1, linewidth = 1
  ) +
  labs(
    title = "Usefulness for Patients (H5)",
    y = "Usefulness Rating (1-5)",
    x = NULL # Remove x-axis title
  ) +
  scale_color_manual(values = c("Genetic" = "#0072B2", "Evolutionary" = "#D55E00")) +
  coord_cartesian(ylim = c(1, 5)) +
  theme_minimal(base_size = 15)

# -- Plot for H6 --
plot_h6 <- ggplot(
  # Filter the data for H6
  data = filter(plot_data_h5h6_long, grepl("H6", hypothesis)), 
  aes(x = EduType, y = score, color = EduType)
) +
  # Add jittered raw data points
  geom_jitter(width = 0.2, alpha = 0.3, size = 2) +
  # Add points and error bars for the means
  geom_point(data = filter(summary_data_h5h6, grepl("H6", hypothesis)), aes(y = mean_score), size = 6) +
  geom_errorbar(
    data = filter(summary_data_h5h6, grepl("H6", hypothesis)),
    aes(ymin = ci_lower, ymax = ci_upper, y = mean_score),
    width = 0.1, linewidth = 1
  ) +
  labs(
    title = "Usefulness for Clinicians (H6)",
    y = NULL, # Remove y-axis title for cleaner look
    x = NULL
  ) +
  scale_color_manual(values = c("Genetic" = "#0072B2", "Evolutionary" = "#D55E00")) +
  coord_cartesian(ylim = c(1, 5)) +
  theme_minimal(base_size = 15)

# --- STEP 4: COMBINE PLOTS WITH PATCHWORK ---

# Place plots side-by-side, add a shared legend and an overall title
plot_h5 + plot_h6 +
  plot_layout(guides = 'collect') &
  plot_annotation(
    title = 'Perceived Usefulness of Education Content',
    theme = theme(plot.title = element_text(size = 20, hjust = 0.5))
  ) &
  theme(legend.position = 'bottom')


