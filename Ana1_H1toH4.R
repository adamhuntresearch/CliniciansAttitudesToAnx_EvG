# THIS IS THE FIRST PRIMARY ANALYSIS R Script for the Clinician reactions to psychoeducation of Anxiety Study
# by Hunt, Carpenter et al. Pre-registration available at https://doi.org/10.17605/OSF.IO/7X98Y
#
# This script runs Bayesian cumulative-link ordinal regression models
# using the 'brms' package to test hypotheses H1, H2, H3, and H4.
# The models predict the post-treatment score while controlling for
# the standardized pre-treatment score.
#
# There are three types of models; logit and probit (which are essentially identical
# but show odds ratio or effect size on the latent variable scale, respectively)
# and lastly, monotonic, which is a principled version of using the pre-score predictor as an ordinal variable, used as a robustness check
#
# It then proceeds to draw figures of the posterior probability and raw scores.

#-----------------------------------------------------------------------------
# SETUP - LOAD LIBRARIES AND DATA
#-----------------------------------------------------------------------------

# Load necessary R packages.
# 'brms' is for Bayesian regression modeling.
# 'tidyverse' is a collection of packages for data manipulation and visualization.
# 'ggplot2' is for data visualisation


library(brms)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)
library(ggridges)

# Set a seed for reproducibility.
set.seed(1234)

# Load the dataset.
file_path <- "FINALLY_analysis_ready_data.csv"
psy_data <- read_csv(file_path)



#-----------------------------------------------------------------------------
# STEP 1: DEMOGRAPHIC ANALYSIS - AGE AND SEX DISTRIBUTION
#-----------------------------------------------------------------------------
# This section analyzes the distribution of participant age and sex across
# the two experimental arms (Evo vs. Gen) to check for baseline differences.

# ---
# 2.1 Age Analysis
# ---

# Step 2.1.1: Clean the Age data
# The 'Q4_Age' column contains non-numeric values. We convert it to numeric,
# which will coerce any non-numeric text into NA.
psy_data_cleaned <- psy_data %>%
  mutate(Q4_Age_Clean = suppressWarnings(as.numeric(Q4_Age)))

# Report how many participants had invalid age entries
age_removed_count <- sum(is.na(psy_data_cleaned$Q4_Age_Clean) & !is.na(psy_data_cleaned$Q4_Age))
cat(sprintf("INFO: Removed %d participants with non-numeric age entries.\n", age_removed_count))

# Filter out NA values for the analysis and visualization
age_analysis_data <- psy_data_cleaned %>%
  filter(!is.na(Q4_Age_Clean))

# Step 2.1.2: Visualize Age Distribution
age_dist_plot <- ggplot(age_analysis_data, aes(x = Q4_Age_Clean, fill = EduType)) +
  geom_density(alpha = 0.6, color = "black", linewidth = 0.2) +
  scale_fill_brewer(palette = "Set1", name = "Arm") +
  labs(
    title = "Distribution of Participant Age by Arm",
    x = "Age",
    y = "Density"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(age_dist_plot)
cat("SUCCESS: Age distribution plot generated. Check the 'Plots' pane.\n\n")

# Step 2.1.3: Run Bayesian Estimation Model (T-test equivalent)
cat("Running Bayesian model for age... (This may take a moment)\n")
age_model <- brm(
  formula = Q4_Age_Clean ~ EduType,
  data = age_analysis_data,
  family = gaussian(),
  seed = 1234,
  silent = 2, refresh = 0 # Suppress verbose Stan compilation output
)

# Step 2.1.4: Print and Interpret Age Model Results
cat("\n--- Age Model Summary ---\n")
print(summary(age_model))
age_summary <- summary(age_model)
age_diff_estimate <- age_summary$fixed["EduTypeGen", "Estimate"]
age_diff_lower <- age_summary$fixed["EduTypeGen", "l-95% CI"]
age_diff_upper <- age_summary$fixed["EduTypeGen", "u-95% CI"]

cat("\n--- Age Interpretation ---\n")
cat(sprintf(
  "The estimated mean difference in age (Genetics - Evolutionary) is %.2f years.\n",
  age_diff_estimate
))
cat(sprintf(
  "The 95%% Credible Interval for this difference is [%.2f, %.2f].\n\n",
  age_diff_lower,
  age_diff_upper
))


# ---
# 2.2 Sex/Gender Analysis (Odds Ratio Approach)
# ---

cat("--- 2.2 Analyzing Sex Distribution (Odds Ratio Approach) ---\n")

# Step 2.2.1: Prepare data for logistic regression
# To get a single Odds Ratio, we must create a binary outcome variable.
# We will compare the largest group (Female) against all other groups combined.
sex_model_data <- psy_data_cleaned %>%
  mutate(Is_Female = if_else(Q7_Gender == "2", 1, 0)) %>%
  filter(!is.na(Q7_Gender)) # Exclude participants with missing gender data

# Step 2.2.2: Run Bayesian Logistic Regression Model
cat("Running Bayesian logistic regression for sex... (This may take a moment)\n")
sex_model <- brm(
  formula = Is_Female ~ EduType,
  data = sex_model_data,
  family = bernoulli(link = "logit"), # Use bernoulli for 0/1 outcomes
  seed = 1234,
  silent = 2, refresh = 0
)

# Step 2.2.3: Extract, Convert, and Interpret the Odds Ratio
cat("\n--- Sex Model Summary (Log-Odds Scale) ---\n")
print(summary(sex_model))

# The model output is on the log-odds scale. We need to exponentiate it
# to get the more interpretable Odds Ratio.
sex_summary <- summary(sex_model)
log_odds_estimate <- sex_summary$fixed["EduTypeGen", "Estimate"]
log_odds_lower <- sex_summary$fixed["EduTypeGen", "l-95% CI"]
log_odds_upper <- sex_summary$fixed["EduTypeGen", "u-95% CI"]

# Convert to Odds Ratio
or_estimate <- exp(log_odds_estimate)
or_lower <- exp(log_odds_lower)
or_upper <- exp(log_odds_upper)

cat("\n--- Sex Interpretation (Odds Ratio Scale) ---\n")
cat(sprintf(
  "The Odds Ratio for a participant being female (Genetics vs. Evolutionary) is %.2f.\n",
  or_estimate
))
cat(sprintf(
  "The 95%% Credible Interval for this Odds Ratio is [%.2f, %.2f].\n",
  or_lower,
  or_upper
))
cat("For an Odds Ratio, the 'no effect' value is 1.0.\n\n")




#-----------------------------------------------------------------------------
# STEP 2: DATA PREPARATION
#-----------------------------------------------------------------------------
# This step involves standardizing the pre-scores, converting post-scores to
# ordered factors, and ensuring all predictors are correctly formatted.

analysis_data <- psy_data %>%
  mutate(
    # --- Standardize (z-score) the pre-score predictors ---
    # The scale() function centers the data (subtracts the mean) and scales it
    # (divides by the standard deviation), creating the z-score.
    H1_pre_z = scale(`ExAnSt1_Given_your_current_understanding_how_optimistic_are_you_about_improvement_in_symptoms_in_patients_presenting_with_anxiety_disorder`),
    H2_pre_z = scale(`ExAnSt2_Given_your_current_understanding_how_willing_do_you_think_patients_would_be_to_share_anxiety_disorder_diagnoses`),
    H3_pre_z = scale(`ExAnSt8_Given_current_public_knowledge_of_anxiety_how_willing_do_you_think_people_are_to_seek_psychiatric_help_for_anxiety`),
    H4_pre_z = scale(`ExAnSt4_Given_your_current_understanding_how_effective_do_you_think_psychosocial_interventions_will_be_in_improving_outcomes_in_anxiety_disorder`),
    
    
    # --- Create ordered factor versions of pre-scores for monotonic models ---
    H1_pre_ord = factor(`ExAnSt1_Given_your_current_understanding_how_optimistic_are_you_about_improvement_in_symptoms_in_patients_presenting_with_anxiety_disorder`, levels = 1:7, ordered = TRUE),
    H2_pre_ord = factor(`ExAnSt2_Given_your_current_understanding_how_willing_do_you_think_patients_would_be_to_share_anxiety_disorder_diagnoses`, levels = 1:7, ordered = TRUE),
    H3_pre_ord = factor(`ExAnSt8_Given_current_public_knowledge_of_anxiety_how_willing_do_you_think_people_are_to_seek_psychiatric_help_for_anxiety`, levels = 1:7, ordered = TRUE),
    H4_pre_ord = factor(`ExAnSt4_Given_your_current_understanding_how_effective_do_you_think_psychosocial_interventions_will_be_in_improving_outcomes_in_anxiety_disorder`, levels = 1:7, ordered = TRUE),
    
    
    # --- Renaming Post Score Variables as the column names ar super long ---
    H1_post_ord = ExAnSt21_Given_the_information_provided_in_this_education_session_how_optimistic_are_you_about_improvement_in_symptoms_in_patients_presenting_with_anxiety_disorder,
    H2_post_ord = ExAnSt22_Given_the_information_provided_in_this_education_session_how_willing_do_you_think_patients_would_be_to_share_anxiety_disorder_diagnoses,
    H3_post_ord = ExAnSt28_If_the_information_provided_in_this_education_session_was_publicly_known_how_willing_do_you_think_people_would_be_to_seek_psychiatric_help_for_anxiety,
    H4_post_ord = ExAnSt24_Given_the_information_provided_in_this_education_session_how_effective_do_you_think_psychosocial_interventions_will_be_in_improving_outcomes_in_anxiety_disorder,
    
    
    # --- Convert other predictors to factors ---
    EduType = factor(EduType, levels = c("Gen", "Evo")),
   
     # The 'Delivery' column corresponds to the 'livesession' variable in the pre-registered plan.
    ## i.e. 1 signifies livesession, 0 signifies virtual session
    Delivery = as.factor(Delivery),
    SessionID = as.factor(SessionID)
  )

## Some quick data checks here
## Check if all our main variables of interest (i.e. the pre- and the post-scores) are whole numbers and are below/equal to 7

## Doing this first for our pre-scores
## All variables here seem to respect the range of values and are whole numbers
analysis_data %>% 
  summarise(all_valid = all(H1_pre_ord %% 1 == 0 & H1_pre_ord <= 7, na.rm = TRUE)) 
  
analysis_data %>% 
  summarise(all_valid = all(H2_pre_ord %% 1 == 0 & H2_pre_ord <= 7, na.rm = TRUE))

analysis_data %>% 
  summarise(all_valid = all(H3_pre_ord %% 1 == 0 & H3_pre_ord <= 7, na.rm = TRUE)) 

analysis_data %>% 
  summarise(all_valid = all(H4_pre_ord %% 1 == 0 & H4_pre_ord <= 7, na.rm = TRUE)) 
  

## Doing this now for our post-scores
## All post-scores seem solid as well!
analysis_data %>% 
  summarise(all_valid = all(H1_post_ord %% 1 == 0 & H1_post_ord <= 7, na.rm = TRUE)) 

analysis_data %>% 
  summarise(all_valid = all(H2_post_ord %% 1 == 0 & H2_post_ord <= 7, na.rm = TRUE)) 

analysis_data %>% 
  summarise(all_valid = all(H3_post_ord %% 1 == 0 & H3_post_ord <= 7, na.rm = TRUE)) 

analysis_data %>% 
  summarise(all_valid = all(H4_post_ord %% 1 == 0 & H4_post_ord <= 7, na.rm = TRUE)) 


### Exporting selected variables from this cleaned dataset for our other robustness check (i.e. multiple comparisions check)
## We will use this dataset in another script. Don't focus on this now.
export_data <- analysis_data %>%
  mutate(H1_pre = ExAnSt1_Given_your_current_understanding_how_optimistic_are_you_about_improvement_in_symptoms_in_patients_presenting_with_anxiety_disorder,
         H2_pre = ExAnSt2_Given_your_current_understanding_how_willing_do_you_think_patients_would_be_to_share_anxiety_disorder_diagnoses,
         H3_pre = ExAnSt8_Given_current_public_knowledge_of_anxiety_how_willing_do_you_think_people_are_to_seek_psychiatric_help_for_anxiety,
         H4_pre = ExAnSt4_Given_your_current_understanding_how_effective_do_you_think_psychosocial_interventions_will_be_in_improving_outcomes_in_anxiety_disorder) %>%
  
  select(EduType, SessionID, Delivery, ParticipantID, H1_post_ord, H1_pre, H2_post_ord, H2_pre, H3_post_ord, H3_pre, H4_post_ord, H4_pre)

# Export this now
write_csv(export_data, "MComparision_Check_h1_h4.csv")



#-----------------------------------------------------------------------------
# STEP 3: RUN MODELS FOR EACH HYPOTHESIS
#-----------------------------------------------------------------------------
# For each hypothesis, we will run the three pre-registered specified models:
# 1. Primary Logit Model for log-odds to be converted to odds ratio
# 2. Probit Model for Effect Size in standard deviations
# 3. Logit Model with Monotonic Pre-Score Predictor for Robustness
#
# We will also add some weakly regularizing priors to test the models with and without
# such priors.
#

# Define the set of regularizing priors
WR_priors <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 1.5), class = "Intercept")
)


# --- HYPOTHESIS 1: Optimism about recovery of patients with anxiety disorder ---
########################################################################################

###### MAIN MODEL ###############

h1_logit_model <- brm(
  formula = bf(H1_post_ord ~ 1 + EduType + H1_pre_z + Delivery + (1 | SessionID)),
  data = analysis_data, family = cumulative(link = "logit"),
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

## Run only the second line below for directly loading the model that we have already run
saveRDS(h1_logit_model, file = "h1_logit_model.rds")
h1_logit_model <- readRDS(file = "h1_logit_model.rds")

## Model summary
summary(h1_logit_model)

## Drawing from the posterior and plotting the main effect
mod1_draws_h1_logit_model <- as_draws_df(h1_logit_model)
hist(mod1_draws_h1_logit_model$b_EduTypeEvo)

# Calculate and print the percentage of posterior draws > 0, rounded to 2 decimal places.
print(paste0("Percentage of posterior distribution above 0: ",
             round(mean(mod1_draws_h1_logit_model$b_EduTypeEvo > 0) * 100, 2), "%"))

## Double-checking this estimate using the hypothesis function in brms
## To check how much of the posterior probability is above zero, look at the post. prob column in the output
## Seems to check out
hypothesis(h1_logit_model, "EduTypeEvo>0")


### H1 logit model with prior ###
##########################################
h1_logit_model_WRprior <- brm(
  formula = bf(H1_post_ord ~ 1 + EduType + H1_pre_z + Delivery + (1 | SessionID)),
  data = analysis_data, family = cumulative(link = "logit"),
  prior = WR_priors,
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

## Run only the second line below for directly loading the model that we have already run
saveRDS(h1_logit_model_WRprior, file = "h1_logit_model_WRprior.rds")
h1_logit_model_WRprior <- readRDS(file = "h1_logit_model_WRprior.rds")


## Model summary
summary(h1_logit_model_WRprior)

## Drawing from the posterior and plotting the main effect
mod1_draws_h1_logit_model_WRprior <- as_draws_df(h1_logit_model_WRprior)
hist(mod1_draws_h1_logit_model_WRprior$b_EduTypeEvo)

# Calculate and print the percentage of posterior draws > 0, rounded to 2 decimal places.
print(paste0("Percentage of posterior distribution above 0: ",
             round(mean(mod1_draws_h1_logit_model_WRprior$b_EduTypeEvo > 0) * 100, 2), "%"))

## Double-checking this estimate using the hypothesis function in brms
## To check how much of the posterior probability is above zero, look at the post. prob column in the output
## Seems to check out
hypothesis(h1_logit_model_WRprior, "EduTypeEvo>0")


### H1 probit model ###
##########################################
h1_probit_model <- brm(
  formula = bf(H1_post_ord ~ 1 + EduType + H1_pre_z + Delivery + (1 | SessionID)),
  data = analysis_data, family = cumulative(link = "probit"),
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

## Run only the second line below for directly loading the model that we have already run
saveRDS(h1_probit_model, file = "h1_probit_model.rds")
h1_probit_model <- readRDS(file = "h1_probit_model.rds")

## Model summary
summary(h1_probit_model)

## Drawing from the posterior and plotting the main effect
mod1_draws_h1_probit_model <- as_draws_df(h1_probit_model)
hist(mod1_draws_h1_probit_model$b_EduTypeEvo)

# Calculate and print the percentage of posterior draws > 0, rounded to 2 decimal places.
print(paste0("Percentage of posterior distribution above 0: ",
             round(mean(mod1_draws_h1_probit_model$b_EduTypeEvo > 0) * 100, 2), "%"))

## Double-checking this estimate using the hypothesis function in brms
## To check how much of the posterior probability is above zero, look at the post. prob column in the output
## Seems to check out
hypothesis(h1_probit_model, "EduTypeEvo>0")


### H1 probit model with prior ###
##########################################
h1_probit_model_WRprior <- brm(
  formula = bf(H1_post_ord ~ 1 + EduType + H1_pre_z + Delivery + (1 | SessionID)),
  data = analysis_data, family = cumulative(link = "probit"),
  prior = WR_priors,
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

## Run only the second line below for directly loading the model that we have already run
saveRDS(h1_probit_model_WRprior, file = "h1_probit_model_WRprior.rds")
h1_probit_model_WRprior <- readRDS(file = "h1_probit_model_WRprior.rds")

## Model summary
summary(h1_probit_model_WRprior)

## Drawing from the posterior and plotting the main effect
mod1_draws_h1_probit_model_WRprior <- as_draws_df(h1_probit_model_WRprior)
hist(mod1_draws_h1_probit_model_WRprior$b_EduTypeEvo)

# Calculate and print the percentage of posterior draws > 0, rounded to 2 decimal places.
print(paste0("Percentage of posterior distribution above 0: ",
             round(mean(mod1_draws_h1_probit_model_WRprior$b_EduTypeEvo > 0) * 100, 2), "%"))

## Double-checking this estimate using the hypothesis function in brms
## To check how much of the posterior probability is above zero, look at the post. prob column in the output
## Seems to check out
hypothesis(h1_probit_model_WRprior, "EduTypeEvo>0")


### H1 robustness check model ###
##########################################
h1_robust_model <- brm(
  formula = bf(H1_post_ord ~ 1 + EduType + Delivery + mo(H1_pre_ord) + (1 | SessionID)),
  data = analysis_data, family = cumulative(link = "logit"),
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

## Run only the second line below for directly loading the model that we have already run
saveRDS(h1_robust_model, file = "h1_robust_model.rds")
h1_robust_model <- readRDS(file = "h1_robust_model.rds")

## Model summary
summary(h1_robust_model)

## Drawing from the posterior and plotting the main effect
mod1_draws_h1_robust_model <- as_draws_df(h1_robust_model)
hist(mod1_draws_h1_robust_model$b_EduTypeEvo)

# Calculate and print the percentage of posterior draws > 0, rounded to 2 decimal places.
print(paste0("Percentage of posterior distribution above 0: ",
             round(mean(mod1_draws_h1_robust_model$b_EduTypeEvo > 0) * 100, 2), "%"))

## Double-checking this estimate using the hypothesis function in brms
## To check how much of the posterior probability is above zero, look at the post. prob column in the output
## Seems to check out
hypothesis(h1_robust_model, "EduTypeEvo>0")



### H1 robustness check model with prior ###
##########################################
h1_robust_model_WRprior <- brm(
  formula = bf(H1_post_ord ~ 1 + EduType + Delivery + mo(H1_pre_ord) + (1 | SessionID)),
  data = analysis_data, family = cumulative(link = "logit"),
  prior = WR_priors,
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

## Run only the second line below for directly loading the model that we have already run
saveRDS(h1_robust_model_WRprior, file = "h1_robust_model_WRprior.rds")
h1_robust_model_WRprior <- readRDS(file = "h1_robust_model_WRprior.rds")

## Model summary
summary(h1_robust_model_WRprior)

## Drawing from the posterior and plotting the main effect
mod1_draws_h1_robust_model_WRprior <- as_draws_df(h1_robust_model_WRprior)
hist(mod1_draws_h1_robust_model_WRprior$b_EduTypeEvo)

# Calculate and print the percentage of posterior draws > 0, rounded to 2 decimal places.
print(paste0("Percentage of posterior distribution above 0: ",
             round(mean(mod1_draws_h1_robust_model_WRprior$b_EduTypeEvo > 0) * 100, 2), "%"))

## Double-checking this estimate using the hypothesis function in brms
## To check how much of the posterior probability is above zero, look at the post. prob column in the output
## Seems to check out
hypothesis(h1_robust_model_WRprior, "EduTypeEvo>0")


# --- HYPOTHESIS 2: Expected willingness for patients to share diagnosis ---
########################################################################################


###### MAIN MODEL ###############
h2_logit_model <- brm(
  formula = bf(H2_post_ord ~ 1 + EduType + H2_pre_z + Delivery + (1 | SessionID)),
  data = analysis_data, family = cumulative(link = "logit"),
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

## Run only the second line below for directly loading the model that we have already run
saveRDS(h2_logit_model, file = "h2_logit_model.rds")
h2_logit_model <- readRDS(file = "h2_logit_model.rds")

## Model summary
summary(h2_logit_model)

## Drawing from the posterior and plotting the main effect
mod1_draws_h2_logit_model <- as_draws_df(h2_logit_model)
hist(mod1_draws_h2_logit_model$b_EduTypeEvo)

# Calculate and print the percentage of posterior draws > 0, rounded to 2 decimal places.
print(paste0("Percentage of posterior distribution above 0: ",
             round(mean(mod1_draws_h2_logit_model$b_EduTypeEvo > 0) * 100, 2), "%"))

## Double-checking this estimate using the hypothesis function in brms
## To check how much of the posterior probability is above zero, look at the post. prob column in the output
## Seems to check out
hypothesis(h2_logit_model, "EduTypeEvo>0")


### H2 logit model with prior ###
##########################################
h2_logit_model_WRprior <- brm(
  formula = bf(H2_post_ord ~ 1 + EduType + H2_pre_z + Delivery + (1 | SessionID)),
  data = analysis_data, family = cumulative(link = "logit"),
  prior = WR_priors,
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

## Run only the second line below for directly loading the model that we have already run
saveRDS(h2_logit_model_WRprior, file = "h2_logit_model_WRprior.rds")
h2_logit_model_WRprior <- readRDS(file = "h2_logit_model_WRprior.rds")

## Model summary
summary(h2_logit_model_WRprior)

## Drawing from the posterior and plotting the main effect
mod1_draws_h2_logit_model_WRprior <- as_draws_df(h2_logit_model_WRprior)
hist(mod1_draws_h2_logit_model_WRprior$b_EduTypeEvo)

# Calculate and print the percentage of posterior draws > 0, rounded to 2 decimal places.
print(paste0("Percentage of posterior distribution above 0: ",
             round(mean(mod1_draws_h2_logit_model_WRprior$b_EduTypeEvo > 0) * 100, 2), "%"))

## Double-checking this estimate using the hypothesis function in brms
## To check how much of the posterior probability is above zero, look at the post. prob column in the output
## Seems to check out
hypothesis(h2_logit_model_WRprior, "EduTypeEvo>0")



### H2 probit model ###
##########################################
h2_probit_model <- brm(
  formula = bf(H2_post_ord ~ 1 + EduType + H2_pre_z + Delivery + (1 | SessionID)),
  data = analysis_data, family = cumulative(link = "probit"),
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

## Run only the second line below for directly loading the model that we have already run
saveRDS(h2_probit_model, file = "h2_probit_model.rds")
h2_probit_model <- readRDS(file = "h2_probit_model.rds")

## Model summary
summary(h2_probit_model)

## Drawing from the posterior and plotting the main effect
mod1_draws_h2_probit_model <- as_draws_df(h2_probit_model)
hist(mod1_draws_h2_probit_model$b_EduTypeEvo)

# Calculate and print the percentage of posterior draws > 0, rounded to 2 decimal places.
print(paste0("Percentage of posterior distribution above 0: ",
             round(mean(mod1_draws_h2_probit_model$b_EduTypeEvo > 0) * 100, 2), "%"))

## Double-checking this estimate using the hypothesis function in brms
## To check how much of the posterior probability is above zero, look at the post. prob column in the output
## Seems to check out
hypothesis(h2_probit_model, "EduTypeEvo>0")


### H2 probit model with prior ###
##########################################
h2_probit_model_WRprior <- brm(
  formula = bf(H2_post_ord ~ 1 + EduType + H2_pre_z + Delivery + (1 | SessionID)),
  data = analysis_data, family = cumulative(link = "probit"),
  prior = WR_priors,
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

## Run only the second line below for directly loading the model that we have already run
saveRDS(h2_probit_model_WRprior, file = "h2_probit_model_WRprior.rds")
h2_probit_model_WRprior <- readRDS(file = "h2_probit_model_WRprior.rds")

## Model summary
summary(h2_probit_model_WRprior)

## Drawing from the posterior and plotting the main effect
mod1_draws_h2_probit_model_WRprior <- as_draws_df(h2_probit_model_WRprior)
hist(mod1_draws_h2_probit_model_WRprior$b_EduTypeEvo)

# Calculate and print the percentage of posterior draws > 0, rounded to 2 decimal places.
print(paste0("Percentage of posterior distribution above 0: ",
             round(mean(mod1_draws_h2_probit_model_WRprior$b_EduTypeEvo > 0) * 100, 2), "%"))

## Double-checking this estimate using the hypothesis function in brms
## To check how much of the posterior probability is above zero, look at the post. prob column in the output
## Seems to check out
hypothesis(h2_probit_model_WRprior, "EduTypeEvo>0")


### H2 robustness check model ###
##########################################
h2_robust_model <- brm(
  formula = bf(H2_post_ord ~ 1 + EduType + Delivery + mo(H2_pre_ord) + (1 | SessionID)),
  data = analysis_data, family = cumulative(link = "logit"),
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

## Run only the second line below for directly loading the model that we have already run
saveRDS(h2_robust_model, file = "h2_robust_model.rds")
h2_robust_model <- readRDS(file = "h2_robust_model.rds")

## Model summary
summary(h2_robust_model)

## Drawing from the posterior and plotting the main effect
mod1_draws_h2_robust_model <- as_draws_df(h2_robust_model)
hist(mod1_draws_h2_robust_model$b_EduTypeEvo)

# Calculate and print the percentage of posterior draws > 0, rounded to 2 decimal places.
print(paste0("Percentage of posterior distribution above 0: ",
             round(mean(mod1_draws_h2_robust_model$b_EduTypeEvo > 0) * 100, 2), "%"))

## Double-checking this estimate using the hypothesis function in brms
## To check how much of the posterior probability is above zero, look at the post. prob column in the output
## Seems to check out
hypothesis(h2_robust_model, "EduTypeEvo>0")



### H2 robustness check with prior ###
##########################################
h2_robust_model_WRprior <- brm(
  formula = bf(H2_post_ord ~ 1 + EduType + Delivery + mo(H2_pre_ord) + (1 | SessionID)),
  data = analysis_data, family = cumulative(link = "logit"),
  prior = WR_priors,
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

## Run only the second line below for directly loading the model that we have already run
saveRDS(h2_robust_model_WRprior, file = "h2_robust_model_WRprior.rds")
h2_robust_model_WRprior <- readRDS(file = "h2_robust_model_WRprior.rds")

## Model summary
summary(h2_robust_model_WRprior)

## Drawing from the posterior and plotting the main effect
mod1_draws_h2_robust_model_WRprior <- as_draws_df(h2_robust_model_WRprior)
hist(mod1_draws_h2_robust_model_WRprior$b_EduTypeEvo)

# Calculate and print the percentage of posterior draws > 0, rounded to 2 decimal places.
print(paste0("Percentage of posterior distribution above 0: ",
             round(mean(mod1_draws_h2_robust_model_WRprior$b_EduTypeEvo > 0) * 100, 2), "%"))

## Double-checking this estimate using the hypothesis function in brms
## To check how much of the posterior probability is above zero, look at the post. prob column in the output
## Seems to check out
hypothesis(h2_robust_model_WRprior, "EduTypeEvo>0")


# --- HYPOTHESIS 3: Predicted willingness for the public to seek help ---
########################################################################################

###### MAIN MODEL ###############

h3_logit_model <- brm(
  formula = bf(H3_post_ord ~ 1 + EduType + H3_pre_z + Delivery + (1 | SessionID)),
  data = analysis_data, family = cumulative(link = "logit"),
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

## Run only the second line below for directly loading the model that we have already run
saveRDS(h3_logit_model, file = "h3_logit_model.rds")
h3_logit_model <- readRDS(file = "h3_logit_model.rds")

## Model summary
summary(h3_logit_model)

## Drawing from the posterior and plotting the main effect
mod1_draws_h3_logit_model <- as_draws_df(h3_logit_model)
hist(mod1_draws_h3_logit_model$b_EduTypeEvo)

# Calculate and print the percentage of posterior draws > 0, rounded to 2 decimal places.
print(paste0("Percentage of posterior distribution above 0: ",
             round(mean(mod1_draws_h3_logit_model$b_EduTypeEvo > 0) * 100, 2), "%"))

## Double-checking this estimate using the hypothesis function in brms
## To check how much of the posterior probability is above zero, look at the post. prob column in the output
## Seems to check out
hypothesis(h3_logit_model, "EduTypeEvo>0")


### H3 logit model with prior ###
##########################################
h3_logit_model_WRprior <- brm(
  formula = bf(H3_post_ord ~ 1 + EduType + H3_pre_z + Delivery + (1 | SessionID)),
  data = analysis_data, family = cumulative(link = "logit"),
  prior = WR_priors,
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

## Run only the second line below for directly loading the model that we have already run
saveRDS(h3_logit_model_WRprior, file = "h3_logit_model_WRprior.rds")
h3_logit_model_WRprior <- readRDS(file = "h3_logit_model_WRprior.rds")

## Model summary
summary(h3_logit_model_WRprior)

## Drawing from the posterior and plotting the main effect
mod1_draws_h3_logit_model_WRprior <- as_draws_df(h3_logit_model_WRprior)
hist(mod1_draws_h3_logit_model_WRprior$b_EduTypeEvo)

# Calculate and print the percentage of posterior draws > 0, rounded to 2 decimal places.
print(paste0("Percentage of posterior distribution above 0: ",
             round(mean(mod1_draws_h3_logit_model_WRprior$b_EduTypeEvo > 0) * 100, 2), "%"))

## Double-checking this estimate using the hypothesis function in brms
## To check how much of the posterior probability is above zero, look at the post. prob column in the output
## Seems to check out
hypothesis(h3_logit_model_WRprior, "EduTypeEvo>0")


### H3 probit model ###
##########################################
h3_probit_model <- brm(
  formula = bf(H3_post_ord ~ 1 + EduType + H3_pre_z + Delivery + (1 | SessionID)),
  data = analysis_data, family = cumulative(link = "probit"),
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

## Run only the second line below for directly loading the model that we have already run
saveRDS(h3_probit_model, file = "h3_probit_model.rds")
h3_probit_model <- readRDS(file = "h3_probit_model.rds")

## Model summary
summary(h3_probit_model)

## Drawing from the posterior and plotting the main effect
mod1_draws_h3_probit_model <- as_draws_df(h3_probit_model)
hist(mod1_draws_h3_probit_model$b_EduTypeEvo)

# Calculate and print the percentage of posterior draws > 0, rounded to 2 decimal places.
print(paste0("Percentage of posterior distribution above 0: ",
             round(mean(mod1_draws_h3_probit_model$b_EduTypeEvo > 0) * 100, 2), "%"))

## Double-checking this estimate using the hypothesis function in brms
## To check how much of the posterior probability is above zero, look at the post. prob column in the output
## Seems to check out
hypothesis(h3_probit_model, "EduTypeEvo>0")


### H3 probit model with prior ###
##########################################
h3_probit_model_WRprior <- brm(
  formula = bf(H3_post_ord ~ 1 + EduType + H3_pre_z + Delivery + (1 | SessionID)),
  data = analysis_data, family = cumulative(link = "probit"),
  prior = WR_priors,
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

## Run only the second line below for directly loading the model that we have already run
saveRDS(h3_probit_model_WRprior, file = "h3_probit_model_WRprior.rds")
h3_probit_model_WRprior <- readRDS(file = "h3_probit_model_WRprior.rds")

## Model summary
summary(h3_probit_model_WRprior)

## Drawing from the posterior and plotting the main effect
mod1_draws_h3_probit_model_WRprior <- as_draws_df(h3_probit_model_WRprior)
hist(mod1_draws_h3_probit_model_WRprior$b_EduTypeEvo)

# Calculate and print the percentage of posterior draws > 0, rounded to 2 decimal places.
print(paste0("Percentage of posterior distribution above 0: ",
             round(mean(mod1_draws_h3_probit_model_WRprior$b_EduTypeEvo > 0) * 100, 2), "%"))

## Double-checking this estimate using the hypothesis function in brms
## To check how much of the posterior probability is above zero, look at the post. prob column in the output
## Seems to check out
hypothesis(h3_probit_model_WRprior, "EduTypeEvo>0")


### H3 robustness check model ###
##########################################
h3_robust_model <- brm(
  formula = bf(H3_post_ord ~ 1 + EduType + Delivery + mo(H3_pre_ord) + (1 | SessionID)),
  data = analysis_data, family = cumulative(link = "logit"),
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

## Run only the second line below for directly loading the model that we have already run
saveRDS(h3_robust_model, file = "h3_robust_model.rds")
h3_robust_model <- readRDS(file = "h3_robust_model.rds")

## Model summary
summary(h3_robust_model)

## Drawing from the posterior and plotting the main effect
mod1_draws_h3_robust_model <- as_draws_df(h3_robust_model)
hist(mod1_draws_h3_robust_model$b_EduTypeEvo)

# Calculate and print the percentage of posterior draws > 0, rounded to 2 decimal places.
print(paste0("Percentage of posterior distribution above 0: ",
             round(mean(mod1_draws_h3_robust_model$b_EduTypeEvo > 0) * 100, 2), "%"))

## Double-checking this estimate using the hypothesis function in brms
## To check how much of the posterior probability is above zero, look at the post. prob column in the output
## Seems to check out
hypothesis(h3_robust_model, "EduTypeEvo>0")


### H3 robustness check with prior ###
##########################################
h3_robust_model_WRprior <- brm(
  formula = bf(H3_post_ord ~ 1 + EduType + Delivery + mo(H3_pre_ord) + (1 | SessionID)),
  data = analysis_data, family = cumulative(link = "logit"),
  prior = WR_priors,
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

## Run only the second line below for directly loading the model that we have already run
saveRDS(h3_robust_model_WRprior, file = "h3_robust_model_WRprior.rds")
h3_robust_model_WRprior <- readRDS(file = "h3_robust_model_WRprior.rds")

## Model summary
summary(h3_robust_model_WRprior)

## Drawing from the posterior and plotting the main effect
mod1_draws_h3_robust_model_WRprior <- as_draws_df(h3_robust_model_WRprior)
hist(mod1_draws_h3_robust_model_WRprior$b_EduTypeEvo)

# Calculate and print the percentage of posterior draws > 0, rounded to 2 decimal places.
print(paste0("Percentage of posterior distribution above 0: ",
             round(mean(mod1_draws_h3_robust_model_WRprior$b_EduTypeEvo > 0) * 100, 2), "%"))

## Double-checking this estimate using the hypothesis function in brms
## To check how much of the posterior probability is above zero, look at the post. prob column in the output
## Seems to check out
hypothesis(h3_robust_model_WRprior, "EduTypeEvo>0")


# --- HYPOTHESIS 4: Belief in efficacy of psychosocial interventions ---
########################################################################################

###### MAIN MODEL ###############

h4_logit_model <- brm(
  formula = bf(H4_post_ord ~ 1 + EduType + H4_pre_z + Delivery + (1 | SessionID)),
  data = analysis_data, family = cumulative(link = "logit"),
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

## Run only the second line below for directly loading the model that we have already run
saveRDS(h4_logit_model, file = "h4_logit_model.rds")
h4_logit_model <- readRDS(file = "h4_logit_model.rds")


## Model summary
summary(h4_logit_model)

## Drawing from the posterior and plotting the main effect
mod1_draws_h4_logit_model <- as_draws_df(h4_logit_model)
hist(mod1_draws_h4_logit_model$b_EduTypeEvo)

# Calculate and print the percentage of posterior draws > 0, rounded to 2 decimal places.
print(paste0("Percentage of posterior distribution above 0: ",
             round(mean(mod1_draws_h4_logit_model$b_EduTypeEvo > 0) * 100, 2), "%"))

## Double-checking this estimate using the hypothesis function in brms
## To check how much of the posterior probability is above zero, look at the post. prob column in the output
## Seems to check out
hypothesis(h4_logit_model, "EduTypeEvo>0")



### H4 logit model with prior ###
########################################################################################
h4_logit_model_WRprior <- brm(
  formula = bf(H4_post_ord ~ 1 + EduType + H4_pre_z + Delivery + (1 | SessionID)),
  data = analysis_data, family = cumulative(link = "logit"),
  prior = WR_priors,
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

## Run only the second line below for directly loading the model that we have already run
saveRDS(h4_logit_model_WRprior, file = "h4_logit_model_WRprior.rds")
h4_logit_model_WRprior <- readRDS(file = "h4_logit_model_WRprior.rds")

## Model summary
summary(h4_logit_model_WRprior)

## Drawing from the posterior and plotting the main effect
mod1_draws_h4_logit_model_WRprior <- as_draws_df(h4_logit_model_WRprior)
hist(mod1_draws_h4_logit_model_WRprior$b_EduTypeEvo)

# Calculate and print the percentage of posterior draws > 0, rounded to 2 decimal places.
print(paste0("Percentage of posterior distribution above 0: ",
             round(mean(mod1_draws_h4_logit_model_WRprior$b_EduTypeEvo > 0) * 100, 2), "%"))

## Double-checking this estimate using the hypothesis function in brms
## To check how much of the posterior probability is above zero, look at the post. prob column in the output
## Seems to check out
hypothesis(h4_logit_model_WRprior, "EduTypeEvo>0")


### H4 probit model ###
########################################################################################
h4_probit_model <- brm(
  formula = bf(H4_post_ord ~ 1 + EduType + H4_pre_z + Delivery + (1 | SessionID)),
  data = analysis_data, family = cumulative(link = "probit"),
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

## Run only the second line below for directly loading the model that we have already run
saveRDS(h4_probit_model, file = "h4_probit_model.rds")
h4_probit_model <- readRDS(file = "h4_probit_model.rds")

## Model summary
summary(h4_probit_model)

## Drawing from the posterior and plotting the main effect
mod1_draws_h4_probit_model <- as_draws_df(h4_probit_model)
hist(mod1_draws_h4_probit_model$b_EduTypeEvo)

# Calculate and print the percentage of posterior draws > 0, rounded to 2 decimal places.
print(paste0("Percentage of posterior distribution above 0: ",
             round(mean(mod1_draws_h4_probit_model$b_EduTypeEvo > 0) * 100, 2), "%"))

## Double-checking this estimate using the hypothesis function in brms
## To check how much of the posterior probability is above zero, look at the post. prob column in the output
## Seems to check out
hypothesis(h4_probit_model, "EduTypeEvo>0")


### H4 probit model with prior ###
########################################################################################
h4_probit_model_WRprior <- brm(
  formula = bf(H4_post_ord ~ 1 + EduType + H4_pre_z + Delivery + (1 | SessionID)),
  data = analysis_data, family = cumulative(link = "probit"),
  prior = WR_priors,
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

## Run only the second line below for directly loading the model that we have already run
saveRDS(h4_probit_model_WRprior, file = "h4_probit_model_WRprior.rds")
h4_probit_model_WRprior <- readRDS(file = "h4_probit_model_WRprior.rds")

## Model summary
summary(h4_probit_model_WRprior)

## Drawing from the posterior and plotting the main effect
mod1_draws_h4_probit_model_WRprior <- as_draws_df(h4_probit_model_WRprior)
hist(mod1_draws_h4_probit_model_WRprior$b_EduTypeEvo)

# Calculate and print the percentage of posterior draws > 0, rounded to 2 decimal places.
print(paste0("Percentage of posterior distribution above 0: ",
             round(mean(mod1_draws_h4_probit_model_WRprior$b_EduTypeEvo > 0) * 100, 2), "%"))

## Double-checking this estimate using the hypothesis function in brms
## To check how much of the posterior probability is above zero, look at the post. prob column in the output
## Seems to check out
hypothesis(h4_probit_model_WRprior, "EduTypeEvo>0")


### H4 robustness check model ###
########################################################################################
h4_robust_model <- brm(
  formula = bf(H4_post_ord ~ 1 + EduType + Delivery + mo(H4_pre_ord) + (1 | SessionID)),
  data = analysis_data, family = cumulative(link = "logit"),
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

## Run only the second line below for directly loading the model that we have already run
saveRDS(h4_robust_model, file = "h4_robust_model.rds")
h4_robust_model <- readRDS(file = "h4_robust_model.rds")

## Model summary
summary(h4_robust_model)

## Drawing from the posterior and plotting the main effect
mod1_draws_h4_robust_model <- as_draws_df(h4_robust_model)
hist(mod1_draws_h4_robust_model$b_EduTypeEvo)

# Calculate and print the percentage of posterior draws > 0, rounded to 2 decimal places.
print(paste0("Percentage of posterior distribution above 0: ",
             round(mean(mod1_draws_h4_robust_model$b_EduTypeEvo > 0) * 100, 2), "%"))

## Double-checking this estimate using the hypothesis function in brms
## To check how much of the posterior probability is above zero, look at the post. prob column in the output
## Seems to check out
hypothesis(h4_robust_model, "EduTypeEvo>0")


### H4 robustness check with prior ###
########################################################################################
h4_robust_model_WRprior <- brm(
  formula = bf(H4_post_ord ~ 1 + EduType + Delivery + mo(H4_pre_ord) + (1 | SessionID)),
  data = analysis_data, family = cumulative(link = "logit"),
  prior = WR_priors,
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

## Run only the second line below for directly loading the model that we have already run
saveRDS(h4_robust_model_WRprior, file = "h4_robust_model_WRprior.rds")
h4_robust_model_WRprior <- readRDS(file = "h4_robust_model_WRprior.rds")

## Model summary
summary(h4_robust_model_WRprior)

## Drawing from the posterior and plotting the main effect
mod1_draws_h4_robust_model_WRprior <- as_draws_df(h4_robust_model_WRprior)
hist(mod1_draws_h4_robust_model_WRprior$b_EduTypeEvo)

# Calculate and print the percentage of posterior draws > 0, rounded to 2 decimal places.
print(paste0("Percentage of posterior distribution above 0: ",
             round(mean(mod1_draws_h4_robust_model_WRprior$b_EduTypeEvo > 0) * 100, 2), "%"))

## Double-checking this estimate using the hypothesis function in brms
## To check how much of the posterior probability is above zero, look at the post. prob column in the output
## Seems to check out
hypothesis(h4_robust_model, "EduTypeEvo>0")

#-----------------------------------------------------------------------------
# STEP 4: VIEW MODEL SUMMARIES
#-----------------------------------------------------------------------------
# Print the summaries for each model to review the results. Find the precise credible interval which doesn't 
# overlap with 0. Make graphs displaying these results.
#
# For new sessions, load the models first (presuming they have been run and saved already)
#

# --- LOAD ALL SAVED MODELS ---

# HYPOTHESIS 1 MODELS
h1_logit_model <- readRDS("h1_logit_model.rds")
h1_logit_model_WRprior <- readRDS("h1_logit_model_WRprior.rds")
h1_probit_model <- readRDS("h1_probit_model.rds")
h1_probit_model_WRprior <- readRDS("h1_probit_model_WRprior.rds")
h1_robust_model <- readRDS("h1_robust_model.rds")
h1_robust_model_WRprior <- readRDS("h1_robust_model_WRprior.rds")

# HYPOTHESIS 2 MODELS
h2_logit_model <- readRDS("h2_logit_model.rds")
h2_logit_model_WRprior <- readRDS("h2_logit_model_WRprior.rds")
h2_probit_model <- readRDS("h2_probit_model.rds")
h2_probit_model_WRprior <- readRDS("h2_probit_model_WRprior.rds")
h2_robust_model <- readRDS("h2_robust_model.rds")
h2_robust_model_WRprior <- readRDS("h2_robust_model_WRprior.rds")

# HYPOTHESIS 3 MODELS
h3_logit_model <- readRDS("h3_logit_model.rds")
h3_logit_model_WRprior <- readRDS("h3_logit_model_WRprior.rds")
h3_probit_model <- readRDS("h3_probit_model.rds")
h3_probit_model_WRprior <- readRDS("h3_probit_model_WRprior.rds")
h3_robust_model <- readRDS("h3_robust_model.rds")
h3_robust_model_WRprior <- readRDS("h3_robust_model_WRprior.rds")

# HYPOTHESIS 4 MODELS
h4_logit_model <- readRDS("h4_logit_model.rds")
h4_logit_model_WRprior <- readRDS("h4_logit_model_WRprior.rds")
h4_probit_model <- readRDS("h4_probit_model.rds")
h4_probit_model_WRprior <- readRDS("h4_probit_model_WRprior.rds")
h4_robust_model <- readRDS("h4_robust_model.rds")
h4_robust_model_WRprior <- readRDS("h4_robust_model_WRprior.rds")


# NOW Printing Summaries of the models

###
# HYPOTHESIS 1 OPTIMISM ABOUT RECOVERY
###

print("--- H1 (Optimism) Logit Model Summary ---")
summary(h1_logit_model)
print("--- H1 (Optimism) Logit Model Summary With Regularizing Priors ---")
summary(h1_logit_model_WRprior)
print("--- H1 (Optimism) Probit Model Summary ---")
summary(h1_probit_model)
print("--- H1 (Optimism) Probit Model Summary With Regularizing Priors ---")
summary(h1_probit_model_WRprior)
print("--- H1 (Optimism) Robustness Check Summary ---")
summary(h1_robust_model)
print("--- H1 (Optimism) Robustness Check Summary With Regularizing Priors ---")
summary(h1_robust_model_WRprior)



#HYPOTHESIS 2 WILLINGNESS TO SHARE DIAGNOSIS

print("--- H2 (Willingness to Share) Logit Model Summary ---")
summary(h2_logit_model)
print("--- H2 (Willingness to Share) Logit Model Summary With Regularizing Priors ---")
summary(h2_logit_model_WRprior)
print("--- H2 (Willingness to Share) Probit Model Summary ---")
summary(h2_probit_model)
print("--- H2 (Willingness to Share) Probit Model Summary With Regularizing Priors ---")
summary(h2_probit_model_WRprior)
print("--- H2 (Willingness to Share) Robustness Check Summary ---")
summary(h2_robust_model)
print("--- H2 (Willingness to Share) Robustness Check Summary With Regularizing Priors ---")
summary(h2_robust_model_WRprior)


#HYPOTHESIS 3 WILLINGNESS TO SEEK HELP IF PUBLICLY KNOWN

print("--- H3 (Willingness to Seek Help) Logit Model Summary ---")
summary(h3_logit_model)
print("--- H3 (Willingness to Seek Help) Logit Model Summary With Regularizing Priors ---")
summary(h3_logit_model_WRprior)
print("--- H3 (Willingness to Seek Help) Probit Model Summary ---")
summary(h3_probit_model)
print("--- H3 (Willingness to Seek Help) Probit Model Summary With Regularizing Priors ---")
summary(h3_probit_model_WRprior)
print("--- H3 (Willingness to Seek Help) Robustness Check Summary ---")
summary(h3_robust_model)
print("--- H3 (Willingness to Seek Help) Robustness Check Summary With Regularizing Priors ---")
summary(h3_robust_model_WRprior)


#### HYPOTHESIS 4 PSYCHOSOCIAL INTERVENTIONS SUMMARIES AND VISUALISATION
####
####
print("--- H4 (Psychosocial Interventions) Logit Model Summary ---")
summary(h4_logit_model)
print("--- H4 (Psychosocial Interventions) Logit Model Summary With Regularizing Priors ---")
summary(h4_logit_model_WRprior)
print("--- H4 (Psychosocial Interventions) Probit Model Summary ---")
summary(h4_probit_model)
print("--- H4 (Psychosocial Interventions) Probit Model Summary With Regularizing Priors ---")
summary(h4_probit_model_WRprior)
print("--- H4 (Psychosocial Interventions) Robustness Check Summary ---")
summary(h4_robust_model)
print("--- H4 (Psychosocial Interventions) Robustness Check Summary With Regularizing Priors ---")
summary(h4_robust_model_WRprior)


########                      #########
######## GRAPH MAKING SECTION #########
########                      #########



################################################################################
# TEST SINGLE POSTERIOR DISTRIBUTION PLOTS (ODDS RATIO SCALE)
################################################################################
#
# This test section creates posterior distribution plots for the main logit models.
# Crucially, it converts the log-odds coefficients from the model into odds ratios
# (by using exp()) before plotting. This makes the x-axis more interpretable.
# An odds ratio of 1 indicates no effect.
#
# We will build the plot for H1 as a detailed example.
#

# --- H1: Posterior Distribution Plot for Odds Ratio ---

# Step 1: Extract the posterior draws for the log-odds coefficient.
# We use the main logit model for H1 (`h1_logit_model`).
log_odds_draws_h1 <- as_draws_df(h1_logit_model_WRprior)$b_EduTypeEvo

# Step 2: Convert the log-odds draws to odds ratio draws.
# The odds ratio is simply the exponent of the log-odds.
odds_ratio_draws_h1 <- exp(log_odds_draws_h1)

# Step 3: Manually calculate the density data for the odds ratios.
# This approach is more stable for creating conditional shading in ggplot2.
density_data_h1 <- density(odds_ratio_draws_h1)
density_df_h1 <- data.frame(x = density_data_h1$x, y = density_data_h1$y)

# Step 4: Calculate the probability of the effect being positive (Odds Ratio > 1).
prob_positive_h1 <- mean(odds_ratio_draws_h1 > 1) * 100

# Step 5: Build the plot using the pre-calculated density data.
ggplot(density_df_h1, aes(x = x, y = y)) +
  
  # Layer 1: Draw the outline of the density curve.
  geom_line(linewidth = 1) +
  
  # Layer 2: Shade the area where the odds ratio is GREATER THAN 1 (positive effect).
  # This indicates that the Evolutionary psychoeducation group had higher odds
  # of increasing their score compared to the Genetic group.
  geom_ribbon(data = subset(density_df_h1, x > 1), 
              aes(ymin = 0, ymax = y), 
              fill = "#f1a226", alpha = 0.8) + # Using the "Evo" color
  
  # Layer 3: Shade the area where the odds ratio is LESS THAN 1 (negative effect).
  geom_ribbon(data = subset(density_df_h1, x < 1), 
              aes(ymin = 0, ymax = y), 
              fill = "#f1a226", alpha = 0.3) +
  
  # Layer 4: Add a vertical dashed line at x = 1.
  # An odds ratio of 1 signifies no difference between the groups.
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +
  
  # Layer 5: Add informative labels and titles.
  labs(
    title = "Optimism about recovery from anxiety disorder (H1)",
    subtitle = paste0("The probability of a positive effect from evolutionary vs. genetic education (Odds Ratio > 1) is ", round(prob_positive_h1, 2), "%"),
    x = "Odds Ratio",
    y = "Density"
  ) +
  
  # Layer 6: Use a clean theme for a professional look.
  theme_minimal(base_size = 15)



################################################################################
# STEP 5: FINAL POSTERIOR DISTRIBUTION PLOTS (ODDS RATIO SCALE) & COMBINED FIGURE
################################################################################
#
# This section creates posterior distribution plots for the main logit models
# for H1, H2, H3, and H4. It converts the log-odds coefficients into odds ratios
# (by using exp()) to make the x-axis more interpretable.
#
# Finally, it uses the 'patchwork' package to combine the four plots into a
# single 2x2 figure for a comprehensive overview.
#


# --- H1: Posterior Distribution Plot for Odds Ratio ---

# Step 1: Extract and convert posterior draws for H1
log_odds_draws_h1 <- as_draws_df(h1_logit_model_WRprior)$b_EduTypeEvo
odds_ratio_draws_h1 <- exp(log_odds_draws_h1)

# Step 2: Calculate density data for H1
density_data_h1 <- density(odds_ratio_draws_h1)
density_df_h1 <- data.frame(x = density_data_h1$x, y = density_data_h1$y)

# Step 3: Calculate probability for H1
prob_positive_h1 <- mean(odds_ratio_draws_h1 > 1) * 100

# Step 4: Build and SAVE the plot for H1 to a variable
plot_h1 <- ggplot(density_df_h1, aes(x = x, y = y)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = subset(density_df_h1, x > 1), 
              aes(ymin = 0, ymax = y), 
              fill = "#D55E00", alpha = 0.8) +
  geom_ribbon(data = subset(density_df_h1, x < 1), 
              aes(ymin = 0, ymax = y), 
              fill = "#D55E00", alpha = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +
  labs(
    title = "Optimism about recovery from anxiety disorder (H1)",
    subtitle = paste0("Probability of positive effect: ", round(prob_positive_h1, 2), "%"),
    x = "Odds Ratio",
    y = "Density"
  ) +
  theme_minimal(base_size = 15)


# --- H2: Posterior Distribution Plot for Odds Ratio ---

# Step 1: Extract and convert posterior draws for H2
log_odds_draws_h2 <- as_draws_df(h2_logit_model_WRprior)$b_EduTypeEvo
odds_ratio_draws_h2 <- exp(log_odds_draws_h2)

# Step 2: Calculate density data for H2
density_data_h2 <- density(odds_ratio_draws_h2)
density_df_h2 <- data.frame(x = density_data_h2$x, y = density_data_h2$y)

# Step 3: Calculate probability for H2
prob_positive_h2 <- mean(odds_ratio_draws_h2 > 1) * 100

# Step 4: Build and SAVE the plot for H2 to a variable
plot_h2 <- ggplot(density_df_h2, aes(x = x, y = y)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = subset(density_df_h2, x > 1), 
              aes(ymin = 0, ymax = y), 
              fill = "#D55E00", alpha = 0.8) +
  geom_ribbon(data = subset(density_df_h2, x < 1), 
              aes(ymin = 0, ymax = y), 
              fill = "#D55E00", alpha = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +
  labs(
    title = "Willingness to share diagnosis (H2)",
    subtitle = paste0("Probability of positive effect: ", round(prob_positive_h2, 2), "%"),
    x = "Odds Ratio",
    y = "Density"
  ) +
  theme_minimal(base_size = 15)


# --- H3: Posterior Distribution Plot for Odds Ratio ---

# Step 1: Extract and convert posterior draws for H3
log_odds_draws_h3 <- as_draws_df(h3_logit_model_WRprior)$b_EduTypeEvo
odds_ratio_draws_h3 <- exp(log_odds_draws_h3)

# Step 2: Calculate density data for H3
density_data_h3 <- density(odds_ratio_draws_h3)
density_df_h3 <- data.frame(x = density_data_h3$x, y = density_data_h3$y)

# Step 3: Calculate probability for H3
prob_positive_h3 <- mean(odds_ratio_draws_h3 > 1) * 100

# Step 4: Build and SAVE the plot for H3 to a variable
plot_h3 <- ggplot(density_df_h3, aes(x = x, y = y)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = subset(density_df_h3, x > 1), 
              aes(ymin = 0, ymax = y), 
              fill = "#D55E00", alpha = 0.8) +
  geom_ribbon(data = subset(density_df_h3, x < 1), 
              aes(ymin = 0, ymax = y), 
              fill = "#D55E00", alpha = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +
  labs(
    title = "Willingness to seek help (H3)",
    subtitle = paste0("Probability of positive effect: ", round(prob_positive_h3, 2), "%"),
    x = "Odds Ratio",
    y = "Density"
  ) +
  theme_minimal(base_size = 15)

# --- H4: Posterior Distribution Plot for Odds Ratio ---

# Step 1: Extract and convert posterior draws for H4
log_odds_draws_h4 <- as_draws_df(h4_logit_model_WRprior)$b_EduTypeEvo
odds_ratio_draws_h4 <- exp(log_odds_draws_h4)

# Step 2: Calculate density data for H4
density_data_h4 <- density(odds_ratio_draws_h4)
density_df_h4 <- data.frame(x = density_data_h4$x, y = density_data_h4$y)

# Step 3: Calculate probability for H4
prob_positive_h4 <- mean(odds_ratio_draws_h4 > 1) * 100

# Step 4: Build and SAVE the plot for H4 to a variable
plot_h4 <- ggplot(density_df_h4, aes(x = x, y = y)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = subset(density_df_h4, x > 1), 
              aes(ymin = 0, ymax = y), 
              fill = "#D55E00", alpha = 0.8) +
  geom_ribbon(data = subset(density_df_h4, x < 1), 
              aes(ymin = 0, ymax = y), 
              fill = "#D55E00", alpha = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +
  labs(
    title = "Efficacy of psychosocial interventions (H4)",
    subtitle = paste0("Probability of positive effect: ", round(prob_positive_h4, 2), "%"),
    x = "Odds Ratio",
    y = "Density"
  ) +
  theme_minimal(base_size = 15)


# --- COMBINE ALL FOUR PLOTS ---

# Arrange the four saved plots into a 2x2 grid and add an overall title.
(plot_h1 + plot_h2) / (plot_h3 + plot_h4) + 
  plot_annotation(
    title = 'Posterior Distributions for the Effect of Evolutionary vs. Genetic Psychoeducation',
    theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))
  )



######## FROM HERE ARE THE RAW DATA FIGURES WITH SCATTERED POINTS AND PRE-POST CHANGES


########                     #########
########   GRAPH FOR H1      #########
########                     #########

# --- STEP 1: PREPARE THE DATA FOR PLOTTING H1 ---

# We need to reshape the data from wide to long format for ggplot.
# We'll start with  original 'psy_data' to use the raw numerical scores.
# Adding a participant ID for good practice.


# Select the relevant columns, rename them for easier use, and pivot to a long format.
plot_data_h1_long <- psy_data %>%
  # Add a unique ID for each participant
  mutate(participant_id = row_number()) %>%
  # Select only the columns we need for this specific plot
  select(
    participant_id,
    EduType,
    H1_pre = `ExAnSt1_Given_your_current_understanding_how_optimistic_are_you_about_improvement_in_symptoms_in_patients_presenting_with_anxiety_disorder`,
    H1_post = `ExAnSt21_Given_the_information_provided_in_this_education_session_how_optimistic_are_you_about_improvement_in_symptoms_in_patients_presenting_with_anxiety_disorder`
  ) %>%
  # Reshape the data from wide (pre/post in different columns) to long
  tidyr::pivot_longer(
    cols = c(H1_pre, H1_post),
    names_to = "time",
    values_to = "score"
  ) %>%
  # Make the 'time' variable a factor to ensure correct ordering on the x-axis
  mutate(time = factor(time, levels = c("H1_pre", "H1_post"), labels = c("Pre-Workshop", "Post-Workshop")))

# --- STEP 2: CALCULATE SUMMARY STATISTICS (MEAN AND 95% CI) ---

summary_data_h1 <- plot_data_h1_long %>%
  # Group by the experimental arm and the time point
  group_by(EduType, time) %>%
  # Calculate the mean, count, standard deviation, and 95% CI
  summarise(
    n = n(),
    mean_score = mean(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    se_score = sd_score / sqrt(n),
    ci_lower = mean_score - 1.96 * se_score,
    ci_upper = mean_score + 1.96 * se_score,
    .groups = 'drop' # Drop the grouping
  )

# --- STEP 3: BUILD THE PLOT WITH GGPLOT2 ---

# We build the plot in layers, starting with the raw data points and adding the summary stats on top.
ggplot(summary_data_h1, aes(x = time, y = mean_score, group = EduType, color = EduType)) +
  
  # 1. Add individual raw data points with jitter and transparency
  #    We specify the 'data' argument here to use the long-format raw data.
  geom_jitter(data = plot_data_h1_long, aes(y = score), width = 0.2, alpha = 0.25, size = 1.8) +
  
  # 2. Add the error bars representing the 95% CI
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.08, linewidth = 0.8) +
  
  # 3. Add lines connecting the mean scores
  geom_line(linewidth = 1.2) +
  
  # 4. Add large dots for the mean scores
  geom_point(size = 4) +
  
  # 5. Customize labels, colors, and theme for a professional look
  labs(
    title = "Optimism About Patient Recovery (H1)",
    subtitle = "Change in clinician optimism from pre- to post-workshop",
    y = "Optimism Score (1-7)",
    x = "Time Point",
    color = "Education Type" # Legend title
  ) +
  scale_color_manual(
    values = c("Gen" = "#0072B2", "Evo" = "#D55E00"), # Blue and Orange/Red
    labels = c("Gen" = "Genetic", "Evo" = "Evolutionary")
  ) +
  coord_cartesian(ylim = c(1, 7)) + # Set y-axis limits to match your 1-7 scale
  theme_minimal(base_size = 15) # Use a clean theme


########                     #########
########   GRAPH FOR H2      #########
########                     #########

# --- STEP 1: PREPARE THE DATA FOR PLOTTING H2 ---

# We select the relevant H2 columns, rename them for clarity, and pivot to a long format.
plot_data_h2_long <- psy_data %>%
  mutate(participant_id = row_number()) %>%
  # ***Select H2 pre and post columns***
  select(
    participant_id,
    EduType,
    H2_pre = `ExAnSt2_Given_your_current_understanding_how_willing_do_you_think_patients_would_be_to_share_anxiety_disorder_diagnoses`,
    H2_post = `ExAnSt22_Given_the_information_provided_in_this_education_session_how_willing_do_you_think_patients_would_be_to_share_anxiety_disorder_diagnoses`
  ) %>%
  # Reshape the data from wide to long
  tidyr::pivot_longer(
    cols = c(H2_pre, H2_post),
    names_to = "time",
    values_to = "score"
  ) %>%
  # Set the correct order and labels for the x-axis
  mutate(time = factor(time, levels = c("H2_pre", "H2_post"), labels = c("Pre-Education", "Post-Education")))

# --- STEP 2: CALCULATE SUMMARY STATISTICS (MEAN AND 95% CI) ---

# This logic remains the same, but now it operates on the H2 data.
summary_data_h2 <- plot_data_h2_long %>%
  group_by(EduType, time) %>%
  summarise(
    n = n(),
    mean_score = mean(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    se_score = sd_score / sqrt(n),
    ci_lower = mean_score - 1.96 * se_score,
    ci_upper = mean_score + 1.96 * se_score,
    .groups = 'drop'
  )

# --- STEP 3: BUILD THE PLOT WITH GGPLOT2 ---

# We use this data to build the plot for H2.
ggplot(summary_data_h2, aes(x = time, y = mean_score, group = EduType, color = EduType)) +
  
  # 1. Add individual raw data points with jitter and transparency
  geom_jitter(data = plot_data_h2_long, aes(y = score), width = 0.2, alpha = 0.25, size = 1.8) +
  
  # 2. Add the error bars representing the 95% CI
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.08, linewidth = 0.8) +
  
  # 3. Add lines connecting the mean scores
  geom_line(linewidth = 1.2) +
  
  # 4. Add large dots for the mean scores
  geom_point(size = 4) +
  
  # 5. ***Customize labels for H2***
  labs(
    title = "Willingness to Share Diagnosis (H2)",
    subtitle = "Change in clinician perception from pre- to post-education",
    y = "Willingness Score (1-7)",
    x = "Time Point",
    color = "Education Type"
  ) +
  scale_color_manual(
    values = c("Gen" = "#0072B2", "Evo" = "#D55E00"),
    labels = c("Gen" = "Genetic", "Evo" = "Evolutionary")
  ) +
  coord_cartesian(ylim = c(1, 7)) + 
  theme_minimal(base_size = 15)



########                     #########
########   GRAPH FOR H3      #########
########                     #########

# --- STEP 1: PREPARE THE DATA FOR PLOTTING H3 ---

# We select the relevant H3 columns, rename them for clarity, and pivot to a long format.
plot_data_h3_long <- psy_data %>%
  mutate(participant_id = row_number()) %>%
  # ***Select H3 pre and post columns***
  select(
    participant_id,
    EduType,
    H3_pre = `ExAnSt8_Given_current_public_knowledge_of_anxiety_how_willing_do_you_think_people_are_to_seek_psychiatric_help_for_anxiety`,
    H3_post = `ExAnSt28_If_the_information_provided_in_this_education_session_was_publicly_known_how_willing_do_you_think_people_would_be_to_seek_psychiatric_help_for_anxiety`
  ) %>%
  # Reshape the data from wide to long
  tidyr::pivot_longer(
    cols = c(H3_pre, H3_post),
    names_to = "time",
    values_to = "score"
  ) %>%
  # Set the correct order and labels for the x-axis
  mutate(time = factor(time, levels = c("H3_pre", "H3_post"), labels = c("Pre-Education", "Post-Education")))

# --- STEP 2: CALCULATE SUMMARY STATISTICS (MEAN AND 95% CI) ---

# This logic remains the same, but now it operates on the H3 data.
summary_data_h3 <- plot_data_h3_long %>%
  group_by(EduType, time) %>%
  summarise(
    n = n(),
    mean_score = mean(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    se_score = sd_score / sqrt(n),
    ci_lower = mean_score - 1.96 * se_score,
    ci_upper = mean_score + 1.96 * se_score,
    .groups = 'drop'
  )

# --- STEP 3: BUILD THE PLOT WITH GGPLOT2 ---

# We build the plot for H3.
ggplot(summary_data_h3, aes(x = time, y = mean_score, group = EduType, color = EduType)) +
  
  # 1. Add individual raw data points with jitter and transparency
  geom_jitter(data = plot_data_h3_long, aes(y = score), width = 0.2, alpha = 0.25, size = 1.8) +
  
  # 2. Add the error bars representing the 95% CI
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.08, linewidth = 0.8) +
  
  # 3. Add lines connecting the mean scores
  geom_line(linewidth = 1.2) +
  
  # 4. Add large dots for the mean scores
  geom_point(size = 4) +
  
  # 5. ***Customize labels for H3***
  labs(
    title = "Willingness to Seek Help (H3)",
    subtitle = "Change in clinician perception from pre- to post-education",
    y = "Willingness Score (1-7)",
    x = "Time Point",
    color = "Education Type"
  ) +
  scale_color_manual(
    values = c("Gen" = "#0072B2", "Evo" = "#D55E00"),
    labels = c("Gen" = "Genetic", "Evo" = "Evolutionary")
  ) +
  coord_cartesian(ylim = c(1, 7)) + 
  theme_minimal(base_size = 15)




########                     #########
########   GRAPH FOR H4      #########
########                     #########

# --- STEP 1: PREPARE THE DATA FOR PLOTTING H4 ---

# We select the relevant H4 columns, rename them for clarity, and pivot to a long format.
plot_data_h4_long <- psy_data %>%
  mutate(participant_id = row_number()) %>%
  # ***Select H4 pre and post columns***
  select(
    participant_id,
    EduType,
    H4_pre = `ExAnSt4_Given_your_current_understanding_how_effective_do_you_think_psychosocial_interventions_will_be_in_improving_outcomes_in_anxiety_disorder`,
    H4_post = `ExAnSt24_Given_the_information_provided_in_this_education_session_how_effective_do_you_think_psychosocial_interventions_will_be_in_improving_outcomes_in_anxiety_disorder`
  ) %>%
  # Reshape the data from wide to long
  tidyr::pivot_longer(
    cols = c(H4_pre, H4_post),
    names_to = "time",
    values_to = "score"
  ) %>%
  # Set the correct order and labels for the x-axis
  mutate(time = factor(time, levels = c("H4_pre", "H4_post"), labels = c("Pre-Education", "Post-Education")))

# --- STEP 2: CALCULATE SUMMARY STATISTICS (MEAN AND 95% CI) ---

# This logic remains the same, but now it operates on the H4 data.
summary_data_h4 <- plot_data_h4_long %>%
  group_by(EduType, time) %>%
  summarise(
    n = n(),
    mean_score = mean(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    se_score = sd_score / sqrt(n),
    ci_lower = mean_score - 1.96 * se_score,
    ci_upper = mean_score + 1.96 * se_score,
    .groups = 'drop'
  )

# --- STEP 3: BUILD THE PLOT WITH GGPLOT2 ---

# We use the finalized specifications to build the plot for H4.
ggplot(summary_data_h4, aes(x = time, y = mean_score, group = EduType, color = EduType)) +
  
  # 1. Add individual raw data points with jitter and transparency
  geom_jitter(data = plot_data_h4_long, aes(y = score), width = 0.2, alpha = 0.25, size = 1.8) +
  
  # 2. Add the error bars representing the 95% CI
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.08, linewidth = 0.8) +
  
  # 3. Add lines connecting the mean scores
  geom_line(linewidth = 1.2) +
  
  # 4. Add large dots for the mean scores
  geom_point(size = 4) +
  
  # 5. ***Customize labels for H4***
  labs(
    title = "Efficacy of Psychosocial Interventions (H4)",
    subtitle = "Change in clinician belief from pre- to post-education",
    y = "Efficacy Score (1-7)",
    x = "Time Point",
    color = "Education Type"
  ) +
  scale_color_manual(
    values = c("Gen" = "#0072B2", "Evo" = "#D55E00"),
    labels = c("Gen" = "Genetic", "Evo" = "Evolutionary")
  ) +
  coord_cartesian(ylim = c(1, 7)) + 
  theme_minimal(base_size = 15)




########                                 #########
########   COMBINING ALL FOUR GRAPHS     #########
########                                 #########

######## Uses the same code as above but edits some Y axis and other aspects to combine ######

# --- SETUP: LOAD NECESSARY LIBRARIES ---
# Make sure you have patchwork installed: install.packages("patchwork")
library(dplyr)
library(ggplot2)
library(patchwork)

# --- GRAPH FOR H1: Optimism About Patient Recovery ---

# 1. Prepare data for H1
plot_data_h1_long <- psy_data %>%
  mutate(participant_id = row_number()) %>%
  select(
    participant_id, EduType,
    H1_pre = `ExAnSt1_Given_your_current_understanding_how_optimistic_are_you_about_improvement_in_symptoms_in_patients_presenting_with_anxiety_disorder`,
    H1_post = `ExAnSt21_Given_the_information_provided_in_this_education_session_how_optimistic_are_you_about_improvement_in_symptoms_in_patients_presenting_with_anxiety_disorder`
  ) %>%
  tidyr::pivot_longer(cols = c(H1_pre, H1_post), names_to = "time", values_to = "score") %>%
  mutate(time = factor(time, levels = c("H1_pre", "H1_post"), labels = c("Pre-Education", "Post-Education")))

# 2. Calculate summary stats for H1
summary_data_h1 <- plot_data_h1_long %>%
  group_by(EduType, time) %>%
  summarise(
    n = n(), mean_score = mean(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE), se_score = sd_score / sqrt(n),
    ci_lower = mean_score - 1.96 * se_score, ci_upper = mean_score + 1.96 * se_score,
    .groups = 'drop'
  )

# 3. Build and SAVE the plot for H1 to a variable
plot_h1 <- ggplot(summary_data_h1, aes(x = time, y = mean_score, group = EduType, color = EduType)) +
  geom_jitter(data = plot_data_h1_long, aes(y = score), width = 0.2, alpha = 0.25, size = 1.8) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.08, linewidth = 0.8) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  labs(
    title = "Optimism About Patient Recovery (H1)",
    # ***Standardized Y-axis label***
    y = "Clinician Rating (1-7)", 
    x = NULL 
  ) +
  scale_color_manual(values = c("Gen" = "#0072B2", "Evo" = "#D55E00"), labels = c("Gen" = "Genetic", "Evo" = "Evolutionary")) +
  coord_cartesian(ylim = c(1, 7)) +
  theme_minimal(base_size = 15)


# --- GRAPH FOR H2: Willingness to Share Diagnosis ---

# 1. Prepare data for H2
plot_data_h2_long <- psy_data %>%
  mutate(participant_id = row_number()) %>%
  select(
    participant_id, EduType,
    H2_pre = `ExAnSt2_Given_your_current_understanding_how_willing_do_you_think_patients_would_be_to_share_anxiety_disorder_diagnoses`,
    H2_post = `ExAnSt22_Given_the_information_provided_in_this_education_session_how_willing_do_you_think_patients_would_be_to_share_anxiety_disorder_diagnoses`
  ) %>%
  tidyr::pivot_longer(cols = c(H2_pre, H2_post), names_to = "time", values_to = "score") %>%
  mutate(time = factor(time, levels = c("H2_pre", "H2_post"), labels = c("Pre-Education", "Post-Education")))

# 2. Calculate summary stats for H2
summary_data_h2 <- plot_data_h2_long %>%
  group_by(EduType, time) %>%
  summarise(
    n = n(), mean_score = mean(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE), se_score = sd_score / sqrt(n),
    ci_lower = mean_score - 1.96 * se_score, ci_upper = mean_score + 1.96 * se_score,
    .groups = 'drop'
  )

# 3. Build and SAVE the plot for H2 to a variable
plot_h2 <- ggplot(summary_data_h2, aes(x = time, y = mean_score, group = EduType, color = EduType)) +
  geom_jitter(data = plot_data_h2_long, aes(y = score), width = 0.2, alpha = 0.25, size = 1.8) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.08, linewidth = 0.8) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  labs(
    title = "Willingness to Share Diagnosis (H2)",
    # ***Removed Y-axis label***
    y = NULL, 
    x = NULL 
  ) +
  scale_color_manual(values = c("Gen" = "#0072B2", "Evo" = "#D55E00"), labels = c("Gen" = "Genetic", "Evo" = "Evolutionary")) +
  coord_cartesian(ylim = c(1, 7)) +
  theme_minimal(base_size = 15)


# --- GRAPH FOR H3: Willingness to Seek Help ---

# 1. Prepare data for H3
plot_data_h3_long <- psy_data %>%
  mutate(participant_id = row_number()) %>%
  select(
    participant_id, EduType,
    H3_pre = `ExAnSt8_Given_current_public_knowledge_of_anxiety_how_willing_do_you_think_people_are_to_seek_psychiatric_help_for_anxiety`,
    H3_post = `ExAnSt28_If_the_information_provided_in_this_education_session_was_publicly_known_how_willing_do_you_think_people_would_be_to_seek_psychiatric_help_for_anxiety`
  ) %>%
  tidyr::pivot_longer(cols = c(H3_pre, H3_post), names_to = "time", values_to = "score") %>%
  mutate(time = factor(time, levels = c("H3_pre", "H3_post"), labels = c("Pre-Education", "Post-Education")))

# 2. Calculate summary stats for H3
summary_data_h3 <- plot_data_h3_long %>%
  group_by(EduType, time) %>%
  summarise(
    n = n(), mean_score = mean(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE), se_score = sd_score / sqrt(n),
    ci_lower = mean_score - 1.96 * se_score, ci_upper = mean_score + 1.96 * se_score,
    .groups = 'drop'
  )

# 3. Build and SAVE the plot for H3 to a variable
plot_h3 <- ggplot(summary_data_h3, aes(x = time, y = mean_score, group = EduType, color = EduType)) +
  geom_jitter(data = plot_data_h3_long, aes(y = score), width = 0.2, alpha = 0.25, size = 1.8) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.08, linewidth = 0.8) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  labs(
    title = "Willingness to Seek Help (H3)",
    # ***Standardized Y-axis label***
    y = "Clinician Rating (1-7)", 
    x = "Time Point"
  ) +
  scale_color_manual(values = c("Gen" = "#0072B2", "Evo" = "#D55E00"), labels = c("Gen" = "Genetic", "Evo" = "Evolutionary")) +
  coord_cartesian(ylim = c(1, 7)) +
  theme_minimal(base_size = 15)


# --- GRAPH FOR H4: Efficacy of Psychosocial Interventions ---

# 1. Prepare data for H4
plot_data_h4_long <- psy_data %>%
  mutate(participant_id = row_number()) %>%
  select(
    participant_id, EduType,
    H4_pre = `ExAnSt4_Given_your_current_understanding_how_effective_do_you_think_psychosocial_interventions_will_be_in_improving_outcomes_in_anxiety_disorder`,
    H4_post = `ExAnSt24_Given_the_information_provided_in_this_education_session_how_effective_do_you_think_psychosocial_interventions_will_be_in_improving_outcomes_in_anxiety_disorder`
  ) %>%
  tidyr::pivot_longer(cols = c(H4_pre, H4_post), names_to = "time", values_to = "score") %>%
  mutate(time = factor(time, levels = c("H4_pre", "H4_post"), labels = c("Pre-Education", "Post-Education")))

# 2. Calculate summary stats for H4
summary_data_h4 <- plot_data_h4_long %>%
  group_by(EduType, time) %>%
  summarise(
    n = n(), mean_score = mean(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE), se_score = sd_score / sqrt(n),
    ci_lower = mean_score - 1.96 * se_score, ci_upper = mean_score + 1.96 * se_score,
    .groups = 'drop'
  )

# 3. Build and SAVE the plot for H4 to a variable
plot_h4 <- ggplot(summary_data_h4, aes(x = time, y = mean_score, group = EduType, color = EduType)) +
  geom_jitter(data = plot_data_h4_long, aes(y = score), width = 0.2, alpha = 0.25, size = 1.8) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.08, linewidth = 0.8) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  labs(
    title = "Efficacy of Psychosocial Interventions (H4)",
    # ***Removed Y-axis label***
    y = NULL, 
    x = "Time Point"
  ) +
  scale_color_manual(values = c("Gen" = "#0072B2", "Evo" = "#D55E00"), labels = c("Gen" = "Genetic", "Evo" = "Evolutionary")) +
  coord_cartesian(ylim = c(1, 7)) +
  theme_minimal(base_size = 15)


# --- STEP 4: COMBINE THE PLOTS WITH PATCHWORK ---

# Arrange the four plots into a 2x2 grid and add an overall title.
(plot_h1 + plot_h2) / (plot_h3 + plot_h4) + 
  plot_layout(guides = 'collect') & 
  plot_annotation(
    title = 'Effect of Psychoeducation Type on Clinician Attitudes Towards Anxiety',
    theme = theme(plot.title = element_text(size = 20, hjust = 0.5))
  ) &
  theme(legend.position = 'bottom')




#-----------------------------------------------------------------------------
# STEP 6: WITHIN-GROUP PRE-POST ANALYSIS FOR EVOLUTIONARY & GENETICS ARMS
#-----------------------------------------------------------------------------
# This section runs a repeated measures ordinal regression to specifically
# assess the pre-post change *within* the evolutionary education group and the genetics arm separately

# It will also assess whether the effect of the treatment differs across each session across both interventions

##################################### H1 ##############################
#######################################################################
#######################################################################


# --- 1. Prepare Data for H1 ---


# Reshape the data from wide to long format for the repeated measures model.
long_data_h1 <- analysis_data %>%
  # Select the participant identifier, covariates, and the raw scores for H1
  select(
    EduType,
    ParticipantID, 
    SessionID, 
    Delivery,
    Pre_Score = `ExAnSt1_Given_your_current_understanding_how_optimistic_are_you_about_improvement_in_symptoms_in_patients_presenting_with_anxiety_disorder`,
    Post_Score = H1_post_ord
  ) %>%
  # Pivot the pre and post scores into a long format
  pivot_longer(
    cols = c(Pre_Score, Post_Score),
    names_to = "Time",
    values_to = "Score"
  ) %>%
  # Convert the Time and Score columns into ordered factors for the model
  mutate(
    Time = factor(Time, levels = c("Pre_Score", "Post_Score"))
  )

# --- 2. Run the Bayesian Repeated Measures Models for H1 ---

# Run the logit model.
prepost_h1_model <- brm(
  formula = bf(Score ~ 1 + Time + Delivery + (Time + Delivery | EduType/SessionID) + (1 | EduType:SessionID:ParticipantID)),
  data = long_data_h1,
  family = cumulative(link = "logit"),
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(normal(0, 1.5), class = "Intercept"),
    prior(exponential(1), class = "sd")
  ),
  cores = 4, chains = 4, iter = 4000, 
  control = list(adapt_delta = 0.99)
)

# --- Save and Load the logit Model ---
saveRDS(prepost_h1_model, file = "prepost_h1_model.rds")
prepost_h1_model <- readRDS("prepost_h1_model.rds")

## Inspecting the model now
summary(prepost_h1_model)

## Drawing posteriors from the model
prepost_h1_model_draws <- as_draws_df(prepost_h1_model)


### Okay, now our goal is to extract the effect of moving from pre-score to post-score for each arm separately
### and plot its posterior
##############################################################################################################

## First we select all the random deviations of the time effect for each arm
post_slope_h1_arm <- prepost_h1_model_draws %>%
  select(starts_with("r_EduType[")) %>%
  select(contains("TimePost_Score"))

## Select the averaged out fixed effect across both arms
fixed_slope_h1 <- prepost_h1_model_draws$b_TimePost_Score

## Now we calculate the effect of moving from pre- to post-score for each arm separately
post_slope_h1_arm <- post_slope_h1_arm %>%
  mutate(
    Evo = fixed_slope_h1 + `r_EduType[Evo,TimePost_Score]`,
    Gen = fixed_slope_h1 + `r_EduType[Gen,TimePost_Score]`
  )

## Now we convert each to an Odds Ratio
post_slope_h1_arm <- post_slope_h1_arm %>%
  mutate(
    Evo_OR = exp(Evo),
    Gen_OR = exp(Gen)
  )

summary(post_slope_h1_arm)

### Okay, now we plot and compare the effects across each arm

## First, we compute all the summary stats
post_slope_H1_summary <- post_slope_h1_arm %>%
  summarise(
    Evo_OR_mean = mean(Evo_OR),
    Gen_OR_mean = mean(Gen_OR),
    Evo_OR_lower = quantile(Evo_OR, probs = 0.025),  # 2.5th percentile
    Evo_OR_upper = quantile(Evo_OR, probs = 0.975),  # 97.5th percentile
    Gen_OR_lower = quantile(Gen_OR, probs = 0.025),  # 2.5th percentile
    Gen_OR_upper = quantile(Gen_OR, probs = 0.975)   # 97.5th percentile
  )

## Then we format the data to a format suitable for plotting
post_slope_h1_long <- data.frame(
  Arm = rep(c("Evo", "Gen"), each = 1),
  Mean = c(post_slope_H1_summary$Evo_OR_mean, post_slope_H1_summary$Gen_OR_mean),
  Lower = c(post_slope_H1_summary$Evo_OR_lower, post_slope_H1_summary$Gen_OR_lower),
  Upper = c(post_slope_H1_summary$Evo_OR_upper, post_slope_H1_summary$Gen_OR_upper)
)


# Plotting the 95% credible intervals now
ggplot(post_slope_h1_long, aes(x = Arm, y = Mean)) +
  geom_point(size = 4) +  # Add points
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +  # Add vertical error bars
  theme_minimal() +
  labs(x = "Category", y = "Odds Ratio", title = "Odds Ratios with Error Bars") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


### Okay, now our goal is to extract the effect of moving from pre-score to post-score for each session
##############################################################################################################
post_slope_H1_session <- prepost_h1_model_draws %>%
  select(starts_with("r_EduType:SessionID[")) %>%
  select(contains("TimePost_Score"))

## Now let's calculate the fixed effect for the evolutionary arm and the genetics arm
fixed_slope_h1_evo <- prepost_h1_model_draws$b_TimePost_Score + prepost_h1_model_draws$`r_EduType[Evo,TimePost_Score]`
fixed_slope_h1_gen <- prepost_h1_model_draws$b_TimePost_Score + prepost_h1_model_draws$`r_EduType[Gen,TimePost_Score]`

## Now we calculate differing effects across each session in the evolutionary arm
post_slope_H1_session <- post_slope_H1_session %>%
  mutate(
    Session_1 = fixed_slope_h1_evo + `r_EduType:SessionID[Evo_Evo_1,TimePost_Score]`,
    Session_2 = fixed_slope_h1_evo + `r_EduType:SessionID[Evo_Evo_2,TimePost_Score]`,
    Session_3 = fixed_slope_h1_evo + `r_EduType:SessionID[Evo_Evo_3,TimePost_Score]`,
    Session_4 = fixed_slope_h1_evo + `r_EduType:SessionID[Evo_Evo_4,TimePost_Score]`,
    Session_5 = fixed_slope_h1_evo + `r_EduType:SessionID[Evo_Evo_5,TimePost_Score]`,
    Session_6 = fixed_slope_h1_evo + `r_EduType:SessionID[Evo_Evo_6,TimePost_Score]`,
    Session_7 = fixed_slope_h1_evo + `r_EduType:SessionID[Evo_Evo_7,TimePost_Score]`,
    Session_8 = fixed_slope_h1_evo + `r_EduType:SessionID[Evo_Evo_8,TimePost_Score]`,
    Session_9 = fixed_slope_h1_evo + `r_EduType:SessionID[Evo_Evo_9,TimePost_Score]`,
    Session_10 = fixed_slope_h1_evo + `r_EduType:SessionID[Evo_Evo_10,TimePost_Score]`,
  )

## Now we calculate differing effects across each session in the genetics arm
post_slope_H1_session <- post_slope_H1_session %>%
  mutate(
    Gen_Session_1 = fixed_slope_h1_gen + `r_EduType:SessionID[Gen_Gen_10,TimePost_Score]`,
    Gen_Session_2 = fixed_slope_h1_gen + `r_EduType:SessionID[Gen_Gen_2,TimePost_Score]`,
    Gen_Session_3 = fixed_slope_h1_gen + `r_EduType:SessionID[Gen_Gen_3,TimePost_Score]`,
    Gen_Session_4 = fixed_slope_h1_gen + `r_EduType:SessionID[Gen_Gen_4,TimePost_Score]`,
    Gen_Session_5 = fixed_slope_h1_gen + `r_EduType:SessionID[Gen_Gen_5,TimePost_Score]`,
    Gen_Session_6 = fixed_slope_h1_gen + `r_EduType:SessionID[Gen_Gen_6,TimePost_Score]`,
    Gen_Session_7 = fixed_slope_h1_gen + `r_EduType:SessionID[Gen_Gen_7,TimePost_Score]`,
    Gen_Session_8 = fixed_slope_h1_gen + `r_EduType:SessionID[Gen_Gen_8,TimePost_Score]`,
    Gen_Session_9 = fixed_slope_h1_gen + `r_EduType:SessionID[Gen_Gen_9,TimePost_Score]`,
  )

## Now we convert all of them to Odds Ratios
post_slope_H1_session <- post_slope_H1_session %>%
  mutate(
    Session_1_OR = exp(Session_1),
    Session_2_OR = exp(Session_2),
    Session_3_OR = exp(Session_3),
    Session_4_OR = exp(Session_4),
    Session_5_OR = exp(Session_5),
    Session_6_OR = exp(Session_6),
    Session_7_OR = exp(Session_7),
    Session_8_OR = exp(Session_8),
    Session_9_OR = exp(Session_9),
    Session_10_OR = exp(Session_10),
    Gen_Session_1_OR = exp(Gen_Session_1),
    Gen_Session_2_OR = exp(Gen_Session_2),
    Gen_Session_3_OR = exp(Gen_Session_3),
    Gen_Session_4_OR = exp(Gen_Session_4),
    Gen_Session_5_OR = exp(Gen_Session_5),
    Gen_Session_6_OR = exp(Gen_Session_6),
    Gen_Session_7_OR = exp(Gen_Session_7),
    Gen_Session_8_OR = exp(Gen_Session_8),
    Gen_Session_9_OR = exp(Gen_Session_9)
  )


### Okay, now we plot and compare the effects across each sessions

## Names of columnns for which we need to compute summary statistics
session_columns <- c(
  "Session_1_OR", "Session_2_OR", "Session_3_OR", "Session_4_OR", "Session_5_OR",
  "Session_6_OR", "Session_7_OR", "Session_8_OR", "Session_9_OR", "Session_10_OR",
  "Gen_Session_1_OR", "Gen_Session_2_OR", "Gen_Session_3_OR", "Gen_Session_4_OR",
  "Gen_Session_5_OR", "Gen_Session_6_OR", "Gen_Session_7_OR", "Gen_Session_8_OR",
  "Gen_Session_9_OR"
)

# Function to compute the mean and 95% credible interval for each column
compute_ci <- function(column_name) {
  column_data <- post_slope_H1_session[[column_name]]
  
  # Calculate the mean and 95% credible interval (using quantiles)
  mean_value <- mean(column_data, na.rm = TRUE)
  lower_ci <- quantile(column_data, 0.025, na.rm = TRUE)
  upper_ci <- quantile(column_data, 0.975, na.rm = TRUE)
  
  tibble(
    Session = column_name,
    mean = mean_value,
    lower_95_CI = lower_ci,
    upper_95_CI = upper_ci
  )
}

# Apply the function across all session and gen session OR columns
session_H1 <- map_dfr(session_columns, compute_ci)

# Plotting the 95% credible intervals now
ggplot(session_H1, aes(x = Session, y = mean)) +
  geom_point(size = 4) +  # Add points
  geom_errorbar(aes(ymin = lower_95_CI, ymax = upper_95_CI), width = 0.2) +  # Add vertical error bars
  theme_minimal() +
  labs(x = "Session", y = "Odds Ratio", title = "Odds Ratios with Error Bars") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red")


# ---  Run the Bayesian Repeated Measures Models for H1 (Probit) ---


# Run the probit model.
prepost_h1_probit_model <- brm(
  formula = bf(Score ~ 1 + Time + Delivery + (Time + Delivery | EduType/SessionID) + (1 | EduType:SessionID:ParticipantID)),
  data = long_data_h1,
  family = cumulative(link = "probit"),
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(normal(0, 1.5), class = "Intercept"),
    prior(exponential(1), class = "sd")
  ),
  cores = 4, chains = 4, iter = 4000,
  control = list(adapt_delta = 0.99)
)

#
# --- Save and Load the Probit Model ---
saveRDS(prepost_h1_probit_model, file = "prepost_h1_probit_model.rds")
prepost_h1_probit_model <- readRDS("prepost_h1_probit_model.rds")

## Drawing posteriors from the model
prepost_h1_probit_model_draws <- as_draws_df(prepost_h1_probit_model)

# --- 3b. Extract and Plot Arm-Level Effects for H1 (Probit) ---
# For probit models, coefficients are effect sizes (SD units). No exp() needed.

## Calculate the effect for each arm separately
post_slope_h1_probit_arm <- prepost_h1_probit_model_draws %>%
  mutate(
    Evo_SD = b_TimePost_Score + `r_EduType[Evo,TimePost_Score]`,
    Gen_SD = b_TimePost_Score + `r_EduType[Gen,TimePost_Score]`
  )

## Compute summary stats for plotting
post_slope_h1_probit_summary <- post_slope_h1_probit_arm %>%
  summarise(
    Evo_SD_mean = mean(Evo_SD),
    Gen_SD_mean = mean(Gen_SD),
    Evo_SD_lower = quantile(Evo_SD, probs = 0.025),
    Evo_SD_upper = quantile(Evo_SD, probs = 0.975),
    Gen_SD_lower = quantile(Gen_SD, probs = 0.025),
    Gen_SD_upper = quantile(Gen_SD, probs = 0.975)
  )


summary(post_slope_h1_probit_summary)


## Format data for plotting
post_slope_h1_probit_long <- data.frame(
  Arm = rep(c("Evo", "Gen"), each = 1),
  Mean = c(post_slope_h1_probit_summary$Evo_SD_mean, post_slope_h1_probit_summary$Gen_SD_mean),
  Lower = c(post_slope_h1_probit_summary$Evo_SD_lower, post_slope_h1_probit_summary$Gen_SD_lower),
  Upper = c(post_slope_h1_probit_summary$Evo_SD_upper, post_slope_h1_probit_summary$Gen_SD_upper)
)

# Plot the 95% credible intervals for H1 Probit
ggplot(post_slope_h1_probit_long, aes(x = Arm, y = Mean)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  theme_minimal() +
  labs(x = "Arm", y = "Standard Deviation Change", title = "H1 (Probit): Pre-Post Effect Size") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")


# --- 4b. Extract and Plot Session-Level Effects for H1 (Probit) ---

## Select session-level random effects
post_slope_h1_probit_session <- prepost_h1_probit_model_draws %>%
  select(starts_with("r_EduType:SessionID[")) %>%
  select(contains("TimePost_Score"))

## Calculate the arm-level effects on the SD scale
fixed_slope_h1_probit_evo <- prepost_h1_probit_model_draws$b_TimePost_Score + prepost_h1_probit_model_draws$`r_EduType[Evo,TimePost_Score]`
fixed_slope_h1_probit_gen <- prepost_h1_probit_model_draws$b_TimePost_Score + prepost_h1_probit_model_draws$`r_EduType[Gen,TimePost_Score]`

## Calculate differing effects across each session
post_slope_h1_probit_session <- post_slope_h1_probit_session %>%
  mutate(
    # Evolutionary Arm Sessions
    Session_1_SD = fixed_slope_h1_probit_evo + `r_EduType:SessionID[Evo_Evo_1,TimePost_Score]`,
    Session_2_SD = fixed_slope_h1_probit_evo + `r_EduType:SessionID[Evo_Evo_2,TimePost_Score]`,
    Session_3_SD = fixed_slope_h1_probit_evo + `r_EduType:SessionID[Evo_Evo_3,TimePost_Score]`,
    Session_4_SD = fixed_slope_h1_probit_evo + `r_EduType:SessionID[Evo_Evo_4,TimePost_Score]`,
    Session_5_SD = fixed_slope_h1_probit_evo + `r_EduType:SessionID[Evo_Evo_5,TimePost_Score]`,
    Session_6_SD = fixed_slope_h1_probit_evo + `r_EduType:SessionID[Evo_Evo_6,TimePost_Score]`,
    Session_7_SD = fixed_slope_h1_probit_evo + `r_EduType:SessionID[Evo_Evo_7,TimePost_Score]`,
    Session_8_SD = fixed_slope_h1_probit_evo + `r_EduType:SessionID[Evo_Evo_8,TimePost_Score]`,
    Session_9_SD = fixed_slope_h1_probit_evo + `r_EduType:SessionID[Evo_Evo_9,TimePost_Score]`,
    Session_10_SD = fixed_slope_h1_probit_evo + `r_EduType:SessionID[Evo_Evo_10,TimePost_Score]`,
    # Genetics Arm Sessions
    Gen_Session_1_SD = fixed_slope_h1_probit_gen + `r_EduType:SessionID[Gen_Gen_10,TimePost_Score]`,
    Gen_Session_2_SD = fixed_slope_h1_probit_gen + `r_EduType:SessionID[Gen_Gen_2,TimePost_Score]`,
    Gen_Session_3_SD = fixed_slope_h1_probit_gen + `r_EduType:SessionID[Gen_Gen_3,TimePost_Score]`,
    Gen_Session_4_SD = fixed_slope_h1_probit_gen + `r_EduType:SessionID[Gen_Gen_4,TimePost_Score]`,
    Gen_Session_5_SD = fixed_slope_h1_probit_gen + `r_EduType:SessionID[Gen_Gen_5,TimePost_Score]`,
    Gen_Session_6_SD = fixed_slope_h1_probit_gen + `r_EduType:SessionID[Gen_Gen_6,TimePost_Score]`,
    Gen_Session_7_SD = fixed_slope_h1_probit_gen + `r_EduType:SessionID[Gen_Gen_7,TimePost_Score]`,
    Gen_Session_8_SD = fixed_slope_h1_probit_gen + `r_EduType:SessionID[Gen_Gen_8,TimePost_Score]`,
    Gen_Session_9_SD = fixed_slope_h1_probit_gen + `r_EduType:SessionID[Gen_Gen_9,TimePost_Score]`
  )

## Define column names for which summary statistics are needed
session_columns_h1_probit <- names(post_slope_h1_probit_session)[grepl("_SD$", names(post_slope_h1_probit_session))]

# Apply the CI function across all relevant columns for H1 Probit
session_H1_probit <- map_dfr(session_columns_h1_probit, ~compute_ci(.x, source_df = post_slope_h1_probit_session))

# Plot the 95% credible intervals now
ggplot(session_H1_probit, aes(x = Session, y = mean)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = lower_95_CI, ymax = upper_95_CI), width = 0.2) +
  theme_minimal() +
  labs(x = "Session", y = "Standard Deviation Change", title = "H1 (Probit): Session-Level Effect Sizes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")


##################################### H2 ##############################
# Hypothesis 2: Willingness for patients to share diagnosis
#######################################################################

# --- 1. Prepare Data for H2 ---

# Reshape the data from wide to long format for the repeated measures model.
long_data_h2 <- analysis_data %>%
  # Select the participant identifier, covariates, and the raw scores for H2
  select(
    EduType,
    ParticipantID, 
    SessionID, 
    Delivery,
    Pre_Score = `ExAnSt2_Given_your_current_understanding_how_willing_do_you_think_patients_would_be_to_share_anxiety_disorder_diagnoses`,
    Post_Score = H2_post_ord
  ) %>%
  # Pivot the pre and post scores into a long format
  pivot_longer(
    cols = c(Pre_Score, Post_Score),
    names_to = "Time",
    values_to = "Score"
  ) %>%
  # Convert the Time and Score columns into ordered factors for the model
  mutate(
    Time = factor(Time, levels = c("Pre_Score", "Post_Score"))
  )

# --- 2. Run the Bayesian Repeated Measures Models for H2 ---

# Run the logit model.
prepost_h2_model <- brm(
  formula = bf(Score ~ 1 + Time + Delivery + (Time + Delivery | EduType/SessionID) + (1 | EduType:SessionID:ParticipantID)),
  data = long_data_h2,
  family = cumulative(link = "logit"),
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(normal(0, 1.5), class = "Intercept"),
    prior(exponential(1), class = "sd")
  ),
  cores = 4, chains = 4, iter = 4000, 
  control = list(adapt_delta = 0.99)
)

# --- Save and Load the logit Model ---
saveRDS(prepost_h2_model, file = "prepost_h2_model.rds")
prepost_h2_model <- readRDS("prepost_h2_model.rds")

## Inspecting the model now
summary(prepost_h2_model)

## Drawing posteriors from the model
prepost_h2_model_draws <- as_draws_df(prepost_h2_model)

# --- 3. Extract and Plot Arm-Level Effects for H2 ---

## Calculate the effect for each arm separately
post_slope_h2_arm <- prepost_h2_model_draws %>%
  mutate(
    Evo = b_TimePost_Score + `r_EduType[Evo,TimePost_Score]`,
    Gen = b_TimePost_Score + `r_EduType[Gen,TimePost_Score]`
  ) %>%
  mutate(
    Evo_OR = exp(Evo),
    Gen_OR = exp(Gen)
  )

## Compute summary stats for plotting
post_slope_h2_summary <- post_slope_h2_arm %>%
  summarise(
    Evo_OR_mean = mean(Evo_OR),
    Gen_OR_mean = mean(Gen_OR),
    Evo_OR_lower = quantile(Evo_OR, probs = 0.025),
    Evo_OR_upper = quantile(Evo_OR, probs = 0.975),
    Gen_OR_lower = quantile(Gen_OR, probs = 0.025),
    Gen_OR_upper = quantile(Gen_OR, probs = 0.975)
  )

## Format data for plotting
post_slope_h2_long <- data.frame(
  Arm = rep(c("Evo", "Gen"), each = 1),
  Mean = c(post_slope_h2_summary$Evo_OR_mean, post_slope_h2_summary$Gen_OR_mean),
  Lower = c(post_slope_h2_summary$Evo_OR_lower, post_slope_h2_summary$Gen_OR_lower),
  Upper = c(post_slope_h2_summary$Evo_OR_upper, post_slope_h2_summary$Gen_OR_upper)
)

# Plot the 95% credible intervals for H2
ggplot(post_slope_h2_long, aes(x = Arm, y = Mean)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  theme_minimal() +
  labs(x = "Arm", y = "Odds Ratio", title = "H2: Pre-Post Change in Willingness to Share Diagnosis") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# --- 4. Extract and Plot Session-Level Effects for H2 ---

## Select session-level random effects
post_slope_h2_session <- prepost_h2_model_draws %>%
  select(starts_with("r_EduType:SessionID[")) %>%
  select(contains("TimePost_Score"))

## Calculate the fixed effect for the evolutionary arm and the genetics arm
fixed_slope_h2_evo <- prepost_h2_model_draws$b_TimePost_Score + prepost_h2_model_draws$`r_EduType[Evo,TimePost_Score]`
fixed_slope_h2_gen <- prepost_h2_model_draws$b_TimePost_Score + prepost_h2_model_draws$`r_EduType[Gen,TimePost_Score]`

## Calculate differing effects across each session, then convert to Odds Ratios
post_slope_h2_session <- post_slope_h2_session %>%
  mutate(
    # Evolutionary Arm Sessions
    Session_1 = fixed_slope_h2_evo + `r_EduType:SessionID[Evo_Evo_1,TimePost_Score]`,
    Session_2 = fixed_slope_h2_evo + `r_EduType:SessionID[Evo_Evo_2,TimePost_Score]`,
    Session_3 = fixed_slope_h2_evo + `r_EduType:SessionID[Evo_Evo_3,TimePost_Score]`,
    Session_4 = fixed_slope_h2_evo + `r_EduType:SessionID[Evo_Evo_4,TimePost_Score]`,
    Session_5 = fixed_slope_h2_evo + `r_EduType:SessionID[Evo_Evo_5,TimePost_Score]`,
    Session_6 = fixed_slope_h2_evo + `r_EduType:SessionID[Evo_Evo_6,TimePost_Score]`,
    Session_7 = fixed_slope_h2_evo + `r_EduType:SessionID[Evo_Evo_7,TimePost_Score]`,
    Session_8 = fixed_slope_h2_evo + `r_EduType:SessionID[Evo_Evo_8,TimePost_Score]`,
    Session_9 = fixed_slope_h2_evo + `r_EduType:SessionID[Evo_Evo_9,TimePost_Score]`,
    Session_10 = fixed_slope_h2_evo + `r_EduType:SessionID[Evo_Evo_10,TimePost_Score]`,
    # Genetics Arm Sessions (following H1 example naming)
    Gen_Session_1 = fixed_slope_h2_gen + `r_EduType:SessionID[Gen_Gen_10,TimePost_Score]`,
    Gen_Session_2 = fixed_slope_h2_gen + `r_EduType:SessionID[Gen_Gen_2,TimePost_Score]`,
    Gen_Session_3 = fixed_slope_h2_gen + `r_EduType:SessionID[Gen_Gen_3,TimePost_Score]`,
    Gen_Session_4 = fixed_slope_h2_gen + `r_EduType:SessionID[Gen_Gen_4,TimePost_Score]`,
    Gen_Session_5 = fixed_slope_h2_gen + `r_EduType:SessionID[Gen_Gen_5,TimePost_Score]`,
    Gen_Session_6 = fixed_slope_h2_gen + `r_EduType:SessionID[Gen_Gen_6,TimePost_Score]`,
    Gen_Session_7 = fixed_slope_h2_gen + `r_EduType:SessionID[Gen_Gen_7,TimePost_Score]`,
    Gen_Session_8 = fixed_slope_h2_gen + `r_EduType:SessionID[Gen_Gen_8,TimePost_Score]`,
    Gen_Session_9 = fixed_slope_h2_gen + `r_EduType:SessionID[Gen_Gen_9,TimePost_Score]`
  ) %>%
  mutate(
    across(c(starts_with("Session_"), starts_with("Gen_Session_")), exp, .names = "{.col}_OR")
  )

## Define column names for which summary statistics are needed
session_columns_h2 <- names(post_slope_h2_session)[grepl("_OR$", names(post_slope_h2_session))]

# Define function to compute the mean and 95% credible interval
compute_ci <- function(column_name, source_df) {
  column_data <- source_df[[column_name]]
  
  mean_value <- mean(column_data, na.rm = TRUE)
  lower_ci <- quantile(column_data, 0.025, na.rm = TRUE)
  upper_ci <- quantile(column_data, 0.975, na.rm = TRUE)
  
  tibble(
    Session = column_name,
    mean = mean_value,
    lower_95_CI = lower_ci,
    upper_95_CI = upper_ci
  )
}

# Apply the function across all relevant columns for H2
session_H2 <- map_dfr(session_columns_h2, ~compute_ci(.x, source_df = post_slope_h2_session))

# Plot the 95% credible intervals now
ggplot(session_H2, aes(x = Session, y = mean)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = lower_95_CI, ymax = upper_95_CI), width = 0.2) +
  theme_minimal() +
  labs(x = "Session", y = "Odds Ratio", title = "H2: Session-Level Odds Ratios with Error Bars") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red")


# --- 2b. Run the Bayesian Repeated Measures Models for H2 (Probit) ---

# Run the probit model.
prepost_h2_probit_model <- brm(
  formula = bf(Score ~ 1 + Time + Delivery + (Time + Delivery | EduType/SessionID) + (1 | EduType:SessionID:ParticipantID)),
  data = long_data_h2,
  family = cumulative(link = "probit"),
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(normal(0, 1.5), class = "Intercept"),
    prior(exponential(1), class = "sd")
  ),
  cores = 4, chains = 4, iter = 4000,
  control = list(adapt_delta = 0.99)
)

# --- Save and Load the Probit Model ---
saveRDS(prepost_h2_probit_model, file = "prepost_h2_probit_model.rds")
prepost_h2_probit_model <- readRDS("prepost_h2_probit_model.rds")

## Drawing posteriors from the model
prepost_h2_probit_model_draws <- as_draws_df(prepost_h2_probit_model)

# --- 3b. Extract and Plot Arm-Level Effects for H2 (Probit) ---

post_slope_h2_probit_arm <- prepost_h2_probit_model_draws %>%
  mutate(
    Evo_SD = b_TimePost_Score + `r_EduType[Evo,TimePost_Score]`,
    Gen_SD = b_TimePost_Score + `r_EduType[Gen,TimePost_Score]`
  )

post_slope_h2_probit_summary <- post_slope_h2_probit_arm %>%
  summarise(
    Evo_SD_mean = mean(Evo_SD),
    Gen_SD_mean = mean(Gen_SD),
    Evo_SD_lower = quantile(Evo_SD, probs = 0.025),
    Evo_SD_upper = quantile(Evo_SD, probs = 0.975),
    Gen_SD_lower = quantile(Gen_SD, probs = 0.025),
    Gen_SD_upper = quantile(Gen_SD, probs = 0.975)
  )

post_slope_h2_probit_long <- data.frame(
  Arm = rep(c("Evo", "Gen"), each = 1),
  Mean = c(post_slope_h2_probit_summary$Evo_SD_mean, post_slope_h2_probit_summary$Gen_SD_mean),
  Lower = c(post_slope_h2_probit_summary$Evo_SD_lower, post_slope_h2_probit_summary$Gen_SD_lower),
  Upper = c(post_slope_h2_probit_summary$Evo_SD_upper, post_slope_h2_probit_summary$Gen_SD_upper)
)

ggplot(post_slope_h2_probit_long, aes(x = Arm, y = Mean)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  theme_minimal() +
  labs(x = "Arm", y = "Standard Deviation Change", title = "H2 (Probit): Pre-Post Effect Size") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")


# --- 4b. Extract and Plot Session-Level Effects for H2 (Probit) ---

post_slope_h2_probit_session <- prepost_h2_probit_model_draws %>%
  select(starts_with("r_EduType:SessionID[")) %>%
  select(contains("TimePost_Score"))

fixed_slope_h2_probit_evo <- prepost_h2_probit_model_draws$b_TimePost_Score + prepost_h2_probit_model_draws$`r_EduType[Evo,TimePost_Score]`
fixed_slope_h2_probit_gen <- prepost_h2_probit_model_draws$b_TimePost_Score + prepost_h2_probit_model_draws$`r_EduType[Gen,TimePost_Score]`

post_slope_h2_probit_session <- post_slope_h2_probit_session %>%
  mutate(
    Session_1_SD = fixed_slope_h2_probit_evo + `r_EduType:SessionID[Evo_Evo_1,TimePost_Score]`,
    Session_2_SD = fixed_slope_h2_probit_evo + `r_EduType:SessionID[Evo_Evo_2,TimePost_Score]`,
    Session_3_SD = fixed_slope_h2_probit_evo + `r_EduType:SessionID[Evo_Evo_3,TimePost_Score]`,
    Session_4_SD = fixed_slope_h2_probit_evo + `r_EduType:SessionID[Evo_Evo_4,TimePost_Score]`,
    Session_5_SD = fixed_slope_h2_probit_evo + `r_EduType:SessionID[Evo_Evo_5,TimePost_Score]`,
    Session_6_SD = fixed_slope_h2_probit_evo + `r_EduType:SessionID[Evo_Evo_6,TimePost_Score]`,
    Session_7_SD = fixed_slope_h2_probit_evo + `r_EduType:SessionID[Evo_Evo_7,TimePost_Score]`,
    Session_8_SD = fixed_slope_h2_probit_evo + `r_EduType:SessionID[Evo_Evo_8,TimePost_Score]`,
    Session_9_SD = fixed_slope_h2_probit_evo + `r_EduType:SessionID[Evo_Evo_9,TimePost_Score]`,
    Session_10_SD = fixed_slope_h2_probit_evo + `r_EduType:SessionID[Evo_Evo_10,TimePost_Score]`,
    Gen_Session_1_SD = fixed_slope_h2_probit_gen + `r_EduType:SessionID[Gen_Gen_10,TimePost_Score]`,
    Gen_Session_2_SD = fixed_slope_h2_probit_gen + `r_EduType:SessionID[Gen_Gen_2,TimePost_Score]`,
    Gen_Session_3_SD = fixed_slope_h2_probit_gen + `r_EduType:SessionID[Gen_Gen_3,TimePost_Score]`,
    Gen_Session_4_SD = fixed_slope_h2_probit_gen + `r_EduType:SessionID[Gen_Gen_4,TimePost_Score]`,
    Gen_Session_5_SD = fixed_slope_h2_probit_gen + `r_EduType:SessionID[Gen_Gen_5,TimePost_Score]`,
    Gen_Session_6_SD = fixed_slope_h2_probit_gen + `r_EduType:SessionID[Gen_Gen_6,TimePost_Score]`,
    Gen_Session_7_SD = fixed_slope_h2_probit_gen + `r_EduType:SessionID[Gen_Gen_7,TimePost_Score]`,
    Gen_Session_8_SD = fixed_slope_h2_probit_gen + `r_EduType:SessionID[Gen_Gen_8,TimePost_Score]`,
    Gen_Session_9_SD = fixed_slope_h2_probit_gen + `r_EduType:SessionID[Gen_Gen_9,TimePost_Score]`
  )

session_columns_h2_probit <- names(post_slope_h2_probit_session)[grepl("_SD$", names(post_slope_h2_probit_session))]

session_H2_probit <- map_dfr(session_columns_h2_probit, ~compute_ci(.x, source_df = post_slope_h2_probit_session))

ggplot(session_H2_probit, aes(x = Session, y = mean)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = lower_95_CI, ymax = upper_95_CI), width = 0.2) +
  theme_minimal() +
  labs(x = "Session", y = "Standard Deviation Change", title = "H2 (Probit): Session-Level Effect Sizes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")




##################################### H3 ##############################
# Hypothesis 3: Willingness for patients to seek help
#######################################################################

# --- 1. Prepare Data for H3 ---

# Reshape the data from wide to long format for the repeated measures model.
long_data_h3 <- analysis_data %>%
  # Select the participant identifier, covariates, and the raw scores for H3
  select(
    EduType,
    ParticipantID, 
    SessionID, 
    Delivery,
    Pre_Score = `ExAnSt8_Given_current_public_knowledge_of_anxiety_how_willing_do_you_think_people_are_to_seek_psychiatric_help_for_anxiety`,
    Post_Score = H3_post_ord
  ) %>%
  # Pivot the pre and post scores into a long format
  pivot_longer(
    cols = c(Pre_Score, Post_Score),
    names_to = "Time",
    values_to = "Score"
  ) %>%
  # Convert the Time and Score columns into ordered factors for the model
  mutate(
    Time = factor(Time, levels = c("Pre_Score", "Post_Score"))
  )

# --- 2. Run the Bayesian Repeated Measures Models for H3 ---

# Run the logit model.
prepost_h3_model <- brm(
  formula = bf(Score ~ 1 + Time + Delivery + (Time + Delivery | EduType/SessionID) + (1 | EduType:SessionID:ParticipantID)),
  data = long_data_h3,
  family = cumulative(link = "logit"),
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(normal(0, 1.5), class = "Intercept"),
    prior(exponential(1), class = "sd")
  ),
  cores = 4, chains = 4, iter = 4000, 
  control = list(adapt_delta = 0.99)
)

# --- Save and Load the logit Model ---
saveRDS(prepost_h3_model, file = "prepost_h3_model.rds")
# prepost_h3_model <- readRDS("prepost_h3_model.rds")

## Inspecting the model now
summary(prepost_h3_model)

## Drawing posteriors from the model
prepost_h3_model_draws <- as_draws_df(prepost_h3_model)

# --- 3. Extract and Plot Arm-Level Effects for H3 ---

## Calculate the effect for each arm separately
post_slope_h3_arm <- prepost_h3_model_draws %>%
  mutate(
    Evo = b_TimePost_Score + `r_EduType[Evo,TimePost_Score]`,
    Gen = b_TimePost_Score + `r_EduType[Gen,TimePost_Score]`
  ) %>%
  mutate(
    Evo_OR = exp(Evo),
    Gen_OR = exp(Gen)
  )

## Compute summary stats for plotting
post_slope_h3_summary <- post_slope_h3_arm %>%
  summarise(
    Evo_OR_mean = mean(Evo_OR),
    Gen_OR_mean = mean(Gen_OR),
    Evo_OR_lower = quantile(Evo_OR, probs = 0.025),
    Evo_OR_upper = quantile(Evo_OR, probs = 0.975),
    Gen_OR_lower = quantile(Gen_OR, probs = 0.025),
    Gen_OR_upper = quantile(Gen_OR, probs = 0.975)
  )

## Format data for plotting
post_slope_h3_long <- data.frame(
  Arm = rep(c("Evo", "Gen"), each = 1),
  Mean = c(post_slope_h3_summary$Evo_OR_mean, post_slope_h3_summary$Gen_OR_mean),
  Lower = c(post_slope_h3_summary$Evo_OR_lower, post_slope_h3_summary$Gen_OR_lower),
  Upper = c(post_slope_h3_summary$Evo_OR_upper, post_slope_h3_summary$Gen_OR_upper)
)

# Plot the 95% credible intervals for H3
ggplot(post_slope_h3_long, aes(x = Arm, y = Mean)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  theme_minimal() +
  labs(x = "Arm", y = "Odds Ratio", title = "H3: Pre-Post Change in Willingness to Seek Help") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# --- 4. Extract and Plot Session-Level Effects for H3 ---

## Select session-level random effects
post_slope_h3_session <- prepost_h3_model_draws %>%
  select(starts_with("r_EduType:SessionID[")) %>%
  select(contains("TimePost_Score"))

## Calculate the fixed effect for the evolutionary arm and the genetics arm
fixed_slope_h3_evo <- prepost_h3_model_draws$b_TimePost_Score + prepost_h3_model_draws$`r_EduType[Evo,TimePost_Score]`
fixed_slope_h3_gen <- prepost_h3_model_draws$b_TimePost_Score + prepost_h3_model_draws$`r_EduType[Gen,TimePost_Score]`

## Calculate differing effects across each session, then convert to Odds Ratios
post_slope_h3_session <- post_slope_h3_session %>%
  mutate(
    # Evolutionary Arm Sessions
    Session_1 = fixed_slope_h3_evo + `r_EduType:SessionID[Evo_Evo_1,TimePost_Score]`,
    Session_2 = fixed_slope_h3_evo + `r_EduType:SessionID[Evo_Evo_2,TimePost_Score]`,
    Session_3 = fixed_slope_h3_evo + `r_EduType:SessionID[Evo_Evo_3,TimePost_Score]`,
    Session_4 = fixed_slope_h3_evo + `r_EduType:SessionID[Evo_Evo_4,TimePost_Score]`,
    Session_5 = fixed_slope_h3_evo + `r_EduType:SessionID[Evo_Evo_5,TimePost_Score]`,
    Session_6 = fixed_slope_h3_evo + `r_EduType:SessionID[Evo_Evo_6,TimePost_Score]`,
    Session_7 = fixed_slope_h3_evo + `r_EduType:SessionID[Evo_Evo_7,TimePost_Score]`,
    Session_8 = fixed_slope_h3_evo + `r_EduType:SessionID[Evo_Evo_8,TimePost_Score]`,
    Session_9 = fixed_slope_h3_evo + `r_EduType:SessionID[Evo_Evo_9,TimePost_Score]`,
    Session_10 = fixed_slope_h3_evo + `r_EduType:SessionID[Evo_Evo_10,TimePost_Score]`,
    # Genetics Arm Sessions (following H1 example naming)
    Gen_Session_1 = fixed_slope_h3_gen + `r_EduType:SessionID[Gen_Gen_10,TimePost_Score]`,
    Gen_Session_2 = fixed_slope_h3_gen + `r_EduType:SessionID[Gen_Gen_2,TimePost_Score]`,
    Gen_Session_3 = fixed_slope_h3_gen + `r_EduType:SessionID[Gen_Gen_3,TimePost_Score]`,
    Gen_Session_4 = fixed_slope_h3_gen + `r_EduType:SessionID[Gen_Gen_4,TimePost_Score]`,
    Gen_Session_5 = fixed_slope_h3_gen + `r_EduType:SessionID[Gen_Gen_5,TimePost_Score]`,
    Gen_Session_6 = fixed_slope_h3_gen + `r_EduType:SessionID[Gen_Gen_6,TimePost_Score]`,
    Gen_Session_7 = fixed_slope_h3_gen + `r_EduType:SessionID[Gen_Gen_7,TimePost_Score]`,
    Gen_Session_8 = fixed_slope_h3_gen + `r_EduType:SessionID[Gen_Gen_8,TimePost_Score]`,
    Gen_Session_9 = fixed_slope_h3_gen + `r_EduType:SessionID[Gen_Gen_9,TimePost_Score]`
  ) %>%
  mutate(
    across(c(starts_with("Session_"), starts_with("Gen_Session_")), exp, .names = "{.col}_OR")
  )

## Define column names for which summary statistics are needed
session_columns_h3 <- names(post_slope_h3_session)[grepl("_OR$", names(post_slope_h3_session))]

# Apply the function across all relevant columns for H3
session_H3 <- map_dfr(session_columns_h3, ~compute_ci(.x, source_df = post_slope_h3_session))

# Plot the 95% credible intervals now
ggplot(session_H3, aes(x = Session, y = mean)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = lower_95_CI, ymax = upper_95_CI), width = 0.2) +
  theme_minimal() +
  labs(x = "Session", y = "Odds Ratio", title = "H3: Session-Level Odds Ratios with Error Bars") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red")

# --- 2b. Run the Bayesian Repeated Measures Models for H3 (Probit) ---

# Run the probit model.
prepost_h3_probit_model <- brm(
  formula = bf(Score ~ 1 + Time + Delivery + (Time + Delivery | EduType/SessionID) + (1 | EduType:SessionID:ParticipantID)),
  data = long_data_h3,
  family = cumulative(link = "probit"),
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(normal(0, 1.5), class = "Intercept"),
    prior(exponential(1), class = "sd")
  ),
  cores = 4, chains = 4, iter = 4000,
  control = list(adapt_delta = 0.99)
)

# --- Save and Load the Probit Model ---
saveRDS(prepost_h3_probit_model, file = "prepost_h3_probit_model.rds")
prepost_h3_probit_model <- readRDS("prepost_h3_probit_model.rds")

## Drawing posteriors from the model
prepost_h3_probit_model_draws <- as_draws_df(prepost_h3_probit_model)

# --- 3b. Extract and Plot Arm-Level Effects for H3 (Probit) ---

post_slope_h3_probit_arm <- prepost_h3_probit_model_draws %>%
  mutate(
    Evo_SD = b_TimePost_Score + `r_EduType[Evo,TimePost_Score]`,
    Gen_SD = b_TimePost_Score + `r_EduType[Gen,TimePost_Score]`
  )

post_slope_h3_probit_summary <- post_slope_h3_probit_arm %>%
  summarise(
    Evo_SD_mean = mean(Evo_SD),
    Gen_SD_mean = mean(Gen_SD),
    Evo_SD_lower = quantile(Evo_SD, probs = 0.025),
    Evo_SD_upper = quantile(Evo_SD, probs = 0.975),
    Gen_SD_lower = quantile(Gen_SD, probs = 0.025),
    Gen_SD_upper = quantile(Gen_SD, probs = 0.975)
  )

post_slope_h3_probit_long <- data.frame(
  Arm = rep(c("Evo", "Gen"), each = 1),
  Mean = c(post_slope_h3_probit_summary$Evo_SD_mean, post_slope_h3_probit_summary$Gen_SD_mean),
  Lower = c(post_slope_h3_probit_summary$Evo_SD_lower, post_slope_h3_probit_summary$Gen_SD_lower),
  Upper = c(post_slope_h3_probit_summary$Evo_SD_upper, post_slope_h3_probit_summary$Gen_SD_upper)
)

ggplot(post_slope_h3_probit_long, aes(x = Arm, y = Mean)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  theme_minimal() +
  labs(x = "Arm", y = "Standard Deviation Change", title = "H3 (Probit): Pre-Post Effect Size") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")


# --- 4b. Extract and Plot Session-Level Effects for H3 (Probit) ---

post_slope_h3_probit_session <- prepost_h3_probit_model_draws %>%
  select(starts_with("r_EduType:SessionID[")) %>%
  select(contains("TimePost_Score"))

fixed_slope_h3_probit_evo <- prepost_h3_probit_model_draws$b_TimePost_Score + prepost_h3_probit_model_draws$`r_EduType[Evo,TimePost_Score]`
fixed_slope_h3_probit_gen <- prepost_h3_probit_model_draws$b_TimePost_Score + prepost_h3_probit_model_draws$`r_EduType[Gen,TimePost_Score]`

post_slope_h3_probit_session <- post_slope_h3_probit_session %>%
  mutate(
    Session_1_SD = fixed_slope_h3_probit_evo + `r_EduType:SessionID[Evo_Evo_1,TimePost_Score]`,
    Session_2_SD = fixed_slope_h3_probit_evo + `r_EduType:SessionID[Evo_Evo_2,TimePost_Score]`,
    Session_3_SD = fixed_slope_h3_probit_evo + `r_EduType:SessionID[Evo_Evo_3,TimePost_Score]`,
    Session_4_SD = fixed_slope_h3_probit_evo + `r_EduType:SessionID[Evo_Evo_4,TimePost_Score]`,
    Session_5_SD = fixed_slope_h3_probit_evo + `r_EduType:SessionID[Evo_Evo_5,TimePost_Score]`,
    Session_6_SD = fixed_slope_h3_probit_evo + `r_EduType:SessionID[Evo_Evo_6,TimePost_Score]`,
    Session_7_SD = fixed_slope_h3_probit_evo + `r_EduType:SessionID[Evo_Evo_7,TimePost_Score]`,
    Session_8_SD = fixed_slope_h3_probit_evo + `r_EduType:SessionID[Evo_Evo_8,TimePost_Score]`,
    Session_9_SD = fixed_slope_h3_probit_evo + `r_EduType:SessionID[Evo_Evo_9,TimePost_Score]`,
    Session_10_SD = fixed_slope_h3_probit_evo + `r_EduType:SessionID[Evo_Evo_10,TimePost_Score]`,
    Gen_Session_1_SD = fixed_slope_h3_probit_gen + `r_EduType:SessionID[Gen_Gen_10,TimePost_Score]`,
    Gen_Session_2_SD = fixed_slope_h3_probit_gen + `r_EduType:SessionID[Gen_Gen_2,TimePost_Score]`,
    Gen_Session_3_SD = fixed_slope_h3_probit_gen + `r_EduType:SessionID[Gen_Gen_3,TimePost_Score]`,
    Gen_Session_4_SD = fixed_slope_h3_probit_gen + `r_EduType:SessionID[Gen_Gen_4,TimePost_Score]`,
    Gen_Session_5_SD = fixed_slope_h3_probit_gen + `r_EduType:SessionID[Gen_Gen_5,TimePost_Score]`,
    Gen_Session_6_SD = fixed_slope_h3_probit_gen + `r_EduType:SessionID[Gen_Gen_6,TimePost_Score]`,
    Gen_Session_7_SD = fixed_slope_h3_probit_gen + `r_EduType:SessionID[Gen_Gen_7,TimePost_Score]`,
    Gen_Session_8_SD = fixed_slope_h3_probit_gen + `r_EduType:SessionID[Gen_Gen_8,TimePost_Score]`,
    Gen_Session_9_SD = fixed_slope_h3_probit_gen + `r_EduType:SessionID[Gen_Gen_9,TimePost_Score]`
  )

session_columns_h3_probit <- names(post_slope_h3_probit_session)[grepl("_SD$", names(post_slope_h3_probit_session))]

session_H3_probit <- map_dfr(session_columns_h3_probit, ~compute_ci(.x, source_df = post_slope_h3_probit_session))

ggplot(session_H3_probit, aes(x = Session, y = mean)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = lower_95_CI, ymax = upper_95_CI), width = 0.2) +
  theme_minimal() +
  labs(x = "Session", y = "Standard Deviation Change", title = "H3 (Probit): Session-Level Effect Sizes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")

##################################### H4 ##############################
# Hypothesis 4: Expected efficacy of psychosocial interventions
#######################################################################

# --- 1. Prepare Data for H4 ---

# Reshape the data from wide to long format for the repeated measures model.
long_data_h4 <- analysis_data %>%
  # Select the participant identifier, covariates, and the raw scores for H4
  select(
    EduType,
    ParticipantID, 
    SessionID, 
    Delivery,
    Pre_Score = `ExAnSt4_Given_your_current_understanding_how_effective_do_you_think_psychosocial_interventions_will_be_in_improving_outcomes_in_anxiety_disorder`,
    Post_Score = H4_post_ord
  ) %>%
  # Pivot the pre and post scores into a long format
  pivot_longer(
    cols = c(Pre_Score, Post_Score),
    names_to = "Time",
    values_to = "Score"
  ) %>%
  # Convert the Time and Score columns into ordered factors for the model
  mutate(
    Time = factor(Time, levels = c("Pre_Score", "Post_Score"))
  )

# --- 2. Run the Bayesian Repeated Measures Models for H4 ---

# Run the logit model.
prepost_h4_model <- brm(
  formula = bf(Score ~ 1 + Time + Delivery + (Time + Delivery | EduType/SessionID) + (1 | EduType:SessionID:ParticipantID)),
  data = long_data_h4,
  family = cumulative(link = "logit"),
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(normal(0, 1.5), class = "Intercept"),
    prior(exponential(1), class = "sd")
  ),
  cores = 4, chains = 4, iter = 4000, 
  control = list(adapt_delta = 0.99)
)

# --- Save and Load the logit Model ---
saveRDS(prepost_h4_model, file = "prepost_h4_model.rds")
prepost_h4_model <- readRDS("prepost_h4_model.rds")

## Inspecting the model now
summary(prepost_h4_model)

## Drawing posteriors from the model
prepost_h4_model_draws <- as_draws_df(prepost_h4_model)

# --- 3. Extract and Plot Arm-Level Effects for H4 ---

## Calculate the effect for each arm separately
post_slope_h4_arm <- prepost_h4_model_draws %>%
  mutate(
    Evo = b_TimePost_Score + `r_EduType[Evo,TimePost_Score]`,
    Gen = b_TimePost_Score + `r_EduType[Gen,TimePost_Score]`
  ) %>%
  mutate(
    Evo_OR = exp(Evo),
    Gen_OR = exp(Gen)
  )

## Compute summary stats for plotting
post_slope_h4_summary <- post_slope_h4_arm %>%
  summarise(
    Evo_OR_mean = mean(Evo_OR),
    Gen_OR_mean = mean(Gen_OR),
    Evo_OR_lower = quantile(Evo_OR, probs = 0.025),
    Evo_OR_upper = quantile(Evo_OR, probs = 0.975),
    Gen_OR_lower = quantile(Gen_OR, probs = 0.025),
    Gen_OR_upper = quantile(Gen_OR, probs = 0.975)
  )

## Format data for plotting
post_slope_h4_long <- data.frame(
  Arm = rep(c("Evo", "Gen"), each = 1),
  Mean = c(post_slope_h4_summary$Evo_OR_mean, post_slope_h4_summary$Gen_OR_mean),
  Lower = c(post_slope_h4_summary$Evo_OR_lower, post_slope_h4_summary$Gen_OR_lower),
  Upper = c(post_slope_h4_summary$Evo_OR_upper, post_slope_h4_summary$Gen_OR_upper)
)

# Plot the 95% credible intervals for H4
ggplot(post_slope_h4_long, aes(x = Arm, y = Mean)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  theme_minimal() +
  labs(x = "Arm", y = "Odds Ratio", title = "H4: Pre-Post Change in Efficacy Belief") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- 4. Extract and Plot Session-Level Effects for H4 ---

## Select session-level random effects
post_slope_h4_session <- prepost_h4_model_draws %>%
  select(starts_with("r_EduType:SessionID[")) %>%
  select(contains("TimePost_Score"))

## Calculate the fixed effect for the evolutionary arm and the genetics arm
fixed_slope_h4_evo <- prepost_h4_model_draws$b_TimePost_Score + prepost_h4_model_draws$`r_EduType[Evo,TimePost_Score]`
fixed_slope_h4_gen <- prepost_h4_model_draws$b_TimePost_Score + prepost_h4_model_draws$`r_EduType[Gen,TimePost_Score]`

## Calculate differing effects across each session, then convert to Odds Ratios
post_slope_h4_session <- post_slope_h4_session %>%
  mutate(
    # Evolutionary Arm Sessions
    Session_1 = fixed_slope_h4_evo + `r_EduType:SessionID[Evo_Evo_1,TimePost_Score]`,
    Session_2 = fixed_slope_h4_evo + `r_EduType:SessionID[Evo_Evo_2,TimePost_Score]`,
    Session_3 = fixed_slope_h4_evo + `r_EduType:SessionID[Evo_Evo_3,TimePost_Score]`,
    Session_4 = fixed_slope_h4_evo + `r_EduType:SessionID[Evo_Evo_4,TimePost_Score]`,
    Session_5 = fixed_slope_h4_evo + `r_EduType:SessionID[Evo_Evo_5,TimePost_Score]`,
    Session_6 = fixed_slope_h4_evo + `r_EduType:SessionID[Evo_Evo_6,TimePost_Score]`,
    Session_7 = fixed_slope_h4_evo + `r_EduType:SessionID[Evo_Evo_7,TimePost_Score]`,
    Session_8 = fixed_slope_h4_evo + `r_EduType:SessionID[Evo_Evo_8,TimePost_Score]`,
    Session_9 = fixed_slope_h4_evo + `r_EduType:SessionID[Evo_Evo_9,TimePost_Score]`,
    Session_10 = fixed_slope_h4_evo + `r_EduType:SessionID[Evo_Evo_10,TimePost_Score]`,
    # Genetics Arm Sessions (following H1 example naming)
    Gen_Session_1 = fixed_slope_h4_gen + `r_EduType:SessionID[Gen_Gen_10,TimePost_Score]`,
    Gen_Session_2 = fixed_slope_h4_gen + `r_EduType:SessionID[Gen_Gen_2,TimePost_Score]`,
    Gen_Session_3 = fixed_slope_h4_gen + `r_EduType:SessionID[Gen_Gen_3,TimePost_Score]`,
    Gen_Session_4 = fixed_slope_h4_gen + `r_EduType:SessionID[Gen_Gen_4,TimePost_Score]`,
    Gen_Session_5 = fixed_slope_h4_gen + `r_EduType:SessionID[Gen_Gen_5,TimePost_Score]`,
    Gen_Session_6 = fixed_slope_h4_gen + `r_EduType:SessionID[Gen_Gen_6,TimePost_Score]`,
    Gen_Session_7 = fixed_slope_h4_gen + `r_EduType:SessionID[Gen_Gen_7,TimePost_Score]`,
    Gen_Session_8 = fixed_slope_h4_gen + `r_EduType:SessionID[Gen_Gen_8,TimePost_Score]`,
    Gen_Session_9 = fixed_slope_h4_gen + `r_EduType:SessionID[Gen_Gen_9,TimePost_Score]`
  ) %>%
  mutate(
    across(c(starts_with("Session_"), starts_with("Gen_Session_")), exp, .names = "{.col}_OR")
  )

## Define column names for which summary statistics are needed
session_columns_h4 <- names(post_slope_h4_session)[grepl("_OR$", names(post_slope_h4_session))]

# Apply the function across all relevant columns for H4
session_H4 <- map_dfr(session_columns_h4, ~compute_ci(.x, source_df = post_slope_h4_session))

# Plot the 95% credible intervals now
ggplot(session_H4, aes(x = Session, y = mean)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = lower_95_CI, ymax = upper_95_CI), width = 0.2) +
  theme_minimal() +
  labs(x = "Session", y = "Odds Ratio", title = "H4: Session-Level Odds Ratios with Error Bars") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red")


# --- 2b. Run the Bayesian Repeated Measures Models for H4 (Probit) ---

# Run the probit model.
prepost_h4_probit_model <- brm(
  formula = bf(Score ~ 1 + Time + Delivery + (Time + Delivery | EduType/SessionID) + (1 | EduType:SessionID:ParticipantID)),
  data = long_data_h4,
  family = cumulative(link = "probit"),
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(normal(0, 1.5), class = "Intercept"),
    prior(exponential(1), class = "sd")
  ),
  cores = 4, chains = 4, iter = 4000,
  control = list(adapt_delta = 0.99)
)

# --- Save and Load the Probit Model ---
saveRDS(prepost_h4_probit_model, file = "prepost_h4_probit_model.rds")
prepost_h4_probit_model <- readRDS("prepost_h4_probit_model.rds")

## Drawing posteriors from the model
prepost_h4_probit_model_draws <- as_draws_df(prepost_h4_probit_model)

# --- 3b. Extract and Plot Arm-Level Effects for H4 (Probit) ---

post_slope_h4_probit_arm <- prepost_h4_probit_model_draws %>%
  mutate(
    Evo_SD = b_TimePost_Score + `r_EduType[Evo,TimePost_Score]`,
    Gen_SD = b_TimePost_Score + `r_EduType[Gen,TimePost_Score]`
  )

post_slope_h4_probit_summary <- post_slope_h4_probit_arm %>%
  summarise(
    Evo_SD_mean = mean(Evo_SD),
    Gen_SD_mean = mean(Gen_SD),
    Evo_SD_lower = quantile(Evo_SD, probs = 0.025),
    Evo_SD_upper = quantile(Evo_SD, probs = 0.975),
    Gen_SD_lower = quantile(Gen_SD, probs = 0.025),
    Gen_SD_upper = quantile(Gen_SD, probs = 0.975)
  )

post_slope_h4_probit_long <- data.frame(
  Arm = rep(c("Evo", "Gen"), each = 1),
  Mean = c(post_slope_h4_probit_summary$Evo_SD_mean, post_slope_h4_probit_summary$Gen_SD_mean),
  Lower = c(post_slope_h4_probit_summary$Evo_SD_lower, post_slope_h4_probit_summary$Gen_SD_lower),
  Upper = c(post_slope_h4_probit_summary$Evo_SD_upper, post_slope_h4_probit_summary$Gen_SD_upper)
)

ggplot(post_slope_h4_probit_long, aes(x = Arm, y = Mean)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  theme_minimal() +
  labs(x = "Arm", y = "Standard Deviation Change", title = "H4 (Probit): Pre-Post Effect Size") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")


# --- 4b. Extract and Plot Session-Level Effects for H4 (Probit) ---

post_slope_h4_probit_session <- prepost_h4_probit_model_draws %>%
  select(starts_with("r_EduType:SessionID[")) %>%
  select(contains("TimePost_Score"))

fixed_slope_h4_probit_evo <- prepost_h4_probit_model_draws$b_TimePost_Score + prepost_h4_probit_model_draws$`r_EduType[Evo,TimePost_Score]`
fixed_slope_h4_probit_gen <- prepost_h4_probit_model_draws$b_TimePost_Score + prepost_h4_probit_model_draws$`r_EduType[Gen,TimePost_Score]`

post_slope_h4_probit_session <- post_slope_h4_probit_session %>%
  mutate(
    Session_1_SD = fixed_slope_h4_probit_evo + `r_EduType:SessionID[Evo_Evo_1,TimePost_Score]`,
    Session_2_SD = fixed_slope_h4_probit_evo + `r_EduType:SessionID[Evo_Evo_2,TimePost_Score]`,
    Session_3_SD = fixed_slope_h4_probit_evo + `r_EduType:SessionID[Evo_Evo_3,TimePost_Score]`,
    Session_4_SD = fixed_slope_h4_probit_evo + `r_EduType:SessionID[Evo_Evo_4,TimePost_Score]`,
    Session_5_SD = fixed_slope_h4_probit_evo + `r_EduType:SessionID[Evo_Evo_5,TimePost_Score]`,
    Session_6_SD = fixed_slope_h4_probit_evo + `r_EduType:SessionID[Evo_Evo_6,TimePost_Score]`,
    Session_7_SD = fixed_slope_h4_probit_evo + `r_EduType:SessionID[Evo_Evo_7,TimePost_Score]`,
    Session_8_SD = fixed_slope_h4_probit_evo + `r_EduType:SessionID[Evo_Evo_8,TimePost_Score]`,
    Session_9_SD = fixed_slope_h4_probit_evo + `r_EduType:SessionID[Evo_Evo_9,TimePost_Score]`,
    Session_10_SD = fixed_slope_h4_probit_evo + `r_EduType:SessionID[Evo_Evo_10,TimePost_Score]`,
    Gen_Session_1_SD = fixed_slope_h4_probit_gen + `r_EduType:SessionID[Gen_Gen_10,TimePost_Score]`,
    Gen_Session_2_SD = fixed_slope_h4_probit_gen + `r_EduType:SessionID[Gen_Gen_2,TimePost_Score]`,
    Gen_Session_3_SD = fixed_slope_h4_probit_gen + `r_EduType:SessionID[Gen_Gen_3,TimePost_Score]`,
    Gen_Session_4_SD = fixed_slope_h4_probit_gen + `r_EduType:SessionID[Gen_Gen_4,TimePost_Score]`,
    Gen_Session_5_SD = fixed_slope_h4_probit_gen + `r_EduType:SessionID[Gen_Gen_5,TimePost_Score]`,
    Gen_Session_6_SD = fixed_slope_h4_probit_gen + `r_EduType:SessionID[Gen_Gen_6,TimePost_Score]`,
    Gen_Session_7_SD = fixed_slope_h4_probit_gen + `r_EduType:SessionID[Gen_Gen_7,TimePost_Score]`,
    Gen_Session_8_SD = fixed_slope_h4_probit_gen + `r_EduType:SessionID[Gen_Gen_8,TimePost_Score]`,
    Gen_Session_9_SD = fixed_slope_h4_probit_gen + `r_EduType:SessionID[Gen_Gen_9,TimePost_Score]`
  )

session_columns_h4_probit <- names(post_slope_h4_probit_session)[grepl("_SD$", names(post_slope_h4_probit_session))]

session_H4_probit <- map_dfr(session_columns_h4_probit, ~compute_ci(.x, source_df = post_slope_h4_probit_session))

ggplot(session_H4_probit, aes(x = Session, y = mean)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = lower_95_CI, ymax = upper_95_CI), width = 0.2) +
  theme_minimal() +
  labs(x = "Session", y = "Standard Deviation Change", title = "H4 (Probit): Session-Level Effect Sizes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")



##############################################################################
# --- CALCULATE AND PRINT ALL SUMMARY STATISTICS ---
##############################################################################

# Helper function to create a formatted summary string
create_summary_string <- function(mean_val, lower_ci, upper_ci, prob_pos, digits = 2) {
  sprintf("Mean = %.2f, 95%% CI [%.2f, %.2f], P(change > 0) = %.2f%%", 
          mean_val, lower_ci, upper_ci, prob_pos)
}

# --- H1 Summaries ---
h1_logit_summary <- post_slope_h1_arm %>%
  summarise(
    OR_mean = mean(Evo_OR),
    OR_lower = quantile(Evo_OR, 0.025),
    OR_upper = quantile(Evo_OR, 0.975),
    Prob_pos = mean(Evo > 0) * 100
  )

h1_probit_summary <- post_slope_h1_probit_arm %>%
  summarise(
    SD_mean = mean(Evo_SD),
    SD_lower = quantile(Evo_SD, 0.025),
    SD_upper = quantile(Evo_SD, 0.975),
    Prob_pos = mean(Evo_SD > 0) * 100
  )

# --- H2 Summaries ---
h2_logit_summary <- post_slope_h2_arm %>%
  summarise(
    OR_mean = mean(Evo_OR),
    OR_lower = quantile(Evo_OR, 0.025),
    OR_upper = quantile(Evo_OR, 0.975),
    Prob_pos = mean(Evo > 0) * 100
  )

h2_probit_summary <- post_slope_h2_probit_arm %>%
  summarise(
    SD_mean = mean(Evo_SD),
    SD_lower = quantile(Evo_SD, 0.025),
    SD_upper = quantile(Evo_SD, 0.975),
    Prob_pos = mean(Evo_SD > 0) * 100
  )

# --- H3 Summaries ---
h3_logit_summary <- post_slope_h3_arm %>%
  summarise(
    OR_mean = mean(Evo_OR),
    OR_lower = quantile(Evo_OR, 0.025),
    OR_upper = quantile(Evo_OR, 0.975),
    Prob_pos = mean(Evo > 0) * 100
  )

h3_probit_summary <- post_slope_h3_probit_arm %>%
  summarise(
    SD_mean = mean(Evo_SD),
    SD_lower = quantile(Evo_SD, 0.025),
    SD_upper = quantile(Evo_SD, 0.975),
    Prob_pos = mean(Evo_SD > 0) * 100
  )

# --- H4 Summaries ---
h4_logit_summary <- post_slope_h4_arm %>%
  summarise(
    OR_mean = mean(Evo_OR),
    OR_lower = quantile(Evo_OR, 0.025),
    OR_upper = quantile(Evo_OR, 0.975),
    Prob_pos = mean(Evo > 0) * 100
  )

h4_probit_summary <- post_slope_h4_probit_arm %>%
  summarise(
    SD_mean = mean(Evo_SD),
    SD_lower = quantile(Evo_SD, 0.025),
    SD_upper = quantile(Evo_SD, 0.975),
    Prob_pos = mean(Evo_SD > 0) * 100
  )

# --- Print all results neatly ---
cat("--- Summary of Results for Evolutionary (Evo) Arm ---\n\n")

cat("H1 (Optimism about recovery):\n")
cat("  Logit (OR): ", create_summary_string(h1_logit_summary$OR_mean, h1_logit_summary$OR_lower, h1_logit_summary$OR_upper, h1_logit_summary$Prob_pos), "\n")
cat("  Probit (SD):", create_summary_string(h1_probit_summary$SD_mean, h1_probit_summary$SD_lower, h1_probit_summary$SD_upper, h1_probit_summary$Prob_pos), "\n\n")

cat("H2 (Willingness to share diagnosis):\n")
cat("  Logit (OR): ", create_summary_string(h2_logit_summary$OR_mean, h2_logit_summary$OR_lower, h2_logit_summary$OR_upper, h2_logit_summary$Prob_pos), "\n")
cat("  Probit (SD):", create_summary_string(h2_probit_summary$SD_mean, h2_probit_summary$SD_lower, h2_probit_summary$SD_upper, h2_probit_summary$Prob_pos), "\n\n")

cat("H3 (Willingness to seek help):\n")
cat("  Logit (OR): ", create_summary_string(h3_logit_summary$OR_mean, h3_logit_summary$OR_lower, h3_logit_summary$OR_upper, h3_logit_summary$Prob_pos), "\n")
cat("  Probit (SD):", create_summary_string(h3_probit_summary$SD_mean, h3_probit_summary$SD_lower, h3_probit_summary$SD_upper, h3_probit_summary$Prob_pos), "\n\n")

cat("H4 (Efficacy of psychosocial interventions):\n")
cat("  Logit (OR): ", create_summary_string(h4_logit_summary$OR_mean, h4_logit_summary$OR_lower, h4_logit_summary$OR_upper, h4_logit_summary$Prob_pos), "\n")
cat("  Probit (SD):", create_summary_string(h4_probit_summary$SD_mean, h4_probit_summary$SD_lower, h4_probit_summary$SD_upper, h4_probit_summary$Prob_pos), "\n")



########################################################################################
# STEP 7: PREPOST EVO PROBIT MODEL POSTERIOR PLOTS (EFFECT SIZE SCALE) & COMBINED FIGURE
########################################################################################
#
# This section creates posterior distribution plots for the pre-post probit models
# for H1, H2, H3, and H4, specifically for the Evolutionary arm. The x-axis
# represents the change from pre- to post-workshop in standard deviation units
# (i.e., the effect size).
#
# Finally, it uses the 'patchwork' package to combine the four plots into a
# single 2x2 figure for a comprehensive overview.
#

# --- H1: Posterior Distribution Plot for Effect Size ---

# Step 1: Extract posterior draws for H1 and calculate the Evo arm specific effect
draws_h1 <- as_draws_df(prepost_h1_probit_model)
effect_size_draws_h1 <- draws_h1$b_TimePost_Score + draws_h1$`r_EduType[Evo,TimePost_Score]`

# Step 2: Calculate density data for H1
density_data_h1 <- density(effect_size_draws_h1)
density_df_h1 <- data.frame(x = density_data_h1$x, y = density_data_h1$y)

# Step 3: Calculate probability of a positive effect for H1
prob_positive_h1 <- mean(effect_size_draws_h1 > 0) * 100

# Step 4: Build and SAVE the plot for H1 to a variable
plot_h1_probit <- ggplot(density_df_h1, aes(x = x, y = y)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = subset(density_df_h1, x > 0), 
              aes(ymin = 0, ymax = y), 
              fill = "#E69F00", alpha = 0.8) +
  geom_ribbon(data = subset(density_df_h1, x < 0), 
              aes(ymin = 0, ymax = y), 
              fill = "#E69F00", alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  labs(
    title = "Optimism about recovery (H1)",
    subtitle = paste0("Probability of positive change: ", round(prob_positive_h1, 2), "%"),
    x = "Standard Deviation Change",
    y = "Density"
  ) +
  theme_minimal(base_size = 15)


# --- H2: Posterior Distribution Plot for Effect Size ---

# Step 1: Extract posterior draws for H2 and calculate the Evo arm specific effect
draws_h2 <- as_draws_df(prepost_h2_probit_model)
effect_size_draws_h2 <- draws_h2$b_TimePost_Score + draws_h2$`r_EduType[Evo,TimePost_Score]`

# Step 2: Calculate density data for H2
density_data_h2 <- density(effect_size_draws_h2)
density_df_h2 <- data.frame(x = density_data_h2$x, y = density_data_h2$y)

# Step 3: Calculate probability of a positive effect for H2
prob_positive_h2 <- mean(effect_size_draws_h2 > 0) * 100

# Step 4: Build and SAVE the plot for H2 to a variable
plot_h2_probit <- ggplot(density_df_h2, aes(x = x, y = y)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = subset(density_df_h2, x > 0), 
              aes(ymin = 0, ymax = y), 
              fill = "#E69F00", alpha = 0.8) +
  geom_ribbon(data = subset(density_df_h2, x < 0), 
              aes(ymin = 0, ymax = y), 
              fill = "#E69F00", alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  labs(
    title = "Willingness to share diagnosis (H2)",
    subtitle = paste0("Probability of positive change: ", round(prob_positive_h2, 2), "%"),
    x = "Standard Deviation Change",
    y = "Density"
  ) +
  theme_minimal(base_size = 15)


# --- H3: Posterior Distribution Plot for Effect Size ---

# Step 1: Extract posterior draws for H3 and calculate the Evo arm specific effect
draws_h3 <- as_draws_df(prepost_h3_probit_model)
effect_size_draws_h3 <- draws_h3$b_TimePost_Score + draws_h3$`r_EduType[Evo,TimePost_Score]`

# Step 2: Calculate density data for H3
density_data_h3 <- density(effect_size_draws_h3)
density_df_h3 <- data.frame(x = density_data_h3$x, y = density_data_h3$y)

# Step 3: Calculate probability of a positive effect for H3
prob_positive_h3 <- mean(effect_size_draws_h3 > 0) * 100

# Step 4: Build and SAVE the plot for H3 to a variable
plot_h3_probit <- ggplot(density_df_h3, aes(x = x, y = y)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = subset(density_df_h3, x > 0), 
              aes(ymin = 0, ymax = y), 
              fill = "#E69F00", alpha = 0.8) +
  geom_ribbon(data = subset(density_df_h3, x < 0), 
              aes(ymin = 0, ymax = y), 
              fill = "#E69F00", alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  labs(
    title = "Willingness to seek help (H3)",
    subtitle = paste0("Probability of positive change: ", round(prob_positive_h3, 2), "%"),
    x = "Standard Deviation Change",
    y = "Density"
  ) +
  theme_minimal(base_size = 15)


# --- H4: Posterior Distribution Plot for Effect Size ---

# Step 1: Extract posterior draws for H4 and calculate the Evo arm specific effect
draws_h4 <- as_draws_df(prepost_h4_probit_model)
effect_size_draws_h4 <- draws_h4$b_TimePost_Score + draws_h4$`r_EduType[Evo,TimePost_Score]`

# Step 2: Calculate density data for H4
density_data_h4 <- density(effect_size_draws_h4)
density_df_h4 <- data.frame(x = density_data_h4$x, y = density_data_h4$y)

# Step 3: Calculate probability of a positive effect for H4
prob_positive_h4 <- mean(effect_size_draws_h4 > 0) * 100

# Step 4: Build and SAVE the plot for H4 to a variable
plot_h4_probit <- ggplot(density_df_h4, aes(x = x, y = y)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = subset(density_df_h4, x > 0), 
              aes(ymin = 0, ymax = y), 
              fill = "#E69F00", alpha = 0.8) +
  geom_ribbon(data = subset(density_df_h4, x < 0), 
              aes(ymin = 0, ymax = y), 
              fill = "#E69F00", alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  labs(
    title = "Efficacy of psychosocial interventions (H4)",
    subtitle = paste0("Probability of positive change: ", round(prob_positive_h4, 2), "%"),
    x = "Standard Deviation Change",
    y = "Density"
  ) +
  theme_minimal(base_size = 15)


# --- COMBINE ALL FOUR PLOTS ---

# Arrange the four saved plots into a 2x2 grid and add an overall title.
(plot_h1_probit + plot_h2_probit) / (plot_h3_probit + plot_h4_probit) + 
  plot_annotation(
    title = 'Posterior Distributions for Pre-Post Change in the Evolutionary Arm',
    theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))
  )

########################################################################################
# STEP 7: PREPOST EVO PROBIT MODEL POSTERIOR PLOTS (EFFECT SIZE SCALE) & COMBINED FIGURE
########################################################################################
#
# This section creates posterior distribution plots for the pre-post probit models
# for H1, H2, H3, and H4. The x-axis represents the change from pre- to post-
# workshop in standard deviation units (i.e., the effect size).
#
# Finally, it uses the 'patchwork' package to combine the four plots into a
# single 2x2 figure for a comprehensive overview.
#

# --- H1: Posterior Distribution Plot for Effect Size ---

# Step 1: Extract posterior draws for H1
effect_size_draws_h1 <- as_draws_df(evo_prepost_h1_probit_model)$b_TimePost_Score

# Step 2: Calculate density data for H1
density_data_h1 <- density(effect_size_draws_h1)
density_df_h1 <- data.frame(x = density_data_h1$x, y = density_data_h1$y)

# Step 3: Calculate probability for H1
prob_positive_h1 <- mean(effect_size_draws_h1 > 0) * 100

# Step 4: Build and SAVE the plot for H1 to a variable
plot_h1_probit <- ggplot(density_df_h1, aes(x = x, y = y)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = subset(density_df_h1, x > 0), 
              aes(ymin = 0, ymax = y), 
              fill = "#E69F00", alpha = 0.8) +
  geom_ribbon(data = subset(density_df_h1, x < 0), 
              aes(ymin = 0, ymax = y), 
              fill = "#E69F00", alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  labs(
    title = "Optimism about recovery (H1)",
    subtitle = paste0("Probability of positive change: ", round(prob_positive_h1, 2), "%"),
    x = "Standard Deviation Change",
    y = "Density"
  ) +
  theme_minimal(base_size = 15)


# --- H2: Posterior Distribution Plot for Effect Size ---

# Step 1: Extract posterior draws for H2
effect_size_draws_h2 <- as_draws_df(evo_prepost_h2_probit_model)$b_TimePost_Score

# Step 2: Calculate density data for H2
density_data_h2 <- density(effect_size_draws_h2)
density_df_h2 <- data.frame(x = density_data_h2$x, y = density_data_h2$y)

# Step 3: Calculate probability for H2
prob_positive_h2 <- mean(effect_size_draws_h2 > 0) * 100

# Step 4: Build and SAVE the plot for H2 to a variable
plot_h2_probit <- ggplot(density_df_h2, aes(x = x, y = y)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = subset(density_df_h2, x > 0), 
              aes(ymin = 0, ymax = y), 
              fill = "#E69F00", alpha = 0.8) +
  geom_ribbon(data = subset(density_df_h2, x < 0), 
              aes(ymin = 0, ymax = y), 
              fill = "#E69F00", alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  labs(
    title = "Willingness to share diagnosis (H2)",
    subtitle = paste0("Probability of positive change: ", round(prob_positive_h2, 2), "%"),
    x = "Standard Deviation Change",
    y = "Density"
  ) +
  theme_minimal(base_size = 15)


# --- H3: Posterior Distribution Plot for Effect Size ---

# Step 1: Extract posterior draws for H3
effect_size_draws_h3 <- as_draws_df(evo_prepost_h3_probit_model)$b_TimePost_Score

# Step 2: Calculate density data for H3
density_data_h3 <- density(effect_size_draws_h3)
density_df_h3 <- data.frame(x = density_data_h3$x, y = density_data_h3$y)

# Step 3: Calculate probability for H3
prob_positive_h3 <- mean(effect_size_draws_h3 > 0) * 100

# Step 4: Build and SAVE the plot for H3 to a variable
plot_h3_probit <- ggplot(density_df_h3, aes(x = x, y = y)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = subset(density_df_h3, x > 0), 
              aes(ymin = 0, ymax = y), 
              fill = "#E69F00", alpha = 0.8) +
  geom_ribbon(data = subset(density_df_h3, x < 0), 
              aes(ymin = 0, ymax = y), 
              fill = "#E69F00", alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  labs(
    title = "Willingness to seek help (H3)",
    subtitle = paste0("Probability of positive change: ", round(prob_positive_h3, 2), "%"),
    x = "Standard Deviation Change",
    y = "Density"
  ) +
  theme_minimal(base_size = 15)


# --- H4: Posterior Distribution Plot for Effect Size ---

# Step 1: Extract posterior draws for H4
effect_size_draws_h4 <- as_draws_df(evo_prepost_h4_probit_model)$b_TimePost_Score

# Step 2: Calculate density data for H4
density_data_h4 <- density(effect_size_draws_h4)
density_df_h4 <- data.frame(x = density_data_h4$x, y = density_data_h4$y)

# Step 3: Calculate probability for H4
prob_positive_h4 <- mean(effect_size_draws_h4 > 0) * 100

# Step 4: Build and SAVE the plot for H4 to a variable
plot_h4_probit <- ggplot(density_df_h4, aes(x = x, y = y)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = subset(density_df_h4, x > 0), 
              aes(ymin = 0, ymax = y), 
              fill = "#E69F00", alpha = 0.8) +
  geom_ribbon(data = subset(density_df_h4, x < 0), 
              aes(ymin = 0, ymax = y), 
              fill = "#E69F00", alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  labs(
    title = "Efficacy of psychosocial interventions (H4)",
    subtitle = paste0("Probability of positive change: ", round(prob_positive_h4, 2), "%"),
    x = "Standard Deviation Change",
    y = "Density"
  ) +
  theme_minimal(base_size = 15)


# --- COMBINE ALL FOUR PLOTS ---

# Arrange the four saved plots into a 2x2 grid and add an overall title.
(plot_h1_probit + plot_h2_probit) / (plot_h3_probit + plot_h4_probit) + 
  plot_annotation(
    title = 'Posterior Distributions for Pre-Post Change in the Evolutionary Arm',
    theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))
  )
