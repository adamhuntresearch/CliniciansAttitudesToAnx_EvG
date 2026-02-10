# THIS IS THE EXPLORATORY ANALYSIS R Script for the Clinician reactions to psychoeducation of Anxiety Study
# by Hunt, Carpenter et al. Pre-registration available at https://doi.org/10.17605/OSF.IO/7X98Y
#
# This script runs Bayesian cumulative-link ordinal regression models
# using the 'brms' package to test the four remaining ExAnSt questions.
# The models predict the post-treatment score while controlling for
# the standardized pre-treatment score.
# 
# There are two types of models; logit and probit (which are essentially identical
# but show odds ratio or effect size on the latent variable scale, respectively)
#
# The logic is the same as the H1-H4 script

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
# STEP 1: DATA PREPARATION FOR EXPLORATORY ANALYSIS
#-----------------------------------------------------------------------------
# This step prepares the data for the four exploratory hypotheses (ExpH1-ExpH4)
# in the same manner as the primary hypotheses.

exploratory_data <- psy_data %>%
  mutate(
    # --- Standardize (z-score) the pre-score predictors for exploratory hypotheses ---
    ExpH1_pre_z = scale(`ExAnSt3_Given_your_current_understanding_how_effective_do_you_think_a_combination_of_medication_and_psychotherapy_will_be_to_improve_outcomes_in_anxiety_disorder`),
    ExpH2_pre_z = scale(`ExAnSt5_Given_your_current_understanding_how_strong_a_role_do_individuals_have_to_instigate_or_mitigate_their_anxiety`),
    ExpH3_pre_z = scale(`ExAnSt6_Given_your_current_understanding_how_unchangeable_are_an_individuals_anxiety_symptoms`),
    ExpH4_pre_z = scale(`ExAnSt7_Given_your_current_understanding_how_deeply_do_you_feel_you_understand_anxiety_disorders`),
    
    # --- Renaming Post Score Variables for exploratory hypotheses ---
    ExpH1_post_ord = ExAnSt23_Given_the_information_provided_in_this_education_session_how_effective_do_you_think_a_combination_of_medication_and_psychotherapy_will_be_to_improve_outcomes_in_anxiety_disorder,
    ExpH2_post_ord = ExAnSt25_Given_the_information_provided_in_this_education_session_how_strong_a_role_do_individuals_have_to_instigate_or_mitigate_their_anxiety,
    ExpH3_post_ord = ExAnSt26_Given_the_information_provided_in_this_education_session_how_unchangeable_are_an_individuals_anxiety_symptoms,
    ExpH4_post_ord = ExAnSt27_Given_the_information_provided_in_this_education_session_how_deeply_do_you_feel_you_understand_anxiety_disorders,
    
    # --- Convert other predictors to factors ---
    # Set 'Gen' as the reference level to get the effect of 'Evo' in the output
    EduType = factor(EduType, levels = c("Gen", "Evo")),
    Delivery = as.factor(Delivery),
    SessionID = as.factor(SessionID)
  )


# --- LOAD SAVED EXPLORATORY MODELS ---
# This script loads the pre-run brms models for the exploratory hypotheses below (ExpH1-ExpH4).

# --- Exploratory Hypothesis 1 Models ---
exph1_logit_model_WRprior <- readRDS("exph1_logit_model_WRprior.rds")
exph1_probit_model_WRprior <- readRDS("exph1_probit_model_WRprior.rds")

# --- Exploratory Hypothesis 2 Models ---
exph2_logit_model_WRprior <- readRDS("exph2_logit_model_WRprior.rds")
exph2_probit_model_WRprior <- readRDS("exph2_probit_model_WRprior.rds")

# --- Exploratory Hypothesis 3 Models ---
exph3_logit_model_WRprior <- readRDS("exph3_logit_model_WRprior.rds")
exph3_probit_model_WRprior <- readRDS("exph3_probit_model_WRprior.rds")

# --- Exploratory Hypothesis 4 Models ---
exph4_logit_model_WRprior <- readRDS("exph4_logit_model_WRprior.rds")
exph4_probit_model_WRprior <- readRDS("exph4_probit_model_WRprior.rds")

#-----------------------------------------------------------------------------
# STEP 3: RUN MODELS FOR EXPLORATORY HYPOTHESES
#-----------------------------------------------------------------------------
# For each exploratory hypothesis, we will run two models as specified:
# 1. Logit Model with Weakly Regularizing Priors
# 2. Probit Model with Weakly Regularizing Priors

# Define the set of regularizing priors
WR_priors <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 1.5), class = "Intercept")
)

# --- EXPLORATORY HYPOTHESIS 1: Effectiveness of combined medication & psychotherapy ---
########################################################################################

### ExpH1 logit model with prior ###
##########################################
exph1_logit_model_WRprior <- brm(
  formula = bf(ExpH1_post_ord ~ 1 + EduType + ExpH1_pre_z + Delivery + (1 | SessionID)),
  data = exploratory_data, family = cumulative(link = "logit"),
  prior = WR_priors,
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

# Save the model
saveRDS(exph1_logit_model_WRprior, file = "exph1_logit_model_WRprior.rds")

# Print model summary
summary(exph1_logit_model_WRprior)

# Check the posterior probability
hypothesis(exph1_logit_model_WRprior, "EduTypeEvo > 0")

# Convert log-odds to odds ratio
h1_logit_summary <- summary(exph1_logit_model_WRprior)
log_odds_h1 <- h1_logit_summary$fixed["EduTypeEvo", "Estimate"]
log_odds_lower_h1 <- h1_logit_summary$fixed["EduTypeEvo", "l-95% CI"]
log_odds_upper_h1 <- h1_logit_summary$fixed["EduTypeEvo", "u-95% CI"]

or_h1 <- exp(log_odds_h1)
or_lower_h1 <- exp(log_odds_lower_h1)
or_upper_h1 <- exp(log_odds_upper_h1)

cat(sprintf(
  "\n--- ExpH1 Odds Ratio (Evolutionary vs. Genetic) ---\nOdds Ratio: %.2f\n95%% CI: [%.2f, %.2f]\n\n",
  or_h1, or_lower_h1, or_upper_h1
))


### ExpH1 probit model with prior ###
##########################################
exph1_probit_model_WRprior <- brm(
  formula = bf(ExpH1_post_ord ~ 1 + EduType + ExpH1_pre_z + Delivery + (1 | SessionID)),
  data = exploratory_data, family = cumulative(link = "probit"),
  prior = WR_priors,
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

# Save the model
saveRDS(exph1_probit_model_WRprior, file = "exph1_probit_model_WRprior.rds")

# Print model summary
summary(exph1_probit_model_WRprior)

# Check the posterior probability
hypothesis(exph1_probit_model_WRprior, "EduTypeEvo > 0")

# --- EXPLORATORY HYPOTHESIS 2: Role of individuals to instigate or mitigate ---
########################################################################################

### ExpH2 logit model with prior ###
##########################################
exph2_logit_model_WRprior <- brm(
  formula = bf(ExpH2_post_ord ~ 1 + EduType + ExpH2_pre_z + Delivery + (1 | SessionID)),
  data = exploratory_data, family = cumulative(link = "logit"),
  prior = WR_priors,
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

# Save the model
saveRDS(exph2_logit_model_WRprior, file = "exph2_logit_model_WRprior.rds")

# Print model summary
summary(exph2_logit_model_WRprior)

# Check the posterior probability
hypothesis(exph2_logit_model_WRprior, "EduTypeEvo > 0")

# Convert log-odds to odds ratio
h2_logit_summary <- summary(exph2_logit_model_WRprior)
log_odds_h2 <- h2_logit_summary$fixed["EduTypeEvo", "Estimate"]
log_odds_lower_h2 <- h2_logit_summary$fixed["EduTypeEvo", "l-95% CI"]
log_odds_upper_h2 <- h2_logit_summary$fixed["EduTypeEvo", "u-95% CI"]

or_h2 <- exp(log_odds_h2)
or_lower_h2 <- exp(log_odds_lower_h2)
or_upper_h2 <- exp(log_odds_upper_h2)

cat(sprintf(
  "\n--- ExpH2 Odds Ratio (Evolutionary vs. Genetic) ---\nOdds Ratio: %.2f\n95%% CI: [%.2f, %.2f]\n\n",
  or_h2, or_lower_h2, or_upper_h2
))


### ExpH2 probit model with prior ###
##########################################
exph2_probit_model_WRprior <- brm(
  formula = bf(ExpH2_post_ord ~ 1 + EduType + ExpH2_pre_z + Delivery + (1 | SessionID)),
  data = exploratory_data, family = cumulative(link = "probit"),
  prior = WR_priors,
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

# Save the model
saveRDS(exph2_probit_model_WRprior, file = "exph2_probit_model_WRprior.rds")

# Print model summary
summary(exph2_probit_model_WRprior)

# Check the posterior probability
hypothesis(exph2_probit_model_WRprior, "EduTypeEvo > 0")


# --- EXPLORATORY HYPOTHESIS 3: Unchangeability of an individual's symptoms ---
########################################################################################

### ExpH3 logit model with prior ###
##########################################
exph3_logit_model_WRprior <- brm(
  formula = bf(ExpH3_post_ord ~ 1 + EduType + ExpH3_pre_z + Delivery + (1 | SessionID)),
  data = exploratory_data, family = cumulative(link = "logit"),
  prior = WR_priors,
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

# Save the model
saveRDS(exph3_logit_model_WRprior, file = "exph3_logit_model_WRprior.rds")

# Print model summary
summary(exph3_logit_model_WRprior)

# Check the posterior probability
hypothesis(exph3_logit_model_WRprior, "EduTypeEvo > 0")

# Convert log-odds to odds ratio
h3_logit_summary <- summary(exph3_logit_model_WRprior)
log_odds_h3 <- h3_logit_summary$fixed["EduTypeEvo", "Estimate"]
log_odds_lower_h3 <- h3_logit_summary$fixed["EduTypeEvo", "l-95% CI"]
log_odds_upper_h3 <- h3_logit_summary$fixed["EduTypeEvo", "u-95% CI"]

or_h3 <- exp(log_odds_h3)
or_lower_h3 <- exp(log_odds_lower_h3)
or_upper_h3 <- exp(log_odds_upper_h3)

cat(sprintf(
  "\n--- ExpH3 Odds Ratio (Evolutionary vs. Genetic) ---\nOdds Ratio: %.2f\n95%% CI: [%.2f, %.2f]\n\n",
  or_h3, or_lower_h3, or_upper_h3
))


### ExpH3 probit model with prior ###
##########################################
exph3_probit_model_WRprior <- brm(
  formula = bf(ExpH3_post_ord ~ 1 + EduType + ExpH3_pre_z + Delivery + (1 | SessionID)),
  data = exploratory_data, family = cumulative(link = "probit"),
  prior = WR_priors,
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

# Save the model
saveRDS(exph3_probit_model_WRprior, file = "exph3_probit_model_WRprior.rds")

# Print model summary
summary(exph3_probit_model_WRprior)

# Check the posterior probability
hypothesis(exph3_probit_model_WRprior, "EduTypeEvo > 0")


# --- EXPLORATORY HYPOTHESIS 4: Depth of understanding of anxiety disorders ---
########################################################################################

### ExpH4 logit model with prior ###
##########################################
exph4_logit_model_WRprior <- brm(
  formula = bf(ExpH4_post_ord ~ 1 + EduType + ExpH4_pre_z + Delivery + (1 | SessionID)),
  data = exploratory_data, family = cumulative(link = "logit"),
  prior = WR_priors,
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

# Save the model
saveRDS(exph4_logit_model_WRprior, file = "exph4_logit_model_WRprior.rds")

# Print model summary
summary(exph4_logit_model_WRprior)

# Check the posterior probability
hypothesis(exph4_logit_model_WRprior, "EduTypeEvo > 0")

# Convert log-odds to odds ratio
h4_logit_summary <- summary(exph4_logit_model_WRprior)
log_odds_h4 <- h4_logit_summary$fixed["EduTypeEvo", "Estimate"]
log_odds_lower_h4 <- h4_logit_summary$fixed["EduTypeEvo", "l-95% CI"]
log_odds_upper_h4 <- h4_logit_summary$fixed["EduTypeEvo", "u-95% CI"]

or_h4 <- exp(log_odds_h4)
or_lower_h4 <- exp(log_odds_lower_h4)
or_upper_h4 <- exp(log_odds_upper_h4)

cat(sprintf(
  "\n--- ExpH4 Odds Ratio (Evolutionary vs. Genetic) ---\nOdds Ratio: %.2f\n95%% CI: [%.2f, %.2f]\n\n",
  or_h4, or_lower_h4, or_upper_h4
))


### ExpH4 probit model with prior ###
##########################################
exph4_probit_model_WRprior <- brm(
  formula = bf(ExpH4_post_ord ~ 1 + EduType + ExpH4_pre_z + Delivery + (1 | SessionID)),
  data = exploratory_data, family = cumulative(link = "probit"),
  prior = WR_priors,
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

# Save the model
saveRDS(exph4_probit_model_WRprior, file = "exph4_probit_model_WRprior.rds")

# Print model summary
summary(exph4_probit_model_WRprior)

# Check the posterior probability
hypothesis(exph4_probit_model_WRprior, "EduTypeEvo > 0")




#-----------------------------------------------------------------------------
# STEP 4: WITHIN-GROUP PRE-POST ANALYSIS FOR EXPLORATORY HYPOTHESES
#-----------------------------------------------------------------------------
# This section runs a repeated measures ordinal regression to specifically
# assess the pre-post change *within* the evolutionary education group and 
# the genetics arm separately for each of the four exploratory questions.

################################## ExpH1 ##############################
#######################################################################

# --- 1. Prepare Data for ExpH1 ---

# Reshape the data from wide to long format for the repeated measures model.
long_data_exph1 <- exploratory_data %>%
  # Select the participant identifier, covariates, and the raw scores for ExpH1
  select(
    EduType,
    ParticipantID, 
    SessionID, 
    Delivery,
    Pre_Score = `ExAnSt3_Given_your_current_understanding_how_effective_do_you_think_a_combination_of_medication_and_psychotherapy_will_be_to_improve_outcomes_in_anxiety_disorder`,
    Post_Score = ExpH1_post_ord
  ) %>%
  # Pivot the pre and post scores into a long format
  pivot_longer(
    cols = c(Pre_Score, Post_Score),
    names_to = "Time",
    values_to = "Score"
  ) %>%
  # Convert the Time column into an ordered factor for the model
  mutate(
    Time = factor(Time, levels = c("Pre_Score", "Post_Score"))
  )

# --- 2. Run the Bayesian Repeated Measures Models for ExpH1 ---

# Run the logit model.
prepost_exph1_model <- brm(
  formula = bf(Score ~ 1 + Time + Delivery + (Time + Delivery | EduType/SessionID) + (1 | EduType:SessionID:ParticipantID)),
  data = long_data_exph1,
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
saveRDS(prepost_exph1_model, file = "prepost_exph1_model.rds")
prepost_exph1_model <- readRDS("prepost_exph1_model.rds")

# Inspecting the model
summary(prepost_exph1_model)

# --- 3. Extract and Analyze Posteriors for ExpH1 ---

# Drawing posteriors from the model
prepost_exph1_model_draws <- as_draws_df(prepost_exph1_model)

# Select all the random deviations of the time effect for each arm
post_slope_exph1_arm <- prepost_exph1_model_draws %>%
  select(starts_with("r_EduType[")) %>%
  select(contains("TimePost_Score"))

# Select the averaged out fixed effect across both arms
fixed_slope_exph1 <- prepost_exph1_model_draws$b_TimePost_Score

# Calculate the effect of moving from pre- to post-score for each arm separately
post_slope_exph1_arm <- post_slope_exph1_arm %>%
  mutate(
    Evo = fixed_slope_exph1 + `r_EduType[Evo,TimePost_Score]`,
    Gen = fixed_slope_exph1 + `r_EduType[Gen,TimePost_Score]`
  )

# Convert each to an Odds Ratio
post_slope_exph1_arm <- post_slope_exph1_arm %>%
  mutate(
    Evo_OR = exp(Evo),
    Gen_OR = exp(Gen)
  )

# Summarize the posterior for the Evo arm
exph1_logit_summary <- post_slope_exph1_arm %>%
  summarise(
    OR_mean = mean(Evo_OR),
    OR_lower = quantile(Evo_OR, 0.025),
    OR_upper = quantile(Evo_OR, 0.975),
    Prob_pos = mean(Evo > 0) * 100
  )

summary(exph1_logit_summary)


################################## ExpH2 ##############################
#######################################################################

# --- 1. Prepare Data for ExpH2 ---

long_data_exph2 <- exploratory_data %>%
  select(
    EduType,
    ParticipantID, 
    SessionID, 
    Delivery,
    Pre_Score = `ExAnSt5_Given_your_current_understanding_how_strong_a_role_do_individuals_have_to_instigate_or_mitigate_their_anxiety`,
    Post_Score = ExpH2_post_ord
  ) %>%
  pivot_longer(
    cols = c(Pre_Score, Post_Score),
    names_to = "Time",
    values_to = "Score"
  ) %>%
  mutate(
    Time = factor(Time, levels = c("Pre_Score", "Post_Score"))
  )

# --- 2. Run the Bayesian Repeated Measures Models for ExpH2 ---

prepost_exph2_model <- brm(
  formula = bf(Score ~ 1 + Time + Delivery + (Time + Delivery | EduType/SessionID) + (1 | EduType:SessionID:ParticipantID)),
  data = long_data_exph2,
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
saveRDS(prepost_exph2_model, file = "prepost_exph2_model.rds")
prepost_exph2_model <- readRDS("prepost_exph2_model.rds")

summary(prepost_exph2_model)

# --- 3. Extract and Analyze Posteriors for ExpH2 ---

prepost_exph2_model_draws <- as_draws_df(prepost_exph2_model)

post_slope_exph2_arm <- prepost_exph2_model_draws %>%
  select(starts_with("r_EduType[")) %>%
  select(contains("TimePost_Score"))

fixed_slope_exph2 <- prepost_exph2_model_draws$b_TimePost_Score

post_slope_exph2_arm <- post_slope_exph2_arm %>%
  mutate(
    Evo = fixed_slope_exph2 + `r_EduType[Evo,TimePost_Score]`,
    Gen = fixed_slope_exph2 + `r_EduType[Gen,TimePost_Score]`
  ) %>%
  mutate(
    Evo_OR = exp(Evo),
    Gen_OR = exp(Gen)
  )

exph2_logit_summary <- post_slope_exph2_arm %>%
  summarise(
    OR_mean = mean(Evo_OR),
    OR_lower = quantile(Evo_OR, 0.025),
    OR_upper = quantile(Evo_OR, 0.975),
    Prob_pos = mean(Evo > 0) * 100
  )

summary(exph2_logit_summary)


################################## ExpH3 ##############################
#######################################################################

# --- 1. Prepare Data for ExpH3 ---

long_data_exph3 <- exploratory_data %>%
  select(
    EduType,
    ParticipantID, 
    SessionID, 
    Delivery,
    Pre_Score = `ExAnSt6_Given_your_current_understanding_how_unchangeable_are_an_individuals_anxiety_symptoms`,
    Post_Score = ExpH3_post_ord
  ) %>%
  pivot_longer(
    cols = c(Pre_Score, Post_Score),
    names_to = "Time",
    values_to = "Score"
  ) %>%
  mutate(
    Time = factor(Time, levels = c("Pre_Score", "Post_Score"))
  )

# --- 2. Run the Bayesian Repeated Measures Models for ExpH3 ---

prepost_exph3_model <- brm(
  formula = bf(Score ~ 1 + Time + Delivery + (Time + Delivery | EduType/SessionID) + (1 | EduType:SessionID:ParticipantID)),
  data = long_data_exph3,
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
saveRDS(prepost_exph3_model, file = "prepost_exph3_model.rds")
prepost_exph3_model <- readRDS("prepost_exph3_model.rds")

summary(prepost_exph3_model)

# --- 3. Extract and Analyze Posteriors for ExpH3 ---

prepost_exph3_model_draws <- as_draws_df(prepost_exph3_model)

post_slope_exph3_arm <- prepost_exph3_model_draws %>%
  select(starts_with("r_EduType[")) %>%
  select(contains("TimePost_Score"))

fixed_slope_exph3 <- prepost_exph3_model_draws$b_TimePost_Score

post_slope_exph3_arm <- post_slope_exph3_arm %>%
  mutate(
    Evo = fixed_slope_exph3 + `r_EduType[Evo,TimePost_Score]`,
    Gen = fixed_slope_exph3 + `r_EduType[Gen,TimePost_Score]`
  ) %>%
  mutate(
    Evo_OR = exp(Evo),
    Gen_OR = exp(Gen)
  )

exph3_logit_summary <- post_slope_exph3_arm %>%
  summarise(
    OR_mean = mean(Evo_OR),
    OR_lower = quantile(Evo_OR, 0.025),
    OR_upper = quantile(Evo_OR, 0.975),
    Prob_pos = mean(Evo > 0) * 100
  )

summary(exph3_logit_summary)


################################## ExpH4 ##############################
#######################################################################

# --- 1. Prepare Data for ExpH4 ---

long_data_exph4 <- exploratory_data %>%
  select(
    EduType,
    ParticipantID, 
    SessionID, 
    Delivery,
    Pre_Score = `ExAnSt7_Given_your_current_understanding_how_deeply_do_you_feel_you_understand_anxiety_disorders`,
    Post_Score = ExpH4_post_ord
  ) %>%
  pivot_longer(
    cols = c(Pre_Score, Post_Score),
    names_to = "Time",
    values_to = "Score"
  ) %>%
  mutate(
    Time = factor(Time, levels = c("Pre_Score", "Post_Score"))
  )

# --- 2. Run the Bayesian Repeated Measures Models for ExpH4 ---

prepost_exph4_model <- brm(
  formula = bf(Score ~ 1 + Time + Delivery + (Time + Delivery | EduType/SessionID) + (1 | EduType:SessionID:ParticipantID)),
  data = long_data_exph4,
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
saveRDS(prepost_exph4_model, file = "prepost_exph4_model.rds")
prepost_exph4_model <- readRDS("prepost_exph4_model.rds")

summary(prepost_exph4_model)

# --- 3. Extract and Analyze Posteriors for ExpH4 ---

prepost_exph4_model_draws <- as_draws_df(prepost_exph4_model)

post_slope_exph4_arm <- prepost_exph4_model_draws %>%
  select(starts_with("r_EduType[")) %>%
  select(contains("TimePost_Score"))

fixed_slope_exph4 <- prepost_exph4_model_draws$b_TimePost_Score

post_slope_exph4_arm <- post_slope_exph4_arm %>%
  mutate(
    Evo = fixed_slope_exph4 + `r_EduType[Evo,TimePost_Score]`,
    Gen = fixed_slope_exph4 + `r_EduType[Gen,TimePost_Score]`
  ) %>%
  mutate(
    Evo_OR = exp(Evo),
    Gen_OR = exp(Gen)
  )

exph4_logit_summary <- post_slope_exph4_arm %>%
  summarise(
    OR_mean = mean(Evo_OR),
    OR_lower = quantile(Evo_OR, 0.025),
    OR_upper = quantile(Evo_OR, 0.975),
    Prob_pos = mean(Evo > 0) * 100
  )

summary(exph4_logit_summary)




#-----------------------------------------------------------------------------
# STEP 5: VISUALIZE EXPLORATORY HYPOTHESES WITH JITTER PLOTS
#-----------------------------------------------------------------------------
# This section creates a 2x2 grid of plots showing the raw pre-post changes
# for each of the four exploratory hypotheses.

# --- GRAPH FOR ExpH1: Effectiveness of Combined Tx ---

# 1. Prepare data for ExpH1
plot_data_exph1_long <- psy_data %>%
  mutate(participant_id = row_number()) %>%
  select(
    participant_id, EduType,
    ExpH1_pre = `ExAnSt3_Given_your_current_understanding_how_effective_do_you_think_a_combination_of_medication_and_psychotherapy_will_be_to_improve_outcomes_in_anxiety_disorder`,
    ExpH1_post = `ExAnSt23_Given_the_information_provided_in_this_education_session_how_effective_do_you_think_a_combination_of_medication_and_psychotherapy_will_be_to_improve_outcomes_in_anxiety_disorder`
  ) %>%
  tidyr::pivot_longer(cols = c(ExpH1_pre, ExpH1_post), names_to = "time", values_to = "score") %>%
  mutate(time = factor(time, levels = c("ExpH1_pre", "ExpH1_post"), labels = c("Pre-Education", "Post-Education")))

# 2. Calculate summary stats for ExpH1
summary_data_exph1 <- plot_data_exph1_long %>%
  group_by(EduType, time) %>%
  summarise(
    n = n(), mean_score = mean(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE), se_score = sd_score / sqrt(n),
    ci_lower = mean_score - 1.96 * se_score, ci_upper = mean_score + 1.96 * se_score,
    .groups = 'drop'
  )

# 3. Build and SAVE the plot for ExpH1
plot_exph1 <- ggplot(summary_data_exph1, aes(x = time, y = mean_score, group = EduType, color = EduType)) +
  geom_jitter(data = plot_data_exph1_long, aes(y = score), width = 0.2, alpha = 0.25, size = 1.8) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.08, linewidth = 0.8) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  labs(
    title = "Effectiveness of Combined Tx (ExpH1)",
    y = "Clinician Rating (1-7)", 
    x = NULL 
  ) +
  scale_color_manual(values = c("Gen" = "#0072B2", "Evo" = "#D55E00"), labels = c("Gen" = "Genetic", "Evo" = "Evolutionary")) +
  coord_cartesian(ylim = c(1, 7)) +
  theme_minimal(base_size = 15)


# --- GRAPH FOR ExpH2: Individual Role in Anxiety ---

# 1. Prepare data for ExpH2
plot_data_exph2_long <- psy_data %>%
  mutate(participant_id = row_number()) %>%
  select(
    participant_id, EduType,
    ExpH2_pre = `ExAnSt5_Given_your_current_understanding_how_strong_a_role_do_individuals_have_to_instigate_or_mitigate_their_anxiety`,
    ExpH2_post = `ExAnSt25_Given_the_information_provided_in_this_education_session_how_strong_a_role_do_individuals_have_to_instigate_or_mitigate_their_anxiety`
  ) %>%
  tidyr::pivot_longer(cols = c(ExpH2_pre, ExpH2_post), names_to = "time", values_to = "score") %>%
  mutate(time = factor(time, levels = c("ExpH2_pre", "ExpH2_post"), labels = c("Pre-Education", "Post-Education")))

# 2. Calculate summary stats for ExpH2
summary_data_exph2 <- plot_data_exph2_long %>%
  group_by(EduType, time) %>%
  summarise(
    n = n(), mean_score = mean(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE), se_score = sd_score / sqrt(n),
    ci_lower = mean_score - 1.96 * se_score, ci_upper = mean_score + 1.96 * se_score,
    .groups = 'drop'
  )

# 3. Build and SAVE the plot for ExpH2
plot_exph2 <- ggplot(summary_data_exph2, aes(x = time, y = mean_score, group = EduType, color = EduType)) +
  geom_jitter(data = plot_data_exph2_long, aes(y = score), width = 0.2, alpha = 0.25, size = 1.8) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.08, linewidth = 0.8) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  labs(
    title = "Individual Role in Anxiety (ExpH2)",
    y = NULL, 
    x = NULL 
  ) +
  scale_color_manual(values = c("Gen" = "#0072B2", "Evo" = "#D55E00"), labels = c("Gen" = "Genetic", "Evo" = "Evolutionary")) +
  coord_cartesian(ylim = c(1, 7)) +
  theme_minimal(base_size = 15)


# --- GRAPH FOR ExpH3: Unchangeability of Symptoms ---

# 1. Prepare data for ExpH3
plot_data_exph3_long <- psy_data %>%
  mutate(participant_id = row_number()) %>%
  select(
    participant_id, EduType,
    ExpH3_pre = `ExAnSt6_Given_your_current_understanding_how_unchangeable_are_an_individuals_anxiety_symptoms`,
    ExpH3_post = `ExAnSt26_Given_the_information_provided_in_this_education_session_how_unchangeable_are_an_individuals_anxiety_symptoms`
  ) %>%
  tidyr::pivot_longer(cols = c(ExpH3_pre, ExpH3_post), names_to = "time", values_to = "score") %>%
  mutate(time = factor(time, levels = c("ExpH3_pre", "ExpH3_post"), labels = c("Pre-Education", "Post-Education")))

# 2. Calculate summary stats for ExpH3
summary_data_exph3 <- plot_data_exph3_long %>%
  group_by(EduType, time) %>%
  summarise(
    n = n(), mean_score = mean(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE), se_score = sd_score / sqrt(n),
    ci_lower = mean_score - 1.96 * se_score, ci_upper = mean_score + 1.96 * se_score,
    .groups = 'drop'
  )

# 3. Build and SAVE the plot for ExpH3
plot_exph3 <- ggplot(summary_data_exph3, aes(x = time, y = mean_score, group = EduType, color = EduType)) +
  geom_jitter(data = plot_data_exph3_long, aes(y = score), width = 0.2, alpha = 0.25, size = 1.8) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.08, linewidth = 0.8) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  labs(
    title = "Unchangeability of Symptoms (ExpH3)",
    y = "Clinician Rating (1-7)", 
    x = "Time Point"
  ) +
  scale_color_manual(values = c("Gen" = "#0072B2", "Evo" = "#D55E00"), labels = c("Gen" = "Genetic", "Evo" = "Evolutionary")) +
  coord_cartesian(ylim = c(1, 7)) +
  theme_minimal(base_size = 15)


# --- GRAPH FOR ExpH4: Depth of Understanding ---

# 1. Prepare data for ExpH4
plot_data_exph4_long <- psy_data %>%
  mutate(participant_id = row_number()) %>%
  select(
    participant_id, EduType,
    ExpH4_pre = `ExAnSt7_Given_your_current_understanding_how_deeply_do_you_feel_you_understand_anxiety_disorders`,
    ExpH4_post = `ExAnSt27_Given_the_information_provided_in_this_education_session_how_deeply_do_you_feel_you_understand_anxiety_disorders`
  ) %>%
  tidyr::pivot_longer(cols = c(ExpH4_pre, ExpH4_post), names_to = "time", values_to = "score") %>%
  mutate(time = factor(time, levels = c("ExpH4_pre", "ExpH4_post"), labels = c("Pre-Education", "Post-Education")))

# 2. Calculate summary stats for ExpH4
summary_data_exph4 <- plot_data_exph4_long %>%
  group_by(EduType, time) %>%
  summarise(
    n = n(), mean_score = mean(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE), se_score = sd_score / sqrt(n),
    ci_lower = mean_score - 1.96 * se_score, ci_upper = mean_score + 1.96 * se_score,
    .groups = 'drop'
  )

# 3. Build and SAVE the plot for ExpH4
plot_exph4 <- ggplot(summary_data_exph4, aes(x = time, y = mean_score, group = EduType, color = EduType)) +
  geom_jitter(data = plot_data_exph4_long, aes(y = score), width = 0.2, alpha = 0.25, size = 1.8) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.08, linewidth = 0.8) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  labs(
    title = "Depth of Understanding (ExpH4)",
    y = NULL, 
    x = "Time Point"
  ) +
  scale_color_manual(values = c("Gen" = "#0072B2", "Evo" = "#D55E00"), labels = c("Gen" = "Genetic", "Evo" = "Evolutionary")) +
  coord_cartesian(ylim = c(1, 7)) +
  theme_minimal(base_size = 15)


# --- COMBINE THE PLOTS WITH PATCHWORK ---

# Arrange the four plots into a 2x2 grid and add an overall title.
(plot_exph1 + plot_exph2) / (plot_exph3 + plot_exph4) + 
  plot_layout(guides = 'collect') & 
  plot_annotation(
    title = 'Exploratory Analyses: Effect of Education on Clinician Attitudes',
    theme = theme(plot.title = element_text(size = 20, hjust = 0.5))
  ) &
  theme(legend.position = 'bottom')




#-----------------------------------------------------------------------------
# STEP 5: VISUALIZE EXPLORATORY POSTERIOR DISTRIBUTIONS
#-----------------------------------------------------------------------------
# This section creates posterior distribution plots for the logit models
# for the four exploratory hypotheses.

# --- ExpH1: Posterior Distribution Plot for Odds Ratio ---

# Step 1: Extract and convert posterior draws
log_odds_draws_exph1 <- as_draws_df(exph1_logit_model_WRprior)$b_EduTypeEvo
odds_ratio_draws_exph1 <- exp(log_odds_draws_exph1)

# Step 2: Calculate density data
density_data_exph1 <- density(odds_ratio_draws_exph1)
density_df_exph1 <- data.frame(x = density_data_exph1$x, y = density_data_exph1$y)

# Step 3: Calculate probability
prob_positive_exph1 <- mean(odds_ratio_draws_exph1 > 1) * 100

# Step 4: Build and SAVE the plot
plot_exph1_posterior <- ggplot(density_df_exph1, aes(x = x, y = y)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = subset(density_df_exph1, x > 1), 
              aes(ymin = 0, ymax = y), 
              fill = "#D55E00", alpha = 0.8) +
  geom_ribbon(data = subset(density_df_exph1, x < 1), 
              aes(ymin = 0, ymax = y), 
              fill = "#D55E00", alpha = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +
  labs(
    title = "Effectiveness of Combined Tx (ExpH1)",
    subtitle = paste0("Probability of positive effect: ", round(prob_positive_exph1, 2), "%"),
    x = "Odds Ratio",
    y = "Density"
  ) +
  theme_minimal(base_size = 15)


# --- ExpH2: Posterior Distribution Plot for Odds Ratio ---

# Step 1: Extract and convert posterior draws
log_odds_draws_exph2 <- as_draws_df(exph2_logit_model_WRprior)$b_EduTypeEvo
odds_ratio_draws_exph2 <- exp(log_odds_draws_exph2)

# Step 2: Calculate density data
density_data_exph2 <- density(odds_ratio_draws_exph2)
density_df_exph2 <- data.frame(x = density_data_exph2$x, y = density_data_exph2$y)

# Step 3: Calculate probability
prob_positive_exph2 <- mean(odds_ratio_draws_exph2 > 1) * 100

# Step 4: Build and SAVE the plot
plot_exph2_posterior <- ggplot(density_df_exph2, aes(x = x, y = y)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = subset(density_df_exph2, x > 1), 
              aes(ymin = 0, ymax = y), 
              fill = "#D55E00", alpha = 0.8) +
  geom_ribbon(data = subset(density_df_exph2, x < 1), 
              aes(ymin = 0, ymax = y), 
              fill = "#D55E00", alpha = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +
  labs(
    title = "Individual Role in Anxiety (ExpH2)",
    subtitle = paste0("Probability of positive effect: ", round(prob_positive_exph2, 2), "%"),
    x = "Odds Ratio",
    y = "Density"
  ) +
  theme_minimal(base_size = 15)


# --- ExpH3: Posterior Distribution Plot for Odds Ratio ---

# Step 1: Extract and convert posterior draws
log_odds_draws_exph3 <- as_draws_df(exph3_logit_model_WRprior)$b_EduTypeEvo
odds_ratio_draws_exph3 <- exp(log_odds_draws_exph3)

# Step 2: Calculate density data
density_data_exph3 <- density(odds_ratio_draws_exph3)
density_df_exph3 <- data.frame(x = density_data_exph3$x, y = density_data_exph3$y)

# Step 3: Calculate probability
prob_positive_exph3 <- mean(odds_ratio_draws_exph3 > 1) * 100

# Step 4: Build and SAVE the plot
plot_exph3_posterior <- ggplot(density_df_exph3, aes(x = x, y = y)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = subset(density_df_exph3, x > 1), 
              aes(ymin = 0, ymax = y), 
              fill = "#D55E00", alpha = 0.8) +
  geom_ribbon(data = subset(density_df_exph3, x < 1), 
              aes(ymin = 0, ymax = y), 
              fill = "#D55E00", alpha = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +
  labs(
    title = "Unchangeability of Symptoms (ExpH3)",
    subtitle = paste0("Probability of positive effect: ", round(prob_positive_exph3, 2), "%"),
    x = "Odds Ratio",
    y = "Density"
  ) +
  theme_minimal(base_size = 15)


# --- ExpH4: Posterior Distribution Plot for Odds Ratio ---

# Step 1: Extract and convert posterior draws
log_odds_draws_exph4 <- as_draws_df(exph4_logit_model_WRprior)$b_EduTypeEvo
odds_ratio_draws_exph4 <- exp(log_odds_draws_exph4)

# Step 2: Calculate density data
density_data_exph4 <- density(odds_ratio_draws_exph4)
density_df_exph4 <- data.frame(x = density_data_exph4$x, y = density_data_exph4$y)

# Step 3: Calculate probability
prob_positive_exph4 <- mean(odds_ratio_draws_exph4 > 1) * 100

# Step 4: Build and SAVE the plot
plot_exph4_posterior <- ggplot(density_df_exph4, aes(x = x, y = y)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = subset(density_df_exph4, x > 1), 
              aes(ymin = 0, ymax = y), 
              fill = "#D55E00", alpha = 0.8) +
  geom_ribbon(data = subset(density_df_exph4, x < 1), 
              aes(ymin = 0, ymax = y), 
              fill = "#D55E00", alpha = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +
  labs(
    title = "Depth of Understanding (ExpH4)",
    subtitle = paste0("Probability of positive effect: ", round(prob_positive_exph4, 2), "%"),
    x = "Odds Ratio",
    y = "Density"
  ) +
  theme_minimal(base_size = 15)


# --- COMBINE ALL FOUR POSTERIOR PLOTS ---

# Arrange the four saved plots into a 2x2 grid and add an overall title.
(plot_exph1_posterior + plot_exph2_posterior) / (plot_exph3_posterior + plot_exph4_posterior) + 
  plot_annotation(
    title = 'Exploratory Analyses: Posterior Distributions for the Effect of Evolutionary vs. Genetic Psychoeducation',
    theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))
  )

