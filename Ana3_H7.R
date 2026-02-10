# THIS IS THE THIRD PRIMARY ANALYSIS R Script for the Clinician reactions to psychoeducation of Anxiety Study 
# by Hunt, Carpenter et al.  Pre-registration available at https://doi.org/10.17605/OSF.IO/7X98Y
#
# This script tests Hypothesis 7 by creating sum-scores from the
# Emotions Beliefs Questionnaire (EBQ) items. It then runs a 
# Bayesian linear mixed-effects model initially (family = gaussian) to
# predict the standardized post-education sum score. It then checks the fit of alternative families.
# The script then proceeds to draw figures of the posterior probability and raw scores.
#
#-----------------------------------------------------------------------------
# STEP 1: SETUP - LOAD LIBRARIES AND DATA
#-----------------------------------------------------------------------------

# Load necessary R packages.
# 'brms' is for Bayesian regression modeling.
# 'tidyverse' is a collection of packages for data manipulation and visualization.

library(brms)
library(tidyverse)
library(ggplot2)
library(dplyr)

# Set a seed for reproducibility.
set.seed(1234)

# Load the dataset.
file_path <- "FINALLY_analysis_ready_data.csv"
psy_data <- read_csv(file_path)

#-----------------------------------------------------------------------------
# STEP 2: DATA PREPARATION - CREATE AND STANDARDIZE SUM SCORES
#-----------------------------------------------------------------------------
# --- REVERSE CODE THE 8 EBQ ITEMS ---
# The original EBQ was coded to have higher scores meaning less belief 
# in usefulness of negative emotions, but in reporting our results we 
# are generally counting higher scores as less stigma, 
# so we need to define the names of the columns we need to reverse and reverse them.

ebq_items_to_reverse <- c(
  "EBQ3_There_is_very_little_use_for_negative_emotions",
  "EBQ7_People_dont_need_their_negative_emotions",
  "EBQ11_Negative_emotions_are_harmful",
  "EBQ15_The_presence_of_negative_emotions_is_a_bad_thing_for_people",
  "EBQ23_There_is_very_little_use_for_negative_emotions",
  "EBQ27_People_dont_need_their_negative_emotions",
  "EBQ211_Negative_emotions_are_harmful",
  "EBQ215_The_presence_of_negative_emotions_is_a_bad_thing_for_people"
)

# Use mutate() with across() to apply the reverse-coding formula to all specified items.
# The formula for a 7-point scale is 8 - score.
psy_data <- psy_data %>%
  mutate(across(all_of(ebq_items_to_reverse), ~ 8 - .x))


# This step then creates the Pre_SumScore and Post_SumScore by adding the relevant
# EBQ items, and then standardizes them (z-scores).

analysis_data_h7 <- psy_data %>%
  # Use rowwise() to perform calculations on each row individually
  rowwise() %>%
  mutate(
    # --- Create the Pre_SumScore ---
    Pre_SumScore = sum(c(
      `EBQ3_There_is_very_little_use_for_negative_emotions`,
      `EBQ7_People_dont_need_their_negative_emotions`,
      `EBQ11_Negative_emotions_are_harmful`,
      `EBQ15_The_presence_of_negative_emotions_is_a_bad_thing_for_people`
    ), na.rm = TRUE), # na.rm = TRUE will sum even if some items are missing
    
    # --- Create the Post_SumScore ---
    Post_SumScore = sum(c(
      `EBQ23_There_is_very_little_use_for_negative_emotions`,
      `EBQ27_People_dont_need_their_negative_emotions`,
      `EBQ211_Negative_emotions_are_harmful`,
      `EBQ215_The_presence_of_negative_emotions_is_a_bad_thing_for_people`
    ), na.rm = TRUE)
  ) %>%
  # ungroup() is important after rowwise operations
  ungroup() %>%
  # Now, standardize the newly created sum scores
  mutate(
    Pre_SumScore_z = as.numeric(scale(Pre_SumScore)),
    Post_SumScore_z = as.numeric(scale(Post_SumScore)),
    
    # --- Convert other predictors to factors ---
    EduType = factor(EduType, levels = c("Gen", "Evo")),
    Delivery = as.factor(Delivery),
    SessionID = as.factor(SessionID)
  )


# Define the set of regularizing priors
WR_priors <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 1.5), class = "Intercept")
)

#-----------------------------------------------------------------------------
# STEP 3: RUN THE GAUSSIAN (CONTINUOUS) MODEL FOR H7
#-----------------------------------------------------------------------------
# This model predicts the standardized post-education sum score. The coefficients
# can be interpreted directly as changes in standard deviation units.

h7_gaussian_model <- brm(
  formula = bf(Post_SumScore_z ~ 1 + EduType + Pre_SumScore_z + Delivery + (1 | SessionID)),
  data = analysis_data_h7,
  family = gaussian(), # Use gaussian() for a continuous outcome
  cores = 4,
  chains = 4,
  iter = 4000
)

saveRDS(h7_gaussian_model, file = "h7_gaussian_model.rds")



h7_gaussian_model_WRprior <- brm(
  formula = bf(Post_SumScore_z ~ 1 + EduType + Pre_SumScore_z + Delivery + (1 | SessionID)),
  data = analysis_data_h7,
  family = gaussian(), # Use gaussian() for a continuous outcome
  prior = WR_priors,
  cores = 4, chains = 4, iter = 4000, control = list(adapt_delta = 0.99)
)

saveRDS(h7_gaussian_model_WRprior, file = "h7_gaussian_model_WRprior.rds")
h7_gaussian_model_WRprior <- readRDS(file = "h7_gaussian_model_WRprior.rds")


#-----------------------------------------------------------------------------
# STEP 4: VIEW MODEL SUMMARY
#-----------------------------------------------------------------------------
# The summary will show the estimated effect of each predictor on the 
# standardized post-education sum score.

# --- Gaussian Model ---
print("--- H7 Gaussian Model Summary ---")
summary(h7_gaussian_model)

# Calculate and print the percentage of the posterior on either side of zero
draws_gauss <- as_draws_df(h7_gaussian_model)
cat("Percentage of posterior for b_EduTypeEvo > 0:",
    round(mean(draws_gauss$b_EduTypeEvo > 0) * 100, 2), "%\n")
cat("Percentage of posterior for b_EduTypeEvo < 0:",
    round(mean(draws_gauss$b_EduTypeEvo < 0) * 100, 2), "%\n\n")


# --- Gaussian Model with Priors ---
print("--- H7 Gaussian Model with Priors Summary ---")
summary(h7_gaussian_model_WRprior)

# ADDED: Calculate and print the percentage of the posterior on either side of zero
draws_gauss_wr <- as_draws_df(h7_gaussian_model_WRprior)
cat("Percentage of posterior for b_EduTypeEvo > 0:",
    round(mean(draws_gauss_wr$b_EduTypeEvo > 0) * 100, 2), "%\n")
cat("Percentage of posterior for b_EduTypeEvo < 0:",
    round(mean(draws_gauss_wr$b_EduTypeEvo < 0) * 100, 2), "%\n\n")


#-----------------------------------------------------------------------------
# STEP 5: POSTERIOR PREDICTIVE CHECK
#-----------------------------------------------------------------------------
# This is a crucial step to check if the gaussian (normal) assumption is
# reasonable for the data. The pp_check() function plots the density of the
# actual outcome (the dark line) against many simulated outcomes from the model
# (the light lines). If the dark line looks similar to the light lines, the
# model is a good fit.

pp_check(h7_gaussian_model_WRprior, ndraws = 100) +
  ggtitle("Posterior Predictive Check for H7 Gaussian Model")


# POSTERIOR CHECK LOOKS NOT QUITE RIGHT - THERE IS A LONG TAIL
# TRYING A LOG NORMAL AND GAMMA MODEL INSTEAD
# BOTH SIMILARLY MISS THE RIGHT SIDE: EVIDENCE OF CEILING EFFECTS

#
#
# Option 1: Log-normal model
h7_lognormal_model <- brm(
  formula = bf(Post_SumScore ~ 1 + EduType + Pre_SumScore_z + Delivery + (1 | SessionID)),
  data = analysis_data_h7,
  family = lognormal(), 
  cores = 4,
  chains = 4,
  iter = 4000
)

summary(h7_lognormal_model)

pp_check(h7_lognormal_model, ndraws = 100) +
  ggtitle("Posterior Predictive Check for H7 Log-Normal Model")

#
#
# Option 2: Gamma model
h7_gamma_model <- brm(
  formula = bf(Post_SumScore ~ 1 + EduType + Pre_SumScore_z + Delivery + (1 | SessionID)),
  data = analysis_data_h7,
  family = Gamma(link = "log"), 
  cores = 4,
  chains = 4,
  iter = 4000
)

summary(h7_gamma_model)

pp_check(h7_gamma_model, ndraws = 100) +
  ggtitle("Posterior Predictive Check for H7 Gamma Model")





################################################################################
# STEP 5: POSTERIOR DISTRIBUTION PLOT FOR GAUSSIAN H7
################################################################################
#
# This script creates a posterior distribution plot for the main Gaussian model
# for H7. The x-axis represents the coefficient for the effect of EduType, which
# can be interpreted as the difference in standard deviation units between the
# Evolutionary and Genetic groups.
#
# A value of 0 indicates no effect.
#

# --- H7: Posterior Distribution Plot for SD Difference ---

# Step 1: Extract the posterior draws for the coefficient.
# These are already on the standard deviation scale, so no transformation is needed.
sd_diff_draws_h7 <- as_draws_df(h7_gaussian_model_WRprior)$b_EduTypeEvo

# Step 2: Manually calculate the density data.
density_data_h7 <- density(sd_diff_draws_h7)
density_df_h7 <- data.frame(x = density_data_h7$x, y = density_data_h7$y)

# Step 3: Calculate the probability of the effect being positive (SD difference > 0).
# A positive effect supports the hypothesis, as the scale has been reversed

prob_positive_h7 <- mean(sd_diff_draws_h7 > 0) * 100

# Step 4: Build the plot using the pre-calculated density data.
ggplot(density_df_h7, aes(x = x, y = y)) +
  
  # Layer 1: Draw the outline of the density curve.
  geom_line(linewidth = 1) +
  
  # Layer 2: Darkly shade the area where the SD difference is NEGATIVE (the hypothesized effect).
  geom_ribbon(data = subset(density_df_h7, x < 0), 
              aes(ymin = 0, ymax = y), 
              fill = "#D55E00", alpha = 0.3) +
  
  # Layer 3: Lightly shade the area where the SD difference is POSITIVE.
  geom_ribbon(data = subset(density_df_h7, x > 0), 
              aes(ymin = 0, ymax = y), 
              fill = "#D55E00", alpha = 0.8) +
  
  # Layer 4: Add a vertical dashed line at x = 0.
  # A difference of 0 signifies no effect.
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  
  # Layer 5: Add informative labels and titles.
  labs(
    title = "Belief in Usefulness of Negative Emotions after Evo vs. Gen Education (H7)",
    subtitle = paste0("Probability of positive effect (SD difference > 0): ", round(prob_positive_h7, 2), "%"),
    x = "Standard Deviation Difference",
    y = "Density"
  ) +
  
  # Layer 6: Use a clean theme for a professional look.
  theme_minimal(base_size = 15)


#-----------------------------------------------------------------------------
# STEP 6: RE-RUNNING THE ANALYSIS WITH A LOG-NORMAL DISTRIBUTION
#-----------------------------------------------------------------------------
#
# The posterior predictive checks in Step 5 indicated that the Gaussian model
# was not a great fit for the data, which has a right-skewed distribution.
# The log-normal distribution was a much better fit. This section re-runs
# the primary analysis for H7 using this more appropriate model.
#
# We will model the raw Post_SumScore, not the standardized version, as
# log-normal models require strictly positive outcomes. The interpretation of
# the coefficients will change: they are on the log scale, so we exponentiate
# them (e.g., exp(b_EduTypeEvo)) to interpret them as multiplicative changes
# (i.e., a ratio) in the median Post_SumScore.
#

# --- 6.1: FIT THE LOG-NORMAL MODEL WITH REGULARIZING PRIORS ---

# We use the same weakly regularizing priors as before for consistency.
h7_lognormal_model_WRprior <- brm(
  formula = bf(Post_SumScore ~ 1 + EduType + Pre_SumScore_z + Delivery + (1 | SessionID)),
  data = analysis_data_h7,
  family = lognormal(),
  prior = WR_priors,
  cores = 4,
  chains = 4,
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

# Save the final model object for future reference
saveRDS(h7_lognormal_model_WRprior, file = "h7_lognormal_model_WRprior.rds")


# --- 6.2: VIEW MODEL SUMMARY AND INTERPRET RESULTS ---

print("--- H7 Log-Normal Model with Priors Summary ---")
summary(h7_lognormal_model_WRprior)

# The coefficient for EduTypeEvo represents the log-ratio of the median score
# for the Evolutionary group compared to the Genetic group.
# A positive coefficient suggests the Evo group has a higher median score.

# Calculate and print the percentage of the posterior on either side of zero
draws_lognorm <- as_draws_df(h7_lognormal_model_WRprior)
cat("Percentage of posterior for b_EduTypeEvo > 0:",
    round(mean(draws_lognorm$b_EduTypeEvo > 0) * 100, 2), "%\n")
cat("Percentage of posterior for b_EduTypeEvo < 0:",
    round(mean(draws_lognorm$b_EduTypeEvo < 0) * 100, 2), "%\n\n")


# --- 6.3: CONFIRM MODEL FIT WITH A POSTERIOR PREDICTIVE CHECK ---

# This plot should confirm that the log-normal model's predictions (light lines)
# are a good match for the distribution of the actual data (dark line).
pp_check(h7_lognormal_model_WRprior, ndraws = 100) +
  ggtitle("Posterior Predictive Check for H7 Log-Normal Model")


# --- 6.4: PLOT THE POSTERIOR DISTRIBUTION FOR THE LOG-NORMAL MODEL ---

# This plot visualizes the uncertainty in our estimate for the effect of
# education type (EduTypeEvo). The x-axis is the log-ratio, where 0 means
# no difference between the groups.

# Step 1: Extract the posterior draws for the coefficient.
log_ratio_draws_h7 <- as_draws_df(h7_lognormal_model_WRprior)$b_EduTypeEvo

# Step 2: Manually calculate the density data for plotting.
density_data_h7_log <- density(log_ratio_draws_h7)
density_df_h7_log <- data.frame(x = density_data_h7_log$x, y = density_data_h7_log$y)

# Step 3: Calculate the probability of the effect being positive.
prob_positive_h7_log <- mean(log_ratio_draws_h7 > 0) * 100

# Step 4: Build the plot using ggplot2.
ggplot(density_df_h7_log, aes(x = x, y = y)) +
  
  # Layer 1: Draw the outline of the density curve.
  geom_line(linewidth = 1) +
  
  # Layer 2: Shade the area where the log-ratio is NEGATIVE (Evo < Gen).
  geom_ribbon(data = subset(density_df_h7_log, x < 0),
              aes(ymin = 0, ymax = y),
              fill = "#D55E00", alpha = 0.3) +
  
  # Layer 3: Shade the area where the log-ratio is POSITIVE (Evo > Gen).
  geom_ribbon(data = subset(density_df_h7_log, x > 0),
              aes(ymin = 0, ymax = y),
              fill = "#D55E00", alpha = 0.8) +
  
  # Layer 4: Add a vertical dashed line at x = 0 (no effect).
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  
  # Layer 5: Add informative labels and titles.
  labs(
    title = "Belief in Usefulness of Negative Emotions after Evo vs. Gen Education (H7)",
    subtitle = paste0("Posterior from Log-Normal Model. Probability of positive effect (Evo > Gen): ", round(prob_positive_h7_log, 2), "%"),
    x = "Log-Ratio of Post-Education Scores (Evo / Gen)",
    y = "Density"
  ) +
  
  # Layer 6: Use a clean theme.
  theme_minimal(base_size = 15)



########                      #########
########  JITTER RAW FOR H7   #########
########                      #########


# --- STEP 1: PREPARE THE DATA FOR PLOTTING H7 ---

# We will calculate the average sum scores for the EBQ items and then
# reshape the data into a long format suitable for ggplot.
plot_data_h7_long <- psy_data %>%
  # Use rowwise() to perform calculations on each row
  rowwise() %>%
  mutate(
    # Calculate the average Pre_SumScore (divided by 4 items)
    Pre_AvgScore = sum(c(
      `EBQ3_There_is_very_little_use_for_negative_emotions`,
      `EBQ7_People_dont_need_their_negative_emotions`,
      `EBQ11_Negative_emotions_are_harmful`,
      `EBQ15_The_presence_of_negative_emotions_is_a_bad_thing_for_people`
    ), na.rm = TRUE) / 4,
    
    # Calculate the average Post_SumScore (divided by 4 items)
    Post_AvgScore = sum(c(
      `EBQ23_There_is_very_little_use_for_negative_emotions`,
      `EBQ27_People_dont_need_their_negative_emotions`,
      `EBQ211_Negative_emotions_are_harmful`,
      `EBQ215_The_presence_of_negative_emotions_is_a_bad_thing_for_people`
    ), na.rm = TRUE) / 4
  ) %>%
  # It's important to ungroup after rowwise operations
  ungroup() %>%
  # Select only the necessary columns
  select(EduType, Pre_AvgScore, Post_AvgScore) %>%
  # Reshape the data from wide to long format
  tidyr::pivot_longer(
    cols = c(Pre_AvgScore, Post_AvgScore),
    names_to = "time",
    values_to = "score"
  ) %>%
  # Set the correct order and labels for the x-axis
  mutate(time = factor(time, levels = c("Pre_AvgScore", "Post_AvgScore"), labels = c("Pre-Education", "Post-Education")))

# --- STEP 2: CALCULATE SUMMARY STATISTICS (MEAN AND 95% CI) ---

summary_data_h7 <- plot_data_h7_long %>%
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

# We use the same aesthetic specifications from the H1-H4 plots.
ggplot(summary_data_h7, aes(x = time, y = mean_score, group = EduType, color = EduType)) +
  
  # 1. Add individual raw data points (jittered)
  geom_jitter(data = plot_data_h7_long, aes(y = score), width = 0.2, alpha = 0.25, size = 1.8) +
  
  # 2. Add error bars for the 95% CI
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.08, linewidth = 0.8) +
  
  # 3. Add lines connecting the means
  geom_line(linewidth = 1.2) +
  
  # 4. Add large dots for the means
  geom_point(size = 4) +
  
  # 5. Customize labels and theme
  labs(
    title = "Belief in Usefulness of Negative Emotions (H7)",
    subtitle = "Change in clinician beliefs from pre- to post-education",
    y = "Average EBQ Score",
    x = "Time Point",
    color = "Education Type"
  ) +
  scale_color_manual(
    values = c("Gen" = "#0072B2", "Evo" = "#D55E00"),
    labels = c("Gen" = "Genetic", "Evo" = "Evolutionary")
  ) +
  coord_cartesian(ylim = c(1, 7)) +
  theme_minimal(base_size = 15)

