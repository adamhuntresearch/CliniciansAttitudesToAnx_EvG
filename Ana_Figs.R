################################################################################
#
#   DEDICATED FIGURE GENERATION SCRIPT
#
#   This script loads all necessary data and saved model objects to reproduce
#   all figures for the Hunt, Carpenter et al. manuscript.
#
#   It assumes that the following files are in the working directory:
#   - FINALLY_analysis_ready_data.csv
#   - All .rds model files saved from the analysis scripts.
#
################################################################################


#-----------------------------------------------------------------------------
# STEP 1: SETUP - LOAD LIBRARIES, DATA, AND MODELS
#-----------------------------------------------------------------------------

# --- Load Libraries ---

library(brms)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggridges)

# --- Load Raw Data ---
# The raw data is needed for the jitter plots.
psy_data <- read_csv("FINALLY_analysis_ready_data.csv")

# --- Load Saved Model Objects ---
# This avoids re-running any time-consuming analyses.

# Models for H1-H4 (Evo vs. Gen Comparison)
h1_logit_model_WRprior <- readRDS("h1_logit_model_WRprior.rds")
h2_logit_model_WRprior <- readRDS("h2_logit_model_WRprior.rds")
h3_logit_model_WRprior <- readRDS("h3_logit_model_WRprior.rds")
h4_logit_model_WRprior <- readRDS("h4_logit_model_WRprior.rds")

# Models for H5-H6 (Usefulness Comparison)
h5_logit_model_WRprior <- readRDS("h5_logit_model_WRprior.rds")
h6_logit_model_WRprior <- readRDS("h6_logit_model_WRprior.rds")

# Models for Pre-Post Analysis within Evolutionary Arm
evo_prepost_h1_probit_model <- readRDS("evo_prepost_h1_probit_model.rds")
evo_prepost_h2_probit_model <- readRDS("evo_prepost_h2_probit_model.rds")
evo_prepost_h3_probit_model <- readRDS("evo_prepost_h3_probit_model.rds")
evo_prepost_h4_probit_model <- readRDS("evo_prepost_h4_probit_model.rds")



#-----------------------------------------------------------------------------
# FIGURE 1: POSTERIOR DISTRIBUTIONS FOR H1-H6 (EVO VS. GEN)
#-----------------------------------------------------------------------------
# This figure combines the six main posterior distribution plots (H1-H6)
# into a single 3x2 grid, showing the effect of Evolutionary vs. Genetic
# psychoeducation in Odds Ratios.

# --- H1: Posterior Plot ---
log_odds_h1 <- as_draws_df(h1_logit_model_WRprior)$b_EduTypeEvo
odds_ratio_h1 <- exp(log_odds_h1)
density_df_h1 <- data.frame(density(odds_ratio_h1)[c("x", "y")])
prob_pos_h1 <- mean(odds_ratio_h1 > 1) * 100
plot_h1 <- ggplot(density_df_h1, aes(x, y)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = subset(density_df_h1, x > 1), aes(ymin = 0, ymax = y), fill = "#D55E00", alpha = 0.8) +
  geom_ribbon(data = subset(density_df_h1, x < 1), aes(ymin = 0, ymax = y), fill = "#D55E00", alpha = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(title = "Optimism about recovery (H1)", subtitle = paste0("Prob. of positive effect: ", round(prob_pos_h1, 2), "%"), x = "Odds Ratio", y = "Density") +
  theme_minimal(base_size = 12)

# --- H2: Posterior Plot ---
log_odds_h2 <- as_draws_df(h2_logit_model_WRprior)$b_EduTypeEvo
odds_ratio_h2 <- exp(log_odds_h2)
density_df_h2 <- data.frame(density(odds_ratio_h2)[c("x", "y")])
prob_pos_h2 <- mean(odds_ratio_h2 > 1) * 100
plot_h2 <- ggplot(density_df_h2, aes(x, y)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = subset(density_df_h2, x > 1), aes(ymin = 0, ymax = y), fill = "#D55E00", alpha = 0.8) +
  geom_ribbon(data = subset(density_df_h2, x < 1), aes(ymin = 0, ymax = y), fill = "#D55E00", alpha = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(title = "Willingness to share diagnosis (H2)", subtitle = paste0("Prob. of positive effect: ", round(prob_pos_h2, 2), "%"), x = "Odds Ratio", y = NULL) +
  theme_minimal(base_size = 12)

# --- H3: Posterior Plot ---
log_odds_h3 <- as_draws_df(h3_logit_model_WRprior)$b_EduTypeEvo
odds_ratio_h3 <- exp(log_odds_h3)
density_df_h3 <- data.frame(density(odds_ratio_h3)[c("x", "y")])
prob_pos_h3 <- mean(odds_ratio_h3 > 1) * 100
plot_h3 <- ggplot(density_df_h3, aes(x, y)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = subset(density_df_h3, x > 1), aes(ymin = 0, ymax = y), fill = "#D55E00", alpha = 0.8) +
  geom_ribbon(data = subset(density_df_h3, x < 1), aes(ymin = 0, ymax = y), fill = "#D55E00", alpha = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(title = "Willingness to seek help (H3)", subtitle = paste0("Prob. of positive effect: ", round(prob_pos_h3, 2), "%"), x = "Odds Ratio", y = NULL) +
  theme_minimal(base_size = 12)

# --- H4: Posterior Plot ---
log_odds_h4 <- as_draws_df(h4_logit_model_WRprior)$b_EduTypeEvo
odds_ratio_h4 <- exp(log_odds_h4)
density_df_h4 <- data.frame(density(odds_ratio_h4)[c("x", "y")])
prob_pos_h4 <- mean(odds_ratio_h4 > 1) * 100
plot_h4 <- ggplot(density_df_h4, aes(x, y)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = subset(density_df_h4, x > 1), aes(ymin = 0, ymax = y), fill = "#D55E00", alpha = 0.8) +
  geom_ribbon(data = subset(density_df_h4, x < 1), aes(ymin = 0, ymax = y), fill = "#D55E00", alpha = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(title = "Efficacy of psychosocial interventions (H4)", subtitle = paste0("Prob. of positive effect: ", round(prob_pos_h4, 2), "%"), x = "Odds Ratio", , y = "Density") +
  theme_minimal(base_size = 12)

# --- H5: Posterior Plot ---
log_odds_h5 <- as_draws_df(h5_logit_model_WRprior)$b_EduTypeEvo
odds_ratio_h5 <- exp(log_odds_h5)
density_df_h5 <- data.frame(density(odds_ratio_h5)[c("x", "y")])
prob_pos_h5 <- mean(odds_ratio_h5 > 1) * 100
plot_h5 <- ggplot(density_df_h5, aes(x, y)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = subset(density_df_h5, x > 1), aes(ymin = 0, ymax = y), fill = "#D55E00", alpha = 0.8) +
  geom_ribbon(data = subset(density_df_h5, x < 1), aes(ymin = 0, ymax = y), fill = "#D55E00", alpha = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(title = "Usefulness for Patients (H5)", subtitle = paste0("Prob. of positive effect: ", round(prob_pos_h5, 2), "%"), x = "Odds Ratio", y = NULL) +
  coord_cartesian(xlim = c(0, 15)) +
  theme_minimal(base_size = 12)

# --- H6: Posterior Plot ---http://127.0.0.1:35143/graphics/09ce9757-e875-4289-953c-017e78b8bba6.png
log_odds_h6 <- as_draws_df(h6_logit_model_WRprior)$b_EduTypeEvo
odds_ratio_h6 <- exp(log_odds_h6)
density_df_h6 <- data.frame(density(odds_ratio_h6)[c("x", "y")])
prob_pos_h6 <- mean(odds_ratio_h6 > 1) * 100
plot_h6 <- ggplot(density_df_h6, aes(x, y)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = subset(density_df_h6, x > 1), aes(ymin = 0, ymax = y), fill = "#D55E00", alpha = 0.8) +
  geom_ribbon(data = subset(density_df_h6, x < 1), aes(ymin = 0, ymax = y), fill = "#D55E00", alpha = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(title = "Usefulness for Clinicians (H6)", subtitle = paste0("Prob. of positive effect: ", round(prob_pos_h6, 2), "%"), x = "Odds Ratio", y = NULL) +
  theme_minimal(base_size = 12)

# --- Combine H1-H6 ---
combined_plot_H1_H6 <- (plot_h1 + plot_h2 + plot_h3) / (plot_h4 + plot_h5 + plot_h6) +
  plot_annotation(
    title = 'Effect of Evolutionary vs. Genetic Education on Clinician Attitudes',
    theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))
  )

# Print the plot to the viewer
print(combined_plot_H1_H6)




#-----------------------------------------------------------------------------
# FIGURE 2: RAW DATA JITTER PLOTS FOR H1-H6
#-----------------------------------------------------------------------------
#
# This code creates the raw data jitter plots for H1 to H6 and combines them.
# The raw output ends up with two legends because pathwork knows there are two different 
# analyses; it's much easier to edit these in an image editor (MS Paint) than code them, so that's what I did. 

# --- Prepare Data for H1-H4 Jitter Plots (Pre-Post Design) ---
plot_data_h1_h4_long <- psy_data %>%
  select(
    EduType,
    H1_pre = `ExAnSt1_Given_your_current_understanding_how_optimistic_are_you_about_improvement_in_symptoms_in_patients_presenting_with_anxiety_disorder`,
    H1_post = `ExAnSt21_Given_the_information_provided_in_this_education_session_how_optimistic_are_you_about_improvement_in_symptoms_in_patients_presenting_with_anxiety_disorder`,
    H2_pre = `ExAnSt2_Given_your_current_understanding_how_willing_do_you_think_patients_would_be_to_share_anxiety_disorder_diagnoses`,
    H2_post = `ExAnSt22_Given_the_information_provided_in_this_education_session_how_willing_do_you_think_patients_would_be_to_share_anxiety_disorder_diagnoses`,
    H3_pre = `ExAnSt8_Given_current_public_knowledge_of_anxiety_how_willing_do_you_think_people_are_to_seek_psychiatric_help_for_anxiety`,
    H3_post = `ExAnSt28_If_the_information_provided_in_this_education_session_was_publicly_known_how_willing_do_you_think_people_would_be_to_seek_psychiatric_help_for_anxiety`,
    H4_pre = `ExAnSt4_Given_your_current_understanding_how_effective_do_you_think_psychosocial_interventions_will_be_in_improving_outcomes_in_anxiety_disorder`,
    H4_post = `ExAnSt24_Given_the_information_provided_in_this_education_session_how_effective_do_you_think_psychosocial_interventions_will_be_in_improving_outcomes_in_anxiety_disorder`
  ) %>%
  pivot_longer(
    cols = -EduType,
    names_to = c("hypothesis", "time"),
    names_sep = "_",
    values_to = "score"
  ) %>%
  mutate(
    time = factor(time, levels = c("pre", "post"), labels = c("Pre-Education", "Post-Education")),
    EduType = factor(EduType, levels = c("Gen", "Evo"), labels = c("Genetic", "Evolutionary"))
  )

summary_data_h1_h4 <- plot_data_h1_h4_long %>%
  group_by(hypothesis, EduType, time) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    ci_lower = mean_score - 1.96 * (sd(score, na.rm = TRUE) / sqrt(n())),
    ci_upper = mean_score + 1.96 * (sd(score, na.rm = TRUE) / sqrt(n())),
    .groups = 'drop'
  )

# --- Prepare Data for H5-H6 Jitter Plots (Post-Only Design) ---
plot_data_h5_h6_long <- psy_data %>%
  select(
    EduType,
    H5 = `Q88_Do_you_think_this_information_is_useful_for_patients`,
    H6 = `Q90_Do_you_think_this_information_is_useful_for_clinicians`
  ) %>%
  pivot_longer(
    cols = -EduType,
    names_to = "hypothesis",
    values_to = "score"
  ) %>%
  mutate(
    EduType = factor(EduType, levels = c("Gen", "Evo"), labels = c("Genetic", "Evolutionary"))
  )

summary_data_h5_h6 <- plot_data_h5_h6_long %>%
  group_by(hypothesis, EduType) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    ci_lower = mean_score - 1.96 * (sd(score, na.rm = TRUE) / sqrt(n())),
    ci_upper = mean_score + 1.96 * (sd(score, na.rm = TRUE) / sqrt(n())),
    .groups = 'drop'
  )

# --- Build Individual Jitter Plots ---

# Function for H1-H4 (Pre-Post)
create_jitter_plot_prepost <- function(hyp_num, title_text, y_label) {
  hyp_data_long <- filter(plot_data_h1_h4_long, hypothesis == paste0("H", hyp_num))
  hyp_data_summary <- filter(summary_data_h1_h4, hypothesis == paste0("H", hyp_num))
  
  ggplot(hyp_data_summary, aes(x = time, y = mean_score, group = EduType, color = EduType)) +
    geom_jitter(data = hyp_data_long, aes(y = score), width = 0.2, alpha = 0.25, size = 1.8) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.08, linewidth = 0.8) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 4) +
    labs(title = title_text, y = y_label, x = NULL) +
    scale_color_manual(values = c("Genetic" = "#0072B2", "Evolutionary" = "#D55E00")) +
    coord_cartesian(ylim = c(1, 7)) +
    theme_minimal(base_size = 14)
}

# Function for H5-H6 (Post-Only)
create_jitter_plot_postonly <- function(hyp_num, title_text, y_label, y_limit) {
  hyp_data_long <- filter(plot_data_h5_h6_long, hypothesis == hyp_num)
  hyp_data_summary <- filter(summary_data_h5_h6, hypothesis == hyp_num)
  
  ggplot(hyp_data_long, aes(x = EduType, y = score, color = EduType)) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 2) +
    geom_point(data = hyp_data_summary, aes(y = mean_score), size = 6) +
    geom_errorbar(
      data = hyp_data_summary,
      aes(ymin = ci_lower, ymax = ci_upper, y = mean_score),
      width = 0.1, linewidth = 1
    ) +
    labs(title = title_text, y = y_label, x = NULL) +
    scale_color_manual(values = c("Genetic" = "#0072B2", "Evolutionary" = "#D55E00")) +
    coord_cartesian(ylim = y_limit) +
    theme_minimal(base_size = 14)
}

# Create all six plots
plot_h1_jitter <- create_jitter_plot_prepost(1, "Optimism About Patient Recovery (H1)", "Clinician Rating (1-7)")
plot_h2_jitter <- create_jitter_plot_prepost(2, "Willingness to Share Diagnosis (H2)", NULL)
plot_h3_jitter <- create_jitter_plot_prepost(3, "Willingness to Seek Help (H3)", NULL)
plot_h4_jitter <- create_jitter_plot_prepost(4, "Efficacy of Psychosocial Interventions (H4)", "Clinician Rating (1-7)")
plot_h5_jitter <- create_jitter_plot_postonly("H5", "Usefulness for Patients (H5)", "Usefulness Rating (1-5)", c(1, 5))
plot_h6_jitter <- create_jitter_plot_postonly("H6", "Usefulness for Clinicians (H6)", NULL, c(1, 5))

# --- Combine Figure 2 ---
combined_plot_H1_H6_jitter <- (plot_h1_jitter + plot_h2_jitter + plot_h3_jitter) / (plot_h4_jitter + plot_h5_jitter + plot_h6_jitter) +
  plot_layout(guides = 'collect') &
  plot_annotation(
    title = 'Original Likert Scale Responses of Clinician Attitudes',
    theme = theme(plot.title = element_text(size = 20, hjust = 0.5))
  ) &
  theme(legend.position = 'bottom')

print(combined_plot_H1_H6_jitter)

