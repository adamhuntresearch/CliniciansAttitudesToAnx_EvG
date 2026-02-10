### This script uses a single multilevel model as a way to model multiple comparision (i.e. running 4 hypothesis tests) in Hunt et al.
## This is a robustness check for H1-H6. H1-H4 are taken together, as they have the same 
# design (ordinal, single item, pre-post) as are H5 and H6 (ordinal, single item, post only)
library(brms)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(patchwork)

# Set a seed for reproducibility.
set.seed(4743)

## Load the dataset
d1 <- read_csv("MComparision_Check_h1_h4.csv")
str(d1)

## Format the dataset a bit so that our models run appropriately
d1 <- d1 %>%
  mutate(
    EduType = factor(EduType, levels = c("Gen", "Evo")),  ## Reading these variables as factors
    SessionID = as.factor(SessionID),
    Delivery = as.factor(Delivery),
    H1_pre_z = scale(H1_pre),                                 ## Standardising the pre-scores
    H2_pre_z = scale(H2_pre),
    H3_pre_z = scale(H3_pre),
    H4_pre_z = scale(H4_pre)
  )

## Now to run all of the hypotheses in one multilevel model, we first need to format the data further.
## That is, we need a single post and pre score column which contains the respective scores from all hypotheses in them
## This would also mean that we need to create a new column which identifies the specific hypothesis that each row of pre- and post- scores belong to
## Doing this below

d2 <- bind_rows(
  d1 %>%
    transmute(
      PID = ParticipantID,
      EduType,
      SessionID,
      Delivery,
      post_score = H1_post_ord,
      pre_score_z = H1_pre_z,
      hypothesis = "H1"
    ),
  
  d1 %>%
    transmute(
      PID = ParticipantID,
      EduType,
      SessionID,
      Delivery,
      post_score = H2_post_ord,
      pre_score_z = H2_pre_z,
      hypothesis = "H2"
    ),
  
  d1 %>%
    transmute(
      PID = ParticipantID,
      EduType,
      SessionID,
      Delivery,
      post_score = H3_post_ord,
      pre_score_z = H3_pre_z,
      hypothesis = "H3"
    ),
  
  d1 %>%
    transmute(
      PID = ParticipantID,
      EduType,
      SessionID,
      Delivery,
      post_score = H4_post_ord,
      pre_score_z = H4_pre_z,
      hypothesis = "H4"
    ))

## Let's do a sanity check to see whether there is anyone in this dataset who doesn't have 4 observations assigned to them
## 1 participant answered 4 questions for 4 different hypotheses. So each participant should have 4 observations.
## Looks like the transformation worked
d2 %>% count(PID) %>% filter(n != 4)

## Another sanity check would be to see that each hypothesis has 171 observations each
## Looks like it worked
d2 %>% count(hypothesis)

## Okay, let's treat hypothesis as a factor
d2 <- d2 %>% mutate(
  hypothesis = as.factor(hypothesis)
)


### Running the model
#################################

mc_check_mod <- brm(post_score ~ 1 + EduType + pre_score_z + Delivery + (1 + EduType + pre_score_z + Delivery | hypothesis) + (1|SessionID) + (1|PID),
                    data = d2,
                    family = cumulative("logit"),
                    prior = c(prior(normal(0, 1.5), class = "Intercept"),
                              prior(normal(0, 1), class = "b"),
                              prior(exponential(1), class = "sd")),
                    warmup = 3000,
                    iter = 6000,
                    chains = 4,
                    cores = 8,
                    control = list(adapt_delta = 0.99),
                    seed = 27)

## Run only the second line below for directly loading the model that we have already run
saveRDS(mc_check_mod, file = "h1_h4_mc_check.rds")
mc_check_mod <- readRDS(file = "h1_h4_mc_check.rds")

## Inspect the model
summary(mc_check_mod)

## Drawing posteriors from the model and plotting first the average effect across all 4 hypotheses
mod_draws <- as_draws_df(mc_check_mod)
plot(density(mod_draws$b_EduTypeEvo))

## What is the probability of an average effect across the four hypotheses?
hypothesis(mc_check_mod, "EduTypeEvo > 0")


### Okay, now our goal is to extract the effect for each hypothesis and plot the posterior to see whether it corresponds well with the models that were run separately.
##############################################################################################################

## First we select all the random deviations of the intervention for each hypothesis
hypothesis_slopes <- mod_draws %>%
  select(starts_with("r_hypothesis[")) %>%
  select(contains("EduTypeEvo"))

## Then, we get the average fixed effect across all the four hypotheses
fixed_effect_evo <- mod_draws$b_EduTypeEvo

## Now, we add the fixed effect + random deviation (for each hypothesis) to get the effect of the intervention for each hypothesis
hypothesis_slopes_full <- hypothesis_slopes %>%
  mutate(
    Hypothesis1 = `r_hypothesis[H1,EduTypeEvo]` + fixed_effect_evo,
    Hypothesis2 = `r_hypothesis[H2,EduTypeEvo]` + fixed_effect_evo,
    Hypothesis3 = `r_hypothesis[H3,EduTypeEvo]` + fixed_effect_evo,
    Hypothesis4 = `r_hypothesis[H4,EduTypeEvo]` + fixed_effect_evo,
  )

## Okay, now we convert each effect to odds ratio
hypothesis_slopes_full_OR <- hypothesis_slopes_full %>%
  mutate(
    Hypothesis1_OR = exp(Hypothesis1),
    Hypothesis2_OR = exp(Hypothesis2),
    Hypothesis3_OR = exp(Hypothesis3),
    Hypothesis4_OR = exp(Hypothesis4)
  )
  
## On to plotting each effect now
## Hypothesis 1
plot(density(hypothesis_slopes_full_OR$Hypothesis1_OR))
mean(hypothesis_slopes_full_OR$Hypothesis1_OR > 1)

## Hypothesis 2
plot(density(hypothesis_slopes_full_OR$Hypothesis2_OR))
mean(hypothesis_slopes_full_OR$Hypothesis2_OR > 1)

## Hypothesis 3
plot(density(hypothesis_slopes_full_OR$Hypothesis3_OR))
mean(hypothesis_slopes_full_OR$Hypothesis3_OR > 1)

## Hypothesis 4
plot(density(hypothesis_slopes_full_OR$Hypothesis4_OR))
mean(hypothesis_slopes_full_OR$Hypothesis4_OR > 1)



### Now for hypotheses H5 and H6 from Hunt et al.
###
###

# Set a seed for reproducibility.
set.seed(4743)


# STEP 1: LOAD AND PREPARE THE DATA
#################################################################

## Load the dataset used for H5 and H6, using a unique name
d1_h5h6 <- read_csv("FINALLY_analysis_ready_data.csv")

## Format the dataset so that our models run appropriately
d1_h5h6 <- d1_h5h6 %>%
  mutate(
    # Ensure outcome variables are ordered factors
    H5_post_ord = factor(`Q88_Do_you_think_this_information_is_useful_for_patients`, ordered = TRUE, levels = 1:5),
    H6_post_ord = factor(`Q90_Do_you_think_this_information_is_useful_for_clinicians`, ordered = TRUE, levels = 1:5),
    
    # Reading predictor variables as factors
    EduType = factor(EduType, levels = c("Gen", "Evo")),
    SessionID = as.factor(SessionID),
    Delivery = as.factor(Delivery),
    PID = as.factor(ParticipantID) # Assuming ParticipantID is the correct column name
  )


# STEP 2: RESHAPE DATA TO LONG FORMAT
#################################################################
## To run both hypotheses in one model, we need a single post-score column
## and a new column to identify which hypothesis each score belongs to.

d2_h5h6 <- bind_rows(
  d1_h5h6 %>%
    transmute(
      PID,
      EduType,
      SessionID,
      Delivery,
      post_score = H5_post_ord,
      hypothesis = "H5"
    ),
  
  d1_h5h6 %>%
    transmute(
      PID,
      EduType,
      SessionID,
      Delivery,
      post_score = H6_post_ord,
      hypothesis = "H6"
    )
)

## Let's do a sanity check. Each participant should now have 2 observations.
## This will return an empty table if the transformation worked correctly.
d2_h5h6 %>% count(PID) %>% filter(n != 2)

## Another sanity check: each hypothesis should have the same number of observations.
d2_h5h6 %>% count(hypothesis)

## Treat 'hypothesis' as a factor for the model
d2_h5h6 <- d2_h5h6 %>% 
  mutate(hypothesis = as.factor(hypothesis))


# STEP 3: RUN THE MULTILEVEL MODEL
#################################################################
## This model formula is adjusted for H5/H6 by removing the pre_score_z term.

mc_check_mod_h5h6 <- brm(
  post_score ~ 1 + EduType + Delivery + (1 + EduType + Delivery | hypothesis) + (1|SessionID) + (1|PID),
  data = d2_h5h6,
  family = cumulative("logit"),
  prior = c(prior(normal(0, 1.5), class = "Intercept"),
            prior(normal(0, 1), class = "b"),
            prior(exponential(1), class = "sd")),
  warmup = 3000,
  iter = 6000,
  chains = 4,
  cores = 8, # Adjust cores as needed for your machine
  control = list(adapt_delta = 0.99),
  seed = 27
)

## Save the model object to avoid re-running it every time
saveRDS(mc_check_mod_h5h6, file = "h5_h6_mc_check.rds")
mc_check_mod_h5h6 <- readRDS(file = "h5_h6_mc_check.rds") # Use this line to load the saved model


# STEP 4: ANALYZE THE RESULTS
#################################################################

## Inspect the model summary
summary(mc_check_mod_h5h6)

## Extract posterior draws from the model
mod_draws_h5h6 <- as_draws_df(mc_check_mod_h5h6)

## What is the probability of an average positive effect across both hypotheses?
hypothesis(mc_check_mod_h5h6, "EduTypeEvo > 0")

## Now, let's extract the specific effect for each hypothesis
#############################################################

## First, select all the random deviations of the intervention for each hypothesis
hypothesis_slopes_h5h6 <- mod_draws_h5h6 %>%
  select(starts_with("r_hypothesis[")) %>%
  select(contains("EduTypeEvo"))

## Then, get the average fixed effect across both hypotheses
fixed_effect_evo_h5h6 <- mod_draws_h5h6$b_EduTypeEvo

## Now, add the fixed effect + random deviation to get the total intervention effect for each hypothesis
hypothesis_slopes_full_h5h6 <- hypothesis_slopes_h5h6 %>%
  mutate(
    Hypothesis5 = `r_hypothesis[H5,EduTypeEvo]` + fixed_effect_evo_h5h6,
    Hypothesis6 = `r_hypothesis[H6,EduTypeEvo]` + fixed_effect_evo_h5h6
  )

## Convert each effect to an odds ratio (OR)
hypothesis_slopes_full_OR_h5h6 <- hypothesis_slopes_full_h5h6 %>%
  mutate(
    Hypothesis5_OR = exp(Hypothesis5),
    Hypothesis6_OR = exp(Hypothesis6)
  )

## Finally, plot the posterior densities and calculate the probability of a positive effect (OR > 1)

## Hypothesis 5
plot(density(hypothesis_slopes_full_OR_h5h6$Hypothesis5_OR), main = "Posterior Density for H5 Effect (Odds Ratio)")
abline(v = 1, col = "red", lty = 2)
prob_h5_positive <- mean(hypothesis_slopes_full_OR_h5h6$Hypothesis5_OR > 1)
print(paste("Probability of a positive effect for H5:", round(prob_h5_positive, 4)))

## Hypothesis 6
plot(density(hypothesis_slopes_full_OR_h5h6$Hypothesis6_OR), main = "Posterior Density for H6 Effect (Odds Ratio)")
abline(v = 1, col = "red", lty = 2)
prob_h6_positive <- mean(hypothesis_slopes_full_OR_h5h6$Hypothesis6_OR > 1)
print(paste("Probability of a positive effect for H6:", round(prob_h6_positive, 4)))