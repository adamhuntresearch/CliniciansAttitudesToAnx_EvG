[ReadMe_Hunt_Carp.txt](https://github.com/user-attachments/files/24436537/ReadMe_Hunt_Carp.txt)
This is the ReadMe the data and R scripts for Clinician’s attitudes to evolutionary and genetic explanations for anxiety: a cluster-randomised study of stigmatisation
Hunt et al. 2025

Analytical plan pre-registered at https://osf.io/7x98y

These are the Rscript files and overview of their purpose (each file also contains substantial in-script comments with more detailed descriptions).
They need to be run sequentially.

To save time running the models yourself, the .rds files are available to download at: https://drive.google.com/drive/folders/1uWH-woUBGUXtpitg1pcg7J-hy42Ry4zR?usp=sharing

The csv. sheets of data which is output throughout different stages of cleaning are in the cleaned_data folder.

1. CLEANING AND PREPPING

Prep1_AnonymisingCleaningCompiling:
Anonymises and cleans Prolific Raw files and combines them into a single datasheet. CANNOT BE RUN BY NON-EXPERIMENTERS TO PRESERVE ANONYMITY.
Outputs "cleaned_psych_data_FULL.csv"

Prep2_EBQRecoding:
Recode EBQ Questionnaire Values from raw Prolific Data before running the analysis.
Outputs "cleaned_psych_data_EBQRECODED.csv"


Prep3_ExcludeNonParticipants:
Excludes participants who failed attention checks, didn't consent, or didn't complete.
Outputs "final_analysis_dataset.csv"

Prep4_RemoveColumns:
A simple script to remove unnecessary columns of completion, duration, and consent.
Outputs "FINALLY_analysis_ready_data.csv"


2. PRIMARY ANALYSES
All Analysis Scripts run with "FINALLY_analysis_ready_data.csv" 

Ana1_H1toH4:
Runs Bayesian cumulative-link ordinal regression models for the pre-post preregistered hypotheses 1-4. Draws some figures.

Ana2_H5H6:
Runs Bayesian cumulative-link ordinal regression models for the post-only preregistered hypotheses 5 & 6. Draws some figures.

Ana3_H7:
Sum scores the EBQ and then runs the Bayesian linear mixed-effects gaussian model. Draws some figures.

Ana_Combined_H1toH4:
A multiple comparison check which puts H1 to H4 into a single multilevel model.

Ana_Figs:
A dedicated figure generation script for the above primary pre-registered analyses


3. EXPLORATORY ANALYSES

ExAna_Descriptive:
Looks at demographics, including participant age, specialisation, and stage of training/position.

ExAna_ExpH1toExpH4:
Runs Bayesian cumulative-link ordinal regression models for the pre-post not preregistered hypotheses from the ExAnSt questionnaire. Draws all their figures.


Contact ah2422 [at] cam.ac.uk (Adam Hunt) for queries.






