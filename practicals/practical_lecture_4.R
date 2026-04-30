############################################################
#### MPhil PHS: QHIA module (Lecutre 4 - Air Pollution) ####
############################################################

# A M D Navaratnam (amdn2@cam.ac.uk) 2026-04-29

# Clean the environment
rm(list = ls())

# Set whether to print or not
print <- FALSE

##### 1. Load packages and read in file #####

##### 1a. Load packages #####


#List of packages you want
mypackages <- c("tidyverse", "ggplot2")

# For loop to check if packages in library, and install +/- load if not
for(p in mypackages) {
  if(!require(p, character.only = TRUE)) {
    install.packages(p)
    library(p, character.only = TRUE)
  }
}

###### 1b. Read the data ######

df <- read.csv("data/lecture_4/synthetic_population_pm25.csv")

###### 2. Set parameters ######

RR_per_5 <- 1.10
increment <- 5
incidence <- 0.06 #60 per 1000 IHD death per year

##### 3. Convert exposure to RR ####

##### What is the RR at baseline for each person
df$RR_baseline # <-

##### What is the RR in the counterfactual scenario for each person
df$RR_scenario # <-


##### 4. Compute PIF

##### What is the PIF


##### 5. Calculate cases
cases_baseline <- incidence*nrow(df)
cases_prevented <- cases_baseline*PIF

paste0("The PIF is ", round(PIF,2), " so the cases prevented are ", round(cases_prevented, 2))
