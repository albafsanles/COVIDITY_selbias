#########################################################################################################
## 1. Create R environment
#########################################################################################################

rm(list=ls())

library(tidyverse)
library(haven)
library(ResourceSelection)
library(pROC)
library(car)
library(mice)
library(data.table)
library(compareGroups)
library(readr)
library(VennDiagram)
library(RColorBrewer)
library("forestplot")



#########################################################################################################
## 2. Set a working folder to save the output
#########################################################################################################



#########################################################################################################
# 3. Read in the datasets with predictors and outcomes
#########################################################################################################



#########################################################################################################
# 4. Make a copy, so that the raw file is kept 
#########################################################################################################



#########################################################################################################
####### 5. Check/prepare the variables 
#########################################################################################################

#########################################################################################################
#Code predictor variables as factors and assign interpretable names
  
db <- db %>%  
  mutate(categorical_predictor1 = as.factor(categorical_predictor1)) %>%
  mutate(categorical_predictor2 = as.factor(categorical_predictor2)) %>%
  mutate(outcome1 = as.factor(outcome1)) %>%
  mutate(outcome2 = as.factor(outcome2)) 

db <- db %>%  
  mutate(continuous_predictor1 = as.numeric(continuous_predictor1)) %>%
  mutate(continuous_predictor2 = as.numeric(continuous_predictor2))

db <- db %>%  
  mutate(identifier = as.character(identifier)) 


#########################################################################################################
# For univariate analyses

var_list <- list("outcome1", "outcome2")
table_OR_names <- c("Outcome", "Variable", "coef", "N", "Coefficient_OR", "Lower_CI", "Upper_CI", "pvalue")



#########################################################################################################
####### 6a. Distribution of predictor variables in the subsamples of cases and controls for each outcome
#########################################################################################################

#########################################################################################################
# Extract descriptive tables of the samples

db_out1 <- filter(db, outcome1==1)
db_out1_res <- compareGroups( ~ categorical_predictor1 +
                                    categorical_predictor2 +
                                    continuous_predictor1 +
                                    continuous_predictor2,
                                  data = db_out1)
db_out1_res <- createTable(db_out1_res)
export2csv(db_out1_res, file="desc_out1.csv", sep=";")

db_out1_controls <- filter(db, outcome1==0)
db_out1_res_controls <- compareGroups( ~ categorical_predictor1 +
                                    categorical_predictor2 +
                                    continuous_predictor1 +
                                    continuous_predictor2,
                                  data = db_out1_controls)
db_out1_res_controls <- createTable(db_out1_res_controls)
export2csv(db_out1_res_controls, file="desc_out1_controls.csv", sep=";")

db_out2 <- filter(db, outcome2==1)
db_out2_res <- compareGroups( ~ categorical_predictor1 +
                                    categorical_predictor2 +
                                    continuous_predictor1 +
                                    continuous_predictor2,
                                  data = db_out1)
db_out2_res <- createTable(db_out2_res)
export2csv(db_out1_res, file="desc_out2.csv", sep=";")

db_out2_controls <- filter(db, outcome2==0)
db_out2_res_controls <- compareGroups( ~ categorical_predictor1 +
                                    categorical_predictor2 +
                                    continuous_predictor1 +
                                    continuous_predictor2,
                                  data = db_out2_controls)
db_out2_res_controls <- createTable(db_out2_res_controls)
export2csv(db_out2_res_controls, file="desc_out2_controls.csv", sep=";")



#########################################################################################################
####### 6b. Descriptive analyses by all the predictors
#########################################################################################################

#########################################################################################################
#### Categorical predictor 1


## Descriptive stats 
#########################################################################################################

db_long <- db %>%
  select(outcome1, outcome2,
         categorical_predictor1) %>%
  gather(outcomeName, outcomeValue, outcome1:outcome2) %>%
  mutate(outcomeName = as.factor(outcomeName)) %>%
  mutate(outcomeName = fct_relevel(outcomeName, "outcome1", "outcome2")) %>%
  mutate(outcomeValue = as.numeric(outcomeValue))

(table_desc_cat_pred1 <- db_long %>%
  filter(!is.na(outcomeValue)) %>%
  group_by(outcomeName, categorical_predictor1) %>%
  summarise(n = n(), n_Yes = sum(outcomeValue == 1, na.rm = TRUE), 
            per_Yes = mean(outcomeValue, na.rm = TRUE) * 100, 
            n_No = sum(outcomeValue == 0, na.rm = TRUE),
            per_No = 100 - mean(outcomeValue, na.rm = TRUE) * 100) %>%
  mutate(n_total = sum(n), per_total = (n / n_total) * 100) %>%
  mutate(categorical_predictor1 = ifelse(categorical_predictor1 == 1, "Name of the category (1)", "Name of the category (2)")))

write.table(table_desc_cat_pred1, file = "desc_cat_pred1.csv", row.names = FALSE, sep=";")

(table_missingness_cat_pred1 <- db_long %>%
    filter(!is.na(outcomeValue)) %>%
    group_by(outcomeName) %>%
    summarise(var = "categorical_predictor1", n = n(), nMiss = sum(is.na(categorical_predictor1)), 
              perMiss = round((sum(is.na(categorical_predictor1)) / n()) * 100, 2)))

write.table(table_missingness_cat_pred1, file = "missingness_cat_pred1.csv", row.names = FALSE, sep=";")


## Univariate logistic regression 
#########################################################################################################

table_OR_cat_predictor1 <- data.frame(matrix(ncol=8, nrow=0))

for (i in var_list) {
  print(i)
  
  ## Odds ratio results
  model <- glm(paste(i[[1]], "~ categorical_predictor1"), data = db, family = "binomial")
  print(summary(model))
  
  # Add these coefficients, CIs, n and p-values to the table
  outcome <- i[[1]]
  varName <- "categorical_predictor1 (ref = reference category)"
  coef <- names(model$coefficients)
  temp <- data.frame(outcome, varName, coef, nobs(model),
                     round(exp(model$coefficients), 2), 
                     round(exp(confint(model)[, 1]), 2),
                     round(exp(confint(model)[, 2]), 2),
                     round(coef(summary(model))[, 4], 3))
  
  colnames(temp) <- table_OR_names
  table_OR_cat_predictor1 <- rbind(table_OR_cat_predictor1, temp)
  table_OR_cat_predictor1 <- table_OR_cat_predictor1 %>%
    filter(coef != "(Intercept)")
}

write.table(table_OR_cat_predictor1, file = "Univar_OR_cat_pred1.csv", row.names = FALSE, sep=";")



#########################################################################################################
#### Continuous predictor 1


# Convert to Z-scores
#########################################################################################################

db <- db %>%
  mutate(continuous_predictor1 = (continuous_predictor1 - mean(continuous_predictor1, na.rm = TRUE)) / sd(continuous_predictor1, na.rm = TRUE))


## Descriptive stats 
#########################################################################################################

db_long <- db %>%
  select(outcome1, outcome2,
         continuous_predictor1) %>%
  gather(outcomeName, outcomeValue, outcome1:outcome2) %>%
  mutate(outcomeName = as.factor(outcomeName)) %>%
  mutate(outcomeName = fct_relevel(outcomeName, "outcome1", "outcome2")) %>%
  mutate(outcomeValue = as.numeric(outcomeValue))

(table_desc_cont_pred1 <- db_long %>%
  filter(!is.na(outcomeValue)) %>%
  group_by(outcomeName, outcomeValue) %>%
    summarise(n = n(), mean = mean(yp_bmi, na.rm = TRUE), 
              SD_Yes = sd(yp_bmi, na.rm = TRUE)) %>%
    mutate(outcomeValue = ifelse(outcomeValue == 0, "No (0)", "Yes (1)")))
    
write.table(table_desc_cont_pred1, file = "desc_cont_pred1.csv", row.names = FALSE, sep=";")

(table_missingness_cont_pred1 <- db_long %>%
    filter(!is.na(outcomeValue)) %>%
    group_by(outcomeName) %>%
    summarise(var = "continuous_predictor1", n = n(), nMiss = sum(is.na(continuous_predictor1)), 
              perMiss = round((sum(is.na(continuous_predictor1)) / n()) * 100, 2)))

write.table(table_missingness_cont_pred1, file = "missingness_cont_pred1.csv", row.names = FALSE, sep=";")


## Univariate logistic regression 
#########################################################################################################

table_OR_cont_predictor1 <- data.frame(matrix(ncol=8, nrow=0))

for (i in var_list) {
  print(i)
  
  ## Odds ratio results
  model <- glm(paste(i[[1]], "~ continuous_predictor1"), data = db, family = "binomial")
  print(summary(model))
  
  # Add these coefficients, CIs, n and p-values to the table
  outcome <- i[[1]]
  varName <- "continuous_predictor1 (z-score)"
  coef <- names(model$coefficients)
  temp <- data.frame(outcome, varName, coef, nobs(model),
                     round(exp(model$coefficients), 2), 
                     round(exp(confint(model)[, 1]), 2),
                     round(exp(confint(model)[, 2]), 2),
                     round(coef(summary(model))[, 4], 3))
  
  colnames(temp) <- table_OR_names
  table_OR_cont_predictor1 <- rbind(table_OR_cont_predictor1, temp)
  table_OR_cont_predictor1 <- table_OR_cont_predictor1 %>%
    filter(coef != "(Intercept)")
}

write.table(table_OR_cont_predictor1, file = "Univar_OR_cont_pred1.csv", row.names = FALSE, sep=";")



