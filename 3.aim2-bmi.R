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
# 3. Read in the datasets with COVID-19 and pre-pandemic data
#########################################################################################################



#########################################################################################################
# 4. Make a copy, so that the raw file is kept 
#########################################################################################################



#########################################################################################################
# 5. Check/prepare the variables
#########################################################################################################

#########################################################################################################
#Code predictor variables as factors and assign interpretable names
  
db <- db %>%  
  mutate(sex = as.factor(sex)) %>%
  mutate(smoking = as.factor(smoking)) %>%
  mutate(education = as.factor(education)) %>%
  mutate(IMD = as.factor(IMD)) %>%
  mutate(outcome1 = as.factor(outcome1)) %>%
  mutate(outcome2 = as.factor(outcome2)) 

db <- db %>%  
  mutate(age_z = as.numeric(age_z)) %>%
  mutate(bmi_z = as.numeric(bmi_z))


#########################################################################################################
# For multivariate regression models

var_list <- list("outcome1", "outcome2")
table_OR_names <- c("Outcome", "Variable", "coef", "N", "Coefficient_OR", "Lower_CI", "Upper_CI", "pvalue")


#########################################################################################################
#For stratified analyses - Stratify sample by sex

db_men <- filter(db, sex==1)
db_women <- filter(db, sex==2)


#########################################################################################################
#Model 1: adjusted for age

table_OR_bmi_model1 <- data.frame(matrix(ncol=8, nrow=0))

for (i in var_list) {
  print(i)
  
  ## Odds ratio results
  model <- glm(paste(i[[1]], "~ bmi_z + age_z"), data = db, family = "binomial")
  print(summary(model))
  
  # Add these coefficients, CIs, n and p-values to the table
  outcome <- i[[1]]
  varName <- "BMI (z-score)"
  coef <- names(model$coefficients)
  temp <- data.frame(outcome, varName, coef, nobs(model),
                     round(exp(model$coefficients), 2), 
                     round(exp(confint(model)[, 1]), 2),
                     round(exp(confint(model)[, 2]), 2),
                     round(coef(summary(model))[, 4], 3))
  
  colnames(temp) <- table_OR_names
  table_OR_bmi_model1 <- rbind(table_OR_bmi_model1, temp)
  table_OR_bmi_model1 <- table_OR_bmi_model1 %>%
    filter(coef != "(Intercept)")
}

write.table(table_OR_bmi_model1, file = "Model1_OR_bmi.csv", row.names = FALSE, sep=";")


###stratified by sex 

#Model for men

table_OR_bmi_men_model1 <- data.frame(matrix(ncol=8, nrow=0))

for (i in var_list) {
  print(i)
  
  ## Odds ratio results
  model <- glm(paste(i[[1]], "~ bmi_z + age_z"), data = db_men, family = "binomial")
  print(summary(model))
  
  # Add these coefficients, CIs, n and p-values to the table
  outcome <- i[[1]]
  varName <- "BMI (z-score)"
  coef <- names(model$coefficients)
  temp <- data.frame(outcome, varName, coef, nobs(model),
                     round(exp(model$coefficients), 2), 
                     round(exp(confint(model)[, 1]), 2),
                     round(exp(confint(model)[, 2]), 2),
                     round(coef(summary(model))[, 4], 3))
  
  colnames(temp) <- table_OR_names
  table_OR_bmi_men_model1 <- rbind(table_OR_bmi_men_model1, temp)
  table_OR_bmi_men_model1 <- table_OR_bmi_men_model1 %>%
    filter(coef != "(Intercept)")
}

write.table(table_OR_bmi_men_model1, file = "Model1_OR_bmi_men.csv", row.names = FALSE, sep=";")


#Model for women

table_OR_bmi_women_model1 <- data.frame(matrix(ncol=8, nrow=0))

for (i in var_list) {
  print(i)
  
  ## Odds ratio results
  model <- glm(paste(i[[1]], "~ bmi_z + age_z"), data = db_women, family = "binomial")
  print(summary(model))
  
  # Add these coefficients, CIs, n and p-values to the table
  outcome <- i[[1]]
  varName <- "BMI (z-score)"
  coef <- names(model$coefficients)
  temp <- data.frame(outcome, varName, coef, nobs(model),
                     round(exp(model$coefficients), 2), 
                     round(exp(confint(model)[, 1]), 2),
                     round(exp(confint(model)[, 2]), 2),
                     round(coef(summary(model))[, 4], 3))
  
  colnames(temp) <- table_OR_names
  table_OR_bmi_women_model1 <- rbind(table_OR_bmi_women_model1, temp)
  table_OR_bmi_women_model1<- table_OR_bmi_women_model1%>%
    filter(coef != "(Intercept)")
}

write.table(table_OR_bmi_women_model1, file = "Model1_OR_bmi_women.csv", row.names = FALSE, sep=";")


#########################################################################################################
#model 2: adjusted for age and smoking

table_OR_bmi_model2 <- data.frame(matrix(ncol=8, nrow=0))

for (i in var_list) {
  print(i)
  
  ## Odds ratio results
  model <- glm(paste(i[[1]], "~ bmi_z + age_z + smoking"), data = db, family = "binomial")
  print(summary(model))
  
  # Add these coefficients, CIs, n and p-values to the table
  outcome <- i[[1]]
  varName <- "BMI (z-score)"
  coef <- names(model$coefficients)
  temp <- data.frame(outcome, varName, coef, nobs(model),
                     round(exp(model$coefficients), 2), 
                     round(exp(confint(model)[, 1]), 2),
                     round(exp(confint(model)[, 2]), 2),
                     round(coef(summary(model))[, 4], 3))
  
  colnames(temp) <- table_OR_names
  table_OR_bmi_model2 <- rbind(table_OR_bmi_model2, temp)
  table_OR_bmi_model2 <- table_OR_bmi_model2 %>%
    filter(coef != "(Intercept)")
}

write.table(table_OR_bmi_model2, file = "Model2_OR_bmi_z.csv", row.names = FALSE, sep=";")


###stratified by sex 

#Model for men

table_OR_bmi_men_model2 <- data.frame(matrix(ncol=8, nrow=0))

for (i in var_list) {
  print(i)
  
  ## Odds ratio results
  model <- glm(paste(i[[1]], "~ bmi_z + age_z + smoking"), data = db_men, family = "binomial")
  print(summary(model))
  
  # Add these coefficients, CIs, n and p-values to the table
  outcome <- i[[1]]
  varName <- "BMI (z-score)"
  coef <- names(model$coefficients)
  temp <- data.frame(outcome, varName, coef, nobs(model),
                     round(exp(model$coefficients), 2), 
                     round(exp(confint(model)[, 1]), 2),
                     round(exp(confint(model)[, 2]), 2),
                     round(coef(summary(model))[, 4], 3))
  
  colnames(temp) <- table_OR_names
  table_OR_bmi_men_model2 <- rbind(table_OR_bmi_men_model2, temp)
  table_OR_bmi_men_model2 <- table_OR_bmi_men_model2 %>%
    filter(coef != "(Intercept)")
}

write.table(table_OR_bmi_men_model2, file = "Model2_OR_bmi_men.csv", row.names = FALSE, sep=";")


#Model for women

table_OR_bmi_women_model2 <- data.frame(matrix(ncol=8, nrow=0))

for (i in var_list) {
  print(i)
  
  ## Odds ratio results
  model <- glm(paste(i[[1]], "~ bmi_z + age_z + smoking"), data = db_women, family = "binomial")
  print(summary(model))
  
  # Add these coefficients, CIs, n and p-values to the table
  outcome <- i[[1]]
  varName <- "BMI (z-score)"
  coef <- names(model$coefficients)
  temp <- data.frame(outcome, varName, coef, nobs(model),
                     round(exp(model$coefficients), 2), 
                     round(exp(confint(model)[, 1]), 2),
                     round(exp(confint(model)[, 2]), 2),
                     round(coef(summary(model))[, 4], 3))
  
  colnames(temp) <- table_OR_names
  table_OR_bmi_women_model2 <- rbind(table_OR_bmi_women_model2, temp)
  table_OR_bmi_women_model2<- table_OR_bmi_women_model2 %>%
    filter(coef != "(Intercept)")
}

write.table(table_OR_bmi_women_model2, file = "Model2_OR_bmi_women.csv", row.names = FALSE, sep=";")


#########################################################################################################
#model 3: adjusted for age, smoking, education and IMD

table_OR_bmi_model3 <- data.frame(matrix(ncol=8, nrow=0))

for (i in var_list) {
  print(i)
  
  ## Odds ratio results
  model <- glm(paste(i[[1]], "~ bmi_z + age_z + smoking + education + IMD"), data = db, family = "binomial")
  print(summary(model))
  
  # Add these coefficients, CIs, n and p-values to the table
  outcome <- i[[1]]
  varName <- "BMI (z-score)"
  coef <- names(model$coefficients)
  temp <- data.frame(outcome, varName, coef, nobs(model),
                     round(exp(model$coefficients), 2), 
                     round(exp(confint(model)[, 1]), 2),
                     round(exp(confint(model)[, 2]), 2),
                     round(coef(summary(model))[, 4], 3))
  
  colnames(temp) <- table_OR_names
  table_OR_bmi_model3 <- rbind(table_OR_bmi_model3, temp)
  table_OR_bmi_model3 <- table_OR_bmi_model3 %>%
    filter(coef != "(Intercept)")
}

write.table(table_OR_bmi_model3, file = "Model3_OR_bmi_z.csv", row.names = FALSE, sep=";")


###stratified by sex 

#Model for men

table_OR_bmi_men_model3 <- data.frame(matrix(ncol=8, nrow=0))

for (i in var_list) {
  print(i)
  
  ## Odds ratio results
  model <- glm(paste(i[[1]], "~ bmi_z + age_z + smoking + education + IMD"), data = db_men, family = "binomial")
  print(summary(model))
  
  # Add these coefficients, CIs, n and p-values to the table
  outcome <- i[[1]]
  varName <- "BMI (z-score)"
  coef <- names(model$coefficients)
  temp <- data.frame(outcome, varName, coef, nobs(model),
                     round(exp(model$coefficients), 2), 
                     round(exp(confint(model)[, 1]), 2),
                     round(exp(confint(model)[, 2]), 2),
                     round(coef(summary(model))[, 4], 3))
  
  colnames(temp) <- table_OR_names
  table_OR_bmi_men_model3 <- rbind(table_OR_bmi_men_model3, temp)
  table_OR_bmi_men_model3 <- table_OR_bmi_men_model3 %>%
    filter(coef != "(Intercept)")
}

write.table(table_OR_bmi_men_model3, file = "Model3_OR_bmi_men.csv", row.names = FALSE, sep=";")


#Model for women

table_OR_bmi_women_model3 <- data.frame(matrix(ncol=8, nrow=0))

for (i in var_list) {
  print(i)
  
  ## Odds ratio results
  model <- glm(paste(i[[1]], "~ bmi_z + age_z + smoking + education + IMD"), data = db_women, family = "binomial")
  print(summary(model))
  
  # Add these coefficients, CIs, n and p-values to the table
  outcome <- i[[1]]
  varName <- "BMI (z-score)"
  coef <- names(model$coefficients)
  temp <- data.frame(outcome, varName, coef, nobs(model),
                     round(exp(model$coefficients), 2), 
                     round(exp(confint(model)[, 1]), 2),
                     round(exp(confint(model)[, 2]), 2),
                     round(coef(summary(model))[, 4], 3))
  
  colnames(temp) <- table_OR_names
  table_OR_bmi_women_model3 <- rbind(table_OR_bmi_women_model3, temp)
  table_OR_bmi_women_model3 <- table_OR_bmi_women_model3 %>%
    filter(coef != "(Intercept)")
}

write.table(table_OR_bmi_women_model3, file = "Model3_OR_bmi_women.csv", row.names = FALSE, sep=";")
