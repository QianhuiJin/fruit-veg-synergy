-----------------------------------------------------------------------------------------------------
  ## nuMoM2b-HHS: Fruit and Vegetable Synergy and Postpartum CVD Project
  
  # Work flow:
  #### Create a 4-level exposure variable representing combinations of high/low fruit and vegetable density 
  #### Create analysis data (restrict to obs with complete FFQ)
  #### Convert categorical covs as dummy variables
  #### Median/mode impute covs + missing indicators
  #### fit outcome missing model 
  #### fit outcome model
  #### fit 4 seperate propensity score models for each level of the exposure
  #### Generate predictions for AIPW
  #### Compute AIPW scores 
  #### Generate a CSV file of all the predicted values and AIPW scores
  
  ## Created by: Qianhui Jin (Jan 6, 2025)
  ## Note: this code is for outcome metabolic syndrome, which is a composite outcome of five conditions
  -----------------------------------------------------------------------------------------------------
  
#-------------------------------------------------------------------------------------------------------------------------------------
# PREPARATION
#-------------------------------------------------------------------------------------------------------------------------------------
# updateR()

install.packages("pacman")
library(pacman)
pacman::p_load(
  rio, here, skimr, tidyverse, lmtest, sandwich, broom,
  fastDummies, SuperLearner, tmle, caret, haven, dplyr
)

# Additional individual package installs (not available on CRAN via pacman)
install.packages("installr")
install.packages("VIM")
install.packages("remotes")
install.packages("randomForest")
install.packages("igraph")
install.packages("data.table")
remove.packages("data.table")
install.packages("igraph")

# Load individual libraries
library(installr)
library(VIM)
library(randomForest)
library(rlang)


remotes::install_github("tlverse/tmle3")
remotes::install_github("tlverse/tlverse")
remotes::install_github("tlverse/sl3")
library(tmle3)
library(sl3)



# Load external crosstab function
source("http://pcwww.liv.ac.uk/~william/R/crosstab.r")

# Show R version info
R.Version()

#--------------------------------------------------------------------------------------------------------------------
# LOAD IN THE UNIMPUTED DATA AND FINALIZE SAMPLE SIZE ---------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------
#setwd("")
hhs<-read_dta("Dataset/numom+hhs_diet_unimputed_v3.dta")
dim(hhs)
table(hhs$hhs_use)
# Restricted to obs with complete FFQ data, n= 3776


# Look at outcome var: metabolic syndrome at v5 
crosstab(hhs, row.vars = "metsx_v5")
sum(is.na(hhs$metsx_v5))
# 64 missing


# Create dichotomized fruit and veg intake denssity vars (cutoff points based on Dietary Guidelines)
summary(hhs$f_totdens_v1)
summary(hhs$v_totdens_v1)

hhs$v_totdens_v1_dich <- ifelse(hhs$v_totdens_v1 >= 1.25, 1, 0) 
hhs$f_totdens_v1_dich  <- ifelse(hhs$f_totdens_v1 >= 1, 1, 0) 

crosstab(hhs, row.vars = "v_totdens_v1_dich")
crosstab(hhs, row.vars = "f_totdens_v1_dich")


# Create a four level fruit and veg combination exposure var
hhs$fv_dens_cat4 <- with(hhs, 
                         ifelse(f_totdens_v1_dich == 0 & v_totdens_v1_dich == 0, 1,  # Low Fruit, Low Veg
                                ifelse(f_totdens_v1_dich == 0 & v_totdens_v1_dich == 1, 2,  # Low Fruit, High Veg
                                       ifelse(f_totdens_v1_dich == 1 & v_totdens_v1_dich == 0, 3,  # High Fruit, Low Veg
                                              ifelse(f_totdens_v1_dich == 1 & v_totdens_v1_dich == 1, 4,  NA )))))  # High Fruit, High Veg
crosstab(hhs, row.vars = "fv_dens_cat4")
crosstab(hhs, row.vars = "fv_dens_cat4", col.vars = "metsx_v5", type = "row.pct")


# Create a binary indicator for each level of exposure
hhs$fv_cat1 <- ifelse(hhs$fv_dens_cat4 == 1, 1, 0)  #  Low Fruit, Low Veg v.s. all others
hhs$fv_cat2 <- ifelse(hhs$fv_dens_cat4 == 2, 1, 0)  #  Low Fruit, High Veg v.s. all others
hhs$fv_cat3 <- ifelse(hhs$fv_dens_cat4 == 3, 1, 0)  #  High Fruit, Low Veg v.s. all others
hhs$fv_cat4 <- ifelse(hhs$fv_dens_cat4 == 4, 1, 0)  #  High Fruit, High Veg v.s. all others


# Clean covs
## Make study sites into numbers
hhs$publicsite_num <- as.numeric(as.factor(hhs$publicsite))
table(hhs$publicsite_num )

## Combine preexisting diabetes and hypertension
table(hhs$prehtn,useNA = "ifany")
table(hhs$prediab,useNA = "ifany")

hhs$pre_combined <- ifelse(hhs$prehtn == 1 | hhs$prediab == 1, 1, 0)
crosstab(hhs, row.vars = "pre_combined")


# Select variables
exposure <- c("fv_dens_cat4","fv_cat1", "fv_cat2","fv_cat3","fv_cat4")
outcomes <- c("metsx_v5")
covariates <- c("d_totdens_v1","p_totdens_v1", "p_seaplantdens_v1", "g_whldens_v1","g_nwhldens_v1","fatratio_v1","sodium_dens_v1", "pct_addsug_v1","pct_satfat_v1", # HEI components
                "smokerpre","married","insurpub","pre_combined", "momeduc4", 
                "momrace4", "momage","bmiprepreg","accult3", "sleepsat3","epds_tot_v1", 
                "stress_tot_v1",  "anx_tot", "pregplanned","pa_totmetwk_new_v1",
                "puqe_tot","gravcat","povperc_v1", "walk_nat_v1",  "v1_alcbinge_3monthsprior", 
                "publicsite_num",  "employcat","adi_nat_v1","artcat" ) 


# Create analytical dataset
a_final <- hhs %>%
  dplyr::select(all_of(outcomes), all_of(exposure), all_of(covariates))
dim(a_final)
#3776      39
sum(!complete.cases(a_final[, covariates]))


#--------------------------------------------------------------------------------------------------------------------
# CREATING DUMMY VARIABLES  -----------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------
# All continuous should be numeric and all categorical vars should be dummy
## Make a list of variables into factors at all one time
factor_names <- c(  "smokerpre","married","insurpub","pre_combined", "momeduc4", 
                    "momrace4", "accult3", "sleepsat3", 
                    "pregplanned",
                    "gravcat",  "v1_alcbinge_3monthsprior", 
                    "publicsite_num",  "employcat","artcat")
a_final[,factor_names] <- lapply(a_final[,factor_names] , factor)

a_final <- dummy_cols(a_final, ignore_na = TRUE, remove_selected_columns = TRUE,
                      remove_first_dummy = TRUE)
names(a_final)
skim(a_final)


#-------------------------------------------------------------------------------------------------------------------------------------
# MEIDAN/MODE IMPUTATION + MISISNG INDICATORS
#-------------------------------------------------------------------------------------------------------------------------------------
nodes_n <- list(W = c(  "d_totdens_v1","p_totdens_v1", "p_seaplantdens_v1", "g_whldens_v1","g_nwhldens_v1","fatratio_v1","sodium_dens_v1", "pct_addsug_v1","pct_satfat_v1", 
                        "momage", "pa_totmetwk_new_v1","bmiprepreg","epds_tot_v1", "stress_tot_v1", "anx_tot","puqe_tot","povperc_v1", "walk_nat_v1","adi_nat_v1",
                        "insurpub_1",  "married_1","pre_combined_1", "pregplanned_1", "v1_alcbinge_3monthsprior_1","smokerpre_1", 
                        "gravcat_2"  ,"gravcat_3"  ,
                        "momrace4_2", "momrace4_3", "momrace4_4", 
                        "momeduc4_2", "momeduc4_3", "momeduc4_4",
                        "accult3_2","accult3_3", "accult3_4",  
                        "sleepsat3_2", "sleepsat3_3",
                        "employcat_2", "employcat_3",
                        "artcat_2", "artcat_3",
                        "publicsite_num_2", "publicsite_num_3", "publicsite_num_4", "publicsite_num_5", "publicsite_num_6", 
                        "publicsite_num_7", "publicsite_num_8"),  
                A = c("fv_dens_cat4","fv_cat1", "fv_cat2","fv_cat3","fv_cat4"),
                Y = "metsx_v5")

# Check how many variables are in each group
sapply(nodes_n, length)



# Impute missing data 
processed_n <- process_missing(a_final, 
                               nodes_n, 
                               complete_nodes = "A") #Drop any rows that have missing exposure data


# Extract processed data and add ID for cross-validation
processed_hhs <- processed_n$data
processed_hhs$id <- 1:nrow(processed_hhs) 

# Shows variables removed during process missing
setdiff(names(a_final), names(processed_hhs))  

# Check
table(processed_hhs$fv_dens_cat4, processed_hhs$metsx_v5)
table(processed_hhs$fv_dens_cat4)
sum(is.na(processed_hhs$metsx_v5))
sum(is.na(processed_hhs))
str(processed_hhs)  
names(processed_hhs)

# Check distribution of other delta vars to avoid 0 variance 
table(processed_hhs$ delta_sleepsat3_2 ) 
table(processed_hhs$delta_sleepsat3_3) 
table(processed_hhs$delta_employcat_2)
table(processed_hhs$delta_employcat_3)
table(processed_hhs$delta_smokerpre_1) #only 4
table(processed_hhs$delta_v1_alcbinge_3monthsprior_1) 
table(processed_hhs$delta_pregplanned_1) #only 2
table(processed_hhs$delta_bmiprepreg) #only 1
table(processed_hhs$delta_epds_tot_v1) 
table(processed_hhs$delta_stress_tot_v1) #only 8
table(processed_hhs$delta_anx_tot) 
table(processed_hhs$delta_povperc_v1) 
table(processed_hhs$delta_walk_nat_v1) 
table(processed_hhs$delta_adi_nat_v1) 

# Check correlation
cor_matrix <- cor(processed_hhs, use = "pairwise.complete.obs")
cor_matrix
high_corr_vars <- findCorrelation(cor_matrix, cutoff = 0.99, names = TRUE)
print(high_corr_vars)

# Removed near zero variance vars
delta_vars_to_remove <- c("delta_smokerpre_1", "delta_pregplanned_1", "delta_bmiprepreg" ,"delta_stress_tot_v1")
processed_hhs[, (delta_vars_to_remove) := NULL]  


# Save
processed_hhs_co <- processed_hhs[delta_metsx_v5 == 1, ] %>%
  select(-c(delta_metsx_v5))

names(processed_hhs_co)

#-------------------------------------------------------------------------------------------------------------------------------------
# SETTING UP OUR DATA TO FEED INTO THE SUPER LEARNER
#-------------------------------------------------------------------------------------------------------------------------------------
# Define exposure and covs for the PS model (full data)
exp_indicators  <- processed_hhs %>% select(fv_cat1, fv_cat2,fv_cat3,fv_cat4)
exposure <- processed_hhs$fv_dens_cat4
covs_ps <- processed_hhs %>% select(-id, -metsx_v5, -delta_metsx_v5, -fv_dens_cat4, -fv_cat1, -fv_cat2,-fv_cat3,-fv_cat4) 
names(covs_ps)
nrow(covs_ps)
#3776

# Define outcome and covs for the outcome missingness model (full data)
outcome_delta <- processed_hhs$delta_metsx_v5
covs_theta <- processed_hhs %>% select(-id, -metsx_v5, -delta_metsx_v5,-fv_cat1, -fv_cat2,-fv_cat3,-fv_cat4)
names(covs_theta)
nrow(covs_theta)
# 3776

# Define outcome and covs for the outcome model (complete outcome data)
outcome_co <- processed_hhs_co$metsx_v5
covs_mu <-  processed_hhs_co %>%
  select( -id, -metsx_v5,-fv_cat1, -fv_cat2,-fv_cat3,-fv_cat4)
names(covs_mu)
nrow(covs_mu)
# 3712

#-------------------------------------------------------------------------------------------------------------------------------------
# CREATE INDEX
#-------------------------------------------------------------------------------------------------------------------------------------
# Define folds 
n <- nrow(processed_hhs)
num_folds <- 10 
folds <- sort(seq(n) %% num_folds) + 1
fold_dat <- tibble(id = processed_hhs$id, folds) ## This id is the row numbers

## Make sure the same observations stay in the same folds 
# NOTE: id in the subset dataset were the row numbers of the full dataset 
# that was created before subset the dataset with only complete 
# outcome in 01_program
fold_dat_mu <-processed_hhs_co %>% select(id) %>% left_join(fold_dat, by = "id")

## Use the row numbers of the current dataset to replace the old IDs in 
# creating fold index list. This will create sequential numbers for the 
# new IDs. We can ensure that during the cross-validation process, there
# won't be row numbers are outside of the range of subset dataset
fold_dat_mu$id_n <- 1:nrow(fold_dat_mu) 

# Create fold index lists
fold_index <- split(fold_dat$id, fold_dat$folds)
fold_index_mu <- split(fold_dat_mu$id_n, fold_dat_mu$folds)

# Check
table(folds) 

length(unlist(fold_index)) == nrow(processed_hhs)
length(unlist(fold_index_mu)) == nrow(processed_hhs_co)

any(duplicated(unlist(fold_index)))  
any(duplicated(unlist(fold_index_mu)))  

all(unlist(fold_index) %in% seq_len(nrow(processed_hhs)))  
all(unlist(fold_index_mu) %in% seq_len(nrow(processed_hhs_co)))  

sapply(fold_index, length)
sapply(fold_index_mu, length)

head(fold_dat)
head(fold_dat_mu)

nrow(fold_dat_mu) == nrow(processed_hhs_co)  
table(fold_dat_mu$folds)  #  balanced
range(fold_dat_mu$id_n)  # 1 to nrow(fold_dat_mu)

length(fold_index) == num_folds  
length(fold_index_mu) == num_folds 

sapply(fold_index, length)  
sapply(fold_index_mu, length) 

#-------------------------------------------------------------------------------------------------------------------------------------
# SETTING UP THE V-FOLD CROSS-VALIDATION SETTINGS 
#-------------------------------------------------------------------------------------------------------------------------------------
# CREATE LEARNER 
.SL.require <- function(package, 
                        message = 
                          paste('loading required package (', package, ') failed', 
                                sep = '')) {
  if(!requireNamespace(package, quietly = FALSE)) {
    stop(message, call. = FALSE)
  }
  invisible(TRUE)
}

screen.glmnet1 <- function (Y, X, family, alpha = 1, minscreen = 10, nfolds = 10,
                            nlambda = 100, ...)
{
  .SL.require("glmnet")
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + ., X)
  }
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, lambda = NULL, type.measure = "deviance",
                             nfolds = nfolds, family = family$family, alpha = alpha,
                             nlambda = nlambda)
  whichVariable <- (as.numeric(coef(fitCV$glmnet.fit, s = fitCV$lambda.min))[-1] !=
                      0)
  if (sum(whichVariable) < minscreen) {
    warning("fewer than minscreen variables passed the glmnet screen, 
            increased lambda to allow minscreen variables")
    sumCoef <- apply(as.matrix(fitCV$glmnet.fit$beta), 2,
                     function(x) sum((x != 0)))
    newCut <- which.max(sumCoef >= minscreen)
    whichVariable <- (as.matrix(fitCV$glmnet.fit$beta)[,
                                                       newCut] != 0)
  }
  return(whichVariable)
}



# Base learners (screen.glmnet1 was removed becasue of warnings: In (function (Y, X, family, alpha = 1, minscreen = 10,  ... :
#fewer than minscreen variables passed the glmnet screen, 
#increased lambda to allow minscreen variables)

sl.lib <- c(list(c("SL.glm", "screen.glmnet1"),
                 c("SL.ranger" ),
                 c("SL.glmnet" ),
                 c("SL.earth"))) 


#-------------------------------------------------------------------------------------------------------------------------------------
# MODELS AND PREIDCTIONS
#-------------------------------------------------------------------------------------------------------------------------------------
# Outcome missingness model
fit_theta <- CV.SuperLearner(Y = outcome_delta,
                             X = covs_theta,
                             method = "method.NNLS", 
                             family = binomial,
                             SL.library = sl.lib,
                             cvControl = list(V = 10, validRows = fold_index),
                             control = list(saveCVFitLibrary = F),
                             parallel = "seq",
                             verbose = T
)

summary(fit_theta)
coef(fit_theta) #GLMnet dominating 
warnings() 





# Outcome model
fit_mu <- CV.SuperLearner(Y = outcome_co,
                          X = covs_mu, 
                          method = "method.NNLS", 
                          family = binomial,
                          SL.library = sl.lib,
                          cvControl = list(V = 10, validRows = fold_index_mu),
                          control = list(saveCVFitLibrary = F), 
                          parallel = "seq",
                          verbose = T)

summary(fit_mu)
coef(fit_mu) 
warnings() 



# Exposure models
fit_pi1 <- CV.SuperLearner(Y = exp_indicators$fv_cat1,
                           X = covs_ps,
                           method = "method.NNLS",
                           family = binomial(),
                           SL.library = sl.lib,
                           cvControl = list(V = 10, validRows = fold_index),
                           control = list(saveCVFitLibrary = F), 
                           parallel = "seq",
                           verbose = TRUE)

summary(fit_pi1)
coef(fit_pi1) #GLM + GLMNet + Earth are contributing the most
warnings()



fit_pi2 <- CV.SuperLearner(Y = exp_indicators$fv_cat2,
                           X = covs_ps,
                           method = "method.NNLS",
                           family = binomial(),
                           SL.library = sl.lib,
                           cvControl = list(V = 10, validRows = fold_index),
                           control = list(saveCVFitLibrary = F), 
                           parallel = "seq",
                           verbose = TRUE)

summary(fit_pi2)
coef(fit_pi2)
warnings()



fit_pi3 <- CV.SuperLearner(Y = exp_indicators$fv_cat3,
                           X = covs_ps,
                           method = "method.NNLS",
                           family = binomial(),
                           SL.library = sl.lib,
                           cvControl = list(V = 10, validRows = fold_index),
                           control = list(saveCVFitLibrary = F), 
                           parallel = "seq",
                           verbose = TRUE)

summary(fit_pi3)
coef(fit_pi3)
warnings()


fit_pi4 <- CV.SuperLearner(Y = exp_indicators$fv_cat4,
                           X = covs_ps,
                           method = "method.NNLS",
                           family = binomial(),
                           SL.library = sl.lib,
                           cvControl = list(V = 10, validRows = fold_index),
                           control = list(saveCVFitLibrary = F), 
                           parallel = "seq",
                           verbose = TRUE)

summary(fit_pi4)
coef(fit_pi4)


# Generate and Evaluate Predictions Needed to Compute Effects
# Create bounding function
bound_func <- function(probs, lower_bound, upper_bound){
  probs <- if_else(probs < lower_bound, lower_bound, probs)
  probs <- if_else(probs > upper_bound, upper_bound, probs)
  return(probs)
}



# Extracting and bounding predicted probabilities pi_hat from 4 models
# Estimated probability of being in each exposure group, given covariates 
pi1 <- fit_pi1$SL.predict[outcome_delta == 1]
pi2 <- fit_pi2$SL.predict[outcome_delta == 1]
pi3 <- fit_pi3$SL.predict[outcome_delta == 1]
pi4 <- fit_pi4$SL.predict[outcome_delta == 1]

summary(pi1)
summary(pi2)
summary(pi3)
summary(pi4)



pi1_hat_t <- bound_func(pi1, .025, .975)
pi2_hat_t <- bound_func(pi2, .025, .975)
pi3_hat_t <- bound_func(pi3, .025, .975)
pi4_hat_t <- bound_func(pi4, .025, .975)

summary(pi1_hat_t)
summary(pi2_hat_t)
summary(pi3_hat_t)
summary(pi4_hat_t)



# Extracting and bounding predicted probabilities from outcome misssing model
# Estimated probabilities of having an observed outcome
theta_hat <- fit_theta$SL.predict[outcome_delta == 1]
summary(theta_hat)
theta_hat_t <- bound_func(theta_hat, .0125, .975)
summary(theta_hat_t)



# Extracting and bounding predicted probabilities from outcome model
mu_hat <- fit_mu$SL.predict
summary(mu_hat) 
mu_hat <- bound_func(mu_hat, .0125, .975)
summary(mu_hat)




# Make predictions for each level of the exposure
# The predicted MetS probability for everyone (with observed outcomes), if they were all assigned exposure level 
mu_hat1 <- NULL
for(i in 1:10){
  mu_hat1 <- rbind(mu_hat1, 
                   predict(fit_mu$AllSL[[i]],
                           newdata = 
                             base::transform(covs_mu[fold_index_mu[[i]], ], 
                                             fv_dens_cat4 = 1), 
                           onlySL=T)$pred)
}
mu_hat1_t <- bound_func(mu_hat1, .0125, .975)
summary(mu_hat1_t)



mu_hat2 <- NULL
for(i in 1:10){
  mu_hat2 <- rbind(mu_hat2, 
                   predict(fit_mu$AllSL[[i]],
                           newdata = 
                             base::transform(covs_mu[fold_index_mu[[i]], ], 
                                             fv_dens_cat4 = 2), 
                           onlySL=T)$pred)
}
mu_hat2_t <- bound_func(mu_hat2, .0125, .975)
summary(mu_hat2_t)


mu_hat3 <- NULL
for(i in 1:10){
  mu_hat3 <- rbind(mu_hat3, 
                   predict(fit_mu$AllSL[[i]],
                           newdata = 
                             base::transform(covs_mu[fold_index_mu[[i]], ], 
                                             fv_dens_cat4 = 3), 
                           onlySL=T)$pred)
}
mu_hat3_t <- bound_func(mu_hat3, .0125, .975)
summary(mu_hat3_t)


mu_hat4 <- NULL
for(i in 1:10){
  mu_hat4 <- rbind(mu_hat4, 
                   predict(fit_mu$AllSL[[i]],
                           newdata = 
                             base::transform(covs_mu[fold_index_mu[[i]], ], 
                                             fv_dens_cat4 = 4), 
                           onlySL=T)$pred)
}
mu_hat4_t <- bound_func(mu_hat4, .0125, .975)
summary(mu_hat4_t)

#-------------------------------------------------------------------------------------------------------------------------------------
# OBTAINING THE AIPW - 'Augmented Inverse Probability Weighting' - VALUES
#-------------------------------------------------------------------------------------------------------------------------------------
# Compute AIPW estimate for level 1 -4
aipw_mu1_fv1 <- (covs_mu$fv_dens_cat4 == 1) / (pi1_hat_t * theta_hat_t) * (outcome_co - mu_hat1_t) + mu_hat1_t
aipw_mu1_fv2 <- (covs_mu$fv_dens_cat4 == 2) / (pi2_hat_t * theta_hat_t) * (outcome_co - mu_hat2_t) + mu_hat2_t
aipw_mu1_fv3 <- (covs_mu$fv_dens_cat4 == 3) / (pi3_hat_t * theta_hat_t) * (outcome_co - mu_hat3_t) + mu_hat3_t
aipw_mu1_fv4 <- (covs_mu$fv_dens_cat4 == 4) / (pi4_hat_t * theta_hat_t) * (outcome_co - mu_hat4_t) + mu_hat4_t

mean(aipw_mu1_fv1)
mean(aipw_mu1_fv2)
mean(aipw_mu1_fv3)
mean(aipw_mu1_fv4)


# Obtain Risk Differences (ATEs)
aipw_res_4_1 <- data.frame(RiskDifference_4_1 = mean(aipw_mu1_fv4 - aipw_mu1_fv1),
                           StandardError = sd(aipw_mu1_fv4 - aipw_mu1_fv1)/sqrt(nrow(covs_mu)),
                           LCL = mean(aipw_mu1_fv4 - aipw_mu1_fv1) - 1.96*sd(aipw_mu1_fv4 - aipw_mu1_fv1)/sqrt(nrow(covs_mu)),
                           UCL = mean(aipw_mu1_fv4 - aipw_mu1_fv1) + 1.96*sd(aipw_mu1_fv4 - aipw_mu1_fv1)/sqrt(nrow(covs_mu)))

round(aipw_res_4_1, 3)



aipw_res_3_1 <- data.frame(RiskDifference_3_1 = mean(aipw_mu1_fv3 - aipw_mu1_fv1),
                           StandardError = sd(aipw_mu1_fv3 - aipw_mu1_fv1)/sqrt(nrow(covs_mu)),
                           LCL = mean(aipw_mu1_fv3 - aipw_mu1_fv1) - 1.96*sd(aipw_mu1_fv3 - aipw_mu1_fv1)/sqrt(nrow(covs_mu)),
                           UCL = mean(aipw_mu1_fv3 - aipw_mu1_fv1) + 1.96*sd(aipw_mu1_fv3 - aipw_mu1_fv1)/sqrt(nrow(covs_mu)))

round(aipw_res_3_1, 3)



aipw_res_2_1 <- data.frame(RiskDifference_2_1 = mean(aipw_mu1_fv2 - aipw_mu1_fv1),
                           StandardError = sd(aipw_mu1_fv2 - aipw_mu1_fv1)/sqrt(nrow(covs_mu)),
                           LCL = mean(aipw_mu1_fv2 - aipw_mu1_fv1) - 1.96*sd(aipw_mu1_fv2 - aipw_mu1_fv1)/sqrt(nrow(covs_mu)),
                           UCL = mean(aipw_mu1_fv2 - aipw_mu1_fv1) + 1.96*sd(aipw_mu1_fv2 - aipw_mu1_fv1)/sqrt(nrow(covs_mu)))

round(aipw_res_2_1, 3)



# Calculate deviation from additive effect (observed joint effect - expected additive effect)
synergy_df <- data.frame(
  synergy = aipw_res_4_1$RiskDifference_4_1 - (aipw_res_2_1$RiskDifference_2_1 + aipw_res_3_1$RiskDifference_3_1),
  se = sqrt(
    aipw_res_4_1$StandardError^2 + 
      aipw_res_2_1$StandardError^2 + 
      aipw_res_3_1$StandardError^2
  )
)
# 95% CI
synergy_df$LCL <- synergy_df$synergy - 1.96 * synergy_df$se
synergy_df$UCL <- synergy_df$synergy + 1.96 * synergy_df$se
round(synergy_df, 4)















