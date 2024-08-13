
#         Author: Pedro Salas Rojo and Paolo Brunori
#         Name of project: Solutions to Exercises
#         Data: Ecuador 2014

rm(list = ls(all.names = TRUE)) 
library(tidyverse)
library(haven)
library(dineq)
library(Metrics)

# Paste your path and plug your name
name <- "Pedro"
path <- ifelse(name == "Pedro", "C:/Users/SALASROJ/Desktop/ECUADOR/",
               "Your_path/")

# Basic setup (do not change)
source(paste0(path,"/scripts/0_functions_exante.R"))
data <- read.csv(paste0(path, "data/ecu_14.csv"), row.names = 1)  
circum <- c("Mother_Occ", "Mother_Edu", "Father_Edu", "Father_Occ")
data[circum] <- lapply(data[circum], as.factor)
model <- income ~ Mother_Occ + Father_Edu + Father_Occ + Mother_Edu

# EXERCISE 1 a) Mocking tune for minbucket ----

# Imagine that you do not want to tune the alpha/mincriterion as in the former
# example, but the minbucket of the tree.

# First, set the grid of values for different minbucket values. In this case, from 20 to 500, in steps of 25.
minbucket_try <- seq(20, 500, 25)
set.seed(1)

# Create folds to create train and test sample (see below) 
data$flag = sample(1:5, nrow(data), replace = TRUE)
rmse_test <- NA

# Loop to perform k-fold cross validation
for (kf in seq(1, 5, 1)) {
  
  # Create train and test sample
  train <- data[data$flag!=kf,]   
  test  <- data[data$flag==kf,]  
  
  # Loop for each depth you want to try
  for (minbu in minbucket_try) {
    
    # Get the tree. Note that the data comes from the "train" sample. All other
    # parameters are fixed. 
    tree <- get_tree(model = model, 
                     data = train,
                     mincri = 0.99,
                     minbu = minbu,
                     maxd = Inf)
    
    # Predict, now using newdate = test. This is the "out of sample" prediction
    test$pred <- predict(tree, type="response", newdata = test)
    
    # Estimate RMSE or any other loss function in your test sample, between the real
    # income value and the prediction.
    rmse_m <- rmse(test$income, test$pred)
    
    # Store and continue with the next fold/value
    rmse_m <- cbind(kf, minbu, rmse_m)
    rmse_test <- as.data.frame(na.omit(rbind(rmse_test, rmse_m)))
    
  }
}

# Now, rmse_test contains the RMSE obtained for each possible depth and
# fold. To tune it is only a matter of averaging across depths.

rmse_test <- rmse_test %>%
  dplyr::select(minbu, rmse_m) %>%
  group_by(minbu) %>%
  summarise_all(mean) %>%
  ungroup()

# Get the tuned value of minbucket as the one corresponding to the smallest 
# out of sample RMSE.

tuned_minbu <- rmse_test$minbu[which.min(rmse_test$rmse_m)]

# Plot minbucket values (x) and y (RMSE) to visually check the minbu value associated to the smallest
# RMSE. Select that in your tree

plot(x = rmse_test$minbu, y = rmse_test$rmse_m)
abline(v = tuned_minbu, col = "red")

tree <- get_tree(model = model, 
                 data = train,
                 mincri = 0.99,
                 minbu = tuned_minbu,
                 maxd = Inf)

# Continue with your analysis.

# EXERCISE 1 b) Mocking tune for other parameters ----

# Note that the k-fold procedure described in EXERCISE 1 a) can be used to tune any other parameter, 
# also in the random forests and in most ML algorithms. 
# Now you have the tools. If you want to tune several parameters at the same
# time, it is only a matter of "nesting" one tuning loop inside the other.
# Note that depending on the grid of values, this can be very time consuming.

# EXERCISE 4 b) Tree and sample size ----
tr_results <- NA

for(sh in seq(0.1, 1, 0.05)){
  for(i in 1:10){
    
    # Set seed such that it iterates
    set.seed(i)
    
    # Get a subsample of the data, sized "sh*100%" of the original data
    perm <- data[sample(1:nrow(data), nrow(data) * sh, replace = FALSE),]
    perm[circum] <- lapply(perm[circum], as.factor)
    
    # Run the tree with the new dataset and estimate IOp and number of types as always
    exante_tree <- get_tree(model = model, 
                            data = perm,
                            mincri = 0.99,
                            minbu = max(round(nrow(perm)*0.01)))
    
    perm$types <- predict(exante_tree, type="node")
    perm <- perm %>%
      group_by(types) %>%
      mutate(y_tilde = mean(income)) %>%
      ungroup()

    iop <- perm %>%
      summarise(types = length(unique(types)),
                gini = gini.wtd(income),
                gini_iop = gini.wtd(y_tilde),
                share = sh,
                iter = i)
    
    tr_results <- na.omit(rbind(tr_results, iop))
    
  }
}

# Since we have run several iterations, we need to average the results
tr_results <- tr_results %>% group_by(share) %>% summarise_all((mean))

# Plot the IOp level and the number of types vs the share of the sample
plot(x = tr_results$share, y = tr_results$gini_iop, 
     xlab = "Shares", ylab = "Gini IOP")
plot(x = tr_results$share, y = tr_results$types, 
     xlab = "Shares", ylab = "Types")
