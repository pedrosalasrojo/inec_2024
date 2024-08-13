
#         Author: Pedro Salas Rojo and Paolo Brunori
#         Name of project: Estimate Ex-ante IOp with Trees and Random Forests
#         Data: Ecuador 2014

rm(list = ls(all.names = TRUE)) 
library(tidyverse)
library(haven)
library(dineq)

# Paste your path and plug your name
name <- "Pedro"
path <- ifelse(name == "Pedro", "C:/Users/SALASROJ/Desktop/ECUADOR/",
               "Your_path/")

# Get ex ante functions, necessary to estimate exante IOp
source(paste0(path,"/scripts/0_functions_exante.R"))

# Load the data 
data <- read.csv(paste0(path, "data/ecu_14.csv"), row.names = 1)  

###### Exercise 0: EXPLORE YOUR DATA ######
gini.wtd(data$income, data$weights)
summary(data$income) # Outcome is in PPP USD 2017

# Arrange Circumstances. Store names
names(data)
# The data contains more circumstances. For the example we will only use the above
circum <- c("Mother_Occ", "Mother_Edu", "Father_Edu", "Father_Occ")
# circum <- c("Sex", "Ethnicity", "Father_Occ", "Mother_Occ", 
#            "Mother_Edu", "Father_Edu", "Birth_Area")

# Store circumstances as "factors", so they are treated as unordered categories
class(data$Ethnicity)
data[circum] <- lapply(data[circum], as.factor)
class(data$Ethnicity)

# Write the model including all circumstances of interest
model <- income ~ Mother_Occ + Father_Edu + Father_Occ + Mother_Edu

# The data contains more circumstances. For the example we will only use the above
# model <- income ~ Sex + Ethnicity + Birth_Area + 
#                  Mother_Occ + Father_Edu + Father_Occ + Mother_Edu

# Tune CTree ----

# Set seed 
set.seed(1)

# Set cross-validation method and number of folds. 
# Note: See package "caret" for details

cv5 <- trainControl(method = "cv", number = 5, verboseIter = FALSE)        

# Define grid of (1-alpha) used to tune the algorithm. 
# Note: See package "caret" for details

grid_tr <- expand.grid(mincriterion = seq(0.6, 0.99, 0.005))

# Tune the algorithm with the cross validation and grid defined above

tune <- tune_tree(data = data,            # Your data
                  model = model,          # The model
                  cv = cv5,               # The cross-validation setting
                  grid = grid_tr,         # The grid for the cross-validation
                  minbu = 100)            # Minbucket (fixed or tuned)

# Get main results: tuned mincriterion, associated RMSE, and results to check
tune_mincrit <- tune[["mincriterion"]]
tune_rmse <- tune[["RMSE"]]
tune_results <- tune[["results"]]

# Plot cross-validated mincriterion and the associated RMSE
plot(x = tune_results$mincriterion, y = tune_results$RMSE)
View(tune_results)

### EXERCISE 1: The default function from "caret" tunes the parameter "mincriterion",
### which corresponds to (1-alpha). 
### Imagine you have a theoretical reason to set mincriterion to 0.99 (alpha = 0.01)
### and that you want to tune the parameter "minbucket", which corresponds to the 
### minimum number of observations in a terminal node.
### a) Write a code to tune the parameter "minbucket", holding mincriterion constant to 0.99.
### b) Write a code to tune both, mincriterion and minbucket at the same time.

# Run tree ----

# Compute tree with the mincriterion selected in the tuning. 
# Note: We need to detach the package "party".

tree <- get_tree(data = data,             # Your data
                 model = model,           # Your model
                 mincri = tune_mincrit,   # Tuned mincriterion        
                 minbu = 500,             # Tuned minbucket (minimum obs in terminal nodes)
                 maxd = Inf)              # Depth of the tree

print(tree)

# Plot tree.
plot_tree(tree,                          # Tree object
          data = data,                   # Name of the dependent variable
          dep = "income",                # Name of dependent variable
          norm = TRUE,                   # Show 1 = mean outcome
          font = 6)                      # Font size of the plot

# Predict types (predict terminal nodes)
data$types <- predict(tree, type="node", newdata = data)
length(unique(data$types))

### EXERCISE 2: In many cases you only want to explore the data, so tuning is not
### really necessary. Still, it is important to understand how the parameters of the
### tree work.
### Imagine that you do not want to tune the tree, and want to plug parameters yourself.
### a) Set mincriterion to 0.99, 0.5 and 0.2. Observe the tree. How does the structure
### change? Do the upper nodes remain the same? Why?
### b) Set minbucket to 500, 250, 50. Observe the tree. How does the structure
### change? Do the upper nodes remain the same? Why?
### c) Combine mincriterion and minbucket to get the highest possible number of types.
### How many can you get? NOTE: do not plot the tree, it may collapse R.

# Results ctree
# Use terminal nodes to estimate the average outcome by types. 
data <- data %>%
  group_by(types) %>%
  mutate(y_tilde = mean(income)) %>%
  ungroup()

# Obtain absolute IOp with Gini and MLD, as well as the number of types and the RMSE
iop <- data %>%
  summarise(name = "Ctree",
            types = length(unique(types)),
            gini = gini.wtd(income),
            mld = mld.wtd(income),
            gini_iop = gini.wtd(y_tilde),
            mld_iop = mld.wtd(y_tilde),
            RMSE = as.numeric(tune_rmse))

### EXERCISE 3: 
# a) Overfit a tree to obtain the highest possible IOp.
# b) Underfit a tree to obtain the lowest possible IOp.
# c) Run a loop with different parameters to get a plot showing IOp values
# d) Give an intuition: how are sample size and IOp associated when you use trees
# e) Give an intuition: how are the number of categories and IOp associated when
#    you use trees?
# f) Erase Mother's Education corresponding to 0, 1, 2, and 13. Reestimate the
# tree and IOp. How does the tree change? Why?

# Tune CForest ----

# Set seed for replicability
set.seed(1)

# Set cross-validation method and number of folds
# Note: See package "caret" for details
cv3 <- trainControl(method = "cv", number = 3, 
                    verboseIter = FALSE)      

# Define grid of "mtry" used to tune the algorithm
# Note: See package "caret" for details
grid_rf <- expand.grid(mtry = seq(1, length(circum), 1))

# Tune the algorithm with the cross validation and grid setting defined above.
# Now, the tuning parameter is "mtry", the number of regressors (X) chosen at
# random in every split in every tree. The default is sqrt(X)
# A better tuning implies a higher number of trees, but it takes time!

tune <- tune_forest(data = data,                 # Your data
                    model = model,               # Your model
                    cv = cv3,                    # The cross-validation setting
                    grid = grid_rf,              # The grid for the cross-validation
                    ntree = 5,                   # Number of trees to run
                    mincri = 0)                  # Mincriterion

tune_mtry <- tune[["mtry"]]
tune_rmse <- tune[["RMSE"]]
tune_results <- tune[["results"]]

plot(x = tune_results$mtry, y = tune_results$RMSE)
View(tune_results)

# Of course, you can also tune other parameters. But that is time consuming. 
# The procedure is analogous to the tuning of the tree.

# Run random forest ----

# Run the random forest with the mtry selected in the tuning.
forest <- get_forest(data = data,                 # Your data
                     model = model,               # Your model
                     ntree = 5,                   # Number of trees to run
                     mtry = tune[["mtry"]],       # Tuned mtry
                     mincri = 0)                  # Mincriterion

# Predict income. We use indexes to make it faster:

data$index<-sample(1:10, dim(data)[1], replace = TRUE)
data$y_tilde_rf_exante <- NA

for (ind in 1:10){
  data$y_tilde_rf_exante[data$index==ind]<-predict(forest, newdata=data[data$index==ind,])
  print(paste0("Random Forest prediction: ", ind," of ", 10, " folds."))
}

# Results from Random Forest. Note that in Random Forest we do not generate types
iop2 <- data %>%
  summarise(name = "Random Forest",
            types = "No Types",
            gini = gini.wtd(income),
            mld = mld.wtd(income),
            gini_iop = gini.wtd(y_tilde_rf_exante),
            mld_iop = mld.wtd(y_tilde_rf_exante),
            RMSE = as.numeric(tune_rmse))

results <- rbind(iop, iop2)

### EXERCISE 4: 
### a) Scheme a code to tune the parameter "mincri", holding "mtry" constant to sqrt(X),
### where X is the number of circumstances in your data.
### b) Advanced: explore IOp with trees and random forest with different sample sizes.
### What algorithm is more stable to changes in the sample size? What is the critical sample size
### or sample share that starts affecting IOp in one algorithm or the other? 
### c) Very Advanced: There are many parameters driving trees and random forests. Explore the 
### algorithm documentation, modify the functions yourself and explore how other parameters
### affect your estimates.

# Compare IOp results from Ctree and Random Forest
results

# Variable importance in the Random Forest ----
importance <- varimp(forest)
print(importance)
importance <- importance*100/max(importance)
print(importance)

# Sometimes one gets "negative" importance. In practice, it means that the contribution
# of that variable to the prediction problem is virtually zero.

# Variable importance in IOp (Shapley)----

# Ctree Shapley

ctreeshap <- shapley(data = data,                   # Your data
                     vars = circum,                 # Vector of circumstances 
                     model = model,                 # Your model
                     depname = "income",            # Name of dependent variable  
                     ntree = 5,                     # Number of trees
                     resample = 0.7,                # Share of observations in resample
                     mincri = 0.5)                  # Mincriterion

# Relative contribution (in %) of Circumstances to IOp
ctreeshap[["rel_shapval"]]

# Marginal contribution of Circumstances to IOp
ctreeshap[["shapval"]]

# Contributions to Circumstances to IOp (max = 100)
# This can be compared with the results from "varimp"
ctreeshap[["rel_shap_max"]]

### EXERCISE 5: 
### As said, "varimp" obtains the importance of a regressor in the prediction problem,
### this is, how good or bad a regressor is to predict the outcome in the random forest.
### The importance has no dimension, is simply a number that have sense when compared to
### the other values.
### However, the Shapley value decomposition measures the marginal and relative importance
### of the regressor to measure IOp. It does have a meaning.
### a) Explore varimp and the Shapley value Decomposition for different sets of circumstances
### How does multicolinearity affect the results? Do they always align?
