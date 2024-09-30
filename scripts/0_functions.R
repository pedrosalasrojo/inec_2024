
#         Author: Pedro Salas Rojo and Paolo Brunori
#         Date: 08/2024
#         Name of project: Ex-Ante IOp Functions

library(caret)
library(partykit)
library(glmnet)
library(grDevices)
library(stringr)
library(dineq)
library(gtools)
library(resample)
library(utils)
library(data.table)
library(stats)

# Function to tune the tree ----

tune_tree <- function(data, model, cv, grid, minbu = 50) {
  
  print(paste0("This function temporarily loads the package 'party'."))
  
  library(party)
  
  set.seed(1)
  
  tr_train <- train(model,
                    data = data, 
                    method = "ctree", 
                    trControl = cv,  
                    tuneGrid = grid,
                    controls = ctree_control(testtype = "Bonferroni", 
                                             minbucket = minbu))
  
  results <- tr_train[["results"]]
  mincri <- tr_train[["bestTune"]][["mincriterion"]]    
  rmse <- round(mean(tr_train[["resample"]][["RMSE"]]), 2)
  
  detach("package:party", unload = TRUE)
  
  print(paste0("The selected alpha is: ", (1-mincri)))
  print(paste0("The selected mincriterion is: ", (mincri)))
  print(paste0("The selected alpha delivers an RMSE of: ",rmse))
  
  return(list(`mincriterion` = mincri, `RMSE` = rmse, `results` = results)) 
}

# Function to get the tree ----

get_tree <- function(data, model, mincri = 0.99, minbu = 50, maxd = Inf) {
  
  set.seed(1)
  
  tree <- partykit::ctree(model,
                          data = data, 
                          control = ctree_control(testtype = "Bonferroni", 
                                                  mincriterion = mincri,
                                                  minbucket = minbu,
                                                  maxdepth = maxd))
  return(tree)
}

# Function to plot trees ----

plot_tree <- function(data, tree, dep, norm = FALSE, font = 6) {
  
  ct_node <- as.list(tree$node)
  data$types <- predict(tree, type = "node")
  
  data$dep <- data[[dep]]

  pred <- data %>%
    group_by(types) %>%
    mutate(x = mean(x = dep)) %>%
    dplyr::select(types, x)  %>%
    summarise_all(funs(mean), na.rm = TRUE) %>%
    ungroup()
  
  a <- data %>%
    mutate(m = mean(x = dep))
  
  mean_pop <- round(mean(a$m),3)
  
  pred <- as.data.frame(pred)
  qi <- pred
  
  for (t in 1:length(qi[,1])){
    typ<-as.numeric(names(table(data$types)[t]))
    qi[t,2]<-length(data$types[data$types==typ])/length(data$types) 
  }
  
  if (norm == TRUE) {
    print(paste0("1 corresponds to the mean: ", mean_pop))
    pred$x <- pred$x/mean_pop
    dig <- 3
  } else {
    pred$x  <-pred$x
    dig <- 0
  }
  
  for(t in 1:nrow(pred)) {
    ct_node[[pred[t,1]]]$info$prediction <- as.numeric(paste(format(round(pred[t, -1], 
                                                                          digits = dig), nsmall = 2)))
    ct_node[[pred[t,1]]]$info$nobs       <- as.numeric(paste(format(round(100*qi[t, -1]  , 
                                                                          digits = 2), nsmall = 2)))
  }
  
  tree$node <- as.partynode(ct_node)
  
  plot((tree), terminal_panel=node_terminal, gp = gpar(fontsize = font),
       tp_args = list(FUN = function(node) 
         c("Rel. Type Mean",node$prediction, "Pop. Share (%)", node$nobs)))
}

# Function to tune the random forest ----

tune_forest <- function(data, model, cv, grid, 
                        ntree = 50, mincri = 0, minbu = 10) {
  
  set.seed(1)
  
  print(paste0("This function temporarily loads the package 'party'."))
  library(party)
  
  rf_train <- caret::train(model,
                           data,
                           method = "cforest",
                           trControl = cv,
                           tuneGrid = grid,
                           controls =  cforest_control(ntree = ntree,
                                                       mincriterion = mincri,
                                                       testtype = "Bonferroni",
                                                       minbucket = minbu))
  
  results <- rf_train[["results"]]
  mtry  <- rf_train[["bestTune"]][["mtry"]]        
  rmse <- round(mean(rf_train[["resample"]][["RMSE"]]), 2)  
  detach("package:party", unload = TRUE)
  
  print(paste0("The selected mtry is: ", mtry))
  print(paste0("The selected mtry delivers an RMSE of: ",rmse))
  return(list(`mtry` = mtry, `RMSE` = rmse, `results` = results))
}

# Function to get the random forest ----

get_forest <- function(data, model, ntree = 50, 
                       mincri = 0, minbu = 10,
                       mtry = "default") {
  
  set.seed(1)
  
  if(mtry == "default"){
    mvar <- ceiling((stringi::stri_count(as.vector(as.character(mod[3])), fixed = "+"))^0.5)
  } else {
    mvar <- mtry
  }
  
  forest <- partykit::cforest(model,
                              data = data,
                              ntree = ntree,
                              mtry = mvar,
                              trace = TRUE,
                              control = ctree_control(testtype = "Bonferroni",
                                                      mincriterion = mincri,
                                                      minbucket = minbu))
  
  return(forest)
}

# Function to run Shapley value decomposition "exante"  ----

shapley <- function(data, model, vars, ntree = 1,  
                    mincri = 0, minbu = 100, resample = 0.632,
                    depname = NA, rel.ineq = TRUE){
  
  time_1 <- Sys.time()
  total <- 2^length(vars)
  num <- 0
  mod_1 <- model
  mod_2 <- update(mod_1, . ~ 1)
  type = "ctree"
  wts = NA
  
  if (is.na(wts)) {
    data$weights <- 1
  } else {
    data$weights <- data[[wts]]
  }
  
  data$income <- data[[depname]]
  
  # Use tree to estimate ineq_base ----
  
  ineq_base <- NA
  
  for(iter in seq(1, ntree, 1)) {
    
    #  print(paste0("Baseline Inequality, Iter : ", iter, " of ", ntree, " trees"))
    
    set.seed(iter)
    perm <- NULL
    perm <- data[sample(1:nrow(data), nrow(data)*resample, replace = FALSE),]  # Bag (Necessary for boostrapping)
    
    if(type=="ctree"){
      
      tree <- get_tree(model = mod_1, 
                       data = perm,
                       mincri = as.numeric(mincri), 
                       minbu = as.numeric(minbu))
      
      perm$types <- predict(tree, type="node")
      
      perm <- perm %>%
        group_by(types) %>%
        mutate(y_tilde = weighted.mean(income, weights)) %>%
        ungroup()
      
      if (rel.ineq==TRUE){
        res_ineq <- gini.wtd(perm$y_tilde, perm$weights)
      } else {
        res_ineq <- modi::weighted.var(perm$y_tilde, perm$weights)
      }
      
    } else if (type=="ols") {
      
    # PLUG CODES TO ESTIMATE IOP WITH OLS      
      
    }  else if (type=="trafotree") {
      
    # PLUG CODES TO ESTIMATE IOP WITH TRANSFORMATION TREE      
      
    } else {
      stop("The method in 'type' does not exist in the Shapley Function")
    }
    
    ineq_base <- rbind(ineq_base, res_ineq)
    
  }
  
  ineq_base <- na.omit(as.data.frame(ineq_base))
  
  ineq_base <- ineq_base %>%                                    
    summarise_all(funs(mean), na.rm = TRUE) 
  
  ineq_base <- as.numeric(ineq_base)
  
  # Loop to get values for each permutation ----
  
  for(value_comb in seq(1, length(circum) , 1)) {           # Values of combinations
    
    comb <- combn(x = vars,                      # Variables to combine (circumstances in this case)
                  m = value_comb,                # Number of combinations. Iterate this to get all possible combinations
                  simplify = TRUE)               # Eliminate combinations with the same meaning (ab = ba)
    
    mat <- matrix(NA, nrow = ncol(comb), ncol = 2)  # Define matrix of results. First column: name of the combination, 
    # in the second goes the ineq value.
    
    for(roww in seq(1, ncol(comb), 1)) {         # For each one of the combinations in "comb"
      
      gap <- NA
      
      for(iter in seq(1, ntree, 1)) {
        
        if(iter == 1){
          num <- num + 1
        } else {
          num <- num 
        }
        
        print(paste0("Combination: ", num," out of: ", total, ", iter : ", iter, " of ", ntree, " trees"))
        
        set.seed(iter)
        nom <- ""
        perm <- NULL
        perm <- data[sample(1:nrow(data), nrow(data)*resample, replace = FALSE),]  # Bag (Necessary for boostrapping)
        
        for (name_var in seq (1, nrow(comb), 1)) {
          
          nom <- paste0(nom,"", comb[name_var, roww])  #Generate name with all variables involved in the combination
          
          circ <- comb[name_var, roww]                 #Select variable that is permuted
          
          perm[[circ]] <- 0                         #Substitute values by 0, so they are not used. Loop by rows if more than one circumstance has to be permuted.
          
        }
        
        # Use tree to estimate the especific contributions to IOp after permutations
        
        if(type=="ctree"){
          
          tree <- get_tree(model = mod_1, 
                           data = perm,
                           mincri = as.numeric(mincri), 
                           minbu = as.numeric(minbu))
          
          perm$types <- predict(tree, type="node")
          
          perm <- perm %>%
            group_by(types) %>%
            mutate(y_tilde = weighted.mean(income, weights)) %>%
            ungroup()
          
          if (rel.ineq==TRUE){
            
            if(value_comb != length(circum)){
              ineq_perm <- round(gini.wtd(perm$y_tilde, perm$weights), 4)
            } else {
              ineq_perm <- 0
            }
            
          } else {
            
            if(value_comb != length(circum)){
              ineq_perm <- round(modi::weighted.var(perm$y_tilde, perm$weights) , 4)
            } else {
              ineq_perm <- 0
            }
            
          }
        
                } else if (type=="ols") {
      
    # PLUG CODES TO ESTIMATE IOP WITH OLS      
      
    }  else if (type=="trafotree") {
      
    # PLUG CODES TO ESTIMATE IOP WITH TRANSFORMATION TREE      
          
        } else{
          stop("The method in 'type' does not exist in the Shapley Function")
        }        
        
        gap <- rbind(gap, ineq_perm)
        
      }
      
      gap <- as.data.frame(gap)
      
      gap <- gap %>%                                    
        summarise_all(funs(mean), na.rm = TRUE) 
      
      mat[roww, 1] <- nom                              
      mat[roww, 2] <- round(as.numeric(gap), 3)        
      
      assign(paste0("matrix_", value_comb), mat)
      
    }
  }
  
  # Get Contributions and Shapley Value  ----
  long <- length(circum) - 1                    
  dem <- factorial(length(circum))
  
  # Shapley weights
  shap_we <- matrix(NA, nrow = long, ncol = 1)   
  for(i in seq(1, long, 1)) {
    num <- factorial(i)*factorial(length(circum)-i-1) 
    shap_we[i, 1] <- round(num/dem, 4)               
  }
  
  # Define matrix of results
  
  all_mat <- NA
  for(k in seq(1,length(circum),1)){
    nam <- (paste0("matrix_",k))
    all_mat <- cbind(all_mat, nam)
  }
  all_mat <- all_mat[,-1]
  shapval <- matrix(NA, ncol = 2, nrow = length(circum))
  
  for(i in circum) {
    contributions <- NA
    
    # Contribution baseline - ineq when only 1 is iterated
    vec <- ifelse(grepl(i, matrix_1[,1]), TRUE, FALSE )   
    num <- factorial(0)*factorial(length(circum)-0-1) 
    shap_we <- num/dem
    cont_1 <- shap_we*(ineq_base - as.numeric(matrix_1[match(TRUE,vec), 2])) 
    contributions <- rbind(contributions, cont_1)
    
    # Get remaining contributions (this is, not the baseline)
    for(k in seq(2, length(circum), 1)){
      
      orig = get(paste0("matrix_", k))                #Get the matrix with the objective (origin) matrix
      prev = get(paste0("matrix_", k-1))              #Get the matrix over which we compare inequality
      
      if(k<=long) {
        l <- k-1                                     # k selects the number of circumstances, but the coalition is always one less!
        num <- factorial(l)*factorial(length(circum)-l-1) #Numerator of weights in Shapley depending on the coalition
        shap_we <- num/dem
      } else {
        num <- factorial(long)*factorial(length(circum)-long-1) # Last contribution, when only the final coalition (all permuted)
        shap_we <- num/dem                               # is considered
      }
      
      vec <- ifelse(grepl(i, orig[,1]), TRUE, FALSE)    #Which value in the objective matrix (2, 3, ..., m-1) when (2, 3, ..., m-1) 
      #variables are permuted, include the circumstance in the first column (names)
      
      for(j in seq(1, length(vec), 1)) {                    #For all values combinations (that we have stored in vec)
        
        if(vec[j]==TRUE){                                 # If the name of the circumstance is TRUE in vec
          y <- ifelse(vec[j] == TRUE, orig[j,2], NA)        # Store the correspondent value of inequality
          
          h <- ifelse(vec[j] == TRUE, orig[j,1], NA)        # Store the name of the combination
          h <- stringr::str_remove(h, i)                                 # Remove from the name in h, the letters corresponding to the circumstance
          
          vec2 <- ifelse(grepl(h, prev[,1]), TRUE, FALSE )  # Search in the matrix including one permutation less, the combination 
          # corresponding to the name stored in h
          cont <- shap_we*(as.numeric(prev[match(TRUE,vec2), 2]) - as.numeric(y)) # Estimate inequality as the inequality in the matrix
          # with one permutation less - inequality in the objective matrix. Multiply by weight.
          
          # assign(paste0("cont_", k,"_",j), cont)                  # Store the contribution. 
          contributions <- rbind(contributions, cont)
          
        } else {                                              # If vec does not contain the objective circumstance, ignore.
          NULL
        }
      }
    }
    
    contributions <- na.omit(contributions)
    contributions <- sum(contributions)
    
    row_mat <- match(i,circum)
    
    shapval[row_mat, 1] <- i                           
    shapval[row_mat, 2] <- contributions   
    
  }
  
  time_2 <- Sys.time()
  
  # Show how much time has passed between time_1 and time_2, to see how long does the loop takes
  print(round(time_2 - time_1, 2))
  
  # Check that the residual is zero (complete decomposition)
  (residual <- sum(as.numeric(shapval[,2])) + as.numeric(orig[,2]) - ineq_base)
  
  print(paste0("Residual: ", residual))
  
  # Get Marginal contribution
  (shapval)
  
  # Get Relative contribution
  rel_shapval <- shapval
  
  for(r in seq(1, nrow(shapval), 1)) {
    rel_shapval[r,2] <- 100*as.numeric(shapval[r,2])/ineq_base
  }
  
  (rel_shapval)
  
  # Check that the sum is 100
  print(paste0("Relative sum rel_shapval: ", sum(as.numeric(rel_shapval[,2]))))
  
  print(paste0("Relative sum resid: ", sum(100*as.numeric(orig[,2])/ineq_base)))
  
  print(paste0("Relative sum both: ", sum(as.numeric(rel_shapval[,2])) + 100*as.numeric(orig[,2])/ineq_base))
  
  # Get maximum importance = 100 and index accordingly the other variables
  
  maximp <- max(as.numeric(shapval[,2]))
  rel_shap_max <- shapval
  
  for(r in seq(1, nrow(shapval), 1)) {
    rel_shap_max[r,2] <- 100*as.numeric(shapval[r,2])/as.numeric(maximp)
  }
  
  (rel_shap_max)
  
  return(list(`shapval` = shapval, `rel_shapval` = rel_shapval, `rel_shap_max` = rel_shap_max)) 
  
}
