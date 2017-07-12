
#######################################################################################################################-
# Initialize ----
#######################################################################################################################-

load("1_explore.rdata")
source("./code/0_init.R")


## Initialize parallel processing
Sys.getenv("NUMBER_OF_PROCESSORS") 
cl = makeCluster(4)
registerDoParallel(cl) 
# stopCluster(cl) #stop cluster


# Do not bootstrap
l.boot = NULL




#######################################################################################################################-
# Performance Plots ----
#######################################################################################################################-

## Data for interpretation
df.interpret = as.data.frame(df.samp)
df.interpret$INT = 1 #Add intercept variable (needed for partial dependence plots)



## Create predictions for holdout data (several simulation runs)
nsim = 3
yhat_holdout = c()
y_holdout = c()
df.holdout_sim = c()
for (sim in 1:nsim) {
  
  # Hold out a k*100% set
  set.seed(sim*999)
  k = 0.2  
  i.holdout = sample(1:nrow(df.interpret), floor(k*nrow(df.interpret)))
  df.holdout = df.interpret[i.holdout,]
  df.train = df.interpret[-i.holdout,]    
  
  # Model and predict holdout
  ctrl = trainControl(method = "none")

  fit = train( df.train[predictors], df.train$target, trControl = ctrl, metric = "spearman",
               method = "gbm", 
               tuneGrid = expand.grid(n.trees = 700, interaction.depth = 6, 
                                      shrinkage = 0.01, n.minobsinnode = 10), 
               verbose = FALSE )
  yhat_holdout = c(yhat_holdout, predict(fit, df.holdout[predictors])) 
  y_holdout = c(y_holdout, df.holdout$target)
  print(cor(yhat_holdout, y_holdout, method = "spearman"))
  df.holdout_sim = bind_rows(df.holdout_sim, df.holdout)
}


## Plot performance & Diagnosis
plot_performance("./output/performance.pdf", yhat_holdout, y_holdout)
plot_diagnosis("./output/diagnosis.pdf", df.holdout_sim, res = yhat_holdout - y_holdout, ylim = c(-2,2))



#######################################################################################################################-
# Importance + Partial Depedence plot ----
#######################################################################################################################-

## Get model for all data
ctrl = trainControl(method = "none")

fit.gbm = train( df.interpret[c("INT",predictors)], df.interpret$target, trControl = ctrl, metric = "spearman",
             method = "gbm", 
             tuneGrid = expand.grid(n.trees = 700, interaction.depth = 6, 
                                    shrinkage = 0.01, n.minobsinnode = 10), 
             verbose = FALSE )



## Get model for bootstrapped data and collect all models into one list
nboot = 20
add_boot = TRUE #for following plots
l.boot = foreach(i = 1:nboot, .combine = c, .packages = c("caret")) %dopar% { 
  
  # Bootstrap
  set.seed(i*9999)
  i.boot = sample(1:nrow(df.interpret), replace = TRUE)
  df.boot = df.interpret[i.boot,]
  
  # Model
  df.boot$INT = 1
  fit.gbm = train( df.boot[c("INT",predictors)], df.boot$target, trControl = ctrl, metric = "spearman",
               method = "gbm", 
               tuneGrid = expand.grid(n.trees = 700, interaction.depth = 6, 
                                      shrinkage = 0.01, n.minobsinnode = 10), 
               verbose = FALSE )
  return(list(fit.gbm))
}



## Variable Importance (works for all caret models)
plot(varImp(fit.gbm)) #Default plot

# Variable importance only for topn important variables 
topn = 10
(topn_varimp = rownames(varImp(fit.gbm)$importance)[order(varImp(fit.gbm)$importance, decreasing = TRUE)][1:topn])
plot_variableimportance("./output/variableimportance_gbm.pdf", vars = topn_varimp, l.boot = l.boot, 
                        ncols = 5, nrows = 2, w = 18, h = 12)



## Partial dependence (only for gbm)
# Default plots
plot(fit.gbm$finalModel, i.var = "age") 
plot(fit.gbm$finalModel, i.var = 7) 

# Partial dependence only for topn important variables 
topn = 10
(topn_varimp = rownames(varImp(fit.gbm)$importance)[order(varImp(fit.gbm)$importance, decreasing = TRUE)][1:topn])
plot_partialdependence("./output/partialdependance_gbm.pdf", vars = topn_varimp, l.boot = l.boot, 
                       ylim = c(1,3), ncols = 5, nrows = 2, w = 18, h = 12)



#######################################################################################################################-
# Interactions
#######################################################################################################################-


## Test for interaction (only for gbm)
(intervars = rownames(varImp(fit.gbm)$importance)[order(varImp(fit.gbm)$importance, decreasing = TRUE)][1:10])
plot_interactiontest("./output/interactiontest_gbm.pdf", vars = intervars, l.boot = NULL)


## -> Relvant interactions
inter1 = c("TT4","TSH_LOG_")
inter2 = c("sex","referral_source_OTHER_") 
inter3 = c("TT4","referral_source_OTHER_")
plot_inter("./output/inter1.pdf", inter1, ylim = c(1,3))
plot_inter("./output/inter2.pdf", inter2, ylim = c(1,3))
plot_inter("./output/inter3.pdf", inter3, ylim = c(1,3))

plot_inter_active("./output/anim", vars = inter1, duration = 3)
