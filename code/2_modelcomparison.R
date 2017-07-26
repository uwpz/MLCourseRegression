
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


#######################################################################################################################-
# Test an algorithm (and determine parameter grid) ----
#######################################################################################################################-

# Tibble to data frame
df.train = as.data.frame(df.samp)
 


## Possible validation controls
ctrl_cv = trainControl(method = "repeatedcv", number = 4, repeats = 1, 
                       summaryFunction = mysummary, returnResamp = "final")
ctrl_cv_fff = trainControl(method = "repeatedcv", number = 4, repeats = 1, 
                    summaryFunction = mysummary, returnResamp = "final", 
                    indexFinal = sample(1:nrow(df.train), 100))  #"Fast" final fit!!! 
set.seed(999)
l.index = list(i = sample(1:nrow(df.train), floor(0.8*nrow(df.train))))
ctrl_index = trainControl(method = "cv", number = 1, index = l.index, 
                          summaryFunction = mysummary, returnResamp = "final")
ctrl_index_fff = trainControl(method = "cv", number = 1, index = l.index, 
                          summaryFunction = mysummary, returnResamp = "final", 
                          indexFinal = sample(1:nrow(df.train), 100))  #"Fast" final fit!!! 
ctrl_none = trainControl(method = "none")



## Fit

fit = train( formula, data = df.train[c("target",predictors)], trControl = ctrl_index_fff, metric = "MAE",
             method = "glmnet", 
             tuneGrid = expand.grid(alpha = c(0,0.2,0.4,0.6,0.8,1), lambda = 2^(seq(-1, -10, -1))),
             #tuneLength = 20, 
             preProc = c("center","scale") ) 
# -> keep alpha=1 to have a full Lasso


fit = train( formula, data = df.train[c("target",predictors)], trControl = ctrl_index_fff, metric = "MAE", 
             method = "glm", 
             tuneLength = 1,
             preProc = c("center","scale") )
# -> no tuning as it is a glm


fit = train( df.train[predictors], df.train$target, trControl = ctrl_index_fff, metric = "MAE", 
             method = "rf", 
             tuneGrid = expand.grid(mtry = seq(1,11,2)), 
             #tuneLength = 2,
             ntree = 500 ) #use the Dots (...) by explicitly specifiying randomForest parameter
# -> keep to the recommended values: mtry = (length(predictors))/3

fit = train( df.train[,predictors], df.train$target, trControl = ctrl_index_fff, metric = "MAE",
             method = "gbm", 
             tuneGrid = expand.grid(n.trees = seq(100,1100,200), interaction.depth = c(3,6,9), 
                                    shrinkage = c(0.1,0.01), n.minobsinnode = c(5,10)), 
             #tuneLength = 6,
             verbose = FALSE )
# -> keep to the recommended values: interaction.depth = 6, shrinkage = 0.01, n.minobsinnode = 10


fit
plot(fit)
varImp(fit) 


skip = function() {
  # Plot tuning result with ggplot
  fit$results %>% 
    ggplot(aes(x = n.trees, y = MdAE, color = as.factor(shrinkage))) +
    geom_line(mapping = aes(linetype = as.factor(interaction.depth))) +
    #geom_errorbar(mapping = aes(ymin = RMSE - RMSESD, ymax = RMSE + RMSESD, linetype = as.factor(numLeaves))) +
    facet_grid(~n.minobsinnode) 
}

# unique(fit$results$lambda)
# plot(fit$finalModel, i.var = 1, type="response", col = twocol)




########################################################################################################################
# Compare algorithms
########################################################################################################################

## Simulation function

df.cv = as.data.frame(df.samp)

perfcomp = function(method, nsim = 5) { 
  
  result = NULL

  for (sim in 1:nsim) {

    # Hold out a k*100% set
    set.seed(sim*999)
    k = 0.2
    i.holdout = sample(1:nrow(df.cv), floor(k*nrow(df.cv)))
    df.holdout = df.cv[i.holdout,]
    df.train = df.cv[-i.holdout,]    
    
    # Control and seed for train
    set.seed(sim*1000) 
    l.index = list(i = sample(1:nrow(df.train), floor(0.8*nrow(df.train))))
    ctrl = trainControl(method = "cv", number = 1, index = l.index, 
                        summaryFunction = mysummary, returnResamp = "final")
    metric = "spearman"


    ## fit data
    fit = NULL
    if (method == "glm") {      
      fit = train( formula, data = df.train[c("target",predictors)], trControl = ctrl, metric = metric, 
                   method = "glm", 
                   tuneLength = 1,
                   preProc = c("center","scale") )
    }

    
    if (method == "glmnet") {      
      fit = train( formula, data = df.train[c("target",predictors)], trControl = ctrl, metric = metric, 
                   method = "glmnet", 
                   tuneGrid = expand.grid(alpha = 1, lambda = 2^(seq(-2, -10, -1))),
                   preProc = c("center","scale") ) 
    }     
    

    if (method == "rpart") {      
      fit = train( df.train[predictors], df.train$target, trControl = ctrl, metric = metric, 
                   method = "rpart",
                   tuneGrid = expand.grid(cp = 2^(seq(-20, -1, 2))) )
    }
    
    
    if (method == "rf") {      
      fit = train( df.train[predictors], df.train$target, trControl = ctrl, metric = metric, 
                   method = "rf", 
                   tuneGrid = expand.grid(mtry = 7), 
                   ntree = 300 )
    }
    
    
    if (method == "gbm") { 
      fit = train( df.train[predictors], df.train$target, trControl = ctrl, metric = metric,
                   method = "gbm", 
                   tuneGrid = expand.grid(n.trees = seq(100,1100,100), interaction.depth = 6, 
                                          shrinkage = 0.01, n.minobsinnode = 10), 
                   verbose = FALSE )
    }


    ## Get metrics

    # Calculate holdout performance
    yhat_holdout = predict(fit, df.holdout[predictors]) 
    perf_holdout = mysummary(data.frame(pred = yhat_holdout, obs=df.holdout$target))
    
    # Put all together
    result = rbind(result, data.frame(sim = sim, method = method, t(perf_holdout)))
  }   
  result
}

df.result = as.data.frame(c())
nsim = 20
df.result = rbind.fill(df.result, perfcomp(method = "glm", nsim = nsim))     
df.result = rbind.fill(df.result, perfcomp(method = "glmnet", nsim = nsim))   
df.result = rbind.fill(df.result, perfcomp(method = "rpart", nsim = nsim))      
df.result = rbind.fill(df.result, perfcomp(method = "rf", nsim = nsim))        
df.result = rbind.fill(df.result, perfcomp(method = "gbm", nsim = nsim))       
df.result$sim = as.factor(df.result$sim)


## Plot results
p = ggplot(df.result, aes(method, spearman)) + 
  geom_boxplot() + 
  geom_point(aes(color = sim), shape = 15) +
  geom_line(aes(group = sim, color = sim), linetype = 2, size = 0.5) +
  coord_flip() +
  labs(title = "Model Comparison") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))
p  
ggsave("./output/model_comparison.pdf", p, width = 8, height = 6)





