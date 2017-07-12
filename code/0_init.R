
#######################################################################################################################-
# Libraries + Parallel Processing Start ----
#######################################################################################################################-

skip = function() {
  install.packages(c("corrplot","vcd","doSNOW","lattice","ggplot2","plyr","reshape2","ROCR","Hmisc",
                     "e1071","pROC","caret","glmnet","randomForest","gbm","rpart"))
}


library(doParallel)
library(plyr)
library(caret)
library(tidyverse)
library(forcats)

library(corrplot)
library(vcd)
library(gridExtra)
library(Hmisc)

library(d3heatmap)
library(htmlwidgets)
library(ROCR)
library(stringr)
library(rgl)

# 
# library(lattice)
# library(ggplot2)
# library(plyr)
# library(reshape2)
# 
# 
# 
# library(pROC)
# library(glmnet)
# library(randomForest)
# library(gbm)
# library(rpart)


# 
# library(pROC)
# library(mboost)
# library(kernlab)
# library(xgboost)
# library(latticeExtra)
# library(partykit)
# library(ICEbox)
# library(gridExtra)



#######################################################################################################################-
# Parameter ----
#######################################################################################################################-

plotloc = ""

theme_my = theme_bw() +  theme(plot.title = element_text(hjust = 0.5))

twocol = c("blue","red")

col3d = colorRampPalette(c("cyan", "white", "magenta"))(100)
colhex = colorRampPalette(c("white", "blue", "yellow", "red"))(100)

manycol = c('#00FF00','#0000FF','#FF0000','#01FFFE','#FFA6FE','#FFDB66','#006401','#010067','#95003A',
            '#007DB5','#FF00F6','#FFEEE8','#774D00','#90FB92','#0076FF','#D5FF00','#FF937E','#6A826C','#FF029D',
            '#FE8900','#7A4782','#7E2DD2','#85A900','#FF0056','#A42400','#00AE7E','#683D3B','#BDC6FF','#263400',
            '#BDD393','#00B917','#9E008E','#001544','#C28C9F','#FF74A3','#01D0FF','#004754','#E56FFE','#788231',
            '#0E4CA1','#91D0CB','#BE9970','#968AE8','#BB8800','#43002C','#DEFF74','#00FFC6','#FFE502','#620E00',
            '#008F9C','#98FF52','#7544B1','#B500FF','#00FF78','#FF6E41','#005F39','#6B6882','#5FAD4E','#A75740',
            '#A5FFD2','#FFB167','#009BFF','#E85EBE')




#######################################################################################################################-
# My Functions ----
#######################################################################################################################-

## Custom summary function for caret training
mysummary = function(data, lev = NULL, model = NULL)
{
  #browser()
  concord = function(obs, pred, n=100000) {
    i.samp1 = sample(1:length(obs), n, replace = TRUE)
    i.samp2 = sample(1:length(obs), n, replace = TRUE)
    obs1 = obs[i.samp1]
    obs2 = obs[i.samp2]
    pred1 = pred[i.samp1]
    pred2 = pred[i.samp2]
    sum((obs1 > obs2) * (pred1 > pred2) + (obs1 < obs2) * (pred1 < pred2) + 0.5*(obs1 == obs2)) / sum(obs1 != obs2)
  }
  if (is.character(data$obs)) 
    data$obs = factor(data$obs, levels = lev)
  
  isNA = is.na(data[, "pred"])
  
  pred = data[, "pred"][!isNA]
  obs = data[, "obs"][!isNA]
  
  spear = cor(pred, obs, method = "spearman")
  pear = cor(pred, obs, method = "pearson")
  AUC = concord(pred, obs)
  
  out = c(spear, pear, AUC)
  names(out) = c("spearman","pearson","AUC")
  out
}


## Calculate probabilty on all data from probabilt from sample data and the corresponding (prior) base probabilities 
prob_samp2full = function(p_sample, b_sample, b_all) {
  p_all = b_all * ((p_sample - p_sample*b_sample) / 
                   (b_sample - p_sample*b_sample + b_all*p_sample - b_sample*b_all))
  p_all
}



## Workaround for ggsave and marrangeGrob not to create first page blank
grid.draw.arrangelist <- function(x, ...) {
  for (ii in seq_along(x)) {
    if (ii > 1) grid.newpage()  # skips grid.newpage() call the first time around
    grid.draw(x[[ii]])
  }
}



## Plot distribution of metric variables per stratum
plot_distr_metr = function(outpdf, df = df.plot, vars = metr, misspct, nbins = 50, color = colhex, ylim = NULL,
                           ncols = 3, nrows = 2, w = 12, h = 8) {
  # Univariate variable importance
  varimp = sqrt(filterVarImp(df[vars], df$target, nonpara = TRUE)[[1]])
  names(varimp) = vars

  # Loop over vars
  plots = map(vars, ~ {
    #. = vars[2]
    print(.)
    
    # Scatterplot
    p = ggplot(data = df, aes_string(x = ., y = "target")) +
      geom_hex() + 
      scale_fill_gradientn(colours = color) +
      geom_smooth(color = "black", method = "gam", level = 0.95, size = 0.5) +
      labs(title = paste0(.," (Imp.: ", round(varimp[.],2),")"),
            x = paste0(.," (NA: ", misspct[.] * 100,"%)"))
    if (length(ylim)) p = p + ylim(ylim)
    
    
    # Inner Histogram
    p.inner = ggplot(data = df, aes_string(x = .)) +
      geom_histogram(aes(y = ..density..), bins = nbins, position = "identity", color = "lightgrey") +
      geom_density(color = "black") +
      theme_void()
    
    # Get underlying data for max of y-value and range of x-value
    tmp = ggplot_build(p)
    yrange = tmp$layout$panel_ranges[[1]]$y.range
    
    # Put all together
    p = p + 
      scale_y_continuous(limits = c(yrange[1] - 0.2*(yrange[2] - yrange[1]), ifelse(length(ylim), ylim[2], NA))) +
      theme(plot.title = element_text(hjust = 0.5)) +
      annotation_custom(ggplotGrob(p.inner), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = yrange[1]) 
    #if (. != vars[1]) p = p + theme(legend.position = "none") 
    p 
  })
  ggsave(outpdf, marrangeGrob(plots, ncol = ncols, nrow = nrows, top = NULL), width = w, height = h)
}



## Plot distribution of nominal variables per stratum
plot_distr_nomi = function(outpdf, df, vars = nomi, ylim = NULL,
                           ncols = 4, nrows = 2, w = 12, h = 8) {
  # Univariate variable importance
  varimp = sqrt(filterVarImp(df[vars], df$target, nonpara = TRUE)[[1]])
  names(varimp) = vars
  
  # Loop over vars
  plots = map(vars, ~ {
    #. = vars[1]
    print(.)
    
    # Main Boxplot
    p = ggplot(df, aes_string(x = ., y = "target")) +
      #geom_violin(scale = "count", fill = "lightgrey") +
      #geom_boxplot(width=0.1, fill="white") +
      geom_boxplot(varwidth = TRUE) +
      coord_flip() +
      scale_x_discrete(labels = paste0(levels(df[[.]]), " (", round(100 * table(df[[.]])/nrow(df), 1), "%)")) +
      labs(title = paste0(.," (Imp.: ", round(varimp[.],2),")"), x = "") +
      theme_my +
      theme(legend.position = "none") 
    if (length(ylim)) p = p + ylim(ylim)

    # Get underlying data for max of y-value and range of x-value
    tmp = ggplot_build(p)
    yrange = tmp$layout$panel_ranges[[1]]$x.range
      
    # Inner Barplot
    p.inner = ggplot(data = df, aes_string(x = .)) +
      geom_bar(fill = "grey", colour = "black", width = 0.9) +
      coord_flip() +
      theme_void()
    
    # Put all together
    p = p + 
      scale_y_continuous(limits = c(yrange[1] - 0.2*(yrange[2] - yrange[1]), ifelse(length(ylim), ylim[2], NA))) +
      theme_my +
      annotation_custom(ggplotGrob(p.inner), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = yrange[1]) 
    p   
  })
  ggsave(outpdf, marrangeGrob(plots, ncol = ncols, nrow = nrows, top = NULL), width = w, height = h)
}

## Plot correlation of metric variables
plot_corr_metr <- function(outpdf, df, vars=metr, misspct, method = "Spearman", 
                           w = 8, h = 8) {

  # Correlation matrix
  if (tolower(method) %in% c("spearman","pearson")) {
    m.corr = abs(cor(df[vars], method = tolower(method), use = "pairwise.complete.obs"))
  } 
  if (tolower(method) ==  "hoeffding") {
    m.corr = hoeffd(as.matrix(df[vars]))$D
  }
  
  # Adapt labels
  rownames(m.corr) = colnames(m.corr) =
    paste0(rownames(m.corr)," (", round(100*misspct[rownames(m.corr)],1),"% imp.)")
  
  # Put in clustered order
  m.corr[which(is.na(m.corr))] = 0 #set NAs to 0
  ord = corrMatOrder(m.corr , order = "hclust")
  m.corr = m.corr[ord, ord]
  d3heatmap(m.corr, colors = "Blues", sym = TRUE, xaxis_font_size = paste0(h, "pt"), xaxis_height = 120, yaxis_width = 160)
  
  # Output as widget (clusters again and might therefore create different order)
  saveWidget(d3heatmap(m.corr, colors = "Blues", sym = TRUE, xaxis_font_size = paste0(h, "pt"),
                         xaxis_height = h*20, yaxis_width = w*40), 
             file = normalizePath(paste0(str_split(outpdf,".pdf", simplify = TRUE)[1,1],".html"), mustWork = FALSE))

  # Output as ggplot
  df.corr = as.data.frame(m.corr) %>% 
    mutate(rowvar = rownames(m.corr)) %>% 
    gather(key = colvar, value = corr, -rowvar) 
  p = ggplot(df.corr, aes(rowvar, colvar)) +
    geom_tile(aes(fill = corr)) + 
    geom_text(aes(label = round(corr, 2))) +
    scale_fill_gradient(low = "white", high = "blue") +
    scale_x_discrete(limits = rev(rownames(m.corr))) +
    scale_y_discrete(limits = rownames(m.corr)) +
    labs(title = paste0(method," Correlation"), fill = "", x = "", y = "") +
    theme_my + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  ggsave(outpdf, p, width = w, height = h)
}



## Plot correlation of nominal variables
plot_corr_nomi <- function(outpdf, df, vars = nomi, 
                           w = 8, h = 8) {
  
  df[vars] = droplevels(df[vars])
  p = length(vars)
  m.chi = matrix(NA, p, p)
  rownames(m.chi) = colnames(m.chi) = paste0(vars," (#Levs: ", map_int(df[vars], ~ length(levels(.))), ")")
  for (i in 1:p) {
    for (j in 1:p) {
      print(paste0(vars[i], "...", vars[j]))
      
      # Corrected contingency coefficient
      tab = table(df[vars[c(i,j)]])
      M = min(dim(tab))
      m.chi[i,j] = assocstats(tab)$cont * sqrt(M / (M - 1))
    }
  }

  # Put in clustered order
  m.chi[which(is.na(m.chi))] = 0 #set NAs to 0
  ord = corrMatOrder(m.chi , order = "hclust")
  m.chi = m.chi[ord, ord]
  
  # Output as widget (clusters again and might therefore create different order)
  saveWidget(d3heatmap(m.chi, colors = "Blues", sym = TRUE, xaxis_font_size = paste0(h, "pt"),
                       xaxis_height = h*20, yaxis_width = w*40), 
             file = normalizePath(paste0(str_split(outpdf,".pdf", simplify = TRUE)[1,1],".html"), mustWork = FALSE))
    
  # Output as ggplot
  df.chi = as.data.frame(m.chi) %>% 
    mutate(rowvar = rownames(m.chi)) %>% 
    gather(key = colvar, value = chi, -rowvar)
  p = ggplot(df.chi, aes(rowvar, colvar)) +
    geom_tile(aes(fill = chi)) + 
    geom_text(aes(label = round(chi, 2))) +
    scale_fill_gradient(low = "white", high = "blue") +
    scale_x_discrete(limits = rev(rownames(m.chi))) +
    scale_y_discrete(limits = rownames(m.chi)) +
    labs(title = "Contig.Coef", fill = "", x = "", y = "") +
    theme_my + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  ggsave(outpdf, p, width = w, height = h)
}



## ROC, Calibration, Gain, Lift, Confusion
plot_performance = function(outpdf, yhat_holdout, y_holdout, color = "blue", colors = twocol, colgrad = colhex, ylim = NULL,
                            ncols = 2, nrows = 2, w = 12, h = 8) {
  # Prepare
  pred_obj = mysummary(data.frame(obs = y_holdout, pred = yhat_holdout))
  spearman = round(pred_obj["spearman"], 2)
  df.perf = data.frame(y = y_holdout, yhat = yhat_holdout, res = y_holdout - yhat_holdout, 
                       midpoint = cut(yhat_holdout, quantile(yhat_holdout, seq(0,1,0.2)), include.lowest = TRUE))
  df.preds = data.frame(type = c(rep("y", length(y_holdout)), rep("yhat", length(y_holdout))),
                        value = c(y_holdout, yhat_holdout))
  df.calib = df.perf %>% group_by(midpoint) %>% summarise(y = mean(y), yhat = mean(yhat))
  
  # Performance
  p_perf = ggplot(data = df.perf, aes_string("yhat", "y")) +
    geom_hex() + 
    scale_fill_gradientn(colors = colgrad, name = "count") +
    geom_smooth(color = "black", method = "gam", level = 0.95, size = 0.5) +
    geom_abline(intercept = 0, slope = 1, color = "grey") + 
    labs(title = bquote(paste("Observed vs. Fitted (", rho[spearman], " = ", .(spearman), ")", sep = "")),
    x = expression(hat(y))) +
    theme(plot.title = element_text(hjust = 0.5))
  if (length(ylim)) p_perf = p_perf + xlim(ylim) + ylim(ylim)
  
  # Residuals
  p_res = ggplot(data = df.perf, aes_string("yhat", "res")) +
    geom_hex() + 
    scale_fill_gradientn(colors = colgrad, name = "count") +
    geom_smooth(color = "black", method = "gam", level = 0.95, size = 0.5) +
    labs(title = "Residuals vs. Fitted", x = expression(hat(y)), y = expression(paste(hat(y) - y))) +
    theme(plot.title = element_text(hjust = 0.5))
  if (length(ylim)) p_res = p_res + xlim(ylim)

  # Distribution of predictions and target (plot similar to plot_distr_metr)
  p_pred = ggplot(data = df.preds, aes_string("value")) +
    geom_histogram(aes(y = ..density.., fill = type), bins = 40, position = "identity") +
    geom_density(aes(color = type)) +
    scale_fill_manual(values = alpha(colors, .2), labels = c("y", expression(paste(hat(y)))), name = " ") + 
    scale_color_manual(values = colors, labels = c("y", expression(paste(hat(y)))), name = " ") +
    labs(title = "Distribution", x = " ") +
    guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE))
  tmp = ggplot_build(p_pred)
  p.inner = ggplot(data = df.preds, aes_string("type", "value")) +
    geom_boxplot(aes_string(color = "type")) +
    coord_flip() +
    scale_y_continuous(limits = c(min(tmp$data[[1]]$xmin), max(tmp$data[[1]]$xmax))) +
    scale_color_manual(values = colors, name = " ") +
    theme_void() +
    theme(legend.position = "none")
  p_pred = p_pred + 
    scale_y_continuous(limits = c(-tmp$layout$panel_ranges[[1]]$y.range[2]/10, NA)) +
    theme_my +
    annotation_custom(ggplotGrob(p.inner), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0) 

  # Calibration
  p_calib = ggplot(df.calib, aes(yhat, y)) +
    geom_line(color = "blue") +
    geom_point(color = "blue") +  
    xlim(range(c(df.calib$y,df.calib$yhat))) +
    ylim(range(c(df.calib$y,df.calib$yhat))) +
    geom_abline(intercept = 0, slope = 1, color = "grey") + 
    labs(title = "Calibration", x = "Prediction Average (in quantile bin)", y = "Observation Average") +
    theme_my 
  
  # Plot
  plots = list(p_perf, p_res, p_calib, p_pred)
    labs(title = "Calibration", x = "Midpoint Predictions", y = "Observed Average") +
    theme_my 
  ggsave(outpdf, marrangeGrob(plots, ncol = ncols, nrow = nrows, top = NULL), width = w, height = h)
}



## Diagnosis: Residuals per predictor
plot_diagnosis = function(outpdf, df.diagnosis, res, vars = predictors, ylim = NULL,
                          ncols = 3, nrows = 3, w = 12, h = 8) {
  plots = map(vars, ~ {
    #. = vars[1]
    print(.)
    
    # Prepare
    df.res = df.diagnosis[.]
    df.res$res = res
    
    # Plot 
    p = ggplot(data = df.res, aes_string(., "res")) +
      geom_hex() + 
      scale_fill_gradientn(colors = colgrad, name = "count")
    if (is.factor(df.res[[.]])) {
      p = p + geom_smooth(color = "black", method = "gam", level = 0.95, size = 0.5, aes(group = 1)) 
    } else {
      p = p + geom_smooth(color = "black", method = "gam", level = 0.95, size = 0.5)
    } 
    p = p + labs(title = "Residuals", y = expression(paste(hat(y) - y))) +
      theme(plot.title = element_text(hjust = 0.5))
    if (length(ylim)) p = p + ylim(ylim)
    p  
  })
  ggsave(outpdf, marrangeGrob(plots, ncol = ncols, nrow = nrows, top = NULL), width = w, height = h)
}



## Variable Importance
plot_variableimportance = function(outpdf, vars, fit = fit.gbm, l.boot = NULL, 
                                   ncols = 5, nrows = 2, w = 18, h = 12) {
  # Group importances
  df.tmp = varImp(fit)$importance %>% 
    mutate(variable = rownames(varImp(fit)$importance)) %>% 
    filter(variable %in% vars) %>% 
    arrange(desc(Overall)) %>% 
    mutate(color = cut(Overall, c(-1,10,50,100), labels = c("low","middle","high"))) 
  
  # Plot
  p = ggplot(df.tmp) +
    geom_bar(aes(variable, Overall, fill = color), stat = "identity") +
    scale_x_discrete(limits = rev(df.tmp$variable)) +
    scale_fill_manual(values = c("blue","orange","red")) +
    labs(title = paste0("Top ", min(topn, length(predictors))," Important Variables (of ", length(predictors), ")"), 
         x = "", y = "") +
    coord_flip() +
    #geom_hline(yintercept = c(10,50), color = "grey", linetype = 2) +
    theme_my + theme(legend.position = "none")
  
  # Bootstrap lines
  if (!is.null(l.boot)) {
    df.tmpboot = map_df(l.boot, ~ {
      df = varImp(.)$importance
      df$variable = rownames(df)
      df
    } , .id = "run")
    p = p + geom_line(aes(variable, Overall, group = run), df.tmpboot, color = "grey", size = 0.1) +
      geom_point(aes(variable, Overall, group = run), df.tmpboot, color = "black", size = 0.3)
  }
  ggsave(outpdf, p, width = 8, height = 6)
} 



## Partial Depdendence
plot_partialdependence = function(outpdf, vars, df = df.interpret, fit = fit.gbm, l.boot = NULL, 
                                  ylim = c(0,1), ncols = 5, nrows = 2, w = 18, h = 12) {

  # Final model
  model = fit$finalModel

  # Plot
  plots = map(vars, ~ {
    #. = vars[1]
    print(.)
    
    # Plot data 
    df.plot = plot(model, i.var = ., return.grid = TRUE) #get plot data 

    if (is.factor(df[[.]])) {
      # Width of bars correspond to freqencies
      tmp = table(df[,.])
      df.plot$width = as.numeric(tmp[df.plot[[.]]])/max(tmp)
      
      # Plot for a nominal variable
      p = ggplot(df.plot, aes_string(., "y")) +
        geom_bar(stat = "identity",  width = df.plot$width, fill = "grey", color = "black") +
        labs(title = ., x = "", y = expression(paste(hat(y)))) +
        #scale_y_continuous(limits = ylim) +
        coord_cartesian(ylim = ylim) +
        theme_my         
    } else {
      # Plot for a metric variable
      df.rug = data.frame(q = quantile(df[,.], prob = seq(.05, .95, .1)), y = 0)
      p = ggplot(df.plot, aes_string(., "y")) +
        geom_line(stat = "identity", color = "black") +
        geom_rug(aes(q, y), df.rug, sides = "b", col = "red") +
        labs(title = ., x = "", y = expression(paste(hat(y)))) +
        scale_y_continuous(limits = ylim) +
        theme_my      
    }
    
    # Add Bootstrap lines and dots
    if (!is.null(l.boot)) {
      
      # Do the same as above for each bootstrapped model
      varactual = .
      df.tmpboot = map_df(l.boot, ~ {
        model_boot = .$finalModel
        df.plot_boot = plot(model_boot, i.var = varactual, return.grid = TRUE)
        df.plot_boot
      } , .id = "run")
      
      if (is.factor(df[,.])) {
        p = p + 
          geom_line(aes_string(., "y", group = "run"), df.tmpboot, color = "lightgrey", size = 0.1) +
          geom_point(aes_string(., "y", group = "run"), df.tmpboot, color = "black", size = 0.3)
      } else {
        p = p + 
          geom_line(aes_string(., "y", group = "run"), df.tmpboot, color = "lightgrey") +
          geom_line(aes_string(., "y"), df.plot, stat = "identity", color = "black") #plot black line again
      }
    }
    # Add overall average
    p + geom_hline(yintercept = mean(df.interpret$target), linetype = 3)
  })
  ggsave(outpdf, marrangeGrob(plots, ncol = ncols, nrow = nrows, top = NULL), width = w, height = h)
}



## Plot Interactiontest
plot_interactiontest = function(outpdf, vars, df = df.interpret, fit = fit.gbm, l.boot = NULL, 
                                ncols = 4, w = 18, h = 12) {
  
  # Derive interaction matrix for topn important variables
  pred_inter = setdiff(vars,"INT") #remove INT from testing variables
  k = length(pred_inter)
  m.inter = matrix(0, k, k)
  for (i in 1:(k - 1)) {
    for (j in (i + 1):k) {
      # Interaction Test
      m.inter[i,j] = interact.gbm(fit$finalModel, df[pred_inter], pred_inter[c(i,j)], fit$finalModel$tuneValue$n.trees)
      m.inter[j,i] = m.inter[i,j]
    }
  }
  colnames(m.inter) = pred_inter
  rownames(m.inter) = pred_inter
  m.inter[is.na(m.inter)] = 0
  
  
  ## Plot in correlation matrix style
  df.inter = as.data.frame(m.inter) %>% 
    mutate(rowvar = rownames(m.inter)) %>% 
    gather(key = colvar, value = inter, -rowvar)
  p = ggplot(df.inter, aes(rowvar, colvar)) +
    geom_tile(aes(fill = inter)) + 
    geom_text(aes(label = round(inter, 2))) +
    scale_fill_gradient(low = "white", high = "blue") +
    scale_x_discrete(limits = rownames(m.inter)) + 
    scale_y_discrete(limits = rev(rownames(m.inter))) +
    labs(title = "Interaction", fill = "", x = "", y = "") +
    theme_my +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  ggsave(outpdf, p, width = w, height = h)
  
  
  
  if (!is.null(l.boot)) {
    # Do the same as above for each bootstrapped model and collect
    df.inter_boot = map_df(l.boot, ~ {
      for (i in 1:(k - 1)) {
        for (j in (i + 1):k) {
          # Interaction Test
          m.inter[i,j] = interact.gbm(.$finalModel, df[pred_inter], pred_inter[c(i,j)], .$finalModel$tuneValue$n.trees)
          m.inter[j,i] = m.inter[i,j]
        }
      }
      m.inter[is.na(m.inter)] = 0
      
      df.inter = as.data.frame(m.inter) %>% 
        mutate(rowvar = rownames(m.inter)) %>% 
        gather(key = colvar, value = inter, -rowvar)
      df.inter
    }, .id = "run")
    
    
    # Same plot but now facetting
    p_boot = ggplot(df.inter_boot, aes(rowvar, colvar)) +
      geom_tile(aes(fill = inter)) + 
      geom_text(aes(label = round(inter, 2))) +
      scale_fill_gradient(low = "white", high = "blue") +
      scale_x_discrete(limits = rownames(m.inter)) + 
      scale_y_discrete(limits = rev(rownames(m.inter))) +
      labs(title = "Interaction per Bootstrap Run", fill = "", x = "", y = "") +
      theme_my +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
      facet_wrap( ~ run, ncol = ncols)
    ggsave(paste0(str_split(outpdf,".pdf", simplify = TRUE)[1,1],"_boot.pdf"), p_boot, width = w, height = h)
  }
}




## Plot interactions of m.gbm
plot_inter = function(outpdf, vars = inter, df = df.interpret, fit = fit.gbm, 
                      ylim = c(0,1), w = 12, h = 8) {
  # outpdf="./output/interaction1.pdf"; vars=inter1; df=df.interpret; fit=fit; w=12; h=8; ylim = c(0,.3)
  
  # Final model
  model = fit$finalModel

  # Marginal plots for digonal
  plots_marginal = map(vars, ~ {
    #.=vars
    # Get interaction data
    df.plot = plot(model, i.var = ., return.grid = TRUE) #get plot data

    if (is.factor(df[[.]])) {
      # Width of bars correspond to freqencies
      tmp = table(df[,.])
      df.plot$width = as.numeric(tmp[df.plot[[.]]])/max(tmp)
      
      # Plot for a nominal variable
      p = ggplot(df.plot, aes_string(., "y", fill = .)) +
        geom_bar(stat = "identity",  width = df.plot$width, color = "black") +
        labs(title = ., x = "", y = expression(paste(hat(y)))) +
        scale_fill_manual(values = manycol) +
        #scale_y_continuous(limits = ylim) +
        coord_cartesian(ylim = ylim) +
        theme_my + 
        theme(legend.position = "none")
    } else {
      # Plot for a metric variable
      df.rug = data.frame(q = quantile(df[,.], prob = seq(.05, .95, .1)), y = 0)
      p = ggplot(df.plot, aes_string(., "y")) +
        geom_line(stat = "identity", color = "black") +
        geom_rug(aes(q, y), df.rug, sides = "b", col = "red") +
        labs(title = ., x = "", y = expression(paste(hat(y)))) +
        scale_y_continuous(limits = ylim) +
        theme_my      
    }
    # Add overall average
    p + geom_hline(yintercept = mean(df.interpret$target), linetype = 3)
  })
  
  # Interaction plots 
  df.plot = plot(model, i.var = vars, return.grid = TRUE) #get plot data

  if (sum(map_lgl(df.plot[vars], ~ is.factor(.))) == 2) {
    # Mosaicplot for nominal-nominal interaction
    plots_inter = map(1:2, ~ { 
      if (.==2) vars = rev(vars)
      #tmp = table(df[vars[2]])
      #df.plot[[vars[1]]] = factor(df.plot[[vars[1]]], levels = rev(levels(df.plot[[vars[1]]])))
      ggplot(df.plot, aes_string(vars[2], "y", fill = vars[1])) +
        geom_bar(stat = "identity", position = "fill") + #, width = rep(tmp/max(tmp), 5)) +
        scale_fill_manual(values = manycol) +
        labs(y = "", x = "") +
        theme_my      
    })
  }  
  if (sum(map_lgl(df.plot[vars], ~ is.factor(.))) == 1) {
    # Grouped line chart for metric-nominal interaction
    p = ggplot(df.plot, aes_string(vars[1], "y", color = vars[2])) +
      geom_line(stat = "identity") +
      labs(y = expression(paste(hat(y)))) +
      scale_color_manual(values = manycol) +
      scale_y_continuous(limits = ylim) +
      guides(color = guide_legend(reverse = TRUE)) +
      theme_my      
    plots_inter = list(p,p)
  }   
  if (sum(map_lgl(df.plot[vars], ~ is.factor(.))) == 0) {
    # Grouped (by quantiles) line chart for metric-metric interaction
    plots_inter = map(1:2, ~ { 
      if (.==2) vars = rev(vars)
      val_near_quant =   map_dbl(quantile(df[[vars[2]]], seq(.05,.95,.1)), ~ {
        df.plot[[vars[2]]][which.min(abs(df.plot[[vars[2]]] - .))]})
      i.tmp = df.plot[[vars[2]]] %in% val_near_quant
      df.tmp = df.plot[i.tmp,]
      df.tmp[vars[2]] = factor( round(df.tmp[[vars[2]]],2) )
      
      ggplot(df.tmp, aes_string(vars[1], "y", color = vars[2])) +
        geom_line(stat = "identity") +
        labs(y = expression(paste(hat(y)))) +
        scale_color_manual(values = manycol) +
        scale_y_continuous(limits = ylim) +
        guides(color = guide_legend(reverse = TRUE)) +
        theme_my   
    })
  } 
  
  # Arrange plots
  plots = list(plots_marginal[[1]], plots_inter[[1]], plots_inter[[2]], plots_marginal[[2]])
  ggsave(outpdf, marrangeGrob(plots, ncol = 2, nrow = 2, top = NULL), width = w, height = h)
}




# Animated Interaction of 2 metric variables
plot_inter_active = function(outfile, vars = inter, df = df.interpret, fit = fit.gbm, duration = 15) {
  
  # Final model
  model = fit$finalModel

  # Get interaction data
  df.plot = plot(model, i.var = vars, return.grid = TRUE) #get plot data 

  # Prepare 3d plot
  x = unique(df.plot[[vars[1]]])
  y = unique(df.plot[[vars[2]]])
  z = matrix(df.plot$y, length(x), length(y), byrow = FALSE)
  nrz = nrow(z)
  colcut = cut((z[-1, -1] + z[-1, -nrz] + z[-nrz, -1] + z[-nrz, -ncol(z)])/4, 100)
  
  # html Widget
  persp3d(x, y, z, col = col3d[colcut], phi = 30, theta = 50, axes = T, ticktype = 'detailed',
          xlab = vars[1], ylab = vars[2], zlab = "")
  writeWebGL(dir = file.path(paste0(outfile)), width = 1000)
  rgl.close()
  
  # animated gif
  open3d("windowRect" = 2*c(20,20,400,400))
  persp3d(x, y, z, col = col3d[colcut], phi = 30, theta = 50, axes = T, ticktype = 'detailed',
          xlab = vars[1], ylab = vars[2], zlab = "")
  movie3d(spin3d(axis = c(0,0,1), rpm = 2), duration = duration, convert = NULL, clean = TRUE, movie = "test",
          dir = paste0(outfile))
  rgl.close()
  
}



#######################################################################################################################-
# Caret definition of MicrosofMl algorithms ----
#######################################################################################################################-

## rxFastTrees (boosted trees)

ms_boosttree = list()
ms_boosttree$label = "MicrosoftML rxFastTrees"
ms_boosttree$library = "MicrosoftML"
ms_boosttree$type = "Classification"
ms_boosttree$parameters = 
  read.table(header = TRUE, sep = ",", strip.white = TRUE, 
             text = "parameter,class,label
             numTrees,numeric,Boosting Interations
             numLeaves,numeric,Number of Leaves
             learningRate,numeric,Shrinkage"                             
  )

ms_boosttree$grid = function(x, y, len = NULL, search = "grid") {
  if (search == "grid") {
    out <- expand.grid(numTrees = floor((1:len) * 50),
                       numLeaves = 2^seq(1, len),
                       learningRate = .1)
  } else {
    out <- data.frame(numTrees = floor(runif(len, min = 1, max = 5000)),
                      numLeaves = 2^sample(1:10, replace = TRUE, size = len),         
                      learningRate = runif(len, min = .001, max = .6)) 
    out <- out[!duplicated(out),]
  }
  out
}

ms_boosttree$fit = function(x, y, wts, param, lev, last, classProbs, ...) { 
  #browser()
  theDots = list(...)
  #if (is.factor(y) && length(lev) == 2) {y = ifelse(y == lev[1], 1, 0)}
  #y = factor(y, levels = c(1,0))
  #x = as.matrix(x)
  modArgs <- list(formula = paste("y~", paste0(names(x), collapse = "+")),
                  data = cbind(x, y),
                  numTrees = param$numTrees,
                  numLeaves = param$numLeaves,
                  learningRate = param$learningRate,
                  type = "binary")
  if (length(theDots) > 0) modArgs <- c(modArgs, theDots)
  do.call("rxFastTrees", modArgs)
}

ms_boosttree$predict = function(modelFit, newdata, submodels = NULL) {
  #browser()
  out = rxPredict(modelFit, newdata)[,"Probability.Y"]
  if (length(modelFit$obsLevels) == 2) {
    out <- ifelse(out >= 0.5, "Y", "N")
  }
  out
}

ms_boosttree$prob = function(modelFit, newdata, submodels = NULL) {
  #browser()
  out = rxPredict(modelFit, newdata)[,"Probability.Y"]
  if (length(modelFit$obsLevels) == 2) {
    out <- cbind(out, 1-out)
    colnames(out) <- c("Y","N")
  }
  out
}

ms_boosttree$levels = function(x) {c("N","Y")}

ms_boosttree$sort = function(x) {
  x[order(x$numTrees, x$numLeaves, x$learningRate), ]
}



## rxForest (random Forest)

ms_forest = list()
ms_forest$label = "MicrosoftML rxFastForest"
ms_forest$library = "MicrosoftML"
ms_forest$type = "Classification"
ms_forest$parameters = 
  read.table(header = TRUE, sep = ",", strip.white = TRUE,
             text = "parameter,class,label
             numTrees,numeric,Number of Trees
             splitFraction,numeric,Fraction of features in split"
  )

ms_forest$grid = function(x, y, len = NULL, search = "grid") {
  if (search == "grid") {
    out <- expand.grid(numTrees = floor((1:len) * 50),
                       splitFraction = seq(0.01, 1, length.out = len))
  } else {
    out <- data.frame(numTrees = floor(runif(len, min = 1, max = 5000)),
                      splitFraction = runif(len, min = 0.01, max = 1))
    out <- out[!duplicated(out),]
  }
  out
}

ms_forest$fit = function(x, y, wts, param, lev, last, classProbs, ...) { 
  theDots = list(...)
  modArgs <- list(formula = paste("y~", paste0(names(x), collapse = "+")),
                  data = cbind(x, y),
                  numTrees = param$numTrees,
                  splitFraction = param$splitFraction,
                  type = "binary")
  if (length(theDots) > 0) modArgs <- c(modArgs, theDots)
  do.call("rxFastForest", modArgs)
}

ms_forest$predict = function(modelFit, newdata, submodels = NULL) {
  out = rxPredict(modelFit, newdata)[,"Probability.Y"]
  if (length(modelFit$obsLevels) == 2) {
    out <- ifelse(out >= 0.5, "Y", "N")
  }
  out
}

ms_forest$prob = function(modelFit, newdata, submodels = NULL) {
  out = rxPredict(modelFit, newdata)[,"Probability.Y"]
  if (length(modelFit$obsLevels) == 2) {
    out <- cbind(out, 1-out)
    colnames(out) <- c("Y","N")
  }
  out
}

ms_forest$levels = function(x) {c("N","Y")}

ms_forest$sort = function(x) {
  x[order(x$numTrees, x$splitFraction), ]
}

