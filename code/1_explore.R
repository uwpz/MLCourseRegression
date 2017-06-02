
#######################################################################################################################-
# Initialize ----
#######################################################################################################################-

source("./code/0_init.R")


#######################################################################################################################-
# ETL ----
#######################################################################################################################-

## Read data
df.orig = read_csv("./data/thyroid.csv", col_names = TRUE)
skip = function() {
  # Check some stuff
  summary(mutate_if(df.orig, is.character, as.factor))
  hist(df.orig$T3, breaks = 20); hist(log(df.orig$T3), breaks = 20)
  summary(df.orig)
}
df = df.orig



## Define target (and impute just for educational reasons)
df = mutate(df, target = T3)
i.na = which(is.na(df$target))
df$target[i.na] = sample(df$target[-i.na], length(i.na))
summary(df$target)



#######################################################################################################################-
# Metric variables: Explore and adapt ----
#######################################################################################################################-

## Metric covariates
metr = c("age","TSH","TT4","T4U","FTI") # NOT USED: "TBG" -> only missings
summary(df[metr]) 



## Handling missings
# Remove covariates with too many missings
(misspct = map_dbl(df[metr], ~ round(sum(is.na(.)/nrow(df)), 3))) # misssing percentage
names(misspct[misspct > 0.9]) # none
metr = setdiff(metr, names(misspct[misspct > 0.9]))
summary(df[metr]) 

# Create mising indicators
(miss = metr[map_lgl(df[metr], ~ any(is.na(.)))])
df[paste0("MISS_",miss)] = map(df[miss], ~ as.factor(ifelse(is.na(.), "miss", "no_miss")))
summary(df[,paste0("MISS_",miss)])

# Impute missings with randomly sampled value
df[miss] = map(df[miss], ~ {
  i.na = which(is.na(.))
  .[i.na] = sample(.[-i.na], length(i.na) , replace = TRUE)
  . }
)
summary(df[metr]) 



## Check for outliers and skewness
plot_distr_metr("./output/distr_metr.pdf", df, vars = metr, misspct = misspct, ylim = c(0,4))



## Outliers + Skewness
# Winsorize
df[,metr] = map(df[metr], ~ {
  .[. > quantile(., 0.99, na.rm = TRUE)] = quantile(., 0.99, na.rm = TRUE)
  .[. < quantile(., 0.01, na.rm = TRUE)] = quantile(., 0.01, na.rm = TRUE)
  . }
)

# Log-Transform
tolog = c("TSH")
df[paste0(tolog,"_LOG_")] = map(df[tolog], ~ {if(min(., na.rm=TRUE) == 0) log(.+1) else log(.)})
metr = map_chr(metr, ~ ifelse(. %in% tolog, paste0(.,"_LOG_"), .)) #adapt metr and keep order



## Removing variables
# Remove Self predictors
metr = setdiff(metr, "xx")

# Remove highly/perfectly (>=98%) correlated (the ones with less NA!)
summary(df[metr])
plot_corr_metr("./output/corr_metr.pdf", df, metr, misspct = misspct) 
metr = setdiff(metr, c("xxx","xxx")) #Put at xxx the variables to remove



## Check
plot_distr_metr("./output/distr_metr_final.pdf", df, vars = metr, misspct = misspct)



#######################################################################################################################-
# Nominal variables: Explore and adapt ----
#######################################################################################################################-


## Nominal covariates
nomi = c("sex","on_thyroxine","query_on_thyroxine","on_antithyroid_medication","sick","pregnant","thyroid_surgery",
         "I131_treatment","query_hypothyroid","query_hyperthyroid","lithium","goitre","tumor","hypopituitary",
         "psych","referral_source","Class")  # NOT USED: "Class" -> its the target
nomi = union(nomi, paste0("MISS_",miss)) #Add missing indicators
summary(df[,nomi])



## Convert missings to own level ("_BLANK_")
df[nomi] = map(df[nomi], ~ as.factor(ifelse(is.na(.), "_BLANK_", as.character(.))))
summary(df[nomi])



## Create compact covariates for "too many members" columns 
topn_toomany = 4
(tmp = map_int(df[nomi], ~ length(levels(.)))) 
(toomany = names(tmp)[which(tmp > topn_toomany)])
(toomany = setdiff(toomany, c("BERUF_GRP_SL"))) #Define exception for important variables
df[paste0(toomany,"_OTHER_")] = map(df[toomany], ~ {
  tmp = table(.)
  as.factor(ifelse(. %in% names(tmp[order(tmp, decreasing=TRUE)][1:topn_toomany]), as.character(.), "_OTHER_")) 
})
nomi = map_chr(nomi, ~ ifelse(. %in% toomany, paste0(.,"_OTHER_"), .)) #Exchange name
summary(df[nomi], topn_toomany + 2)



## Check
plot_distr_nomi("./output/distr_nomi.pdf", df, vars = nomi, ylim = c(0,5), nrows = 4, w = 18, h = 12)



## Removing variables
# Remove Self-predictors
nomi = setdiff(nomi, "Class")

# Remove highly/perfectly (>=99%) correlated (the ones with less levels!)
plot_corr_nomi("./output/corr_nomi.pdf", df, nomi, w = 12, h = 12) 
nomi = setdiff(nomi, "MISS_FTI")



## Check
plot_distr_nomi("./output/distr_nomi_final.pdf", df, vars = nomi, nrows = 4, w = 18, h = 12)



#######################################################################################################################-
# Prepare final data ----
#######################################################################################################################-

##Final predictors
predictors = c(metr, nomi)



## Sample
i.samp = sample(1:nrow(df), 2000)
df.samp = df[i.samp,]



## Remove factors with less than 2 levels
df.samp = droplevels(df.samp) # Drop levels
onelev = which(sapply(df.samp[,predictors], function(x) ifelse(is.factor(x), length(levels(x)), NA)) <= 1)
predictors = setdiff(predictors, predictors[onelev])



## Final formula
predictors
formula = as.formula(paste("target", "~", paste(predictors, collapse = " + ")))


## Check
setdiff(colnames(df.samp), predictors)
setdiff(predictors,colnames(df.samp))
summary(df.samp[,predictors])


## Save
save.image("1_explore.rdata")







