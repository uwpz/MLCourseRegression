
#######################################################################################################################-
# Initialize ----
#######################################################################################################################-

source("0_init.R")


#######################################################################################################################-
# ETL ----
#######################################################################################################################-

## Read data
df.orig = read_tsv("AmesHousing.txt", col_names = TRUE)
colnames(df.orig) = str_replace_all(colnames(df.orig), "[ /]","")
colnames(df.orig)[45:46] = c("firstFlrSF","secondFlrSF")
skip = function() {
  # Check some stuff
  hist(df.orig$SalePrice); hist(log(df.orig$SalePrice))
  summary(df.orig)
  summary(mutate_if(df.orig, is.character, as.factor))
}
df = df.orig



## Define target
df = mutate(df, target = log(SalePrice))
summary(df$target)



#######################################################################################################################-
# Metric variables: Explore and adapt ----
#######################################################################################################################-

## Metric covariates
metr = setdiff(colnames(df)[which(map_lgl(df, ~ is.numeric(.)))], c("SalePrice","target"))[1:30]
summary(df[metr]) 



## Handling missings
# Remove covariates with too many missings
tmp = map_dbl(df[metr], ~ round(sum(is.na(.)/nrow(df)), 3)) # misssing percentage
tmp[order(tmp, decreasing = TRUE)]
names(tmp[tmp > 0.9]) # none
metr = setdiff(metr, names(tmp[tmp > 0.9]))
summary(df[metr]) 

# Create mising indicators
(miss = metr[map_lgl(df[metr], ~ any(is.na(.)))])
df[paste0("MISS_",miss)] = map(df[miss], ~ as.factor(ifelse(is.na(.), "miss", "no_miss")))
summary(df[paste0("MISS_",miss)])

# Impute missings with randomly sampled value
imputepct = map_dbl(df[metr], ~ sum(is.na(.))/nrow(df))
df[miss] = map(df[miss], ~ {
  i.na = which(is.na(.))
  .[i.na] = sample(.[-i.na], length(i.na) , replace = TRUE)
  . }
)
summary(df[metr]) 



## Check for outliers and skewness
plot_distr_metr("./output/distr_metr.pdf", df, vars = metr, imputepct = imputepct) 



## Outliers + Skewness
# Winsorize
df[,metr] = map(df[metr], ~ {
  upper = quantile(., 0.99, na.rm = TRUE)
  lower = quantile(., 0.01, na.rm = TRUE)
  .[. > upper] = upper
  .[. < lower] = lower
  . }
)

# Log-Transform
tolog = c("TSH","age")
df[paste0(tolog,"_LOG_")] = map(df[tolog], ~ {if(min(., na.rm=TRUE) == 0) log(.+1) else log(.)})
metr = map_chr(metr, ~ ifelse(. %in% tolog, paste0(.,"_LOG_"), .)) #adapt metr and keep order



## Removing variables
# Remove Self predictors
metr = setdiff(metr, "T3")

# Remove highly/perfectly (>=98%) correlated (the ones with less NA!)
summary(df[metr])
plot_corr_metr("./output/corr_metr.pdf", df, metr, imputepct, w = 8, h = 8) 
metr = setdiff(metr, c("xxx","xxx")) #Put at xxx the variables to remove

## Check for outliers and skewness
plot_distr_metr("./output/distr_metr.pdf", df, vars = metr, imputepct = imputepct)



#######################################################################################################################-
# Nominal variables: Explore and adapt ----
#######################################################################################################################-

## Nominal covariates
nomi = setdiff(colnames(df), c(metr,"PID","SalePrice","target"))  # NOT USED: "Class" -> its the target
nomi = union(nomi, paste0("MISS_",miss)) #Add missing indicators
df[nomi] = map(df[nomi], ~ as.factor(as.character(.)))
summary(df[,nomi])




## Convert missings to own level ("_BLANK_")
df[nomi] = map(df[nomi], ~ fct_explicit_na(., "_BLANK_"))
summary(df[nomi])



## Create compact covariates for "too many members" columns 
topn_toomany = 20
tmp = map_int(df[nomi], ~ length(levels(.)))
tmp[order(tmp, decreasing = TRUE)]
(toomany = names(tmp)[which(tmp > topn_toomany)])
(toomany = setdiff(toomany, c("BERUF_GRP_SL"))) #Define exception for important variables
df[paste0(toomany,"_OTHER_")] = map(df[toomany], ~ fct_lump(., topn_toomany, other_level = "_OTHER_"))
nomi = map_chr(nomi, ~ ifelse(. %in% toomany, paste0(.,"_OTHER_"), .)) #Exchange name
summary(df[nomi], topn_toomany + 2)



## Check
plot_distr_nomi("./output/distr_nomi.pdf", df, nomi[1:37])


## Removing variables
# Remove Self-predictors
nomi = setdiff(nomi, "xxx")

# Remove highly/perfectly (>=99%) correlated (the ones with less levels!)
plot_corr_nomi("./output/corr_nomi.pdf", df, nomi, w = 8, h = 8) 
nomi = setdiff(nomi, "MISS_FTI")



#######################################################################################################################-
# Prepare final data ----
#######################################################################################################################-

##Final predictors
predictors = c(metr, nomi)


## Undersample (DO NOT OVERSAMPLE!!! as this prevents an honest test error)
summary(df$target)
df.samp = bind_rows(map(0:1, ~ {
  i.samp = which(df$target_num == .)
  set.seed((i+1)*999)
  sample_n(df[i.samp,], min(1000, length(i.samp))) #take all but 1000 at most
}))
summary(df.samp$target)



## Remove factors with just 1 level
#df.samp = droplevels(df.samp) # Drop levels
# predictors = setdiff(predictors, 
#                      predictors[which(sapply(df.samp[,predictors], 
#                                              function(x) ifelse(is.factor(x), length(levels(x)), NA)) == 1)])


## Final formula
predictors
formula = as.formula(paste("target", "~", paste(predictors, collapse = " + ")))


## Check
setdiff(colnames(df.samp), predictors)
setdiff(predictors,colnames(df.samp))
summary(df.samp[,predictors])


## Save
save.image("1_explore.rdata")







