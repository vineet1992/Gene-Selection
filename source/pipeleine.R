if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

if(!require(pipeliner)) devtools::install_github('alexioannides/pipeliner')
if(!require(pipeliner)) devtools::install_github('tidyverse/tidyverse')
if(!require(pipeliner)) devtools::install_github('tidyverse/modelr')

library(pipeliner)
library(modelr)
library(tidyverse)


data <- read.csv(paste0(getwd(),'/../data/GSE11121_2.csv'),header=FALSE)

col = read.csv(paste0(getwd(),"/../data/headers.csv"),header=TRUE)[1:ncol(data)]
col <- colnames(col)
## remove the for loop ##




####Cannot use this for-loop, it should be abandoned
####A legal R variable cannot start with a number, so this causes problems when using a formula object
for(i in 1:length(col)) ## removes any leading X in the header names
{
  if(substring(col[i], 1, 1) == 'X')
    col[i] <- substr(col[i], 2, nchar(col[i]))
}
colnames(data) <- col

pipeline_func <- function(df,FUN = lm) {
  
  ###This formula assumes that the first variable is the Target variable of interest (you probably need to include target as a param)
  form = as.formula(paste(colnames(df)[1]," ~ .",sep=""))
  
  ###Correctly runs the linear model (or any other model)
  mdl = FUN(form,df)
  
  ##Return the model
  return(mdl)
}

###Linear model doesn't work if  | features | > | samples | so take a small chunk of features
temp = data[,1:100]


# 5-fold cross-validation using machine learning pipelines


###Again here the target variable name should not be hardcoded
cv_rmse <- crossv_kfold(temp, 5) %>% 
  mutate(model = map(train, ~ pipeline_func(as.data.frame(.x),lm)),
         predictions = map2(model, test, ~ predict(.x, as.data.frame(.y))),
         residuals = map2(predictions, test, ~ .x - as.data.frame(.y)$'X1007_s_at'),
         rmse = map_dbl(residuals, ~ sqrt(mean(.x ^ 2)))) %>% 
  summarise(mean_rmse = mean(rmse), sd_rmse = sd(rmse))

cv_rmse