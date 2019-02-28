if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

if(!require(pipeliner)) devtools::install_github('alexioannides/pipeliner')
if(!require(pipeliner)) devtools::install_github('tidyverse/tidyverse')
if(!require(pipeliner)) devtools::install_github('tidyverse/modelr')

library(pipeliner)
library(modelr)
library(tidyverse)

data <- read.csv(paste0(getwd(),'/data/GSE11121_2.csv'),header=FALSE)

col = read.csv(paste0(getwd(),"/data/headers.csv"),header=TRUE)[1:22284]
col <- colnames(col)
## remove the for loop ##
for(i in 1:length(col)) ## removes any leading X in the header names
{
  if(substring(col[i], 1, 1) == 'X')
    col[i] <- substr(col[i], 2, nchar(col[i]))
}
colnames(data) <- col


pipeline_func <- function(df,y,algo) {
  pipeline(
    df,
    #transform_features(function(df) {
    ## CHANGES for all columns
    #  as.data.frame(apply(df, MARGIN=2 , scale))
    #}),
    
    estimate_model(function(df) {
      algo(y ~ ., df)
    })
  )
}

test_func <- function(model){
  print.Date(model)
  return(25)
}
.x <- data
.y <- data$`202472_at`
# 5-fold cross-validation using machine learning pipelines
cv_rmse <- crossv_kfold(data, 5) %>% 
  mutate(model = map(train, ~ pipeline_func(as.data.frame(.x[,1:100]),.y,lm)),
         test <- test_func(model),
         predictions = map2(model, test, ~ predict(.x, as.data.frame(.y))),
         residuals = map2(predictions, test, ~ .x - as.data.frame(.y)),
         rmse = map_dbl(residuals, ~ sqrt(mean(.x ^ 2)))) %>% 
  summarise(mean_rmse = mean(rmse), sd_rmse = sd(rmse))

cv_rmse
