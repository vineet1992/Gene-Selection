if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

if(!require(pipeliner)) devtools::install_github('alexioannides/pipeliner')
if(!require(pipeliner)) devtools::install_github('tidyverse/tidyverse')
if(!require(pipeliner)) devtools::install_github('tidyverse/modelr')

library(pipeliner)
library(modelr)
library(tidyverse)

data <- faithful

pipeline_func <- function(df) {
  pipeline(
    df,
    transform_features(function(df) {
      transmute(df, x1 = (waiting - mean(waiting)) / sd(waiting))
    }),
    
    transform_response(function(df) {
      transmute(df, y = (eruptions - mean(eruptions)) / sd(eruptions))
    }),
    
    estimate_model(function(df) {
      lm(y ~ 1 + x1, df)
    }),
    
    inv_transform_response(function(df) {
      transmute(df, pred_eruptions = pred_model * sd(eruptions) + mean(eruptions))
    })
  )
}

test_func <- function(model){
   print.Date(model)
   return(25)
}

# 5-fold cross-validation using machine learning pipelines
cv_rmse <- crossv_kfold(data, 5) %>% 
  mutate(model = map(train, ~ pipeline_func(as.data.frame(.x))),
         test <- test_func(model),
         predictions = map2(model, test, ~ predict(.x, as.data.frame(.y))),
         residuals = map2(predictions, test, ~ .x - as.data.frame(.y)$eruptions),
         rmse = map_dbl(residuals, ~ sqrt(mean(.x ^ 2)))) %>% 
  summarise(mean_rmse = mean(rmse), sd_rmse = sd(rmse))

cv_rmse