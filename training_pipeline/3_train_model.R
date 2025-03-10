# Loading libraries
library(dplyr)
library(caret)
library(Metrics)
library(beepr)

# Loading scaled data set
train_data_3_scaled_xgbt_01 <- read.csv("descriptors_dataframe_desc_center_0to1_maccs.csv")

# Removing ordinal numbers and SMILES notations
train_data_3_scaled_xgbt_01 <- select(train_data_3_scaled_xgbt_01,
                                      -X,
                                      -SMILES)

# Preparing control parameters which include three times repeated ten-fold Cross-Validation, tune grid and enabling verbose iter messages to control the model generating process
train_control_scaled_xgbt_01 <- trainControl(method = "repeatedcv",
                                             number = 10,
                                             repeats = 3,
                                             search = "grid",
                                             verboseIter = TRUE)

# Preparing xgbTree tune grid with viable hyperparameters
tune_grid_3_scaled_xgbt_01 <- expand.grid(nrounds = seq(250, 1500, by = 250),
                                          max_depth = seq(1, 9, by = 1),
                                          colsample_bytree = seq(0.33, 0.52, by = 0.02),
                                          min_child_weight = seq(1, 9, by = 1),
                                          eta = seq(0.01, 0.05, by = 0.01),
                                          gamma = seq(0.01, 0.15, by = 0.02),
                                          subsample = seq(0.3, 1.2, by = 0.3))

model_test_3_scaled_xgbt_01 <- train(HC50.exp ~ .,
                                     data = train_data_3_scaled_xgbt_01,
                                     method = "xgbTree",
                                     trControl = train_control_scaled_xgbt_01,
                                     #tuneLength = 250, # tuneLength is only used when 'random' search mode is enabled in control parameters
                                     tuneGrid = tune_grid_3_scaled_xgbt_01,
                                     verbosity = 0)

# Saving results
print(model_test_3_scaled_xgbt_01)
write.csv((as.data.frame(model_test_3_scaled_xgbt_01$results)), 
          file = "xgbt_scaled_3_0to1.csv")
print(model_test_3_scaled_xgbt_01_$results)
View(model_test_3_scaled_xgbt_01_$results)
saveRDS(model_test_3_scaled_xgbt_01,
        file = "xgbt_scaled_3_0to1.rds")