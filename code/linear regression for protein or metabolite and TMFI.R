#protein
data_protein <- fread('/path_to_data/data.csv')

protein_columns <- colnames(data_protein)[2:2924]
data_protein <-data.frame(data_protein)

vars_to_factor <- c("sex", "education_class", "Asian", "Mixed", "Black", "Chinese", "others")
data_meta[vars_to_factor] <- lapply(data_meta[vars_to_factor], as.factor)

linear_regression <- function(protein_column, data) {
  formula <- as.formula(paste("total ~ age_i0 + sex + 
                              education_class + p53_time + Mixed + Asian + Black + Chinese+
                              others +", protein_column))
  model <- glm(formula, data = data)
  summary_model <- summary(model)
  results <- summary_model$coefficients
  return(results)
}

for (protein_column in protein_columns) {
  
  results <- linear_regression(protein_column, data_protein)
  
  result_df <- as.data.frame(results)
  write.csv(result_df, file = 
              paste0('/path_to_save/',
                     protein_column, "_results.csv"), row.names = TRUE)
}


# metabolite
data_meta <- fread('/path_data/data.csv')

meta_columns <- colnames(data_meta)[2:252]
data_meta <-data.frame(data_meta)

vars_to_factor <- c("sex", "education_class", "Asian", "Mixed", "Black", "Chinese", "others")
data_meta[vars_to_factor] <- lapply(data_meta[vars_to_factor], as.factor)

linear_regression <- function(meta_column, data) {
  formula <- as.formula(paste("total ~ age_i0 + sex + 
                              education_class + p53_time + Mixed + Asian + Black + Chinese+
                              others +", meta_column))
  model <- glm(formula, data = data)
  summary_model <- summary(model)
  results <- summary_model$coefficients
  return(results)
}


for (meta_column in meta_columns) {
  
  results <- linear_regression(meta_column, data_meta)
  
  result_df <- as.data.frame(results)
  write.csv(result_df, file = 
              paste0('/path_to_save/',
                     meta_column, "_results.csv"), row.names = TRUE)
}
