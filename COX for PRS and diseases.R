data <- fread('/path_to_data/data.csv')
data <- data.frame(data)

vars_to_factor <- c("sex", "education_class")
data[vars_to_factor] <- lapply(data[vars_to_factor], as.factor)


diseases <- c('death','cancer', 'CVD', 'digestive', 'respiratory', 'endocrine', 'neurological', 'psychiatric',
              'musculoskeleta', 'genitourinary',
              'eye', 'ear', 'skin')

for (disease in diseases) {
  dis_time <- paste0(disease, '_time')
  
  cox_formula <- as.formula(
    paste0('Surv(', dis_time, ', ', disease, ') ~ PRS + age_i0 + sex + education_class')
  )

  print(cox_formula)
  
  
  cox_model <- coxph(cox_formula, data = data)
  print(summary(cox_model))
  
  model_summary <- summary(cox_model)

  col <- c(colnames(model_summary$coefficients), colnames(model_summary$conf.int))
  results_matrix <- matrix(nrow = nrow(model_summary$coefficients), ncol = 9)
  
  results_matrix[,1:5] <- model_summary$coefficients
  results_matrix[,6:9] <- model_summary$conf.int
  colnames(results_matrix) <- col
  rownames(results_matrix) <- rownames(model_summary$coefficients)
  disease
  
  print(disease)
  print(vif(cox_model))
  
  predictors <- model.matrix(cox_model)[, -1]
  
  cor_matrix <- cor(predictors)
  print(cor_matrix)
  
  write.csv(results_matrix, paste0('path_to_save/cox_', 
                                   disease, '.csv'))
}
