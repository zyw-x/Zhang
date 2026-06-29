data <- fread('/path_to_data/data.csv')
data <- data.frame(data)
data$sex <- as.factor(data$sex)
data$education_class <- as.factor(data$education_class)
data$Asian <- as.factor(data$Asian)
data$Mixed <- as.factor(data$Mixed)
data$Black <- as.factor(data$Black)
data$Chinese <- as.factor(data$Chinese)
data$others <- as.factor(data$others)

diseases <- c('cancer', 'CVD','digestive', 'respiratory','endocrine',
              'neurological',  'psychiatric','musculoskeleta', 'genitourinary', 
              'eye', 'ear', 'skin')

for (disease in diseases) {
  data.lo <- data[, c(disease,'age_i2',
                      'total',"sex",'education_class',
                      'Mixed' , 'Asian', 'Black', 'Chinese', 'others')]
  colnames(data.lo) <- c('status','age_i2',
                         'total',"sex",'education_class',
                         'Mixed' , 'Asian', 'Black', 'Chinese', 'others')
  
  logit_model <- glm(status ~ total + age_i2 + 
                       sex + education_class +Mixed + 
                       Asian + Black + Chinese + others, 
                     data = data.lo, 
                     family = binomial(link = "logit"))

  summary(logit_model)
  result <- data.frame(summary(logit_model)$coefficients)
  
  coefficients <- coef(logit_model)
  std_errors <- sqrt(diag(vcov(logit_model)))
  
  result['OR'] <- exp(coefficients)
  result['lower.95'] <- exp(coefficients - 1.96 * std_errors)
  result['upper.95'] <- exp(coefficients + 1.96 * std_errors)
  
  write.csv(result, paste0('/path_to_save/logit_', 
                           disease, '.csv'))
}
