data <- fread('/path_to_data/data.csv')
data <-data.frame(data)

data <- data %>%
  mutate(
    age = case_when(
      age_i2 < 60 ~ 0,      
      age_i2 >= 60 ~ 1)
  )

data$age <- as.factor(data$age)
data$sex <- as.factor(data$sex)
data$education_class <- as.factor(data$education_class)
data$Asian <- as.factor(data$Asian)
data$Mixed <- as.factor(data$Mixed)
data$Black <- as.factor(data$Black)
data$Chinese <- as.factor(data$Chinese)
data$others <- as.factor(data$others)

lifes <- c("age","sex")
diseases <- c("death","CVD", "cancer", "digestive", "neurological", "psychiatric",
              "endocrine", "genitourinary", "musculoskeleta", "eye", "ear",
              "skin", "respiratory")

for (life in lifes) {
  for (disease in diseases) {
    dis_time <- paste0(disease, '_time')
    
    data.cox <- data[, c(disease, dis_time, life, 'total', 'sex', 'age_i2', 'p53_time', 
                         'education_class', 'Mixed', 'Asian', 'Black', 'Chinese', 'others')]
    data.cox <- na.omit(data.cox)
    
    data.cox[[life]] <- as.factor(data.cox[[life]])
    
    cox_formula <- as.formula(
      paste0('Surv(', dis_time, ', ', disease, ') ~ ', life ,'* total + age_i2 + p53_time + sex + education_class + Mixed + Asian + Black + Chinese + others')
    )
    
    cox_model <- coxph(cox_formula, data = data.cox)
    
    model_summary <- summary(cox_model)
    
    col <- c(colnames(model_summary$coefficients), colnames(model_summary$conf.int))
    results_matrix <- matrix(nrow = nrow(model_summary$coefficients), ncol = 9)
    
    results_matrix[,1:5] <- model_summary$coefficients
    results_matrix[,6:9] <- model_summary$conf.int
    colnames(results_matrix) <- col
    rownames(results_matrix) <- rownames(model_summary$coefficients)
    
    p_values <- model_summary$coefficients[, "Pr(>|z|)"]
    
    print(p_values)
    write.csv(results_matrix, paste0('/path_to_save/', 
                                     life,'_',disease, '_cox.csv'))
    
  }
}
