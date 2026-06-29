data <- fread('/path_to_data/data.csv')
data <- data.frame(data)

data <- data %>%
  mutate(
    longsleep = case_when(
      sleep_duration < 6 ~ NA_real_,                  
      sleep_duration >= 6 & sleep_duration < 8 ~ 0, 
      sleep_duration >= 8 ~ 1                        
    ),
    shortsleep = case_when(
      sleep_duration < 6 ~ 1,                 
      sleep_duration >= 6 & sleep_duration < 8 ~ 0, 
      sleep_duration >= 8 ~ NA_real_                        
    ),
    obesity = case_when(
      bmi < 30 ~ 0,      
      bmi >= 30 ~ 1
    ),
    highactivity = case_when(
      Summed_MET_minutes_per_week_for_all_activity < 1719.5 ~ 0,      
      Summed_MET_minutes_per_week_for_all_activity >= 1719.5 ~ 1
    )
  )

data$longsleep <- as.factor(data$longsleep)
data$shortsleep <- as.factor(data$shortsleep)
data$obesity <- as.factor(data$obesity)
data$highactivity <- as.factor(data$highactivity)
data$sex <- as.factor(data$sex)
data$education_class <- as.factor(data$education_class)
data$Asian <- as.factor(data$Asian)
data$Mixed <- as.factor(data$Mixed)
data$Black <- as.factor(data$Black)
data$Chinese <- as.factor(data$Chinese)
data$others <- as.factor(data$others)

lifes <- c("longsleep","shortsleep","obesity","highactivity",
           "smoke_pack","drink","diet_i0_score")
diseases <- c("death","CVD", "cancer", "digestive", "neurological", "psychiatric",
              "endocrine", "genitourinary", "musculoskeleta", "eye", "ear",
              "skin", "respiratory")

for (life in lifes) {
  for (disease in diseases) {
    dis_time <- paste0(disease, '_time')
    
    cox_formula <- as.formula(
      paste0('Surv(', dis_time, ', ', disease, ') ~ ', life, ' + age_i2 + sex + education_class + p53_time + Mixed + Asian + Black + Chinese + others')
    )
    
    data.cox <- data[,c(disease, dis_time,life,'sex','age_i2', 'p53_time','education_class','Mixed' , 'Asian', 'Black', 'Chinese', 'others')]
    data.cox <- na.omit(data.cox)
    
    cox_model <- coxph(cox_formula, data = data.cox)
    
    model_summary <- summary(cox_model)
    
    col <- c(colnames(model_summary$coefficients), colnames(model_summary$conf.int))
    results_matrix <- matrix(nrow = nrow(model_summary$coefficients), ncol = 9)
    
    results_matrix[,1:5] <- model_summary$coefficients
    results_matrix[,6:9] <- model_summary$conf.int
    colnames(results_matrix) <- col
    rownames(results_matrix) <- rownames(model_summary$coefficients)
    
    p_values <- model_summary$coefficients[, "Pr(>|z|)"]
    
    print(p_values[1])
  
    write.csv(results_matrix, paste0('/path_to_save/', life,'@',disease, '_cox.csv'))
  }
}
