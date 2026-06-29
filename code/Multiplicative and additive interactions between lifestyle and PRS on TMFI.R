data <- fread('/path_to_data/data.csv')
data <-data.frame(data)

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
      Summed_MET_minutes_per_week_for_all_activity < 1719.5 ~ 1,      
      Summed_MET_minutes_per_week_for_all_activity >= 1719.5 ~ 0
    ),
    PRS = ifelse(PRS >= median(PRS,na.rm = TRUE), 1, 0),
    drink = ifelse(drink >= median(drink,na.rm = TRUE), 1, 0),
    smoke_pack = ifelse(smoke_pack > median(smoke_pack, na.rm = TRUE), 1, 0),
    diet_i0_score = ifelse(diet_i0_score < median(diet_i0_score,na.rm = TRUE), 1, 0)
  )

data$longsleep <- as.factor(data$longsleep)
data$shortsleep <- as.factor(data$shortsleep)
data$obesity <- as.factor(data$obesity)
data$highactivity <- as.factor(data$highactivity)
data$PRS <- as.factor(data$PRS)
data$drink <- as.factor(data$drink)
data$smoke_pack <- as.factor(data$smoke_pack)
data$diet_i0_score <- as.factor(data$diet_i0_score)
data$sex <- as.factor(data$sex)
data$education_class <- as.factor(data$education_class)

diseases <- c('total')
lifes <- c("longsleep","shortsleep","obesity","highactivity",
           "smoke_pack","drink","diet_i0_score")

for (life in lifes) {
  for (disease in diseases) {
    data.cox <- data[, c(disease, life, 'PRS', 'sex', 'age_i0','p53_time', 'education_class')]
    data.cox <- na.omit(data.cox)
    
    data.cox$PRS <- as.factor(data.cox$PRS)
    data.cox[[life]] <- as.factor(data.cox[[life]])
    
    
    
    lm_formula <- as.formula(
      paste0('total ~ ', life,'*PRS + age_i0 + p53_time + sex + education_class')
    )
    
    model <- glm(lm_formula,data=data.cox)
    
    result <- data.frame(summary(model)$coefficients)
    
    coefficients <- coef(model)
    std_errors <- sqrt(diag(vcov(model)))
    
    result['OR'] <- exp(coefficients)
    result['lower.95'] <- exp(coefficients - 1.96 * std_errors)
    result['upper.95'] <- exp(coefficients + 1.96 * std_errors)
    write.csv(result, paste0('/path_to_save/', 
                             life,'.csv'))
    
    out <- interactionR(model,
                        exposure_names =c("PRS", life),
                        ci.type ="delta", ci.level = 0.95,em = F, recode = F )
    
    result <- out$dframe
    result$life <- life
    result$disease <- disease 
    write.csv(result, paste0('/path_to_save/', 
                             life, '.csv'), row.names = F)
  }
}
