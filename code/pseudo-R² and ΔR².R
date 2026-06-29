data <- fread('/path_to_data/data.csv')

data <- data.frame(data)

data$sex <- as.factor(data$sex)
data$education_class <- as.factor(data$education_class)

events <- c(
  "death","cancer","CVD","digestive","respiratory",
  "endocrine","neurological","psychiatric",
  "musculoskeleta","genitourinary",
  "eye","ear","skin"
)

results <- data.frame()

for(event in events){
  
  time_var <- paste0(event, "_time")
  
  formula_null <- as.formula(
    paste0(
      "Surv(", time_var, ", ", event,
      ") ~ 1"
    )
  )
  
  formula0 <- as.formula(
    paste0(
      "Surv(", time_var, ", ", event,
      ") ~ age_i0 + sex + education_class"
    )
  )
 
  formula1 <- as.formula(
    paste0(
      "Surv(", time_var, ", ", event,
      ") ~ PRS + age_i0 + sex + education_class"
    )
  )
  
  model_null <- coxph(formula_null, data = data)
  
  model0 <- coxph(formula0, data = data)
  
  model1 <- coxph(formula1, data = data)
  
  n <- model1$n
  
  LL_null <- as.numeric(logLik(model_null))
  LL0 <- as.numeric(logLik(model0))
  LL1 <- as.numeric(logLik(model1))
  
  R2_model0 <- 1 - exp(
    -(2 / n) * (LL0 - LL_null)
  )
  
  R2_model1 <- 1 - exp(
    -(2 / n) * (LL1 - LL_null)
  )
  
  Delta_R2 <- 1 - exp(
    -(2 / n) * (LL1 - LL0)
  )
  
  results <- rbind(
    results,
    data.frame(
      Outcome = event,
      
      R2_model0 = R2_model0,
      R2_model1 = R2_model1,
      Delta_R2 = Delta_R2
    )
  )
}
