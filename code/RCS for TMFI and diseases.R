data <- fread('/path_to_data/data.csv')
data <- data.frame(data)

vars_to_factor <- c("sex", "education_class", "Asian", "Mixed", "Black", "Chinese", "others")
data[vars_to_factor] <- lapply(data[vars_to_factor], as.factor)

iron_label <- c(
  "total" = "Thigh muscle fat infiltration(%)"
)

iron_cols <- c("total")
covariates <- c("age_i2" ,"sex","Mixed","Asian",                      
                "Black","Chinese","others","education_class")

for (iron_col in iron_cols) {
  df <- na.omit(data[,c(covariates,iron_col)])
  print(colnames(df))
  names(df)[names(df) == iron_col] <- "iron"
  
  dd <- datadist(df)
  options(datadist='dd') 
  aic_values <- sapply(3:5, function(k) {
    model <- cph(
      Surv(death_time, death) ~ rcs(iron, k) + age_i2 + sex +
        Mixed + Asian + Black + Chinese + others +
        education_class,
      data = df
    )
    
    AIC(model)
  })
  aic_values
  names(aic_values) <- c("3", "4", "5")
  best_k <- as.numeric(names(which.min(aic_values)))
  
  dd$limits$iron[2]<-min(df$iron)
  fit1 <- cph(
    Surv(death_time, death) ~ rcs(iron, best_k) + age_i2 + sex +
      Mixed + Asian + Black + Chinese + others +
      education_class,
    data = df
  )
  an1 <- anova(fit1)
  an1
  HR1 <- rms::Predict(fit1, iron, fun=exp, ref.zero = TRUE)
}
