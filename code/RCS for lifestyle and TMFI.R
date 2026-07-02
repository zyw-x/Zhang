data <- fread('/path_to_data/data.csv')
data <- data.frame(data)
data$sex <- as.factor(data$sex)
data$education_class <- as.factor(data$education_class)
data$Asian <- as.factor(data$Asian)
data$Mixed <- as.factor(data$Mixed)
data$Black <- as.factor(data$Black)
data$Chinese <- as.factor(data$Chinese)
data$others <- as.factor(data$others)
colnames(data)

covariates <- c("age_i2" ,"total","sex","Mixed","Asian",                      
                "Black","Chinese","others","education_class" ,"p53_time")

data.sleep <- na.omit(data[,c(covariates,'sleep_duration')])
data.smoke_pack <-  na.omit(data[,c(covariates,'smoke_pack')])
data.bmi <-  na.omit(data[,c(covariates,'bmi')])
data.drink <-  na.omit(data[,c(covariates,'drink')])
data.diet <-  na.omit(data[,c(covariates,'diet_i0_score')])
data.Summed <-  na.omit(data[,c(covariates,'Summed_MET_minutes_per_week_for_all_activity')])

data.Summed$Summed_MET_minutes_per_week_for_all_activity <- as.numeric(data.Summed$Summed_MET_minutes_per_week_for_all_activity)
data.Summed <- na.omit(data.Summed)
dd <- datadist(data.Summed)
options(datadist='dd') 
colnames(data.Summed)

aic_values <- sapply(3:5, function(k) {
  model <- ols(total ~ rcs(Summed_MET_minutes_per_week_for_all_activity, k) + age_i0 + 
                 sex + education_class + p53_time+ Mixed + Asian + Black + Chinese + others,
               data = data.Summed)
  AIC(model)
})
aic_values
names(aic_values) <- c("3", "4", "5")
best_k <- as.numeric(names(which.min(aic_values)))

dd$limits$Summed_MET_minutes_per_week_for_all_activity[2]<-median(data.Summed$Summed_MET_minutes_per_week_for_all_activity)
fit1<-ols(total ~ rcs(Summed_MET_minutes_per_week_for_all_activity, best_k) + age_i0 + 
            sex + education_class + p53_time+ Mixed + Asian + Black + Chinese + others,
          data = data.Summed)
an1 <- anova(fit1)
an1
HR1 <- rms::Predict(fit1,Summed_MET_minutes_per_week_for_all_activity, fun=exp, ref.zero = TRUE)
