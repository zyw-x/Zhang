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

boot.med <- function(data, mediators, expose){
  data <- data[sample(1:nrow(data), replace = T), ]
  
  lm_formula <- as.formula(paste(mediators, "~ ", expose, " + sex + age_i0 + p53_time + Mixed + Asian + Black + Chinese + others + education_class"))
  alpha.temp <- coefficients(lm(lm_formula, data))[2]
  
  cox_formula1 <- as.formula(paste("Surv(time, status == 1) ~", mediators, "+ ", expose, " + 
                                    sex + age_i2 + p53_time + Mixed + Asian + Black + Chinese + others + education_class"))
  beta.temp <- coefficients(coxph(cox_formula1, data, method = "breslow"))[1]
  
  cox_formula2 <- as.formula(paste("Surv(time, status == 1) ~ ", expose, " +", mediators, "+  
                                   sex + age_i2 + p53_time + Mixed + Asian + Black + Chinese + others + education_class"))
  c_prime.temp <- coefficients(coxph(cox_formula2, data, method = "breslow"))[1]
  
  cox_formula3 <- as.formula(paste("Surv(time, status == 1) ~ ", expose, " +  
                                   sex + age_i2 + p53_time + Mixed + Asian + Black + Chinese + others + education_class"))
  c.temp <- coefficients(coxph(cox_formula3, data, method = "breslow"))[1]
  
  IE1.l <- alpha.temp * beta.temp
  IE2.l <- c.temp - c_prime.temp
  DE.l <- c_prime.temp
  TOT.l <- c.temp
  
  results <- c(IE1.l, IE2.l, DE.l, TOT.l)
  return(results)
}

num_cores <- 120

for (life in lifes) {

  result.m <- matrix(nrow = length(diseases), ncol = 18)
  result.m <- data.frame(result.m)
  rownames(result.m) <- diseases
  
  for (disease in diseases) {
    time <- paste0(disease, '_time')
    data.md <- data[,c(disease, time, 'total',life,'sex','age_i0','age_i2', 'p53_time','education_class','Mixed' , 'Asian', 'Black', 'Chinese', 'others')]
    colnames(data.md) <- c('status','time', 'total',life,'sex','age_i0','age_i2', 'p53_time','education_class','Mixed' , 'Asian', 'Black', 'Chinese', 'others')
    data.md <- na.omit(data.md)
    
    G <- 1000
    
    med.boot.cox <- mclapply(1:G, FUN = function(i) boot.med(data.md, 'total', life), mc.cores = num_cores)
    saveRDS(med.boot.cox, file = paste0('/path_to_save/',life,'_',disease,'_md.rds'))

    IE1.cox <- unlist(lapply(med.boot.cox, '[[', 1)) 
    IE2.cox <- unlist(lapply(med.boot.cox, '[[', 2))
    DE.cox <- unlist(lapply(med.boot.cox, '[[', 3)) 
    TOT.cox <- unlist(lapply(med.boot.cox, '[[', 4)) 
    
    calc_pval <- function(x) {
      pval <- mean(x < 0)
      if (pval > 0.5) {
        pval <- 1 - pval
      }
      return(2 * pval)
    }
    
    dir.cox <- mean(DE.cox)
    ci.dir.cox <- as.table(quantile(DE.cox, c(0.025, 0.975)))
    
    tot.cox <- mean(TOT.cox)
    ci.tot.cox <- as.table(quantile(TOT.cox, c(0.025, 0.975)))
    
    ind.ab.cox <- mean(IE1.cox)
    ind.ccp.cox <- mean(IE2.cox)
    ci.ab.cox <- as.table(quantile(IE1.cox, c(0.025, 0.975))) 
    ci.ccp.cox <- as.table(quantile(IE2.cox, c(0.025, 0.975)))
    
    PM_1 <- ind.ab.cox / tot.cox  
    PM_2 <- ind.ccp.cox / tot.cox
    
    pval_IE1 <- round(calc_pval(IE1.cox), 6)
    pval_IE2 <- round(calc_pval(IE2.cox), 6)
    pval_DE <- round(calc_pval(DE.cox), 6)
    pval_TOT <- round(calc_pval(TOT.cox), 6)
    
    result.m[which(rownames(result.m) == disease), ] <- c(dir.cox, ci.dir.cox[1], ci.dir.cox[2], pval_DE,
                                                          tot.cox, ci.tot.cox[1], ci.tot.cox[2],pval_TOT,
                                                          ind.ab.cox, ci.ab.cox[1], ci.ab.cox[2], pval_IE1,
                                                          ind.ccp.cox, ci.ccp.cox[1], ci.ccp.cox[2], pval_IE2,
                                                          PM_1, PM_2)
  }
  
  colnames(result.m) <- c('dir.cox', 'ci.dir.cox_2.5', 'ci.dir.cox_97.5','dir.p',
                          'tot.cox', 'ci.tot.cox_2.5', 'ci.tot.cox_97.5', 'tot.p',
                          'ind.ab.cox', 'ci.ab.cox_2.5', 'ci.ab.cox_97.5', 'ind.ab.p',
                          'ind.ccp.cox', 'ci.ccp.cox_2.5', 'ci.ccp.cox_97.5','ind.ccp.p',
                          'PM_1','PM_2')
  
  write.csv(result.m, paste0('/path_to_save/', life, "_mediation.csv"), row.names = TRUE)
}
