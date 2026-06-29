data <- fread('/path_to_data/data.csv')

data <- data.frame(data)
data$sex <- as.factor(data$sex)
data$education_class <- as.factor(data$education_class)

boot.med <- function(data, mediator, exposure, outcome, covariates) {
  data <- data[sample(1:nrow(data), replace = T), ]
  
  lm_formula_mediator <- as.formula(paste(mediator, "~", exposure, "+", "sex + age_i0 + education_class"))
  alpha.temp <- coefficients(lm(lm_formula_mediator, data))[2]
  
  lm_formula_outcome_mediator <- as.formula(paste(outcome, "~", mediator, "+", exposure, "+", covariates))
  beta.temp <- coefficients(lm(lm_formula_outcome_mediator, data))[2]
  
  lm_formula_outcome_total <- as.formula(paste(outcome, "~", exposure, "+", mediator, "+", covariates))
  c_prime.temp <- coefficients(lm(lm_formula_outcome_total, data))[2] 
  
  lm_formula_outcome_exposure <- as.formula(paste(outcome, "~", exposure, "+", covariates))
  c.temp <- coefficients(lm(lm_formula_outcome_exposure, data))[2]
  
  IE1.l <- alpha.temp * beta.temp        
  IE2.l <- c.temp - c_prime.temp     

  TOT.l <- c.temp

  DE.l <- c_prime.temp
  
  results <- c(IE1.l, IE2.l, DE.l, TOT.l)
  return(results)
}

outcome <- 'total'
mediators <- colnames(data[,16:266])
exposure <- 'prs'
covariates = "sex + age_i0 + p53_time + education_class"
set.seed(123)


num_cores <- 80

result.m <- matrix(nrow = length(mediators), ncol = 18)
result.m <- data.frame(result.m)
rownames(result.m) <- mediators

for (mediator in mediators) {
  
  G <- 1000
  
  med.boot.cox <- mclapply(1:G, FUN = function(i) boot.med(data, mediator, exposure, outcome, covariates), mc.cores = num_cores)
  
  saveRDS(med.boot.cox, file = paste0('/path_to_save/', mediator,'_meta.rds'))
  
  IE1.cox <- unlist(lapply(med.boot.cox, '[[', 1)) 
  IE2.cox <- unlist(lapply(med.boot.cox, '[[', 2))
  DE.cox  <- unlist(lapply(med.boot.cox, '[[', 3)) 
  TOT.cox <- unlist(lapply(med.boot.cox, '[[', 4)) 
  
  calc_pval <- function(x) {
    pval <- mean(x < 0)
    if (pval > 0.5) {
      pval <- 1 - pval
    }
    return(2 * pval)
  }
  
  dir.cox <- mean(DE.cox)
  ci.dir.cox <- quantile(DE.cox, c(0.025, 0.975))
  
  tot.cox <- mean(TOT.cox)
  ci.tot.cox <- quantile(TOT.cox, c(0.025, 0.975))
  
  ind.ab.cox <- mean(IE1.cox)
  ind.ccp.cox <- mean(IE2.cox)
  ci.ab.cox <- quantile(IE1.cox, c(0.025, 0.975))
  ci.ccp.cox <- quantile(IE2.cox, c(0.025, 0.975))
  
  PM_1 <- ind.ab.cox / tot.cox  
  PM_2 <- ind.ccp.cox / tot.cox
  
  pval_IE1 <- calc_pval(IE1.cox)
  pval_IE2 <- calc_pval(IE2.cox)
  pval_DE   <- calc_pval(DE.cox)
  pval_TOT  <- calc_pval(TOT.cox)
  
  result.m[which(rownames(result.m) == mediator), ] <- c(dir.cox, ci.dir.cox[1], ci.dir.cox[2], pval_DE,
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
