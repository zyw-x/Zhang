data <- fread('/path_to_data/data.csv')
data <- data.frame(data)

vars_to_factor <- c("sex", "education_class", "Asian", "Mixed", "Black", "Chinese", "others")
data[vars_to_factor] <- lapply(data[vars_to_factor], as.factor)

muscle_vars <- c("p24353_i2", "p24354_i2", "p23355_i2", "p23356_i2", "total")
muscle_titles <- c("Left anterior TMFI",
                   "Right anterior TMFI",
                   "Left posterior TMFI",
                   "Right posterior TMFI",
                   "TMFI")

cors_i0 <- c("tl", "reaction_time", "FEV1", "FVC", 'BMD', "grip_min", "grip_max", "frailty_score", "p30120_i0_Lymphocyte", "p30130_i0_Monocyte",
             "p30140_i0_Neutrophill", "p30710_i0_CRP", "SII")
cors_i2 <- c("father_death_age", "mother_death_age")

all_results <- list()

for (i in 1:length(muscle_vars)) {
  muscle_var <- muscle_vars[i]
  muscle_title <- muscle_titles[i]
  
  pcor_results <- data.frame()
  
  
  for (cor in cors_i0) {
    data.cor <- data[, c("sex", "age_i2", muscle_var, cor,
                         "education_class", 'Mixed', 'Asian', 'Black', 'Chinese', 
                         'others', 'p53_time')]
    data.cor <- na.omit(data.cor)
    
    pcor_result <- pcor.test(
      x = data.cor[[muscle_var]], 
      y = data.cor[[cor]],
      z = data.cor[, c("sex", "age_i2", "education_class",
                       'Mixed', 'Asian', 'Black', 'Chinese', 
                       'others', 'p53_time')],
      method = "pearson"  
    )
    pcor_result['R2'] <- pcor_result$estimate * pcor_result$estimate
    row.names(pcor_result) <- cor
    pcor_results <- rbind(pcor_results, pcor_result)
  }
  
  for (cor in cors_i2) {
    data.cor <- data[, c("sex", "age_i2", muscle_var, cor,
                         "education_class", 'Mixed', 'Asian', 'Black', 'Chinese',
                         'others')]
    data.cor <- na.omit(data.cor)
    
    pcor_result <- pcor.test(
      x = data.cor[[muscle_var]],
      y = data.cor[[cor]],
      z = data.cor[, c("sex", "age_i2", "education_class",
                       'Mixed', 'Asian', 'Black', 'Chinese',
                       'others')],
      method = "pearson"
    )
    pcor_result['R2'] <- pcor_result$estimate * pcor_result$estimate
    row.names(pcor_result) <- cor
    pcor_results <- rbind(pcor_results, pcor_result)
  }
  
  acc_r2 <- data.frame(
    Variable = rownames(pcor_results),
    Muscle_var = muscle_var,
    Muscle_title = muscle_title,
    R2 = pcor_results$R2,
    R_value = pcor_results$estimate,
    p_value = pcor_results$p.value,
    full_name = c("Telomere length", "Reaction time", 
                  'FEV1', 'FVC', 
                  'Heel BMD',       
                  "Weaker HGS", "Stronger HGS", "Frailty", "Lymphocyte count", "Monocyte count", "Neutrophil count", "C-reactive protein", "SII",
                  "Father's age at death", "Mother's age at death")
  )
  
  acc_r2$full_name <- factor(acc_r2$full_name, 
                             levels = c("Father's age at death", "Mother's age at death",
                                        "Frailty", "Stronger HGS", "Weaker HGS", "Reaction time", 
                                        'FVC', 'FEV1', 
                                        'Heel BMD',
                                        "Telomere length", "Neutrophil count", "Lymphocyte count", "Monocyte count", "C-reactive protein", "SII"))
  
  all_results[[muscle_var]] <- acc_r2
}


combined_results <- do.call(rbind, all_results)
rownames(combined_results) <- NULL
