data <- fread('/path_to_data/data.csv')
data <- data.frame(data)

diseases <- c('cancer', 'CVD','digestive', 'neurological',  'psychiatric',
              'respiratory','endocrine', 'genitourinary', 'musculoskeleta',
              'eye', 'ear', 'skin')
for (disease in diseases) {
  print(paste0(disease, ': ', threshold_age))
  threshold_age <- quantile(data$age_i2[data[disease] == 1], 0.10, na.rm = TRUE)
  
  data <- data %>%
    mutate(!!paste0(disease, "_group") := case_when(
      !!sym(disease) == 0 ~ "Control",
      !!sym(disease) == 1 & age_i2 <= threshold_age ~ "Early-onset",
      !!sym(disease) == 1 & age_i2 > threshold_age ~ "Other-onset"
    ))
}

all_results <- data.frame()

for (disease in diseases) {
  dis_group <- paste0(disease, '_group')
  data.lm <- data[,c('total', 'age_i2', dis_group)]
  colnames(data.lm) <- c('total', 'age_i2', 'group')
  
  model <- lm(total ~ group + age_i2,
              data=data.lm)
  
  model <- lm(total ~ group + age_i2,data=data.lm)
  summary(model)
  emm <- emmeans(model, ~group)
  
  summary_results <- summary(emm, infer = TRUE)
  
  summary_results <- summary_results %>%
    mutate(significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    disease = disease)
  
  pairwise_results <- data.frame(pairs(emm, adjust = "BH"))
  
  pairwise_results <- pairwise_results %>%
    mutate(
      significance = case_when(
        pairwise_results$p.value < 0.001 ~ "***",
        pairwise_results$p.value < 0.01 ~ "**",
        pairwise_results$p.value < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      disease = disease
    )
  
  summary_results['const'] <- pairwise_results$contrast
  summary_results['const_estimate'] <-  pairwise_results$estimate
  summary_results['const_SE'] <-  pairwise_results$SE
  summary_results['const_p'] <-  pairwise_results$p.value
  summary_results['const_sig'] <- pairwise_results$significance

  all_results <- rbind(all_results, summary_results)
}
