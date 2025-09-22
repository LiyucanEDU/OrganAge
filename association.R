# =============================================================================
# Complete Analysis Pipeline for Organ-specific Biological Age
# =============================================================================
# Author: [Your Name]
# Date: 2025-09-21
# Dataset: Taizhou Imaging Study
# Objective: Construct organ-specific biological age indicators and assess 
#            their associations with lifestyle factors and chronic diseases
# =============================================================================

# =============================================================================
# 0. Environment Setup and Data Loading
# =============================================================================

# Load required R packages
required_packages <- c("compareGroups", "readxl", "dplyr", "tidyr", "ggplot2", 
                      "ggpubr", "rstatix", "correlation", "psych", "mediation",
                      "bruceR", "survival", "openxlsx", "modelsummary")
lapply(required_packages, library, character.only = TRUE)

# Load pre-processed data files
load("./03output/svm_res_0120.Rdata")  # SVM model results
load("./03output/08lifestyle/lifestyle1129.Rdata")  # Lifestyle data
load("./03output/08lifestyle/lifestyle_use0920.Rdata")  # Lifestyle variable labels

# Read disease and baseline data
BaseFollowDiseaseBA <- read.csv("02derived_data/BaseFollowDiseaseBA.csv")
dat2022_1017 <- read.csv("01input_data/dat2022_1017.csv")

# =============================================================================
# 1. Data Preprocessing and Biological Age Calculation
# =============================================================================

# 1.1 Extract SVM results for organ gaps and predictions
extract_svm_results <- function(svm_results) {
  # Extract gap data (4th element) and prediction data (6th element)
  bagaplm <- svm_results[[1]][[4]]  # Gap data including age
  bagap <- svm_results[[1]][[6]]    # Pure gap data without age
  
  return(list(gap_with_age = bagaplm, gap = bagap))
}

# Extract SVM-derived data
svm_data <- extract_svm_results(svm_res_0120)
bagaplm <- svm_data$gap_with_age
bagap <- svm_data$gap

# 1.2 Define organ names for analysis
organ_complete <- c('Cardiovascular', 'Bone', 'Gut', 'Kidney', 'Metabolism', 
                   'Brain', 'Cognition', 'Body')
organ_simple <- c('Cardiac', 'Bone', 'Gut', 'Kidney', 'Metab.', 'Brain', 
                 'Cognition', 'Body')

# 1.3 Calculate predicted biological age (BA) by adding chronological age to gaps
calculate_predicted_ba <- function(bagap_data, organ_names_simple, organ_names_complete) {
  # Predicted BA = Chronological age + Organ-specific gap
  ba_predicted <- bagap_data %>%
    mutate(
      Cardiovascular = Cardiac + age,
      Bone = Bone + age,
      Gut = Gut + age,
      Kidney = Kidney + age,
      Metabolism = Metab. + age,
      Brain = Brain + age,
      Cognition = Cognition + age,
      Body = Body + age
    ) %>%
    select(all_of(c('ID', 'sex', 'age', organ_names_complete)))
  
  return(ba_predicted)
}

ba_predicted <- calculate_predicted_ba(bagap, organ_simple, organ_complete)

# 1.4 Generate descriptive statistics tables
generate_descriptive_tables <- function(data_list, organ_names) {
  # Generate overall descriptive statistics
  ba_descr <- descrTable(~., data = data_list$ba, show.n = TRUE, show.all = TRUE)
  
  # Generate sex-stratified descriptive statistics
  bagap_descr <- descrTable(sex ~ ., data = data_list$bagap, show.n = TRUE, show.all = TRUE)
  bagaplm_descr <- descrTable(sex ~ ., data = data_list$bagaplm, show.n = TRUE, show.all = TRUE)
  
  # Export results to CSV files
  export2csv(ba_descr, file = './06newoutput/16table/ba_descriptive_stats.csv')
  export2csv(bagap_descr, file = './06newoutput/16table/bagap_sex_stratified.csv')
  export2csv(bagaplm_descr, file = './06newoutput/16table/bagaplm_sex_stratified.csv')
  
  return(list(ba = ba_descr, bagap = bagap_descr, bagaplm = bagaplm_descr))
}

# Create data list for descriptive analysis
data_list <- list(ba = ba_predicted, bagap = bagap, bagaplm = bagaplm)
descriptive_tables <- generate_descriptive_tables(data_list, organ_simple)

# =============================================================================
# 2. Biological Age and Chronological Age Correlation Analysis
# =============================================================================

# 2.1 Calculate correlations between biological age and chronological age
calculate_ba_age_correlation <- function(ba_data, organ_names) {
  cor_results <- correlation::correlation(
    ba_data, 
    select = 'age', 
    select2 = organ_names,
    p_adjust = 'none', 
    include_factors = TRUE
  )
  return(cor_results)
}

ba_age_correlation <- calculate_ba_age_correlation(ba_predicted, organ_complete)

# 2.2 Plot biological age-chronological age scatter plots
source("25BA-CA_BioAge.R")  # Load plotting functions

plot_ba_ca_scatter <- function(ba_data, organ_order, filename, width = 13, height = 4) {
  p <- plot_BACA(
    ba_data, 
    organ_order = organ_order,
    nrow = 1, 
    width = width, 
    height = height, 
    legend.position = 'bottom'
  )
  
  ggsave(
    filename = filename, 
    plot = p,
    device = 'pdf',
    width = width, 
    height = height, 
    dpi = 600, 
    units = 'in'
  )
  
  return(p)
}

# Overall scatter plot
filename_total <- './06newoutput/14BACAdot/f2_baca_0311.pdf'
scatter_plot_total <- plot_ba_ca_scatter(ba_predicted, organ_simple, filename_total)

# Sex-stratified scatter plots
ba_man <- ba_predicted %>% filter(sex == 1)    # Male participants
ba_woman <- ba_predicted %>% filter(sex == 2)  # Female participants

filename_man <- './06newoutput/14BACAdot/f2_bacaman_0305.pdf'
scatter_plot_men <- plot_ba_ca_scatter(ba_man, organ_simple, filename_man, height = 3)

filename_woman <- './06newoutput/14BACAdot/f2_bacawoman_0305.pdf'
scatter_plot_women <- plot_ba_ca_scatter(ba_woman, organ_simple, filename_woman, height = 3)

# =============================================================================
# 3. Inter-organ Biological Age Correlation Analysis
# =============================================================================

# 3.1 Load correlation analysis functions
source("D:/work/01衰老/02标志物/20230417/22interplay_func.R")
source("D:/work/01衰老/02标志物/20230417/24corheatmap_Heatmap.R")

# 3.2 Calculate partial correlations between organs
calculate_partial_correlations <- function(data, organ_names, group = NULL) {
  if (!is.null(group)) {
    data <- data %>% filter(sex == group)
  }
  
  # Calculate partial correlations between organ gaps
  gap_pcor <- ba_gap_pcor_func(data, organBA_name = organ_names)
  
  # Calculate BA by adding chronological age to gaps
  ba_data <- data %>%
    mutate(
      Cardiac = Cardiac + age,
      Bone = Bone + age,
      Gut = Gut + age,
      Brain = Brain + age,
      Kidney = Kidney + age,
      Metab. = Metab. + age,
      Cognition = Cognition + age,
      Body = Body + age
    ) %>%
    select(all_of(c('ID', 'sex', 'age', organ_names)))
  
  # Calculate partial correlations between biological ages
  ba_pcor <- ba_pcor_func(ba_data, organ_name = organ_names)
  
  return(list(gap_pcor = gap_pcor, ba_pcor = ba_pcor))
}

# Calculate correlations for different groups
correlations <- list(
  total = calculate_partial_correlations(bagaplm, organ_simple),
  men = calculate_partial_correlations(bagaplm, organ_simple, group = 1),
  women = calculate_partial_correlations(bagaplm, organ_simple, group = 2)
)

# 3.3 Standardize and export correlation results
standardize_correlation_results <- function(pcor_data) {
  pcor_data <- pcor_data %>% 
    select(v1, v2, n, estimate, p.value, p.fdr) %>%
    rename(
      organ1 = v1, 
      organ2 = v2, 
      n = n, 
      r = estimate, 
      p.value = p.value, 
      p.fdr = p.fdr
    )
  
  pcor_data$r <- sprintf("%.3f", pcor_data$r)
  pcor_data$p.value <- format(signif(pcor_data$p.value, digits = 4), scientific = TRUE)
  pcor_data$p.fdr <- format(signif(pcor_data$p.fdr, digits = 4), scientific = TRUE)
  
  return(pcor_data)
}

# Standardize all correlation results
standardized_results <- list()
for (group_name in names(correlations)) {
  standardized_results[[paste0("ba_", group_name)]] <- 
    standardize_correlation_results(correlations[[group_name]]$ba_pcor)
  standardized_results[[paste0("gap_", group_name)]] <- 
    standardize_correlation_results(correlations[[group_name]]$gap_pcor)
}

# Export standardized results to CSV files
output_files <- c(
  "f2_ba_pcor_0120.csv", "f2_ba_pcor_woman_0120.csv", "f2_ba_pcor_man_0120.csv",
  "f2_bagaplm_pcor_0120.csv", "f2_bagaplm_woman_pcor_0120.csv", "f2_bagaplm_man_pcor_0120.csv"
)

for (i in seq_along(standardized_results)) {
  write.csv(standardized_results[[i]], 
            file = paste0('./06newoutput/16table/', names(standardized_results)[i], '.csv'), 
            row.names = FALSE)
}

# 3.4 Plot correlation heatmaps
plot_correlation_heatmap <- function(pcor_data, filename, organ_order) {
  cor_plot_withPvalue_heatmap(
    pcor_data[, c(1, 2, 3, 9)], 
    filename = filename,
    organ_order = organ_order
  )
}

# Create heatmaps for different groups
heatmap_files <- list(
  total = './06newoutput/14BACAdot/baheat0120.pdf',
  men = './06newoutput/14BACAdot/ba_man_heat0120.pdf',
  women = './06newoutput/14BACAdot/ba_woman_heat0120.pdf'
)

for (group in names(correlations)) {
  plot_correlation_heatmap(
    correlations[[group]]$ba_pcor, 
    heatmap_files[[group]], 
    organ_simple
  )
}

# 3.5 Construct correlation network plots
plot_correlation_network <- function(pcor_data, filename) {
  ba_cor_net_func(pcor_data, filename = filename)
}

# Create network plots
network_files <- list(
  total = './06newoutput/14BACAdot/banet0120.pdf',
  women = './06newoutput/14BACAdot/banet_woman_1128.pdf',
  gap = './06newoutput/14BACAdot/baGaplmNet0310.pdf'
)

plot_correlation_network(correlations$total$ba_pcor, network_files$total)
plot_correlation_network(correlations$women$ba_pcor, network_files$women)
plot_correlation_network(correlations$total$gap_pcor, network_files$gap)

# =============================================================================
# 4. Biological Age and Lifestyle Factors Correlation Analysis
# =============================================================================

# 4.1 Prepare lifestyle data for analysis
prepare_lifestyle_data <- function(lifestyle_data, ba_gap_data, svm_results) {
  # Handle family income calculation
  lifestyle_data$family_num[lifestyle_data$family_num == 0] <- NA
  lifestyle_data$income_average <- ifelse(
    is.na(lifestyle_data$total_income) | is.na(lifestyle_data$family_num), 
    NA, 
    lifestyle_data$total_income / lifestyle_data$family_num
  )
  
  # Define lifestyle variables for analysis
  life_vars <- c("edu_y", "family_num", 'income_average', "smoke", 'drink', 'tea',
                'exercise_frequency', 'sleep_cate', 'hPDI', 'uPDI', 'PDI', 'TUG', 
                'Tinetti', 'Identification', 'menophania_age', 'menses_regular', 
                'menses_interval1', 'menopause_age')
  
  # Create display names for lifestyle variables
  life_display_names <- c("Education years", "No. of family members", 'Per capita income',
                         'Smoke', 'Alcohol', 'Tea', 'Exercise frequency', 'Sleep score',
                         'hPDI', 'uPDI', 'PDI', 'TUG', 'Tinetti', 'Olfactory identification',
                         'Age at menarche', 'Menses regular', 'Menstrual cycle length', 
                         'Age at menopause')
  
  # Merge biological age data with lifestyle data
  ba_lifestyle <- left_join(svm_results[[1]][[4]], lifestyle_data)
  
  # Clean and rename columns
  colnames(ba_lifestyle) <- gsub("\nGap", "", colnames(ba_lifestyle), fixed = TRUE)
  colnames(ba_lifestyle) <- multireplace(colnames(ba_lifestyle), organ_simple, organ_complete)
  colnames(ba_lifestyle) <- multireplace(colnames(ba_lifestyle), life_vars, life_display_names)
  
  return(list(ba_lifestyle = ba_lifestyle, life_vars = life_vars, 
              life_display_names = life_display_names))
}

lifestyle_results <- prepare_lifestyle_data(lifestyle, bagaplm, svm_res_0121)
ba_lifestyle <- lifestyle_results$ba_lifestyle
life_display_names <- lifestyle_results$life_display_names

# 4.2 Create sex-stratified datasets
ba_lifestyle_men <- ba_lifestyle %>% filter(sex == '1')    # Male participants
ba_lifestyle_women <- ba_lifestyle %>% filter(sex == '2')  # Female participants

# 4.3 Calculate partial correlations between biological age and lifestyle factors
calculate_ba_lifestyle_correlation <- function(ba_lifestyle_data, organ_names, 
                                             life_names, filename) {
  
  # Load disease data for covariate adjustment
  disease_fac <- read.csv('./06newoutput/12violinplot/disease_fac_1129.csv')
  disease_num <- disease_fac[, c('ID', 'MetabolicHealthsum')]
  disease_num$MetabolicHealthsum <- as.numeric(disease_num$MetabolicHealthsum)
  
  # Calculate partial correlations adjusting for age and sex
  cor_results <- ba_lifestyle_pcor_plot(
    ba_lifestyle_data,
    organBA_name = organ_names,
    life_new_name = life_names,
    filename = filename,
    filter_p = 'p.fdr',
    t_mat = FALSE,
    pcor_var = c('age', 'sex'),
    width = 5.5,
    height = 6
  )
  
  # Standardize output format
  cor_stand <- cor_results %>% 
    select(Parameter1, Parameter2, n, r, p, p.fdr) %>%
    rename(Agegap = Parameter1, Lifestyle = Parameter2) %>%
    mutate(
      r = sprintf("%.3f", r),
      p = format(signif(p, digits = 4), scientific = TRUE),
      p.fdr = format(signif(p.fdr, digits = 4), scientific = TRUE)
    )
  
  return(cor_stand)
}

# Define reproductive-related labels for lifestyle variables
reproductive_labels <- rep(FALSE, length(life_display_names))
reproductive_labels[(length(reproductive_labels)-3):length(reproductive_labels)] <- TRUE
reproductive_labels <- as.logical(reproductive_labels)

# Calculate lifestyle correlations
filename_cor <- './06newoutput/11lifestyle/0323_gaplm_fdr_long.pdf'
ba_lifestyle_correlations <- calculate_ba_lifestyle_correlation(
  ba_lifestyle, organ_complete, life_display_names, filename_cor
)

# Export correlation results
openxlsx::write.xlsx(ba_lifestyle_correlations, 
                    file = './06newoutput/16table/f3_BALifestyleCor0123.xlsx')

# =============================================================================
# 5. Biological Age and Chronic Disease Association Analysis
# =============================================================================

# 5.1 Prepare and categorize chronic disease data
prepare_disease_data <- function(disease_data) {
  # Select major chronic disease variables
  disease_fac <- disease_data %>% 
    select(ID, Hypertension, Hyperlipidemia, Diabetes, Osteoporosis, CVD, Dementia) %>%
    mutate(
      # Calculate metabolic health composite score
      MetabolicHealthsum = select(., Hypertension, Hyperlipidemia, Diabetes) %>% rowSums(),
      
      # Create binary metabolic health classification
      MetabolicHealth_2fac2 = case_when(
        MetabolicHealthsum == 0 ~ 0,
        MetabolicHealthsum %in% c(1, 2, 3) ~ 1,
        TRUE ~ NA_real_
      ),
      
      # Create ternary metabolic health classification
      MetabolicHealth_3fac = case_when(
        MetabolicHealthsum == 0 ~ 0,
        MetabolicHealthsum == 1 ~ 1,
        MetabolicHealthsum %in% c(2, 3) ~ 2,
        TRUE ~ NA_real_
      ),
      
      # Create factor labels for metabolic health
      UnhealthMetabolic = factor(MetabolicHealth_2fac2, levels = c(0, 1), 
                                labels = c('NonUnhealthMetabolic', 'UnhealthMetabolic')),
      MetabolicHealth = factor(MetabolicHealth_3fac, levels = c(0, 1, 2), 
                              labels = c('healthy', 'moderately\nunhealthy', 'highly\nunhealthy')),
      
      # Dichotomize osteoporosis for analysis
      Osteoporosis = ifelse(Osteoporosis > 1, 1, 0),
      Osteoporosis = factor(Osteoporosis, levels = c(0, 1), 
                           labels = c('NonOsteoporosis', 'Osteoporosis')),
      
      # Dichotomize other diseases
      Hyperlipidemia = factor(Hyperlipidemia, levels = c(0, 1), 
                             labels = c('NonHyperlipidemia', 'Hyperlipidemia')),
      Hypertension = factor(Hypertension, levels = c(0, 1), 
                           labels = c('NonHypertension', 'Hypertension')),
      Diabetes = factor(Diabetes, levels = c(0, 1), 
                       labels = c('NonDiabetes', 'Diabetes')),
      CVD = factor(CVD, levels = c(0, 1), labels = c('NonCVD', 'CVD')),
      Dementia = factor(Dementia, levels = c(0, 1), labels = c('NonDementia', 'Dementia'))
    )
  
  # Export processed disease data
  write.csv(disease_fac, file = './06newoutput/12violinplot/disease_fac_1129.csv', row.names = FALSE)
  return(disease_fac)
}

disease_data <- prepare_disease_data(BaseFollowDiseaseBA)

# 5.2 Convert data to long format for statistical testing
convert_to_long_format <- function(ba_gap_data, disease_data, organ_names) {
  # Merge biological age gap data with disease data
  ba_disease <- left_join(ba_gap_data, disease_data)
  
  # Convert to long format for organ-specific analysis
  dat_long <- ba_disease %>% 
    pivot_longer(all_of(organ_names), names_to = "organ", values_to = "measure") %>%
    mutate(organ = factor(organ, levels = organ_names, labels = organ_names))
  
  # Handle outliers in metabolic organ gap (remove values > 5)
  dat_long_test <- dat_long %>%
    mutate(measure = ifelse(organ == 'Metab.' & as.numeric(measure) > 5, NA, measure))
  
  return(dat_long_test)
}

dat_long_test <- convert_to_long_format(bagaplm, disease_data, organ_simple)

# 5.3 Create sex-stratified long format datasets
dat_long_men <- dat_long_test %>% filter(sex == 1)
dat_long_women <- dat_long_test %>% filter(sex == 2)

# 5.4 Statistical testing functions for metabolic health
perform_metabolic_health_tests <- function(dat_long_data) {
  library(ggpubr)
  library(rstatix)
  
  # Set metabolic health grouping levels
  dat_long_data$MetabolicHealth <- factor(
    dat_long_data$MetabolicHealth,
    levels = c('healthy', 'moderately\nunhealthy', 'highly\nunhealthy')
  )
  
  # Wilcoxon rank-sum test across metabolic health groups
  wilcox_test <- compare_means(
    measure ~ MetabolicHealth,
    data = dat_long_data,
    group.by = 'organ',
    method = 'wilcox.test',
    p.adjust.method = 'none',
    alternative = 'greater'
  )
  
  # Calculate effect sizes
  wilcox_effsize_res <- dat_long_data %>%
    group_by(organ) %>%
    wilcox_effsize(
      measure ~ MetabolicHealth,
      p.adjust.method = "fdr",
      alternative = 'less'
    )
  
  # Combine test results
  wilcox_all <- full_join(wilcox_test, wilcox_effsize_res) %>%
    mutate(p.adj.fdr = p.adjust(p, method = 'fdr')) %>%
    select(organ, group1, group2, n1, n2, p, p.adj.fdr, effsize)
  
  # Format output for publication
  wilcox_output <- wilcox_all %>%
    mutate(
      p = format(signif(p, digits = 4), scientific = TRUE),
      p.adj.fdr = format(signif(p.adj.fdr, digits = 4), scientific = TRUE),
      effsize = sprintf("%.3f", effsize)
    )
  
  return(wilcox_output)
}

# Perform metabolic health association tests
metabolic_results <- list(
  total = perform_metabolic_health_tests(dat_long_test),
  men = perform_metabolic_health_tests(dat_long_men),
  women = perform_metabolic_health_tests(dat_long_women)
)

# Export metabolic health results
openxlsx::write.xlsx(metabolic_results, 
                    file = './06newoutput/16table/f4_metabolichealth_results_0120.xlsx')

# 5.5 Multi-disease comparison function
perform_multi_disease_tests <- function(dat_long_data, diseases) {
  wilcox_results <- data.frame()
  
  # Loop through each disease for comparison
  for (disease in diseases) {
    # Create formula for Wilcoxon test
    formula <- as.formula(paste('measure ~', disease))
    
    # Wilcoxon rank-sum test
    wilcox_test <- compare_means(
      formula,
      data = dat_long_data,
      group.by = 'organ',
      method = 'wilcox.test',
      p.adjust.method = 'none'
    )
    
    # Calculate effect sizes
    wilcox_effsize_res <- dat_long_data %>%
      group_by(organ) %>%
      wilcox_effsize(formula, p.adjust.method = "fdr")
    
    # Combine results for current disease
    disease_results <- full_join(wilcox_test, wilcox_effsize_res) %>%
      mutate(p.adj.fdr = p.adjust(p, method = 'fdr')) %>%
      select(organ, group1, group2, n1, n2, p, p.adj.fdr, effsize)
    
    wilcox_results <- rbind(wilcox_results, disease_results)
  }
  
  # Format results for publication
  wilcox_output <- wilcox_results %>%
    mutate(
      p = format(signif(p, digits = 4), scientific = TRUE),
      p.adj.fdr = format(signif(p.adj.fdr, digits = 4), scientific = TRUE),
      effsize = sprintf("%.3f", effsize)
    )
  
  return(wilcox_output)
}

# Define diseases for multi-disease analysis
diseases <- c("UnhealthMetabolic", "Dementia", "Osteoporosis")

# Perform multi-disease association tests
multi_disease_results <- perform_multi_disease_tests(dat_long_test, diseases)

# Export multi-disease results
openxlsx::write.xlsx(multi_disease_results, 
                    file = './06newoutput/16table/f4_multi_disease_results_0120.xlsx')

# 5.6 Create violin plots for disease associations
# Note: Requires violinmerge_func to be sourced from external file
source("violin_plot_functions.R")  # Assumes this contains violinmerge_func

plot_disease_violin <- function(dat_long_data, disease_names, filename, p_label = 'p.adj.label.star') {
  # Create merged violin plots with statistical annotations
  violinmerge_func(dat_long_data, disease_names, filename, plabel = p_label)
}

# Create violin plots for major diseases
disease_names_main <- c("UnhealthMetabolic", 'Osteoporosis', "Dementia")
filename_main <- './06newoutput/12violinplot/violin0120.pdf'
plot_disease_violin(dat_long_test, disease_names_main, filename_main)

# =============================================================================
# 6. Mediation Analysis
# =============================================================================

# 6.1 Prepare mediation analysis dataset
prepare_mediation_data <- function(ba_gap_data, lifestyle_data, disease_data, base_data) {
  # Calculate BMI from height and weight
  base_data$BMI <- base_data$Weight_now / (base_data$Height_now * 0.01)^2
  
  # Select key variables for mediation analysis
  med_data_origin <- base_data %>% 
    select(ID, smoke, Identification, HTD2, DM2, drink, BMI, CVD, sleep_cate, Staple) %>%
    full_join(disease_data[, c('ID', 'CVD')]) %>%
    full_join(lifestyle_data[, c('ID', 'sleep_cate')])
  
  # Merge with biological age gap data
  med_data <- left_join(ba_gap_data, med_data_origin)
  
  # Add metabolic health score
  disease_num <- disease_data[, c('ID', 'MetabolicHealthsum')]
  disease_num$MetabolicHealthsum <- as.numeric(disease_num$MetabolicHealthsum)
  med_data <- full_join(med_data, disease_num)
  
  return(med_data)
}

mediation_data <- prepare_mediation_data(bagaplm, lifestyle, disease_data, dat2022_1017)

# 6.2 General mediation analysis function
perform_mediation_analysis <- function(data, x_var, m_var, y_var, covar_vars) {
  library(mediation)
  library(bruceR)
  
  # Select complete cases for exposure, mediator, and outcome
  complete_cases <- complete.cases(data[, c(x_var, m_var, y_var)])
  complete_data <- data[complete_cases, ] %>% select(all_of(c(x_var, m_var, y_var, covar_vars)))
  
  
  # Construct model formulas
  outcome_formula <- as.formula(paste(y_var, "~", x_var, "+", m_var, "+", paste(covar_vars, collapse = "+")))
  mediator_formula <- as.formula(paste(m_var, "~", x_var, "+", paste(covar_vars, collapse = "+")))
  
  # Fit regression models
  outcome_model <- lm(outcome_formula, data = complete_data)
  mediator_model <- lm(mediator_formula, data = complete_data)
  
  # Perform mediation analysis with bootstrapping
  set.seed(123)  # Ensure reproducibility
  mediation_model <- mediation::mediate(
    mediator_model,
    outcome_model,
    treat = x_var,
    mediator = m_var,
    covariates = covar_vars,
    boot = TRUE
  )
  
  return(mediation_model)
}

# 6.3 Execute mediation analyses for smoking and olfactory identification
mediation_analyses <- list()

# Brain as mediator
brain_mediation <- perform_mediation_analysis(
  mediation_data,
  x_var = "smoke",
  m_var = "Brain",
  y_var = "Identification",
  covar_vars = c('age', 'sex', 'BMI', 'MetabolicHealthsum', 'drink', 'family_num')
)

# Kidney as mediator
kidney_mediation <- perform_mediation_analysis(
  mediation_data,
  x_var = "smoke",
  m_var = "Kidney",
  y_var = "Identification",
  covar_vars = c('age', 'sex', 'BMI', 'MetabolicHealthsum', 'drink', 'family_num', 'sleep_cate')
)

# Save mediation analysis results
save(brain_mediation, kidney_mediation, 
     file = '06newoutput/11lifestyle/mediation_results_0120.Rdata')

# =============================================================================
# 7. Survival Analysis and Prediction Models
# =============================================================================

# 7.1 Load survival analysis functions
source("42predict_HR_func.R")

# 7.2 Prepare survival analysis dataset
prepare_survival_data <- function(ba_data, ba_gap_data, disease_data, covar_data) {
  # Read baseline medication data
  base_med <- read.csv('./02derived_data/base_med.csv')
  
  # Merge biological age data with disease and covariate data
  f20_all_covar <- left_join(round_df_func(ba_gap_data), round_df_func(disease_data))
  f20_all_covar <- left_join(f20_all_covar, covar_data)
  
  # Rename organ columns to full names
  colnames(f20_all_covar) <- multireplace(
    colnames(f20_all_covar), organ_simple, organ_complete
  )
  
  # Define covariates for adjustment
  covar_vars <- c('sex', 'smoke', 'BMI', 'Diabetes', 'Hypertension', 
                 'Hyperlipidemia', 'Hypotensive_Druguse', 'Hypoglycemic_Druguse')
  
  return(list(covar_data = f20_all_covar, covar_vars = covar_vars))
}

survival_data <- prepare_survival_data(ba_predicted, bagaplm, BaseFollowDiseaseBA, cve_covar)
f20_all_covar <- survival_data$covar_data
covar_vars <- survival_data$covar_vars

# 7.3 Cox proportional hazards regression analysis
perform_cox_analysis <- function(data, disease_outcome, organ_names, covar_vars, 
                                complete_cases = TRUE) {
  if (complete_cases) {
    # Complete case analysis
    model_results <- CoxModelFunc_bacovar(
      data = data,
      disease = disease_outcome,
      covar = covar_vars,
      organ_order = organ_names,
      complete_organ = TRUE
    )
  } else {
    # Full cohort analysis
    model_results <- CoxModelFunc_bacovar(
      data = data,
      disease = disease_outcome,
      covar = covar_vars,
      organ_order = organ_names,
      complete_organ = FALSE
    )
  }
  
  return(model_results)
}

# Execute Cox regression analyses
cox_results_complete <- perform_cox_analysis(
  f20_all_covar, 'cve', organ_complete, covar_vars, complete_cases = TRUE
)

cox_results_full <- perform_cox_analysis(
  f20_all_covar, 'cve', organ_complete, covar_vars, complete_cases = FALSE
)

# 7.4 Export Cox regression results
export_cox_results <- function(cox_model, filename_append, include_n = TRUE) {
  CoxOutputFunc(
    cox_model,
    pathname_append = filename_append,
    output = 'default',
    n = include_n
  )
}

# Export results for both complete cases and full cohort
export_cox_results(cox_results_complete, 'cve_complete_cases')
export_cox_results(cox_results_full, 'cve_full_cohort')

# 7.5 Create forest plots for Cox regression results
plot_cox_forest <- function(cox_model) {
  # Extract results table
  cox_table <- CoxOutputFunc(
    cox_model,
    estimate = '{estimate}',
    statistic = c('{conf.low}', '{conf.high}', '{p.value}'),
    n = TRUE,
    output = 'data.frame'
  )
  
  # Select key results for forest plot
  cox_table_selected <- cox_table[, c(1:2, 31:42)]
  
  # Create forest plot
  forest_plot <- bacacox_forest_plot(cox_table_selected)
  
  # Save plot as PDF
  pdf(paste0('06newoutput/17cox/cox_forest_', Sys.Date(), '.pdf'), 
      width = 6.6, height = 7)
  print(forest_plot)
  dev.off()
  
  return(forest_plot)
}

# Generate forest plot
cox_forest_plot <- plot_cox_forest(cox_results_full)

# =============================================================================
# 8. Utility Functions
# =============================================================================

# 8.1 Round numeric variables while preserving integers
round_df_func <- function(df, round_digit = 4) {
  # Identify numeric columns
  numeric_cols <- sapply(df, class) %in% c("numeric")
  round_cols <- numeric_cols
  
  # Determine which numeric columns need rounding (contain decimals)
  for (i in which(numeric_cols)) {
    if (any(df[, i] %% 1 != 0, na.rm = TRUE)) {
      round_cols[i] <- TRUE
    }
  }
  
  # Apply rounding to columns with decimal values
  if (sum(round_cols) > 1) {
    df[, round_cols] <- apply(df[, round_cols, drop = FALSE], 2, 
                             function(x) round(x, round_digit))
  } else if (sum(round_cols) == 1) {
    df[, round_cols] <- round(df[, round_cols], round_digit)
  }
  
  return(df)
}

# 8.2 BMI classification function
classify_bmi <- function(bmi_values) {
  case_when(
    bmi_values < 18.5 ~ 1,    # Underweight
    bmi_values < 24 ~ 2,      # Normal weight
    bmi_values < 28 ~ 3,      # Overweight
    TRUE ~ 4                  # Obese
  )
}

# 8.3 Min-max standardization function
standard_func <- function(x) {
  min_x <- min(x, na.rm = TRUE)
  max_x <- max(x, na.rm = TRUE)
  (x - min_x) / (max_x - min_x)
}

# 8.4 Quantile-based categorization
my_cut <- function(x, n_quantiles = 3) {
  if (n_quantiles == 3) {
    # Tertiles
    q <- quantile(x, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  } else if (n_quantiles == 2) {
    # Median split
    q <- quantile(x, probs = c(0, 1/2, 1), na.rm = TRUE)
  }
  cut(x, breaks = q, include.lowest = TRUE, labels = FALSE)
}

# 8.5 Weighted sum calculation
weight_func <- function(data, select_vars, weights) {
  # Calculate weighted sum of selected variables
  weighted_sum <- as.matrix(data[, select_vars]) %*% weights / sum(weights)
  return(weighted_sum)
}
