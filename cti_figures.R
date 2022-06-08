# ########## cTI Results ################## 

# To create figures from the output of the cTI model (Paper Figures 2-5 - Figure 1 created in matlab)

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(cocor)
library(viridis)

#setwd("C:/Users/Admin/Documents/McGill/PhD_Project/cTI_paper/sensitivity_analysis")
setwd("C:/Users/Admin/Documents/McGill/PhD_Project/cTI_paper/results_may10")

# read in cti table
cti_table <- read.table('cti_table.txt', header = TRUE, sep = ',')

# read in clinical info

clin_path <- 'C:/Users/Admin/Documents/McGill/PhD_Project/genfitable_clinical_neuropsych_bl.csv'
clin_path <- 'C:/Users/Admin/Documents/McGill/PhD_Project/cTI_paper/sensitivity_analysis/genfi_clinical_sensitivity.txt'

clin_table <- read.table(clin_path, header = TRUE, sep = ',')

#extract gene status and age columns from clinical table

#for sensitivity analysis - need to rename columns to match
if (str_detect(clin_path, 'sensitivity')) {
   clin_table <- clin_table %>% 
    rename(Genetic.status.2 = GeneticStatus2,
            Age.at.visit = AgeAtVisit)
  
}
data <- clin_table[ , c("Genetic.status.2", "Age.at.visit")]

#new column for gene status categories
data <- mutate(data, disease.status = case_when(Genetic.status.2 == 'neg' ~ 'non-carrier',
                                                Genetic.status.2 == 'pos - asymp' ~ 'presymptomatic',
                                                Genetic.status.2 == 'pos - symp' ~ 'symptomatic'))


#---- Create Scatterplots of cTI scores vs all clinical/neuropsych/EYO ----

# Scatterplot of cTI pseudotime and disease severity, clinical/neuropsych tests

#cti score columns
pseudotime_names <- names(cti_table)[grepl("cTI_pseudotime" , names(cti_table))]

# get all relevant clinical/neuropsych/eyo columns
index_clinical <- match(c("EYO", "MMSE"), names(cti_table))
clin_columns <- c(index_clinical[1], index_clinical[length(index_clinical)]:(ncol(cti_table)) )

# use only gene carriers

#add gene status to cti table
cti_table$GeneStatus <- data$disease.status

#remove non-carriers (only want carriers for correlations)
cti_table <- cti_table[cti_table$GeneStatus != 'non-carrier', ]

# SCATTERPLOTS

plot_list = list()
counter <- 1

for (clinical_column in clin_columns) {
  
  # want scatterplots for combined modalities
  pseudotime_name <- "cTI_pseudotime_All_factors"
  
  #extract cti scores (all modalities) and clinical/neuropsych/eyo
  table <- cti_table[, c(colnames(cti_table[clinical_column]), pseudotime_name)]
  
  #create scatterplot
  p <- ggplot(table, aes_string(colnames(cti_table[clinical_column]), pseudotime_name)) + 
    geom_point(colour = "black", size = 0.8) +
    labs(x = colnames(cti_table[clinical_column]), y = pseudotime_name) +
    stat_smooth(method = "lm", size = 1, se = FALSE, color = '#0072BD') +
    theme_classic()

  plot_list[[counter]] = p
  counter <- counter + 1
}
#to view plots
plot_list[[14]] #change number for each parameter


#---- Paper Figure 2: associations between cTI scores and neuropsychological tests ----

# Extract and re-label plots from above for combined figure

elements <- theme(axis.title.y = element_text(size = 8, vjust = 0), 
                  axis.title.x = element_text(size = 8), plot.margin =  unit(c(1,0.5,0.5,0.5), "lines"))

dsf_all_factors <- plot_list[[4]] + labs(x = "Digit span forward", y = "cTI disease score") +
  elements

dsb_all_factors <- plot_list[[5]] + labs(x = "Digit span backward", y = "cTI disease score") +
  elements

tmta_all_factors <- plot_list[[6]] + labs(x = "TMTA time", y = "cTI disease score") +
  elements
tmtb_all_factors <- plot_list[[7]] + labs(x = "TMTB time", y = "cTI disease score") +
  elements

digit_symbol_all_factors <- plot_list[[8]] + labs(x = "Digit Symbol", y = "cTI disease score") +
  elements

boston_naming_all_factors <- plot_list[[9]] + labs(x = "Boston Naming", y = "cTI disease score") +
  elements

vf_animals_all_factors <- plot_list[[10]] + labs(x = "VF animals", y = "cTI disease score") +
  elements

block_design_all_factors <- plot_list[[14]] + labs(x = "Block design", y = "cTI disease score") +
  elements


# Make combined figure
figure_all_factors <- ggarrange(dsf_all_factors, dsb_all_factors, tmta_all_factors, 
                                tmtb_all_factors, digit_symbol_all_factors, 
                                boston_naming_all_factors, vf_animals_all_factors,
                                block_design_all_factors, align = "hv",
                                vjust = 1, hjust = 0) #labels="AUTO"
figure_all_factors


#---- Paper Figure 3: cTI scores vs age by gene status ----


#add cti scores (combined modalities) to data table (with gene status and age)
data$pseudotime <- cti_table$cTI_pseudotime_All_factors
data$disease.status <- factor(data$disease.status)
data$disease.status <- relevel(data$disease.status, ref = "non-carrier")

# figure
ggplot(data, aes(x = Age.at.visit, y = pseudotime)) +
  geom_point(aes(color = disease.status), size = 1.1) +
  stat_smooth(method = "lm", size = 1, se = FALSE, color = "black") +
  theme_classic() +
  theme(legend.title = element_blank()) +
  labs(x = 'Age', y = 'cTI disease score') +
  scale_color_viridis(discrete = TRUE, end = 0.7)

#---- Paper Figures 4 and 5: feature contributions ----

#pdf(file = "Figure4_Nov19.pdf", width = 7, height = 7) # defaults to 7 x 7 inches

#read in features table - summed across factors
features_mod <- read.table('features_mod_table_sorted.txt', header = TRUE, sep = ',')
features_mod <- features_mod %>% rename(Factor = Var1) %>% rename(Weight = Var2)

#figure
ggplot(features_mod, aes(reorder(Factor, -Weight), Weight)) +
  geom_col(fill = 'dodgerblue3') + 
  theme(axis.text.x = element_text(size = 10)) +
  theme_classic() +
  theme(axis.title.x = element_blank())

#dev.off()

#read in features table - summed across regions
features_region <- read.table('features_regions_table_sorted.txt', header = TRUE, sep = ',')

#reorganize data
features_region_long <- features_region %>% 
  pivot_longer(cols = starts_with("features_matrix"), names_to = "Factor",
               names_prefix = "features_matrix_", values_to = "Weight") %>% 
  mutate(Factor = case_when(Factor == 1 ~ 'GM density',
                            Factor == 2 ~ 'T1/T2 ratio',
                            Factor == 3 ~ 'fALFF',
                            Factor == 4 ~ 'FA',
                            Factor == 5 ~ 'MD',))
factor_order = c("GM density", "T1/T2 ratio", "fALFF", "FA", "MD")
features_region_long$Factor <- factor(features_region_long$Factor, levels = factor_order)

#pdf(file = "Figure5_Nov19.pdf", width = 11, height = 8) # defaults to 7 x 7 inches

ggplot(features_region_long, aes(reorder(region_names, -Weight), Weight, fill = Factor)) +
  geom_col(width = 0.5) +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(), legend.title = element_blank(),
        legend.position = c(.95, .95), legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))

#dev.off()


#---- Do cTI scores significantly differ by gene status (anovas) ----


# run anova and tukey tests for gene status and cti scores (combined and individual modalities)

#get column names for cTI scores for each factor
pseudotime_names <- names(cti_table)[grepl("cTI_pseudotime" , names(cti_table))]

for (pseudotime_name in pseudotime_names) {
  pseudotime <- cti_table[, pseudotime_name]
  mdl <- aov(pseudotime ~ data$disease.status)
  print(pseudotime_name)
  print(summary(mdl)) # print results of anova
  print(TukeyHSD(mdl))
}

