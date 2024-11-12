
library(tidyr)
library(dplyr)

raw_biomarkers <- read.csv("./data/data-import/biomarkers.csv",header = TRUE, sep = ",")
raw_covariates <- read.csv("./data/data-import/covariates.csv",header = TRUE, sep = ",")
###########################
# 0 data treatment
# split data (Biomarker contains two elements)
biomarkers <-raw_biomarkers %>% separate(Biomarker, into = c("PatientID", "time"), sep = "-")
# clean useless columns or rows
biomarkers <- na.omit(biomarkers)
covariates <- raw_covariates[, colSums(is.na(raw_covariates)) < nrow(raw_covariates)]
covariates <- na.omit(covariates)
# rename the columns
colnames(covariates) <- c("PatientID","Age","Gender", "Smoker", "VAS_inclusion", "VAS_12months")

# 1 data description
summary(biomarkers)
summary(covariates)


-------------------------
#task 1 hypothesis testing
  
# 1 data treatment
# data group by VAS
covariates$VAS_Group <- ifelse(covariates$VAS_inclusion >= 6, "High", "Low")
# Merge the datasets based on the PatientID
merged_data <- merge(biomarkers, covariates)
# Filter the data for the time in inclusion
high_VAS_group <- merged_data[merged_data$time == "0weeks" & merged_data$VAS_Group == "High", ]
low_VAS_group <- merged_data[merged_data$time == "0weeks" & merged_data$VAS_Group == "Low", ]

# 2 check whether normal distribution
biomarkers_list <- c("IL.8", "VEGF.A", "OPG", "TGF.beta.1", "IL.6", "CXCL9", "CXCL1", "IL.18", "CSF.1")
for (biomarker in biomarkers_list) {
  high_VAS_values <- high_VAS_group[[biomarker]]
  low_VAS_values <- low_VAS_group[[biomarker]]
  
  print(paste("Normality test for", biomarker))
  print(shapiro.test(high_VAS_values))
  print(shapiro.test(low_VAS_values))
}

# 3 for normal distributed data (IL.18,VEGF.A), check homogeneity of variance
norm_biolmarkers_list <- c("IL.18", "VEGF.A")
for (biomarker in norm_biolmarkers_list) {
  high_VAS_values <- high_VAS_group[[biomarker]]
  low_VAS_values <- low_VAS_group[[biomarker]]
  
  # do the F-test 
  print(paste("Homogeneity of variance test for", biomarker))
  variance_result <- var.test(high_VAS_values, low_VAS_values)$p.value
  print(variance_result)
}


# 4 do the test
# for non-normal distributed data :  the Mann-Whitney U test 
unnorm_biolmarkers_list <- c("IL.8", "OPG", "TGF.beta.1", "IL.6", "CXCL9", "CXCL1", "CSF.1")

# Create a data frame to store the results
results <- data.frame(Biomarker = character(), 
                      V = numeric(), 
                      P_Value = numeric(), 
                      stringsAsFactors = FALSE)

# Loop through each biomarker
for (biomarker in unnorm_biolmarkers_list) {
  # Subset data for high and low VAS groups
  high_VAS_values <- high_VAS_group[[biomarker]]
  low_VAS_values <- low_VAS_group[[biomarker]]
  # Perform the Mann-Whitney U test
  mw_test_result <- wilcox.test(high_VAS_values, low_VAS_values, paired = FALSE)
  
  # Store the results
  results <- rbind(results, data.frame(Biomarker = biomarker, 
                                       V = mw_test_result$statistic, 
                                       P_Value = mw_test_result$p.value))
}

# Print the results
print(results)

# do the t-test for IL.18, VEGF.A
for (biomarker in norm_biolmarkers_list) {
  high_VAS_values <- high_VAS_group[[biomarker]]
  low_VAS_values <- low_VAS_group[[biomarker]]
  
  # Perform the t-test
  t_test_result <- t.test(high_VAS_values, low_VAS_values, var.equal = TRUE)
  
  # Print the results
  print(paste("T-test for", biomarker))
  print(t_test_result)
}


# 4 redo the test with Bonferroni correction
# Define the significance level
alpha_bonferroni <- 0.05 / 9 
for (biomarker in norm_biolmarkers_list) {
  high_VAS_values <- high_VAS_group[[biomarker]]
  low_VAS_values <- low_VAS_group[[biomarker]]
  
  # do t-test
  test_result <- t.test(high_VAS_values, low_VAS_values, var.equal = TRUE)
  cat("T-test for", biomarker, "\n")
  print(test_result)
  
  # use new  significance level to explain the result
  if (test_result$p.value < alpha_bonferroni) {
    cat(biomarker, "the differences vary betwwen low and high VAS group（p < ", alpha_bonferroni, "）\n")
  } 
}

# for unnormal distributed variables use Mann-Whitney U test 
for (biomarker in unnorm_biolmarkers_list) {
  high_VAS_values <- high_VAS_group[[biomarker]]
  low_VAS_values <- low_VAS_group[[biomarker]]
  
  # do Mann-Whitney U test
  test_result <- wilcox.test(high_VAS_values, low_VAS_values)
  cat("Mann-Whitney U test for", biomarker, "\n")
  print(test_result)
  
  # use new  significance level to explain the result
  if (test_result$p.value < alpha_bonferroni) {
    cat(biomarker, "differences vary bewteen low and high VAS group（p < ", alpha_bonferroni, "）\n")
  } 
  
  cat("\n")  
}


# ········································
# task 2 regression
# 1 data treatment
filtered_data <-merged_data %>% filter(time == "0weeks")
# Set a random seed for reproducibility
set.seed(66)
# create a train set for 80% of the data
train_set <- sample(1:nrow(filtered_data), size = 0.8 * nrow(filtered_data))

# Split the data into train  and test sets
train_data <- filtered_data[train_set, ]
test_data <- filtered_data[-train_set, ]

# 2 feature selection
# check the normal distribution of each features
biomarkers_list <- c("IL.8", "VEGF.A", "OPG", "TGF.beta.1", "IL.6", "CXCL9", "CXCL1", "IL.18", "CSF.1", "VAS_inclusion", "Age","Smoker","Gender")
for (biomarker in biomarkers_list) {
  print(paste("Normality test for", biomarker))
  if (biomarker == "Smoker" || biomarker == "Gender") {
    shapiro_result <- shapiro.test(train_data$VAS_12months[train_data[[biomarker]] == "1"])
    shapiro_result <-shapiro.test(train_data$VAS_12months[train_data[[biomarker]] == "2"])
  } else {
    shapiro_result <-shapiro.test(train_data[[biomarker]])
  }
  print(shapiro_result)
}
# check the correlation
  #for continuous variables
continuous_vars <- c("IL.8", "VEGF.A", "OPG", "TGF.beta.1", "IL.6", "CXCL9", "CXCL1", "IL.18", "CSF.1", "VAS_inclusion", "Age")
for (var in continuous_vars) {
  if (var == "VEGF.A" || var == "IL.18") { # features satisfy the normal distribution
    correlations <- cor(train_data[continuous_vars], train_data$VAS_12months, method = "pearson")
  } else { # features do not satisfy the normal distribution
    correlations <- cor(train_data[continuous_vars], train_data$VAS_12months, method = "spearman")
  }
}
print(correlations)
  # plot Correlation  heatmap
continuous_vars <- train_data[, c("IL.8", "VEGF.A", "OPG", "TGF.beta.1", "IL.6", "CXCL9", "CXCL1", "IL.18", "CSF.1", "VAS_inclusion", "Age")]
cor_matrix <- cor(continuous_vars, use = "complete.obs")
print(cor_matrix)
library(corrplot)
corrplot(cor_matrix, method = "color", type = "lower", order = "hclust", 
         tl.col = "black", tl.srt = 45)

# for discrete variables
cor_test_gender_vas <- cor.test(as.numeric(train_data$Gender), train_data$VAS_12months, method = "spearman")
cor_test_smoker_vas <- cor.test(as.numeric(train_data$Smoker), train_data$VAS_12months, method = "spearman")
print(cor_test_gender_vas)
print(cor_test_smoker_vas)

# 3 Define the regression formula
formula1 <- VAS_12months ~  IL.8+ VEGF.A+ OPG+ TGF.beta.1+ IL.6+ CXCL9+ CXCL1+ IL.18+ CSF.1+ VAS_inclusion+ Age
formula2 <- VAS_12months ~  0+ VAS_inclusion +  IL.6 + IL.8 + VEGF.A + OPG  + TGF.beta.1+ CSF.1+ IL.18
# Fit the model
model1 <- lm(formula1, data = train_data)
model2 <- lm(formula2, data = train_data)
# Summary of the model 
summary(model1)
summary(model2)

# 4 Predict the VAS_12months values using the test data
predictions1 <- predict(model1, newdata = test_data)
predictions2 <- predict(model2, newdata = test_data)
# Create a data frame for comparison
comparison1 <- data.frame(
  Actual_VAS = test_data$VAS_12months,
  Predicted_VAS = predictions1
)
comparison2 <- data.frame(
  Actual_VAS = test_data$VAS_12months,
  Predicted_VAS = predictions2
)
# Calculate error metrics
comparison1$Absolute_Error <- abs(comparison1$Actual_VAS - comparison1$Predicted_VAS)
comparison1$Squared_Error <- (comparison1$Actual_VAS - comparison1$Predicted_VAS)^2
comparison2$Absolute_Error <- abs(comparison2$Actual_VAS - comparison2$Predicted_VAS)
comparison2$Squared_Error <- (comparison2$Actual_VAS - comparison2$Predicted_VAS)^2
# Calculate Mean Absolute Error (MAE) and Root Mean Squared Error (RMSE) and Mean Squared Error (MSE)
MAE1 <- mean(comparison1$Absolute_Error)
RMSE1 <- sqrt(mean(comparison1$Squared_Error))
MSE1 <- mean(comparison1$Squared_Error)
MAE2 <- mean(comparison2$Absolute_Error)
RMSE2 <- sqrt(mean(comparison2$Squared_Error))
MSE2 <- mean(comparison2$Squared_Error)

# Display comparison and metrics
print(comparison1)
cat("Model 1 - MAE:", MAE1, "\nModel 1 - MSE:", MSE1, "\nModel 1 - RMSE:", RMSE1, "\n")
print(comparison2)
cat("Model 2 - MAE:", MAE2, "\nModel 2 - MSE:", MSE2, "\nModel 2 - RMSE:", RMSE2, "\n")


# Residuals for Model 1
residuals1 <- comparison1$Actual_VAS - comparison1$Predicted_VAS
plot(comparison1$Predicted_VAS, residuals1,
     main = "Residuals vs. Predicted Values for Model 1",
     xlab = "Predicted VAS (Model 1)", ylab = "Residuals")
abline(h = 0, col = "red")

# Residuals for Model 2
residuals2 <- comparison2$Actual_VAS - comparison2$Predicted_VAS
plot(comparison2$Predicted_VAS, residuals2,
     main = "Residuals vs. Predicted Values for Model 2",
     xlab = "Predicted VAS (Model 2)", ylab = "Residuals")
abline(h = 0, col = "red")
