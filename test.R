
library(KODAMA)

data("MetRef")
data <- MetRef$data
rm(MetRef)

# Remove missing values from the original data
data_clean <- na.omit(data)
rm(data)
# Scale the cleaned data
data_scaled <- scale(data_clean)
rm(data_clean)

# Calculate variances and filter columns based on variance threshold
variances <- apply(data_scaled, 2, var)
threshold <- 0.01  # Adjust threshold as needed
data_filtered <- data_scaled[, variances > threshold]

# Identify columns with any missing data
cols_with_na <- colSums(is.na(data_scaled)) > 0

# Remove columns with missing data
data_no_na_columns <- data_scaled[, !cols_with_na]
rm(data_scaled)
# Remove missing values from filtered data
data_filtered <- na.omit(data_no_na_columns)
rm(data_no_na_columns)

# Check for missing values in the filtered data
if (sum(is.na(data_filtered)) > 0) {
  stop("Missing values are still present in data_filtered")
}

# Install necessary packages if not already installed
install.packages("cluster")
install.packages("factoextra")
install.packages("scales")  # For rescaling data



# Load libraries
library(cluster)
library(factoextra)
library(scales)

# Assuming 'expression' is your matrix of expression values (rows are genes, columns are samples)
# Normalize the expression values
expression_scaled <- t(scale(t(expression)))

# Define the range of clusters to evaluate
range_n_clusters <- 2:4

# Initialize vector to store silhouette scores
silhouette_avg <- numeric(length(range_n_clusters))

# Evaluate silhouette scores for each number of clusters
for (i in seq_along(range_n_clusters)) {
  n_clusters <- range_n_clusters[i]

  # Perform k-means clustering
  kmeans_result <- kmeans(expression_scaled, centers = n_clusters, nstart = 25)

  # Calculate silhouette score
  silhouette_scores <- silhouette(kmeans_result$cluster, dist(expression_scaled))
  silhouette_avg[i] <- mean(silhouette_scores[, 3])
}

# Determine the optimal number of clusters
optimal_n_clusters <- range_n_clusters[which.max(silhouette_avg)]

# Print the optimal number of clusters
print(paste("Optimal number of clusters:", optimal_n_clusters))

# Plot the silhouette scores to visualize
plot(range_n_clusters, silhouette_avg, type = "b", pch = 19,
     xlab = "Number of clusters (k)", ylab = "Average Silhouette Score",
     main = "Optimal Number of Clusters Based on Silhouette Score")
# Load necessary libraries
library(cluster)   # For silhouette score
library(ggplot2)   # For visualization (optional)


data("MetRef")
data <- MetRef$data
# Remove missing values from the original data
data_clean <- na.omit(data)
# Step 1: Prepare the data
data_scaled <- scale(data_clean)  # Normalize the features

# Identify columns with any missing data
cols_with_na <- colSums(is.na(data_scaled)) > 0

# Remove columns with missing data
data_no_na_columns <- data_scaled[, !cols_with_na]
rm(data_scaled)
# Remove missing values from filtered data
data_filtered <- na.omit(data_no_na_columns)
rm(data_no_na_columns)

# Step 2: Perform K-means clustering with k = 2
set.seed(10)  # For reproducibility
kmeans_model <- kmeans(data_filtered, centers = 2, nstart = 10)

# Step 3: Get the cluster assignments
cluster_labels <- kmeans_model$cluster

# Step 4: Visualize clustering (Optional, using first two principal components)
pca_res <- prcomp(data_scaled, scale = TRUE)
pca_data <- data.frame(pca_res$x[, 1:2], Cluster = factor(cluster_labels))

ggplot(pca_data, aes(PC1, PC2, color = Cluster)) +
  geom_point(size = 3) +
  ggtitle("K-means Clustering (k=2)") +
  theme_minimal()

# Step 5: Evaluate the silhouette score for k = 2
sil_score <- silhouette(cluster_labels, dist(data_scaled))
mean_silhouette <- mean(sil_score[, 3])  # Mean silhouette score
cat("Mean silhouette score for k=2:", mean_silhouette, "\n")

# Optional: View the K-means model results
print(kmeans_model)



### 3. Univariate Hypothesis Testing:

Often, the data you are dealing with is a subset (sample) of the complete data (population). Thus, the common question here is:

  -   *Can the findings of the sample be extrapolated to the population? i.e., Is the sample representative of the population, or has the population changed?*

  Such questions are answered using specific hypothesis tests designed to deal with such univariate data-based problems.

a.  **Z Test:** Used for numerical (quantitative) data where the sample size is greater than 30 and the population’s standard deviation is known.

```{r}
# Z Test: Test if sample mean is significantly different from population mean
library(stats)

# Perform Z Test
# z_score <- (mean(sample_data_large) - population_mean) / (population_sd / sqrt(length(sample_data_large)))
# z_score
# p_value_z <- 2 * pnorm(-abs(z_score))  # Two-tailed test
# p_value_z
```

Interpretation: If the p-value is less than the significance level (commonly 0.05), the sample mean is significantly different from the population mean.

b.  **One-Sample t-Test:** Used for numerical (quantitative) data where the sample size is less than 30 or the population’s standard deviation is unknown.

```{r}
# One-Sample t-Test: Test if sample mean is significantly different from population mean
t_test_result <- t.test(sample_data_small, mu = population_mean)
t_test_result
```

Interpretation: The t-test result provides a p-value and confidence interval for the sample mean. A p-value less than 0.05 indicates a significant difference from the population mean.

c.  **Chi-Square Test:** Used with ordinal categorical data

```{r}
# Chi-Square Test: Test the distribution of categorical data
observed_counts <- table(category_data)
expected_counts <- rep(length(category_data) / length(observed_counts), length(observed_counts))

chi_square_result <- chisq.test(observed_counts, p = expected_counts / sum(expected_counts))
chi_square_result
```

Interpretation: The Chi-Square test assesses whether the observed frequencies differ from the expected frequencies. A p-value less than 0.05 suggests a significant difference.

d.  **Kolmogorov-Smirnov Test:** Used with nominal categorical data

```{r}
# Kolmogorov-Smirnov Test: Compare sample distribution to a normal distribution
ks_test_result <- ks.test(sample_data_large, "pnorm", mean = population_mean, sd = population_sd)
ks_test_result
```

Interpretation: The KS test assesses whether the sample follows the specified distribution. A p-value less than 0.05 indicates a significant deviation from the normal distribution.


#### Step 2: Preliminary Checks

Before applying Pearson's correlation, check the assumptions:

i. Normality Check --> Use the Shapiro-Wilk test to assess if the variables are normally distributed.

- Null hypothesis: the data = normally distributed

- Alternative hypothesis: the data = not normally distributed

- If the p-value is less than 0.05, the null hypothesis is rejected

```{r}
# Shapiro-Wilk test for normality
shapiro.test(mtcars$mpg)
shapiro.test(mtcars$wt)
```

<details>

<summary>*Does from data of each of the 2 variables (mpg, wt) follow a normal distribution?*:</summary>

`mpg`: The p-value is 0.1229, which is greater than 0.05. Therefore, we do not reject the null hypothesis. This suggests that the mpg variable does not significantly deviate from a normal distribution.

`wt`: The p-value is 0.09265, which is also greater than 0.05. Thus, we do not reject the null hypothesis. This indicates that the wt variable does not significantly deviate from a normal distribution.

