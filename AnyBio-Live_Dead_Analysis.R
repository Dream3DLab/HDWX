#This script imports csv files containing Volume and Dead dye intensity data from single cell segmented imaging data.
#Install all necessary packages and adjust the names to your naming where indicated

library(plyr)
library(readr)
library(dplyr)
library(ggplot2)
library(mclust)
library(pracma)
#Function to read csv files and import the filenames. ATTENTION: The first three 
#rows are skipped as for us, row four contained the column headers. Set another 
#number if your csv file is formatted differently.
read_plus <- function(flnm, filename) {
  data <- read_csv(flnm, skip = 3)
  data$filename <- filename  # Add filename as a new column
  return(data)
}

#Set folder path
working_directory <- "/YouFolderPath"
setwd(working_directory)

# import Volume data, replace "Volume" pattern if you named your files differently.
pat <- "Volume"
files <- list.files(path = working_directory, pattern = pat, full.names = TRUE, recursive = TRUE)
Volume <- ldply(Map(read_plus, files, basename(files)))

# import Mean intensity dead dye. Change "Ch=2" pattern if you named your files differently.
pat <- "Ch=2"
files <- list.files(path = working_directory, pattern = pat, full.names = TRUE, recursive = TRUE)
Red <- ldply(Map(read_plus, files, basename(files)))

# Combine data frames using column names. Change column names if named differently.
Live_Dead <- cbind(Volume[, c("Volume", "ID", "filename")], 
                   Red[, c("Intensity Mean")])

# Rename columns for clarity
colnames(Live_Dead) <- c("Volume", "ID", "Filename", "Intensity_Red")

# Combine Filename and ID into Unique_ID
Live_Dead$Unique_ID <- paste(Live_Dead$Filename, Live_Dead$ID, sep = "_")

#Make extra columns for conditions Casted/Printed and the addition of Vitamin C
Live_Dead$Casted_Printed <- Live_Dead$Filename
Live_Dead$VitaminC <- Live_Dead$Filename

#Replace the column content based on whether the Filename contains the indicator "Casted" or "Printed"
Live_Dead$Casted_Printed <- ifelse(grepl("Casted", Live_Dead$Casted_Printed), "Casted", Live_Dead$Casted_Printed)
Live_Dead$Casted_Printed <- ifelse(grepl("Printed", Live_Dead$Casted_Printed), "Printed", Live_Dead$Casted_Printed)

#Replace the column content based on whether the Filename contains the indicator "withVitC" or "withoutVitC"
Live_Dead$VitaminC <- ifelse(grepl("withVitC", Live_Dead$VitaminC), "Yes", Live_Dead$VitaminC)
Live_Dead$VitaminC <- ifelse(grepl("withoutVitC", Live_Dead$VitaminC), "No", Live_Dead$VitaminC)

#Filter out debris below 1500 volume. Check your smallest cell volumes in your imaging data and change number.
Live_Dead <- Live_Dead %>% filter(Volume > 1500)

# Group by Casted_Printed and VitaminC and calculate quartiles and IQR for each group. Filter out outliers.
Live_Dead_Filtered <- Live_Dead %>%
  group_by(Casted_Printed, VitaminC) %>%
  mutate(
    Q1 = quantile(Intensity_Red, 0.25),
    Q3 = quantile(Intensity_Red, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR
  ) %>%
  ungroup() %>%
  filter(Intensity_Red >= lower_bound & Intensity_Red <= upper_bound) %>%
  select(-Q1, -Q3, -IQR, -lower_bound, -upper_bound)

# Create violin plot of the red dye intensity. You should see clearly distinguishable populations.
ggplot(Live_Dead, aes(x = Casted_Printed, y = Intensity_Red, fill = VitaminC)) +
  geom_violin() +
  labs(title = "Intensity Red by Printed/Casted and VitaminC",
       x = "Printed/Casted",
       y = "Intensity_Red") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  ylim(0, 2500)


# Function to calculate threshold per group based on local minimum after the first peak to separate populations
calculate_group_threshold <- function(data) {
  density_est <- density(data$Intensity_Red, na.rm = TRUE)
  peaks <- findpeaks(density_est$y, minpeakheight = 0.001)
  if (length(peaks[, 1]) > 0) {
    first_peak_index <- peaks[1, 2]
    local_min_indices <- which(diff(sign(diff(density_est$y[first_peak_index:length(density_est$y)]))) == 2) + first_peak_index
    if (length(local_min_indices) > 0) {
      threshold_index <- local_min_indices[1]
      threshold <- density_est$x[threshold_index]
      return(threshold)
    }
  }
  return(NA)
}

# Calculate the threshold for each group
thresholds <- Live_Dead %>%
  group_by(Casted_Printed, VitaminC) %>%
  summarise(threshold = calculate_group_threshold(pick(everything())), .groups = 'drop')

# Merge the thresholds back to the original data
Live_Dead <- Live_Dead %>%
  left_join(thresholds, by = c("Casted_Printed", "VitaminC"))

# Calculate the percentage of total IDs above the threshold for each file (percentage dead cells)
result <- Live_Dead %>%
  group_by(Filename, Casted_Printed, VitaminC, threshold) %>%
  summarise(
    total_ids = n(),
    above_threshold = sum(Intensity_Red > threshold, na.rm = TRUE),
    percentage_above_threshold = (above_threshold / total_ids) * 100,
    .groups = 'drop'
  )

# Calculate the percentage of total IDs below the threshold for each file (viability)
result_viability <- Live_Dead %>%
  group_by(Filename, Casted_Printed, VitaminC, threshold) %>%
  summarise(
    total_ids = n(),
    below_threshold = sum(Intensity_Red > threshold, na.rm = TRUE),
    viability = 100-((below_threshold / total_ids) * 100),
    .groups = 'drop'
  )

# Create the boxplot
ggplot(result_viability, aes(x = Casted_Printed, y = viability, fill = VitaminC)) +
  geom_boxplot() +
  labs(title = "Viability",
       x = "Casted/Printed and VitaminC",
       y = "Viability [%]") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")

# Adding new rows containing viability of cells after cell counting prior to printing/casting
result_viability_all <- add_row(result_viability, 
                                     Casted_Printed = "Prior Casting", VitaminC = "No", viability = 80, 
                                     .before = Inf)

result_viability_all <- add_row(result_viability, 
                                     Casted_Printed = "Prior Casting", VitaminC = "No", viability = 81, 
                                     .before = Inf)

# Ensure Casted_Printed is a factor and set the levels in the desired order
result_viability_all$Casted_Printed <- factor(result_viability_all$Casted_Printed, 
                                                   levels = c("Prior Casting", "Casted", "Printed"))

# Create the boxplot
ggplot(result_viability_all, aes(x = Casted_Printed, y = viability, fill = VitaminC)) +
  geom_boxplot() +
  labs(title = "Viability",
       x = "Prior to processing, Casted/Printed and VitaminC",
       y = "Viability [%]") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") 

# Perform ANOVA to test for significant differences between the groups
anova_result <- aov(viability ~ Casted_Printed * VitaminC, data = filtered_result_viability)
summary(anova_result)

# Check assumptions of ANOVA
par(mfrow = c(2, 2))
plot(anova_result)

# Check the p-value of the main effects
anova_summary <- summary(anova_result)
casted_printed_p_value <- anova_summary[[1]]["Casted_Printed", "Pr(>F)"]
vitaminC_p_value <- anova_summary[[1]]["VitaminC", "Pr(>F)"]

# Perform Tukey's HSD post-hoc test if any main effect is significant
if (vitaminC_p_value < 0.05) {
  tukey_result <- TukeyHSD(anova_result, "VitaminC")
  print(tukey_result)
}