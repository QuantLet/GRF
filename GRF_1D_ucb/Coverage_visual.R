install.packages('reshape2')
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(reshape2)

# Data preparation
data <- data.frame(
  Function = rep(c("sin(x)", "sin(x)","sin(8x)", "sin(8x)", "sin(5x)", "sin(3x)", 
                   "max(0, 1-|x_{1}|/0.2)","max(0, 1-|x_{1}|/0.2)","max(0, 1-|x_{1}|/0.2)", "1+x+2x^2+3x^3", "1+x+2x^2+3x^3", "1+x+2x^2+3x^3"), times = 3),
  Bandwidth = rep(c(10, 5, 10, 5, 10, 10, 10, 5, 1, 10, 5, 1), times = 3),
  GridSize = rep(c(200, 100, 50), each = 12),
  Coverage_PW_BS = c(93.39, 92.525, 93.295, 92.425, 93.425, 93.455, 93.35, 92.495, 90.045, 93.4, 92.505, 90.14, 
                     93.46, 92.61, 93.3, 92.37, 93.32, 93.47, 93.36, 92.69, 89.95, 93.42, 92.5, 89.84, 
                     93.32, 92.4, 92.88, 92.2, 92.96, 93.32, 93.32, 92.34, 90.1, 93.34, 92.42, 90.14),
  Coverage_PW_Std_Normal = c(94.665, 94.97, 94.605, 94.885, 94.59, 95.575, 94.375, 94.615, 95.17, 94.595, 94.885, 95.19,
                             94.73, 95.02, 94.57, 94.93, 94.68, 94.53, 94.42, 94.77, 95.1, 94.56, 94.86, 95.2,
                             94.76, 94.92, 94.56, 94.92, 94.74, 94.18, 94.3, 94.48, 95.52, 94.68, 94.84, 95.52),
  Coverage_Uniform_BS = c(97, 87, 92, 88, 99, 93, 96, 88, 86, 98, 87, 86,
                          94, 88, 93, 91, 98, 92, 93, 89, 81, 95, 90, 84,
                          93, 87, 91, 87, 94, 93, 94, 89, 79, 92, 89, 77),
  Coverage_Uniform_Std_Normal = c(87, 74, 73, 72, 84, 73, 85, 72, 79, 80, 79, 77,
                                  82, 79, 81, 79, 82, 80, 84, 71, 76, 84, 79, 80,
                                  88, 76, 85, 85, 85, 87, 88, 79, 80, 82, 83, 79)
)


# Modify the data: Remove GridSize of 1 and exclude sin(3x) and sin(5x)
filtered_data <- subset(data, Bandwidth != 1 & !Function %in% c("sin(3x)", "sin(5x)"))

# Melt the data to long format
long_data <- melt(filtered_data, 
                  id.vars = c("Function", "Bandwidth", "GridSize"), 
                  measure.vars = c("Coverage_Uniform_BS", "Coverage_Uniform_Std_Normal"),
                  variable.name = "Coverage_Type", 
                  value.name = "Coverage")

# Plot the data without offset, with a continuous y-axis and larger axis labels
ggplot(long_data, aes(x = factor(GridSize), y = Coverage, color = as.factor(Bandwidth), shape = Coverage_Type)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "red") +  # Add a dashed line at 95%
  
  facet_wrap(~ Function) +
  labs(title = "Coverage Uniform BS and Standard Normal Across Different Grid Sizes", 
       x = "Grid Size", y = "Coverage (%)", color = "Bandwidth", shape = "Coverage Type") +
  theme_minimal() +
  scale_y_continuous(breaks = seq(floor(min(long_data$Coverage)), ceiling(max(long_data$Coverage)), by = 2)) +  # Adjust y-axis labels
  theme(
    axis.title.x = element_text(size = 14),  # Increase x-axis label size
    axis.title.y = element_text(size = 14),  # Increase y-axis label size
    axis.text.x = element_text(size = 12),   # Increase x-axis tick label size
    axis.text.y = element_text(size = 12)    # Increase y-axis tick label size
  )

