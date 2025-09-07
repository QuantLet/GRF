rm(list = ls())
library(grf)
library(drf)
library(ggplot2)
library(glue)
library(av)

set.seed(42)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))

dir.create("weights_frames_grf", showWarnings = FALSE)

get_x <- function(n) {
  X <- matrix((-(n / 2 - 1):(n / 2)) / (n / 2), nrow = n)
  colnames(X) <- 'X1'
  return(X)
}

theta_triangle <- function(x, width) {
  pmax(1 - abs(x / width), 0)
}

get_y <- function(X, theta, sigma) {
  n <- nrow(X)
  return(theta(X) + rnorm(n, 0, sigma))
}

# Create transparent theme
transparent_theme <- theme_classic() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

frame_counter <- 0

# Train forests for sigma = 0.1 only (or choose one sigma value)
sig <- 0
rfs <- list()
rf_Xs <- list()

# Generate RF fits for different sample sizes
for (n in c(500, 1000, 2000)) {
  width <- 0.2
  X <- get_x(n)
  theta <- function(X) theta_triangle(X[, 1], width)
  Y <- get_y(X, theta, sig)
  rf <- regression_forest(X, Y)
  
  rfs[[as.character(n)]] <- rf
  rf_Xs[[as.character(n)]] <- X[, 'X1']
}

# Create x values from -1 to 1 with step of 0.02 for the red dot animation
x1s <- matrix(seq(-1, 1, by = 0.02), ncol = 1)
colnames(x1s) <- 'X1'
alpha_list <- lapply(rfs, function(x) get_sample_weights(x, newdata = x1s))

# Create effective weights plots with traveling red dot
for (i in 1:nrow(x1s)) {
  df <- data.frame(
    x = rf_Xs[['500']],
    alpha_500 = alpha_list[['500']][i, ],
    alpha_1000 = alpha_list[['1000']][i, ],
    alpha_2000 = alpha_list[['2000']][i, ]
  )
  
  # Create a separate data frame for the red dot
  red_dot_df <- data.frame(x = x1s[i], y = 0)
  
  p <- ggplot(df, aes(x = x)) +
    geom_line(aes(y = alpha_500), color = "black") +
    geom_line(aes(y = alpha_1000), color = "blue") +
    geom_line(aes(y = alpha_2000), color = "red") +
    # Add red dot at the bottom
    geom_point(data = red_dot_df, aes(x = x, y = y), color = "red", size = 3) +
    labs(
      title = glue("Effective Weights for x₀ = {round(x1s[i], 3)}, σ = {sig}"),
      x = "x", y = expression(alpha[i](x))
    ) +
    ylim(0, 0.1) +
    transparent_theme +
    xlim(-1, 1)  # Set x-axis limits from -1 to 1
  
  frame_counter <- frame_counter + 1
  
  ggsave(
    filename = sprintf("weights_frames_grf/frame_%03d.png", frame_counter),
    plot = p,
    width = 6, height = 4, dpi = 150,
    bg = "transparent"
  )
}