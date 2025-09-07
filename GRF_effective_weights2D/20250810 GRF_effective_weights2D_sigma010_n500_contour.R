rm(list=ls())

# ===== Packages =====
library(grf)
library(dplyr)
library(ggplot2)
library(glue)
library(RColorBrewer)

suppressWarnings({
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    pth <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(pth)) setwd(dirname(pth))
  }
})

set.seed(42)

# ===== Helpers =====
get_design <- function(n_x1) {
  expand.grid(
    X1 = seq(-0.5, 0.5, length.out = n_x1),
    X2 = seq(0, 1, by = 0.02)
  )
}
theta_triangle <- function(X, width) pmax(1 - abs(X$X1) / width, 0)
simulate_y <- function(X, width, sigma) theta_triangle(X, width) + rnorm(nrow(X), 0, sigma)
make_dirs <- function(...) { for (d in c(...)) if (!dir.exists(d)) dir.create(d, recursive = TRUE) }

# ===== Plot theme =====
theme_fix <- theme_bw() +
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.background    = element_rect(fill = "transparent", colour = NA),
    panel.background   = element_rect(fill = "transparent", colour = NA),
    legend.background  = element_rect(fill = "transparent", colour = NA),
    legend.key         = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA)
  )

# ===== I/O =====
make_dirs("contour")

# ===== Parameters（sigma=0.10, n=500）=====
n        <- 500
sig      <- 0.10
eta      <- 0.2
x2_star  <- 0.5
x1_path  <- sort(unique(c(seq(-0.5, 0.5, length.out = 61), -0.200, 0.300)), decreasing = FALSE)

X <- get_design(n)
Y <- simulate_y(X, width = eta, sigma = sig)

rf <- regression_forest(
  X, Y,
  num.trees = 2000,
  mtry = 2,
  min.node.size = 5,
  honesty = TRUE,
  sample.fraction = 0.5
)

x1_probe  <- seq(-0.5, 0.5, length.out = 21)
alpha_max <- 0
for (x1p in x1_probe) {
  w_probe <- get_forest_weights(rf, data.frame(X1 = x1p, X2 = x2_star))
  w_probe <- if (is.matrix(w_probe)) as.numeric(w_probe) else as.numeric(w_probe)
  alpha_max <- max(alpha_max, max(w_probe, na.rm = TRUE))
}
brks <- c(0, 0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.06, 0.08, round(alpha_max, 3))

ix <- 0
for (x1 in x1_path) {
  ix <- ix + 1
  
  w <- get_forest_weights(rf, data.frame(X1 = x1, X2 = x2_star))
  w <- if (is.matrix(w)) as.numeric(w) else as.numeric(w)
  
  df_plot <- cbind(as.data.frame(rf$X.orig), alpha = w)
  colnames(df_plot)[1:2] <- c("X1","X2")
  
  p_con <- ggplot(df_plot, aes(X1, X2, z = alpha)) +
    geom_contour_filled(aes(fill = after_stat(level)), breaks = brks) +
    scale_fill_brewer(palette = "Blues", direction = 1,
                      name = expression(alpha[i](x[0])), drop = FALSE) +
    theme_fix +
    annotate("point", x = x1, y = x2_star, pch = 4, colour = "red", size = 2) +
    coord_fixed(xlim = c(-0.5, 0.5), ylim = c(0, 1), expand = FALSE) +
    labs(x = "X1", y = "X2")
  
  ggsave(
    filename = sprintf("contour/frame_%04d.png", ix),
    plot = p_con,
    width = 6, height = 4, dpi = 150,
    bg = "transparent"
  )
}