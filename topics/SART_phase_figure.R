# R code to reconstruct phase-distance diagrams from TomoSAR presentation
# Slide 5: Phase as a relative distance measure

library(ggplot2)
library(grid)
library(gridExtra)
library(patchwork)

# ============================================================================
# LEFT IMAGE: Phasor Diagram + Sinusoidal Wave
# ============================================================================

# --- Part 1: Phasor Diagram (top) ---
create_phasor_diagram <- function() {
  # Set up the plot
  p <- ggplot() + 
    coord_fixed(ratio = 1) +
    xlim(-1.5, 1.5) +
    ylim(-1.5, 1.5) +
    theme_void() +
    theme(plot.margin = margin(10, 10, 10, 10))
  
  # Draw the circle
  circle_data <- data.frame(
    x = cos(seq(0, 2*pi, length.out = 100)),
    y = sin(seq(0, 2*pi, length.out = 100))
  )
  p <- p + geom_path(data = circle_data, aes(x = x, y = y), 
                     linetype = "dotted", size = 0.5)
  
  # Draw axes
  p <- p + 
    geom_segment(aes(x = -1.5, xend = 1.5, y = 0, yend = 0), 
                 arrow = arrow(length = unit(0.2, "cm")), size = 0.3) +
    geom_segment(aes(x = 0, xend = 0, y = -1.5, yend = 1.5), 
                 arrow = arrow(length = unit(0.2, "cm")), size = 0.3)
  
  # Define angle for the phasor (45 degrees = pi/4)
  phi <- pi/4
  
  # Draw E_0 vector (horizontal)
  p <- p + 
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 0), 
                 arrow = arrow(length = unit(0.2, "cm")), 
                 size = 0.8, color = "black")
  
  # Draw rotating vector E (at angle phi)
  E_x <- cos(phi)
  E_y <- sin(phi)
  p <- p + 
    geom_segment(aes(x = 0, xend = E_x, y = 0, yend = E_y), 
                 arrow = arrow(length = unit(0.2, "cm")), 
                 size = 0.8, color = "black")
  
  # Draw angle arc
  arc_angles <- seq(0, phi, length.out = 30)
  arc_data <- data.frame(
    x = 0.3 * cos(arc_angles),
    y = 0.3 * sin(arc_angles)
  )
  p <- p + geom_path(data = arc_data, aes(x = x, y = y), size = 0.4)
  
  # Add labels
  p <- p + 
    annotate("text", x = -1.3, y = -0.1, label = "E[v]", parse = TRUE, size = 4) +
    annotate("text", x = 1.3, y = -0.15, label = "a", size = 4) +
    annotate("text", x = -0.15, y = 1.3, label = "b", size = 4) +
    annotate("text", x = 1.1, y = 0.15, label = "E[0]", parse = TRUE, size = 3.5) +
    annotate("text", x = E_x + 0.15, y = E_y + 0.15, label = "E[0]^{b}", parse = TRUE, size = 3.5) +
    annotate("text", x = 0.4, y = 0.15, label = "phi", parse = TRUE, size = 3.5) +
    annotate("text", x = E_x - 0.2, y = E_y/2 + 0.1, label = "omega*t", parse = TRUE, size = 3)
  
  # Draw vertical dashed line from tip of vector
  p <- p + 
    geom_segment(aes(x = E_x, xend = E_x, y = 0, yend = E_y), 
                 linetype = "dashed", size = 0.3, color = "black")
  
  # Draw horizontal dashed projection
  p <- p + 
    geom_segment(aes(x = 0, xend = E_x, y = E_y, yend = E_y), 
                 linetype = "dashed", size = 0.3, color = "black")
  
  return(p)
}

# --- Part 2: Sinusoidal Wave (bottom) ---
create_sine_wave <- function() {
  # Create wave data
  z <- seq(0, 2.5*pi, length.out = 200)
  wave1 <- sin(z)
  wave2 <- sin(z - pi/4)  # Shifted wave
  
  wave_data <- data.frame(
    z = z,
    wave1 = wave1,
    wave2 = wave2
  )
  
  p <- ggplot(wave_data) + 
    geom_hline(yintercept = 0, size = 0.3) +
    geom_line(aes(x = z, y = wave1), size = 0.8) +
    geom_line(aes(x = z, y = wave2), linetype = "dashed", size = 0.6) +
    scale_x_continuous(
      breaks = c(0, pi/2, pi, 3*pi/2, 2*pi),
      labels = c("", "pi/2", "pi", "3*pi/2", "2*pi"),
      expand = c(0.01, 0)
    ) +
    ylim(-1.5, 1.5) +
    labs(x = "z", y = "") +
    theme_minimal() +
    theme(
      axis.title.x = element_text(hjust = 1, size = 12),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(5, 10, 10, 10)
    )
  
  # Add E_0 label with arrow
  p <- p + 
    annotate("segment", x = -0.3, xend = -0.3, y = -1.2, yend = 1.2,
             arrow = arrow(length = unit(0.15, "cm"), ends = "both"), size = 0.4) +
    annotate("text", x = -0.55, y = 0, label = "E[0]", parse = TRUE, size = 4)
  
  # Add wavelength lambda indicator
  p <- p + 
    annotate("segment", x = 0, xend = 2*pi, y = 1.3, yend = 1.3,
             arrow = arrow(length = unit(0.15, "cm"), ends = "both"), size = 0.4) +
    annotate("text", x = pi, y = 1.45, label = "lambda", parse = TRUE, size = 4)
  
  return(p)
}

# ============================================================================
# RIGHT IMAGE: Phase Difference Between Two Waves
# ============================================================================

create_phase_difference_diagram <- function() {
  # Create wave data
  z <- seq(0, 3*pi, length.out = 300)
  wave_A <- sin(z)
  wave_B <- sin(z - 3*pi/2)  # 270 degree shift (3/4 wavelength)
  
  wave_data <- data.frame(
    z = z,
    wave_A = wave_A,
    wave_B = wave_B
  )
  
  # Find points B1 and B2 for marking
  B1_idx <- which.min(abs(z - pi/2))
  B2_idx <- which.min(abs(z - (pi/2 + 2*pi)))
  
  p <- ggplot(wave_data) + 
    geom_hline(yintercept = 0, size = 0.3) +
    geom_line(aes(x = z, y = wave_A), size = 0.8) +
    geom_line(aes(x = z, y = wave_B), size = 0.8) +
    coord_cartesian(xlim = c(0, 3*pi), ylim = c(-1.3, 1.5)) +
    theme_void() +
    theme(plot.margin = margin(10, 10, 10, 10))
  
  # Mark points A, B1, B2
  p <- p + 
    annotate("point", x = 0, y = 0, size = 2) +
    annotate("text", x = -0.2, y = -0.15, label = "A", size = 5) +
    annotate("point", x = z[B1_idx], y = wave_A[B1_idx], size = 2) +
    annotate("text", x = z[B1_idx] + 0.3, y = wave_A[B1_idx], label = "B[1]", parse = TRUE, size = 4) +
    annotate("point", x = z[B2_idx], y = wave_A[B2_idx], size = 2) +
    annotate("text", x = z[B2_idx] + 0.3, y = wave_A[B2_idx], label = "B[2]", parse = TRUE, size = 4)
  
  # Draw wavelength lambda indicator on top wave
  lambda_start <- 0.5
  lambda_end <- lambda_start + 2*pi
  p <- p + 
    annotate("segment", x = lambda_start, xend = lambda_end, y = 1.3, yend = 1.3,
             arrow = arrow(length = unit(0.15, "cm"), ends = "both"), size = 0.4) +
    annotate("text", x = (lambda_start + lambda_end)/2, y = 1.45, 
             label = "lambda", parse = TRUE, size = 4.5)
  
  # Draw distance d indicator on right side
  d_x <- 2.8*pi
  p <- p + 
    annotate("segment", x = d_x, xend = d_x, 
             y = wave_A[which.min(abs(z - d_x))], 
             yend = wave_B[which.min(abs(z - d_x))],
             arrow = arrow(length = unit(0.15, "cm"), ends = "both"), size = 0.4) +
    annotate("text", x = d_x + 0.25, 
             y = (wave_A[which.min(abs(z - d_x))] + wave_B[which.min(abs(z - d_x))])/2, 
             label = "d", size = 5)
  
  # Add phase difference equations
  eq_x <- 2*pi + 0.5
  eq_y_start <- -0.9
  p <- p + 
    annotate("text", x = eq_x, y = eq_y_start, 
             label = "Delta*phi == (d/lambda)*2*pi", parse = TRUE, 
             size = 3.5, hjust = 0) +
    annotate("text", x = eq_x + 0.3, y = eq_y_start - 0.15, 
             label = "== (3/4)*2*pi", parse = TRUE, 
             size = 3.5, hjust = 0) +
    annotate("text", x = eq_x + 0.3, y = eq_y_start - 0.30, 
             label = "== (3/2)*pi", parse = TRUE, 
             size = 3.5, hjust = 0) +
    annotate("text", x = eq_x + 0.3, y = eq_y_start - 0.45, 
             label = "== 270*degree", parse = TRUE, 
             size = 3.5, hjust = 0)
  
  return(p)
}

# ============================================================================
# Create and Save Both Diagrams
# ============================================================================

# Create left image (phasor + wave combined vertically)
phasor_plot <- create_phasor_diagram()
wave_plot <- create_sine_wave()

left_image <- phasor_plot / wave_plot + 
  plot_layout(heights = c(1, 0.7))

# Create right image (phase difference)
right_image <- create_phase_difference_diagram()

# Save individual plots
ggsave("phase_phasor_and_wave.png", left_image, width = 8, height = 10, dpi = 300, bg = "white")
ggsave("phase_difference.png", right_image, width = 10, height = 6, dpi = 300, bg = "white")

# Display plots
print("Left image (Phasor + Wave):")
print(left_image)

print("\nRight image (Phase Difference):")
print(right_image)

# Optional: Create side-by-side comparison
combined_plot <- left_image | right_image
ggsave("phase_distance_complete.png", combined_plot, width = 16, height = 10, dpi = 300, bg = "white")

cat("\nPlots saved as:\n")
cat("- phase_phasor_and_wave.png (left image)\n")
cat("- phase_difference.png (right image)\n")
cat("- phase_distance_complete.png (both images side by side)\n")