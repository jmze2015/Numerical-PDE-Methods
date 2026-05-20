library(ggplot2)
library(gganimate)
library(dplyr)

# PARAMETERS
Nx <- 200
Nt <- 250

L  <- 2
c  <- 1

dx <- L / Nx
dt <- 0.4 * dx / c   # CFL condition

x <- seq(0, L, length.out = Nx)

# Initial condition: Gaussian pulse
u <- exp(-100 * (x - 0.3)^2)

# Store data
results <- data.frame()

# Time stepping
for (n in 1:Nt) {
  
  # Store current frame
  temp <- data.frame(
    x = x,
    u = u,
    time = n
  )
  
  results <- bind_rows(results, temp)
  
  # Upwind scheme
  u_new <- u
  
  for (i in 2:Nx) {
    u_new[i] <- u[i] - c * dt/dx * (u[i] - u[i-1])
  }
  
  # Periodic BC
  u_new[1] <- u[1] - c * dt/dx * (u[1] - u[Nx])
  
  u <- u_new
}

# Animation
p <- ggplot(results, aes(x, u, group = time)) +
  geom_line(linewidth = 1.2) +
  ylim(0, 1.2) +
  labs(
    title = "Linear Advection Equation",
    subtitle = "Time Step: {closest_state}",
    x = "x",
    y = "u(x,t)"
  ) +
  transition_states(time, transition_length = 1, state_length = 1)

animate(p, fps = 30, width = 800, height = 400)