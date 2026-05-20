
library(plotly)
library(dplyr)

## parameters

M <- 200 ## number of spacial subdivisions
#N <- 300 ## number of temporal subdivisions
L <- 1 ## length of interval
a <- 3 ## speed of wave


h <- L / M
k <- 0.5 * h / abs(a)
lambda <- a * k / h

T_final <- 1 ## seconds
N <- ceiling(T_final / k)

## initial distribution
x <- seq(0, L-h, by = h)
u <- exp(-50*((x-0.25)^2))

results <- data.frame(
  x = x,
  u = u,
  time_step = 1,
  time = k
)

# results

for(n in 2:N){
  
  u_new <- u
  
  for(i in 2:M){
    u_new[i] <- u[i] - lambda* (u[i] - u[i-1])
  }
  
  ## periodic boundary
  u_new[1] <- u[1] - lambda*(u[1] - u[M])
  
  results <- rbind(results, data.frame(
    x = x,
    u = u_new,
    time_step = n,
    time = n*k
  ))
  
  u <- u_new
  
}

fig <- plot_ly(
  data = results,
  x = ~x,
  y = ~u,
  frame = ~time,
  type = "scatter",
  mode = "lines"
) %>%
  layout(
    title = "linear advection: upwind",
    xaxis = list(title = "x"),
    yaxis = list(title = "u(x,t)", range = c(0, 1.2))
  ) %>%
  animation_opts(
    frame = 10,
    transition = 0,
    redraw = FALSE
  ) %>%
  animation_slider(
    currentvalue = list(prefix = "t = ")
  )

fig














