
library(plotly)
library(dplyr)

## parameters

M <- 500 ## number of spacial subdivisions
#N <- 300 ## number of temporal subdivisions

L <- 5 ## length of interval

## x spacing
h <- L / M

## initial distribution
x <- seq(-L/2, (L/2)-h, by = h)
u <- exp(-16*((x)^2)) 

## CFL condition for t spacing
k <- 0.8 * h / max(abs(u))
lambda <- k / h

## Simulation Time
T_max <- 150 ## seconds
N <- ceiling(T_max / k)



Flux <- function(x){
  return(x^2 / 2)
}

results <- data.frame(
  x = x,
  u = u,
  time_step = 0,
  time = 0
)

# results

for(n in 1:N){
  
  u_new <- u
  
  for(i in 2:M){
    u_new[i] <- u[i] - lambda * (Flux(u[i]) - Flux(u[i-1]))
  }
  
  ## boundary condition
  u_new[1] <- u[1] - lambda* (Flux(u[1]) - Flux(u[M]))
  
  u <- u_new
  
  if(n %% 10 == 0){
    results <- rbind(results, data.frame(
      x = x, 
      u = u,
      time_step = n,
      time = n*k
    ))
  }
  
}

#results

fig <- plot_ly(
  data = results,
  x = ~x,
  y = ~u,
  frame = ~time,
  type = "scatter",
  mode = "lines"
) %>%
  layout(
    title = "Inviscid Burgers': Upwind",
    xaxis = list(title = "x", range = c(-5.1, 5.1)),
    yaxis = list(title = "u(x,t)", range = c(0, 1.1)),
    shapes = list(
      
      list(
      type = "line",
      x0 = L/2,
      x1 = L/2,
      y0 = 0,
      y1 = 2.2,
      line = list(
        color = "blue",
        width = 2, 
        dash = "dash"
        )
      ),
      
      list(
        type = "line",
        x0 = -(L/2),
        x1 = -(L/2),
        y0 = 0,
        y1 = 2.2,
        line = list(
          color = "orange",
          width = 2, 
          dash = "dash"
        )
      ),
      
      list(
        type = "line",
        x0 = -(L/2),
        x1 = (L/2),
        y0 = 0.443/5,
        y1 = 0.443/5,
        line = list(
          color = "red",
          width = 2, 
          dash = "dash"
        )
      )
  
      )
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












