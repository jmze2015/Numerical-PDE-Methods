
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
k <- 0.5 * h / max(abs(u))
lambda <- k / h

## Simulation Time
T_max <- 15 ## seconds
N <- ceiling(T_max / k)


f <- function(x){
  return(x^2 / 2)
}

Lax_Flux <- function(u, v){
  f(u) - f(v)
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
  
  for(i in 2:(M-1)){
    u_new[i] <- 0.5*(u[i+1]+u[i-1]) - 0.5* lambda * (Lax_Flux(u[i+1], u[i-1]))
  }
  
  ## boundary condition
  u_new[1] <- 0.5*(u[2]+u[M]) - 0.5 * lambda * (Lax_Flux(u[2],u[M]))
  u_new[M] <- 0.5*(u[1]+u[M-1]) - 0.5 * lambda * (Lax_Flux(u[1],u[M-1]))
  
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


fig <- plot_ly(
  data = results,
  x = ~x,
  y = ~u,
  frame = ~time,
  type = "scatter",
  mode = "lines"
) %>%
  layout(
    title = "Inviscid Burgers': Lax-Friedrichs",
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


























