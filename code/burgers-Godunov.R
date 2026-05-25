library(plotly)
library(dplyr)

## parameters
M <- 500 ## number of spacial subdivisions
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

## Burgers f(u)
f <- function(x){
  return(x^2 / 2)
}

## Godunov Numerical Flux Rule: From Riemann Problem Set up
Gud_Flux <- function(u, v){
  if (u <= v && u > 0){
    f(u)
  } else if(u <= v && v < 0){
    f(v)
  } else if(u <= 0 && 0 <= v){
    0
  } else if(u > v && ((u + v)/2) >= 0){
    f(u)
  } else if(u > v && ((u + v)/2) < 0){
    f(v)
  }
}

## Data frame for storing PDE data
results <- data.frame(
  x = x,
  u = u,
  time_step = 0,
  time = 0
)

## Main scheme update loop
for(n in 1:N){
  
  u_new <- u
  
  for(i in 2:(M-1)){
    u_new[i] <- u[i]-lambda*(Gud_Flux(u[i], u[i+1]) - Gud_Flux(u[i-1], u[i]))
  }
  
  ## boundary condition
  u_new[1] <- u[1] - lambda*(Gud_Flux(u[1], u[2]) - Gud_Flux(u[M], u[1]))
  
  ## boundary condition
  u_new[M] <- u[M] - lambda * (Gud_Flux(u[M], u[1]) - Gud_Flux(u[M-1], u[M]))
  
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
    title = "Inviscid Burgers': Godunov",
    xaxis = list(title = "x", range = c(-2.6, 2.6)),
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

htmlwidgets::saveWidget(fig, "burgers_godunov_animation.html")

## Screenshot Code
times <- c(0, 0.24, 0.80, 1.28, 2.40)

dir.create("Images_Godunov", showWarnings = FALSE)

for (tt in times) {
  
  df_t <- results[results$time == tt, ]
  
  fig_t <- plot_ly(
    data = df_t,
    x = ~x,
    y = ~u,
    type = "scatter",
    mode = "lines"
  ) %>%
    layout(
      title = paste0("Inviscid Burgers': Godunov, t = ", tt),
      xaxis = list(title = "x", range = c(-2.6, 2.6)),
      yaxis = list(title = "u(x,t)", range = c(0, 1.1)),
      shapes = list(
        list(
          type = "line",
          x0 = L/2,
          x1 = L/2,
          y0 = 0,
          y1 = 2.2,
          line = list(color = "blue", width = 2, dash = "dash")
        ),
        list(
          type = "line",
          x0 = -L/2,
          x1 = -L/2,
          y0 = 0,
          y1 = 2.2,
          line = list(color = "orange", width = 2, dash = "dash")
        )
      )
    )
  
  filename <- paste0(
    "Images_Godunov/Godunov_t_",
    gsub("\\.", "p", sprintf("%.2f", tt)),
    ".png"
  )
  
  htmlwidgets::saveWidget(
    fig_t,
    "temp.html",
    selfcontained = TRUE
  )
  
  webshot2::webshot(
    "temp.html",
    file = filename,
    vwidth = 800,
    vheight = 600
  )
  
}

























