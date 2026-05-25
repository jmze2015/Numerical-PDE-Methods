
library(plotly)
library(dplyr)
library(webshot2)
## parameters
M <- 500 ## number of spacial subdivisions
L <- 5 ## length of interval

## x spacing
h <- L / M

## initial distribution
x <- seq(-L/2, (L/2)-h, by = h)
u <- exp(-16*((x)^2)) 

max(abs(u))
## CFL condition for t spacing
k <- 0.8 * h / max(abs(u))
lambda <- k / h

## Simulation Time
T_max <- 150 ## seconds
N <- ceiling(T_max / k)

## Burger's flux function
Flux <- function(x){
  return(x^2 / 2)
}

## Data frame to hold time
results <- data.frame(
  x = x,
  u = u,
  time_step = 0,
  time = 0
)

## General Update For Loop
for(n in 1:N){
  u_new <- u
  
  ## Upwind Flux Rule
  for(i in 2:M){
    u_new[i] <- u[i] - lambda * (Flux(u[i]) - Flux(u[i-1]))
  }
  
  ## Boundary Condition
  u_new[1] <- u[1] - lambda* (Flux(u[1]) - Flux(u[M]))
  
  u <- u_new
  
  ## Rendering Helper
  if(n %% 10 == 0){
    results <- rbind(results, data.frame(
      x = x, 
      u = u,
      time_step = n,
      time = n*k
    ))
  }
}


## full animation
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
      
      # ,list(
      #   type = "line",
      #   x0 = -(L/2),
      #   x1 = (L/2),
      #   y0 = 0.443/5,
      #   y1 = 0.443/5,
      #   line = list(
      #     color = "red",
      #     width = 2, 
      #     dash = "dash"
      #   )
      # )
  
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


## fig ## full animation call


## Screenshot Code
times <- c(0, 0.24, 0.80, 1.28, 2.40)

dir.create("Images", showWarnings = FALSE)

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
      title = paste0("Inviscid Burgers': Upwind, t = ", tt),
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
    "Images/upwind_t_",
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
  
  file.remove(temp.html)
}















