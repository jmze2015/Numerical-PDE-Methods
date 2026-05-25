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

## Lax Flux function
Lax_Flux <- function(u, v){
  f(u) - f(v)
}

## LW rule
A <- function(i){
  if(i == M){
    (u[M]+u[1])/2
  } else{
    (u[i]+u[i+1])/2
  }
}

## Data Frame for storing data
results <- data.frame(
  x = x,
  u = u,
  time_step = 0,
  time = 0
)

## Main Update For Loop for simulation
for(n in 1:N){
  
  u_new <- u
  
  ## Lax-Wendroff Update Rule
  for(i in 2:(M-1)){
    u_new[i] <- (
      u[i] 
      - (lambda / 2)*Lax_Flux(u[i+1], u[i-1]) 
      + (lambda^2 / 2) * (A(i)*Lax_Flux(u[i+1], u[i]) 
                          - A(i-1)*Lax_Flux(u[i], u[i-1]))
    )
  }
  
  ## boundary condition
  u_new[1] <- (
    u[1] 
    - (lambda / 2)*Lax_Flux(u[2], u[M]) 
    + (lambda^2 / 2) * (A(1)*Lax_Flux(u[2], u[1]) 
                        - A(M)*Lax_Flux(u[1], u[M]))
  )
  ## boundary condition
  u_new[M] <- (
    u[M] 
    - (lambda / 2)*Lax_Flux(u[1], u[M-1]) 
    + (lambda^2 / 2) * (A(M)*Lax_Flux(u[1], u[M]) 
                        - A(M-1)*Lax_Flux(u[M], u[M-1]))
  )

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

## animation code
fig <- plot_ly(
  data = results,
  x = ~x,
  y = ~u,
  frame = ~time,
  type = "scatter",
  mode = "lines"
) %>%
  layout(
    title = "Inviscid Burgers': Lax-Wendroff",
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


htmlwidgets::saveWidget(fig, "burgers_lax_wendroff_animation.html")

## Screenshot Code
times <- c(0, 0.24, 0.80, 1.28, 2.40)

dir.create("Images_LW", showWarnings = FALSE)

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
      title = paste0("Inviscid Burgers': Lax-Wendroff, t = ", tt),
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
    "Images_LW/LW_t_",
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




















