library(neuralnet)
library(caret)
library(plotly)
library(dplyr)

## parameter
a <- 1 ## wave speed

## mesh sizes
M <- 300 ## spatial subdivisions
N <- 200 ## temporal subdivisions
L <- 6
Tmax <- 5 ## seconds

## spacings
dx <- L/M
dt <- Tmax/N

## grid points
x <- seq(-L/2, L/2, length.out = M)
t <- seq(0,Tmax, length.out = N)

## We consider the linear advection formula u_t + a*u_x = 0
## with u(0,x) = g(x) = e^{-2x^2}.

## initial distribution
g <- function(x) exp(-0.5*x^2)

## dataframe necessities
tlist <- numeric(N*M)
for(i in 1:N){
  tlist[((i-1)*M+1):(i*M)] <- t[i]
}
xlist <- rep(x, N)

## Data frame we will use for validation.
truth <- data.frame(t = tlist, x = xlist , u = g(xlist - a*tlist))

# ## Time dependence visualization | not necessary
# viz <- truth[truth$t %in% t[c(1,50,100,150,200)],]
# plot(viz$x, viz$u, col = c(rep("red", 300),
#                            rep("blue", 300),
#                            rep("green", 300),
#                            rep("orange", 300),
#                            rep("black", 300))
#      )

## creating a training data set
training_data <- truth[sample(1:nrow(truth), 3000),]

## Creating Neural Net Model (still only simple loss function)
NN <- neuralnet(u~ t + x, data = training_data, hidden = c(10), threshold = 1.0)
## Regression using data frame
NN_prediction <- neuralnet::compute(NN, truth[,1:2])$net.result



## Truth Figure

fig_truth <- plot_ly(
  data = truth,
  x = ~x,
  y = ~u,
  frame = ~t,
  type = "scatter",
  mode = "lines"
) %>%
  layout(
    title = "Linear Advection: Entropy Solution",
    xaxis = list(title = "x", range = c(-3.1, 3.1)),
    yaxis = list(title = "u(x,t)", range = c(0, 1.1))
    ) %>%
  animation_opts(
    frame = 20,
    transition = 0,
    redraw = FALSE
  ) %>%
  animation_slider(
    currentvalue = list(prefix = "t = ")
  )

fig_truth

## NN Figure

fig_NN <- plot_ly(
  data = data.frame(t=tlist, x=xlist, u = NN_prediction),
  x = ~x,
  y = ~u,
  frame = ~t,
  type = "scatter",
  mode = "lines"
) %>%
  layout(
    title = "Linear Advection: NN Approximation",
    xaxis = list(title = "x", range = c(-3.1, 3.1)),
    yaxis = list(title = "u(x,t)", range = c(0, 1.1))
  ) %>%
  animation_opts(
    frame = 20,
    transition = 0,
    redraw = FALSE
  ) %>%
  animation_slider(
    currentvalue = list(prefix = "t = ")
  )

fig_NN


