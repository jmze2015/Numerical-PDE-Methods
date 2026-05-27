library(caret)
library(neuralnet)

data(iris)
head(iris)

## our reference on which we train on
reference <- iris$Species

## Notes on NN Function

# neuralnet( typish function
#   formula, ## Y~x1 + \cdots || df-feature~. & call data
#   data, ## data.frame
#   hidden = 1, ## hidden layer size vector
#   threshold = 0.01, 
#   stepmax = 1e5, ## max number of training steps
#   rep = 1,
#   algorithm = "rprop+", ## backpropagation
#   err.fct = "sse", ## typical SSE loss function | can be custom f(x)
#   act.fct = "logistic", ##activation functions tanh, ReLU?, others
#   linear.output = TRUE, ## regression if true, classificatio else 
#   learningrate = NULL,
#   lifesign = "none"
# )
NN <- neuralnet(Species ~., data = iris, hidden = c(5,3), threshold = 0.1)
plot(NN)

regression <- as.data.frame(neuralnet::compute(NN, iris[,1:4])$net.result)
names(regression) <- unique(iris$Species)

##prediction dataframe
pred <- apply(regression, 1, which.max)
pred <- factor(pred, levels = 1:3, labels= levels(iris$Species))
truth <- iris$Species


confusionMatrix(pred, truth)




