library(glmnet)
library("ggplot2")
library("ggmap")
library("maptools")
library("maps")
library("tictoc")
library(ISLR)
options(rgl.printRglwidget = TRUE)
library(rgl)
#############################################################################################################################################
#
# This file is for the cross-validation. As the calculations take around 30 minutes. we decided pto excluded from the markdown-file to not 
# slow down the knitting process too much. The used code snippets are copies of the markdown file
#
#############################################################################################################################################

# calculation of the distance as the main comparison metric
distances <- function(predicted, actual_value) {
  dif <- predicted-actual_value
  dif <- dif * (40030/360) # scaling coordinates to km by the factor circumference (km) / 360°
  mse <- sqrt(dif[,1]^2 + dif[,2]^2)
  return(mse)
}

# importing of data
data <- read.csv("Data/default_plus_chromatic_features_1059_tracks.txt", header=FALSE)
data <- as.data.frame(data)
colnames(data)[117:118] <- c("Latitude", "Longitude")

# defining of training and test data
set.seed(1)
n <-dim(data)[1]
train <- sample(1:n, 0.8*n)
test <- (1:n)[-train]
x <- model.matrix(cbind(Longitude, Latitude)~., data)[,1:116]
y <- data[, c("Latitude","Longitude")]


#########################################################
# new code below 

# cross validation function similiar to the CV lesson in ML 1
cross.validation <- function(grid) {
  cv.err <- rep(NA, length(grid))
  tic("Cross-Validation")
  for (k in 1:length(grid)) {
    
    cv.loo <- rep(NA, n)
    
    for (i in 1:n) {
      # defining data as LOOCV
      loo.x <- x.train[-i,]
      loo.y <- y.train[-i,]
      # do multivariate regression model
      loo.fit <- glmnet(loo.x, loo.y, family = "mgaussian", alpha=0, lambda=grid[k])
      pred.y <- predict(loo.fit, s=grid[k], newx=t(as.matrix(x.train[i,])))
      cv.loo[i] <- distances(pred.y, y.train[i,])
    }
    
    cv.err[k] <- mean(cv.loo)
  }
  toc()
  return (cv.err)
}


# preparations for cross validation
n<-length(train)
x.train <- x[train,]
y.train <- y[train,]

grid.rough <- 10^seq(10, -2, length=100)
grid.fine <- 10^seq(0, -1, length=100)

# doing cross-validation
cv.err.rough <- cross.validation(grid.rough)
cv.err.fine <- cross.validation(grid.fine)

plot(log(grid.rough), cv.err.rough)
plot(log(grid.fine), cv.err.fine)
bestlambda <- grid.fine[which(min(cv.err.fine))]

