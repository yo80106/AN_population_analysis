install.packages("RCurl")
library(RCurl)
# Add row.names = 1 to the argument list of read.csv. 
# The reason this works is that read.csv is being told that 
# the first column (hence the row.names equals 1) contains the row names.
cereals = read.csv("~/cereals.csv",header=TRUE,row.names = 1)
str(cereals)
head(cereals)

## Centered Data
# "scale" each element by those values by subtracting the mean and dividing by the sd
# If you use  scale(x, scale=FALSE), 
# it will only subtract the mean but not divide by the std deviation.
X = scale(cereals, center = TRUE , scale = FALSE)
dim(X)
colMeans(X)

## Cross-Products among variables
# association matrix of variables 
assoc_variables = t(X) %*% X
dim(assoc_variables)
# first 5 rows and 5 columns 
assoc_variables[1:5,1:5]

## Cross-products between objects
# association matrix of objects 
assoc_objects = X %*% t(X)
dim(assoc_objects)
# first 5 rows and 5 columns 
assoc_objects[1:5, 1:5]

## Inertia
# inertia of variables
sum(diag(assoc_variables))
# inertia of objects
sum(diag(assoc_objects))

## EVD example
EVD = eigen(assoc_variables)
lambda = EVD$values
U = EVD$vectors
# compare to assoc_variable
head(U %*% diag(lambda) %*% t(U))
head(assoc_variables)
