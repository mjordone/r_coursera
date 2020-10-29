# Functions

## Simple sum function
add2 <- function(x,y) {
        x + y
}
add2(3,2) ##example of how it is used.

## Vector function (substract the numbers that are higher than a given number)
above10 <- function(x) {
        use <- x>10 #this return a logical vector of trues and falses, indicating which element of vector x is greater than 10
        x[use]
}

above <- function(x,n) {
        use <- x>n #similar to the previous function, but the n is provided by the user, not fixed in the function
        x[use]
} 

x<-1:20
above(x,12)

above <- function(x,n = 10) { #similar to the previous function, but n is 10 by default (i.e., if not specified to be otherwise)
        use <- x>n 
        x[use]
} 

## dataframe functions (calculating the mean of each column)
columnmean <- function(y, removeNA = TRUE) {
        nc <- ncol(y) # identifies the number of columns in the matrix
        means <- numeric(nc) # creates an empty numeric vector that stores the means of each column
        for(i in 1:nc) {
                means[i] <- mean(y[,i], na.rm = removeNA) #le asinga valor al elemento i en el vector means, de acuerdo al promedio de los elementos en la columna i de la matriz y. Considera o no los NA, de acuerdo a lo especidficado en removeNA
        }
        means # imprime el vector means
}

columnmean(airquality)
columnmean(airquality, FALSE)

### args(function) returns all the possible arguments in the function "function"

##Lazy evaluation
f<-function(a,b) {
        a^2
}
f(2) # function "f" never actually uses the argument b, so calling f(2) will not produce an error, because the 2 gets positionally matched to a.

f <- function(a,b) {
        print(a)
        print(b)
}
f(45) # "45" gets printed before the error was triggered. This is because "b" did not have to be evaluated until after "print(a)". 45 was matched to "a".

#The "..." Argument
myplot <- function(x, y, type = "l", ...) {
        plot(x, y, type = type, ...) # "..." absorbs all the other arguments in the function "plot"
}

#Lexical Scoping
make.power <- function(n) {  ## this function retunrs another funciton as its value (one that creates powers (potencias)). Here, n is a free variable in the "pow" function, but not in the "make.power" function
        pow <- function(x) {
                x^n
        }
        pow
}
cube <- make.power(3)
square <- make.power(2)

##To know what is in a function's environment
ls(environment(cube)) # returns "n" "pow", which are the elements in this function
get("n", envirnoment(cube)) # returns the value of the element "n" (in this case, 3)


# Loops: lapply
x <- list(a = 1:5, b = rnorm(10))
lapply(x,mean)

x <- 1:4
lapply(x, runif, min = 0, max = 10) # the last two arguments are specifications of the "runif" function.

x <- list(a = matrix(1:4, 2, 2), b = matrix(1:6, 3 , 2))
lapply(x, function(lala) lala[,1]) #to extract the first column of each matrix in list "x". This requires creating a function that does not exist outside this application of lapply (that's an anonymous function)

# Loops: sapply
...

##Loops: apply
x <- matrix(rnorm(200), 20, 10)
apply(x, 2, mean) #calculates the mean of each column in the matrix (the "2" indicates the dimension--the margin--over which the function is applied--here, the columns)
apply(x, 1, sum)

###### Alternatives for doing this same thing
rowSums = apply(x, 1, sum)
rowMeans = apply(x, 1, mean)
colSums = ...
colMeans = ...

#### More on apply
x <- matrix(rnorm(200), 20, 10)
apply(x, 1, quantile, probs = c(0.25, 0.75)) #calculates the 25 and 75 percentile of each row. Apply creates a matrix of this.

a <- array(rnorm(2*2*10), c(2,2,10))
apply(a, c(1,2), mean) # calculates the average of all these 2x2 matrices. Here We maintain dimensions 1 and 2, but collapse dim 3.

rowMeans(a, dims = 2) #another way to do teh same


#Loops: mapply
mapply(rep, 1:4, 4:1) #rep is for "repeat", it neads two arguments: what will be repeated and the number of repetitions. This is equivalent to: list (rep(1,4), rep(2,3), rep(3,2), rep(4,1))

noise <- function(n, mean, sd) {
        rnorm(n, mean, sd)
}
mapply(noise, 1:5, 1:5, 2)


# Loops: tapply
x <- c(rnorm(10), runif(10), rnorm(10, 1))
f <- gl(3, 10)
tapply(x,f, mean) # takes the mean of each groups of numbers in x. (the second argument in the formula is an index, identifying the different elements over which the function "mean" must pass over).

tapply(x,f, mean, simplify = FALSE) # Not simplifying results in returning a list instead of an array/matrix.
tapply(x, f, range) # returns two (not one) value per element in x.


#Loops: Split
x <- c(rnorm(10), runif(10), rnorm(10, 1)) 
f <- gl(3, 10)
split (x, f)

lapply(split(x,f), mean)

library(datasets)
head(airquality)

s <- split(airquality, airquality$Month) # Separate the cases in the dataset according to months.
lapply(s, function(x) colMeans(x[,c("Ozone", "Solar.R", "Wind")])) # Again, this requires creating a function that extract the columns.
sapply(s, function(x) colMeans(x[,c("Ozone", "Solar.R", "Wind")])) # Simplify the results of lapply, returning the results in one single matrix (rather than separated vectors)
sapply(s, function(x) colMeans(x[,c("Ozone", "Solar.R", "Wind")], na.rm = TRUE)) # Same as the previous one, but removing NAs 


##### Splitting on more than one level
x <- rnorm(10)
f1 <- gl(2, 5) # E.g., sex (there are two levels and multiple members in each of them--in this case, 5).
f2 <- gl(5, 2) # E.g., race
interaction(f1, f2) # to look at the combination of factors (e.g., sex and race). Combines all the levels of f1 with all the levles of f2 (in this example, generates 10 levels.)

str(split(x, list(f1, f2))) # splitting vector x according to the levels of the two otehr factors (it is not necessary to use the interaction).

str(split(x, list(f1, f2), drop = TRUE)) # The same as before, but dropping empty levels.


# Messages
printmessage <- function(x) {
        if(x > 0) {
                print("x is greater than zero")
        } else {
                print("x is less or equal to zero")
        }
        invisible(x) #Invisible is a function that prevents autoprinting, the object that is returned does not get printed to the console.
}

printmessage(1)
printmessage(NA) # returns an error, because NA is not defined and the function cannot go on.

printmessage2 <- function(x) { # It corrects the error of NAs
        if(is.na(x)) {
                print("x is a missing value!")
        } else if(x > 0) {
                print("x is greater than zero")
        } else {
                print("x is less or equal to zero")
        }
        invisible(x)
}

printmessage2(NA)


# Generating random numbers
dnorm(x, mean = 0, sd = 1, log = FALSE)
pnorm(q, mean = 0, lower.tail = TRUE, log.p = FALSE) # The lower tail is the part of the distribution that goes to the left. If you want to evaluate the upper tail, you have to set lower.tail to FALSE.
qnorm(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
rnorm(n, mean = 0, sd = 1)

x <- rnorm(10)
x <- rnorm(10, 20, 2) #generates 10 numbers with mean 20 and sd 2

set.seed(1) #if I run this and the next function two times in a row I will get the same 5 numbers. An I will get the same 5 numbers anytime y reset the seed to 1 and run rnomr(5). If I don't reset the seed, the next 5 numbers will be different.
rnorm(5)


#### Generating Poisson data
rpois(10, 1)

ppois(4, 2) # Cumulative distribution. Pr(x <= 4)



## Generating Random Numbers from a Linear Model
#Suppose we want to simulate from the following linear model  
#`y=ß0+ß1x+€` 
#Where €~N(0,2^2) (mean=0, sd=2). Assume x~(0,1^2) (mean=0, sd=1), ß0=0.5 and ß1=2

set.seed(20)  
x <- rnorm(100)
e <- rnorm(100, 0, 2)
y <- 0.5 + 2 * x + e
plot(x,y)

#What if x is binary?
set.seed(10)
x <- rbinom(100, 1, 0.5) # the probability of having 1 is 0.5
e <- rnorm(100, 0, 2) #mean = 0 sd = 2
y <- 0.5 + 2 * x + e

summary(y)
plot(x, y)


## Suppose we want to simulate from a poisson model where (output data are cat variables, not continuous variables. Here the error distribution will not be normal, but Poisson):
# Y~ Poisson(µ)
# logµ = ß0 + ß1x
# and ß0 = 0.5 and ß1 = 0.3. We need to use the `rpois` function for this
set.seed(1)
x <- rnorm(100)
log.mu <- 0.5 + 0.3 * x
y <- rpois(100, exp(log.mu)) #to get the mean for the Poisson varaible, we need to exponentiate log.mu

summary(y)
plot(x,y)


# Random Sampling
set.seed(1)
sample(1:10, 4) #sample of 4 numbers between 1 and 10, without replacements (not repeating numbers)
sample(letters, 5) #sample of 5 letters

sample(1:10) #permutation of all the numbers in the vector (because the size of the sample is not specified)
sample(1:10)
sample(1:10, replace = TRUE) #Sample w/replacement (allows repeating numbers)


# Profiling: system.time()
## Elapsed time > user time
system.time(readLines("http://www.jhsph.edu")) # the computer waste time waiting for the website to open, then the elapsed time is greater than the user.

## Elapsed time < user time
hilbert <- function(n) { #the function creates a hilbert type matrix
        i <- 1:n
        1/outer(i - 1, i, "+")
}
x <- hilbert(1000)
system.time(svd(x)) #caluclates the singular value decomposition of this matrix with the "svd" (this function makes use of the acceleerate framwork on the Mac, which is a muti-threaded linear algebra library. So it can take advantage of two different cores of the computer--if it has them.)

## Timing longer expressions
system.time({
        n <- 1000
        r <- numeric(n)
        for (i in 1:n) {
                x <- rnorm(n)
                r[i] <- mean(x)
        }
})


## R Profiler Raw Output
lm(x ~ y)
sample.interval=10000


