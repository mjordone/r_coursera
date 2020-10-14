# See working directory
getwd()

# Instalar paquetes
install.packages("")

#to create objects
x <- 1
x <- 0:6

# to obtain the class of an object
class(x)

# to create vectors of objects (if two objects of different types are put in the same vector, it does not return an error, 
# but coerce the list into the least common denominator--either number or character)
x <- c() ## "c" from "concatenate"?
x <- c(0.5,0.6)     ## numeric
x <- c(TRUE,FALSE)  ## logical
x <- c(T,F)         ## logical
x <- c("a","b","c") ## character
x <- 9:29           ## integer
x <- c(1+0i, 2+4i)  ## complex
x <- vector("numeric", length = 10) ## vectors through specifications

#To coerce a vector into a specific class
as.numeric(x)
as.logical(x)
as.character(x)

#to create lists
x <- list(1,"a",TRUE,1+4i)


#to create matrices
m <- matrix(nrow=2,ncol=3)
m <- matrix(1:6,nrow=2,ncol=3) #gets constructed columnwise (i.e. filling columns first)

# Matrices created vectorwise
m <- 1:10
dim(m) <- c(2,5)

# Matrices by binding columns or binding rows
x <- 1:3
y <- 10:12
cbind(x,y) # creates columns x and y
rbind(x,y) # creates rows x and y

# to create factors
x <- factor(c("yes", "yes", "no", "yes", "no")) # asinga los valores automáticamente
table(x)                                        # para visualizar la matriz donde "yes" y "no" son columnas/variables
unclass(x)                                      # muestra los valores numéricos y las etiquetas por separado 

x <- factor(c("yes", "yes", "no", "yes", "no"),
            levels = c("yes", "no"))            # defino el nivel base (o si no lo define alfabéticamente. Aquí defino que "yes" va antes -es the base level)

# to test whether an object is a missing value
is.na()
is.nan()

# To create data frames
x <- data.frame(foo = 1:4, bar = c(T,T,F,F)) # a data frame with two variables: "foo" and "bar".

# to convert a data frame into a matrix
data.matrix() #it coerce the data frame into a matrix, so some objects can end up being something unexpected.

# to name objects
x <-1:3
names(x) <- c("foo", "bar", "norf")

# to create lists with names
x<- list(a=1,b=2,c=3) #a,g,c are the names

# to name matrices
m <- matrix(1:4, nrow=2, ncol=2)
dimnames(m) <- list(c("a","b"), c("c","d")) # las filas son "a" y "b"; las columnas son "c" y "d". dimnames assign a list, where the first element of the list is the vector of the row names and the second element is the vector of the column names

# To read data
read.table() #e.g., data <- read.table("foo.txt)
read.csv() # identical to read.table but default separator is a comma
readLines() # For reading lines of a text file
source() # For reading in R code files (eg, objects or codes)
dget() # for reading in R code files
load() # for reading in saved workspaces
unserialize() # for reading single R objectis in binary form

# to write data to files
write.table()
writeLines()
dump()
dput()
save()
serialize()

#dput-ting R Objects (basically creating a document that then can be read using dget)
y <- data.frame(a=1,b="a")
dput(y)
Structure(list(a=1,
               b=Structure(1L, .label="a",
                           class = "factor")),
          .Names = c("a","b"), row.names=c(NA,-1L),
          class = "data.frame")
dput(y,file="y.R")
new.y <- dget("y.R")

# Dumping R Objects (difference with dput is that dump can be used on multpl objects. In this example, objects "x" and "y")
x <- "foo"
y <- data.frame(a=1, b="a")
dump(c("x","y"), file="data.R")
rm(x,y)
source("data.R")

#file connections (description is the name of the file. open is a code indicating "r" read only, "w" writing, "a" appending)
str(file)
function (description = "", open = "", blocking = TRUE,
          encoding = getOption("encoding"))
  
  con <- file("foo.txt", "r")
data <- read.csv(con)
close(con) # all this block is the same as data <-read.csv("foo.txt")

#Connections to read parts of a file (here, first, we are telling R to read the first 10 lines of a gz file. writeLines can be used to write lines back to the text)
con <- gzfile("words.gz")
x <- readLines(con,10)

#### same, for url
con <- url("http://www.jhsph.edu", "r") #"r" is for "reading
x <- readLines(con)
head(x)


# Subsetting examples
x <- c("a","b", "c", "c", "d", "a")
x[1] #it returns the first element of the object "x", using a numeric index
x[2]
x[1:4] #it returns elements 1 to 4 of the object "x"

x[x>"a"] #it returns all the elements that are greater than "a", using a logical index (works because in R letters are associated to numbers)
u <- x>"a" #it creates a logical vector saying "True" or "False" for each element in the object
x[u] #it subsets the elements of the "x" vector that are greater than "a"

### Subsetting lists (can use [, [[, or $)
x <- list(foo=1:4,bar=0.6) #two elements of the list: first a sequence from 1 to 4 (named foo), and the second is a number, 0.6 (named bar)
x[1] #returns foo (the sequence from 1 to 4) as a list
x[[1]] # does the same, although returning just the sequence (not as a list)
x$bar # returns the element that is associated with the name bar
x[["bar"]] # same as before
x["bar"] # list with the element bar in it

###to extract multipe elements of a list (use [ )
x <- list(foo=1:4, bar=0.6, baz="hello")
x[c(1,3)]

###example of [[ with computed indices (e.g., when the name of an element is created after the creation of the object)
x <- list(foo=1:4, bar=0.6, baz="hello")
name <- "foo"
x[[name]] # computed index for "foo"
x$name #Error! because the element "name" does not exist
x$foo #No problem. element "foo" does exist

###example of [ with computed indices (e.g., when the name of an element is created after the creation of the object)
x <- list(foo=1:4, bar=0.6, baz="hello")
name<-c("foo","bar")
x[name] #returns foo and bar

#Subsetting nested elements of a list
x <- list(a=list(10,12,14), b=c(3.14,2.81))
x[[c(1,3)]] # substracts the third element of the first element (which is a list)
x[[1]][[3]]
x[[c(2,1)]]

#Subsetting matrices (element is by default returned as a vector of length 1, not as a 1x1 matrix)
x <-matrix(1:6,2,3)
x[1,2] #returns the element in row 1 column 2 
x[1, ] #returns all the elements in row 1 (as a vector)
x[ ,2] #returns all the elements in column 2 (as a vector)
x[1,2, drop = FALSE] # returns element in row 1 column 2 ***as a matrix*** (set to FALSE the option to drop the second dimension)
x[1, , drop = FALSE] # returns all the elements in row 1 ***as a matrix***

#Subsetting with partial matching
x <- list(aardvark = 1:5)
x$a  #With dollar sign
x[["a", exact=FALSE]]  #with double bracket


#Removing NA Values (is.na)
x <- c(1,2,NA,4,NA,5)
bad <-is.na(x) # creats a logical vector where TRUE=NA in x
x[!bad] #substracts all the elements that are not "bad"

####for several vectors (complete.cases)
x <- c(1,2,NA,4,NA,5)
y <- c("a","b",NA,"d",NA,"f")
good <- complete.cases(x,y) #substract all the non-missing values in both vectors **(vectors must be of the same length and have the NAs in the same position!)**
x[good] #returns all the non-missing values in x
y[good] #returns all the non-missing values in y

####for matrices (complete.cases)
airquality[1:6, ]
good <-complete.cases(airquality) #extract the dataframe, discounting the cases where there are one or more NAs.
airquality[good, ][1:6, ] #returns the dataframe without the cases where there are at least one NA.

#Vectorized operations
x<-1:4 ; y<-6:9
x+y #aritmetic operations (returns a vector adding the elements that are in the same position)
x*y #aritmetic operations (returns a vector multiplying the elements that are in the same position)
x/y #aritmetic operations (returns a vector dividing the elements that are in the same position)
x>2 #using logical operations (returns a logical vector showing where the operation is TRUE)
x>=2 #using logical operations (returns a logical vector showing where the operation is TRUE)
y==8 #using logical operations (returns a logical vector showing where the operation is TRUE)

#Vectorized Matrix Operations
x<- matrix(1:4, 2, 2); y <- matrix (rep(10,4), 2,2)
x*y # element-wise multiplication
x/y
x %*% y #true matrix multiplication


# Control Structure *IF-ELSE* (Else is optional)
if(<condition>) {
  ##do something
} else if(condition2) {
  ##do something differnet
} else {
  ## do something else
}

if(x>3) {
  y<-10
} else {
  y<-0
}

y <- if (x>3) {
  10
} else {
  0
}

if (x>3) {
  
} else {
  
}

# Control Structure *FOR LOOPS*
for(i in 1:10) {
  print(i)
}

x <- c("a","b","c","d")
for(i in 1:4) {
  print(x[i])
}

for (i in seq_along(x)) { ##seq-along(x) creates an integer sequence of length x
  print(x[i])
}

for(letter in x) {
  print(letter)
}

for(i in 1:4) print(x[i]) ##omiting the curly braces works if you have one single expression

####Nested FOR loops
x<- matrix(1:6, 2 ,3)
for(i in seq_len(nrow(x))) {
  for(j in seq_len(ncol(x))) {
    print(x[i,j])
  }
}

# Control Structure *WHILE* (caan result in infinite loops if not written properly)
count <-0
while (count < 10) {
  print(count)
  count<-count+1
}

### testing multple condition in a while/if loop
z<-5
while (z>=3 && z<=10) {
  print (z)
  conin <- rbinom(1, 1, 0.5) #creates random number according to the parameters specified
  
  if(coin == 1) { ## random walk
    z<-z+1
  } else {
    z<-z-1
  }
}

# Control Structures *REPEAT*, *NEXT*, *BREAK* (repeat initiates an infinite loop, unleas there is a "break" in the loop)
x0 <- 1
tol <- le-8

repeat {
  x1 <- computeEstimate()
  
  if(abs(x1-x0)<tol) {
    break
  } else {
    x0 <- x1
  }
}

### NEXT (skips the specified iterations)
for(i in 1:100) {
  if(i <= 20) {
    ## Skip the first 20 iterations
    next
  }
  ## Do something here
}

# Dates and Tiems
x <- as.Date("1970-01-01") # coerces dates from a character string
unclass(x) #returns 0, which is the number of dates since 1970-01-01 (the date use as reference to store the dates)

x <- Sys.time() # date and time as of now (on POSIXct. I cannot extract seconds, hours, etc... It is stored just as a long number)
p <- as.POSIXlt(x) # creates date and time in POSIXlt
names(unclass(p)) # Returns the classes of the elements
p$sec #Returns the seconds element

datestring <- c("January 10, 2012 10:40", "December 9, 2011 9:10")
x <- strptime(datestring, "%B %d, %Y %H:%M") #Allows extracting dates and times written in a different format

x <- as.Date("2012-01-01")
y <- strptime("9 Jan 2011 11:34:21", "%d %b %Y %H:%M:%S")
x <- as.POSIXlt(x) #Different classes of date cannot be mixed. They must be taken to the same class.
x-y

y <- as.POIX.ct("2012-10-25 06:00:00", tz = "TMG") #specifies the time zone. Useful to do operations with times (eg add or substract)


# Debugging tools: Traceback
mean(b) #where b doesn't exist, therefore, calling the mean is going to return an error
traceback() #shows all what the last function executed did.

lm(y - x)
traceback()

# Debugging tools: Debug
debug(lm)
lm(y - x) # once here, you can pring n + enter to run the lines. It will stop running lines when the error is found.

# Debugging tools: recover
options(error = recover) #then, whenever there is an error, R will offer some options to dig into the failed function.
read.csv("nosuchfile")
options(error=NULL) # to get back to the normal mode


# str
str(str)
str(lm)

x <- rnorm(100, 2, 4)
str(x)
summary(x)

f <- gl(40, 10) # creating a factor
str(f)

str(mtcars)

m <- matrix(rnorm(100), 10, 10)
str(m)
m[, 1] # the first colum is the same to what str shows

library(datasets)
head(airquality)
s <- split(airquality, airquality$Month)
str(s)
