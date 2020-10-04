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

#dput-ting R Objects (basically greating a document that then can be read using dget)
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



