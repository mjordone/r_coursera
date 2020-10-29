# Get/set your working directory
getwd()
setwd()

setwd("../") #moves the directory one directory up

# Checking for and creating directories
file.exists("directoryName") #Will check to see if the directory exists (as a subdirectory within the directory where we are). Returns TRUE or FALSE
dir.create("directoryName") #Will create a directory if it doesn't exist

if (!file.exist("data")) {
  dir.create("data")
}


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
read.table() #e.g., data <- read.table("foo.txt"). You can turn the "sep" argument to commas and header to TRUE to read csv
read.csv() # identical to read.table but default separator is a comma
readLines() # For reading lines of a text file
source() # For reading in R code files (eg, objects or codes)
dget() # for reading in R code files
load() # for reading in saved workspaces
unserialize() # for reading single R objectis in binary form
read.xlsx() # for excel files (requres xlsx package). requires specifiying the sheet and the header. Can be set to read specific rows and columns (colIndex, rowIndex)

# to write data to files
write.table()
writeLines()
dump()
dput()
save()
serialize()
write.xlsx()

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

x[x$var1 <= 3 & x$var3 >11,] #subsetting based on logical commands. Here "var1" and "var3" are the names of the first and third columns.
x[x$var1 <= 3 | x$var3 >15,] #same idea using an "or" instead of an "and"

#Subsetting with partial matching
x <- list(aardvark = 1:5)
x$a  #With dollar sign
x[["a", exact=FALSE]]  #with double bracket


#Removing NA Values (is.na)
x <- c(1,2,NA,4,NA,5)
bad <-is.na(x) # creats a logical vector where TRUE=NA in x
x[!bad] #substracts all the elements that are not "bad"

#######for several vectors (complete.cases)
x <- c(1,2,NA,4,NA,5)
y <- c("a","b",NA,"d",NA,"f")
good <- complete.cases(x,y) #substract all the non-missing values in both vectors **(vectors must be of the same length and have the NAs in the same position!)**
x[good] #returns all the non-missing values in x
y[good] #returns all the non-missing values in y

#######for matrices
airquality[1:6, ]
good <-complete.cases(airquality) #extract the dataframe, discounting the cases where there are one or more NAs.
airquality[good, ][1:6, ] #returns the dataframe without the cases where there are at least one NA.

x[which(x$var2 > 8),] #returns the indices where var2 is greater than 8. When you do that, it doesn't return the NAs


#Sorting
sort(x$var1)
sort(x$var1, decreasing=TRUE)
sort(x$var2, na.last=TRUE) #puts the NAs at the end

#Ordering
x[order(x$var1),] #ordering a data frame by a particular variable (it is applied to the rows)
x[order(x$var1, x$var3),] #same, using multiple variables. The first variable mention is the primary sorting criteria

library(plyr)
arrange(x, var1) #same, using plyr package
arrange(x,desc(var1)) #same, but decreasing order

#Adding rows and columns
x$var4 <- rnorm(5) #creates var4 in dataframe 4 (here, it creates random numbers)
y <- cbind(x,rnorm(5)) #same, using the cbind (column bind) commnad. Adds the column at the right side of x. If we change the order, the new column is added to the left side.

y <- rbind(x,rnorm()) #adds a row (row bind)


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


###################################################################################################################

## Download data from internet
fileUrl <- "https://data.baltimorecity.gov/api/views/2weq-566u/rows.csv?accessType=DOWNLOAD"
download.file(fileUrl, destfile = "cameras.csv", method = "curl") # Method is "curl" because it is an https website. Note that the name and the extension of the file are included in destfile (Here I am just saving the file in the dw, but I can be more specific about where i want it to be). download.file is agnostic to the type of file (it simply downlad and save the file), but requires that extension in destfile corresponds to the actual extension of the file.

dateDownloaded <- date() #to keep track of the date of download
dateDownloaded


## Reading XML
library(XML)
library(RCurl)
fileUrl <- "https://www.w3schools.com/xml/simple.xml"
xData <- getURL(fileUrl)
doc <- xmlParse(xData, useInternal=TRUE)

rootNode <- xmlRoot(doc) #xmlRoot is a wrapper for the entire document
xmlName(rootNode) #provides the name for the root node
names(rootNode) # provides all the nested elements within that root node
rootNode[1] #Este es similar al siguiente, al parecer un nivel más arriba, pero entrega poca información adicional
rootNode[[1]] # provides the first element (like in a list)
rootNode[[1]][[1]] # provides the first element nested in the first element (which is like a list)

xmlSApply(rootNode, xmlValue) #Programmatically extract parts of the file. It loops through the all the elements of the root node and, in this case, implement the "xmlValue" function (xmlValue returns the values)

#############Get the items on the menu and prices
xpathSApply(rootNode,"//name",xmlValue)
xpathSApply(rootNode,"//price",xmlValue)


########Another example

fileUrl <- "https://www.espn.com/nfl/team/_/name/bal/baltimore-ravens"
xData <- getURL(fileUrl)
doc <- htmlTreeParse(xData,useInternal=TRUE) #Here we are parsing an HTML file (not an XML)

scores <- xpathSApply(doc,"//li[@class='score']",xmlValue) #"li" stands for "list"
teams <- xpathSApply(doc,"//li[@class='team-name']",xmlValue)


##JSON
library(jsonlite)
jsonData <- fromJSON("https://api.github.com/users/jtleek/repos")
names(jsonData) #provide all the top level variables
names(jsonData$owner) #access one column of the dataframe
jsonData$owner$login

myjson <- toJSON(iris, pretty = TRUE) #transforming one dataset into JSON
cat(myjson) #to print it out

iris2 <- fromJSON(myjson) #from JSON to dataframe. It is the same function we used for the URL
head(iris2)

#data.table 
library(data.table)
DF = data.frame(x=rnorm(9),y=rep(c("a","b","c"), each=3),z=rnorm(9))
head(DF,3)
DT = data.table(x=rnorm(9),y=rep(c("a","b","c"), each=3),z=rnorm(9)) #creating a databframe using data.table (very similar to what we do when we use data.frame)
head(DT,3)

tables() #to see all of that data tables in memory ("tables", with final "s")

DT[2,] # to subset rows
DT[DT$y="a",] #subsets all cases where variable "y" is equal to "a"
DT[c(2,3)] #subsets the second and third row of that table
DT[,c(2,3)] #this does not subsets columns, as could be expected

####subsetting columns in data.table
{
  x = 1
  y = 2
}
k = {print(10);5}

DT[, list(mean(x), sum(z))] #we can pass a list of functions that are applied to vierables named by columns (here, variables "x" and "z").
DT[,table(y)] #same idea
DT[,w:=z^2] #adding new columns (here, new column "w" which is equal to the "z" variable squared). When we add a new variable to a dataframe R copy over the entire dataframe and add a new variable to it (we get two copies of the dataframe in memory). OJO: when copying data tables, changes made to the original are transferred to the copy

#####Mulitple operations
DT[,m:={tmp <- (x+z); log2(tmp+5)}] #Creating variable M with two statements. Each statement is separated by a smicolon. First statement assigns the temporary variable by adding the variables X and Z. The second statement takes the log base two of that temporary variable plus five

#####plyr like operations
DT[,a:=x>0] # Creates variable A. Values are TRUE or FALSE according to x>0
DT[,b:= mean(x+w), by=a] #creates variable b taking the mean of x+w, and is going to set that mean for all the cases where a=TRUE and the mean for all the cases where a=FALSE

####Special variables
set.seed(123);
DT <- data.table(x=sample(letters[1:3], 1E5, TRUE))
DT[, .N, by=x] #.N its a containment number of times that a particular group appears. Here, it counts the number of times grouped by the x variable 

####Keys (allows sorting and subsetting the dataframe faster)
DT <- data.table(x=rep(c("a","b","c"), each=100), y=rnorm(300))
setkey(DT, x)
DT['a'] #Since we set the key to be x, R knows that has to look for "a" within "x" to subset the data frame.

DT1 <- data.table(x=c('a', 'a', 'b', 'dt1'), y=1:4)
DT2 <- data.table(x=c('a', 'b', 'dt2'), z=5:7)
setkey(DT1, x); setkey(DT2, x)
merge(DT1, DT2) #using keys to merge two data tables. It only keept the cases where "x" coincide

##### Fast reading
big_df <- data.frame(x=rnorm(1E6), y=rnorm(1E6))
file <- tempfile()
write.table(big_df, file=file, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
system.time(fread(file))

system.time(read.table(file, header=TRUE, sep="\t"))


##MySQL
library(RMySQL)
ucscDb <- dbConnect(MySQL(),user="genome", host="genome-mysql.cse.ucsc.edu") #This opens a connection given by this handle (ucscDb)
result <- dbGetQuery(ucscDb,"show databases;"); dbDisconnect(ucscDb); #"show databses" is actually a MySQL command that we are sending to this database through the dbGetQuery function. "dbDisconnect()" disconnects you from the MySQL server (this is something important to do). It retunrs TRUE, meaning tha we did disconnect from the server.

result # List all the databases that are available wthin the MySQL server that is located at this particular host address ("genome-mysql.cse.ucsc.edu")

hg19 <- dbConnect(MySQL(),user="genome", db="hg19", host="genome-mysql.cse.ucsc.edu") #Adding the database we want to access (the server can have multiple databases)
allTables <- dbListTables(hg19) # List all the tables in the database "hg19" (which was specified in the code above)

length(allTables) #Entrega el número de tablas
allTables[1:5] #Entrega el nombre de las 5 primeras tablas

dbListFields(hg19,"affyU133Plus2") #List the fields for one specific table within "hg19". Fields corresponds to something like the column names in tha dataframe.

dbGetQuery(hg19, "select count(*) from affyU133Plus2") #provides the number of rows in that specific table. It works by passing, in quotes, a MySQL command. This command counts all the records in the table

affyData <- dbReadTable(hg19, "affyU133Plus2") #To extract one of the dables. So, you can extract the data in this databse, one table at the time.
head((affyData))

query <- dbSendQuery(hg19, "select * from affyU133Plus2 where misMatches between 1 and 3") #To select a subset of the data (useful because one single table can be too big to handle in R). Here, we use a MySQL command saying to select * (i.e. select all the observations in this table affy...), where the misMatches variable is between 1 and 3.
affyMis <- fetch(query); quantile(affyMis$misMatches) #fetches the results of query.
affyMisSmall <- fetch(query,n=10); dbClearResult(query); #fetches only the first 10 records. When you do that, you need to clear the query (with dbClearResult)
dim(affyMisSmall)

dbDisconnect(hg19)

##R HDF5 package
BiocManager::install("rhdf5")                     

library(rhdf5)
created = h5createFile("example.h5") #to create a hdf5 file
created

###############Create groups
created = h5createGroup("example.h5","foo")
created = h5createGroup("example.h5","baa")
created = h5createGroup("example.h5","foo/foobaa") #creates a subgroup of foo called foobaa
h5ls("example.h5") #"ls" stands for list

A = matrix(1:10,nr=5,nc=2)
h5write(A,"example.h5","foo/A") #write matrix A to group foo
B = array(seq(0.1,2.0,by=0.1),dim=c(5,2,2))
attr(B,"scale") <- "liter"
h5write(B, "example.h5","foo/foobaa/B")
h5ls("example.h5")

df = data.frame(1L:5L, seq(0,1,length.out=5), c("ab","cde","fghi","a","s"), stringsAsFactors=FALSE)
h5write(df,"example.h5","df") #write the dataset to the top level group.
h5ls("example.h5")

################Reading data

readA = h5read("example.h5","foo/A")
readB = h5read("example.h5","foo/foobaa/B")
readdf = h5read("example.h5","df")
readA

#################Writing and reading chunks
h5write(c(12,13,14),"example.h5","foo/A",index=list(1:3,1)) #writing these numbers to the first 3 rows in the first column of matrix A in group foo in file example.h5. We can pass the index argument to read specific portions of a dataset (using h5read instead of h5write).
h5read("example.h5","foo/A")


##Webscraping (three different approaches)
con = url("http://scholar.google.com/citations?user=HI-I6C0AAAAJ&hl=en") # Open a connection to a url (this is what the url function is for).
htmlCode = readLines(con) #to read out some of the data from that connection
close(con)
htmlCode

library(XML)
library(RCurl)
url <- "https://scholar.google.com/citations?user=HI-I6C0AAAAJ&amp;hl=e"
tabs <- getURL(url)
html <- htmlTreeParse(tabs,useInternalNodes=T) #internal nodes true, so that we can get the internal structure out
xpathSApply(html, "//title",xmlValue)
xpathSApply(html,"//td[@id='col-citedby']", xmlValue)


library(httr)
html2 = GET(url)
content2 = content(html2,as="text") #the content was extracted as text (just one large text string)
parsedHtml = htmlParse(content2, asText=TRUE)
xpathSApply(parsedHtml, "//title", xmlValue)


#################Accessing websites with passwords (httr package is useful for this)
pg1 = GET("http://httpbin.org/basic-auth/user/passwd", authenticate("user","passwd"))
pg1
names(pg1)

###################Handles (allows using authentication across multiple accesses to a website)
google = handle("http://google.com") #Here we created a handle through google, so if I use this handle I can save the cookies for the next access
pg1 = GET(handle=google,path="/")
pg1 = GET(handle=google,path="search")


## Accessing Twitter from R (with httr package)
myapp = oauth_app("twitter", key="YourConsumerKeyHere", secret="YourConsumerSecretHere") #it starts the authorization process for the application. "twitter" because here, it is the twitter app. Consumer key and secret obtained from the application website
sig = sign_oauth1.0(myapp, token="YourTokenHere", token_secret="yourTokenSecretHere") #signing in with the application. Token and Token Secret collected from the website. This (and the previous) are the credentials that will allow to access data that is privately held by Twitter, taht's only available to people with an application.
homeTL = GET("https://api.twitter.com/1.1/statuses/home_timeline.json", sig) #Similar to what we do for reading data off the internet. We use a very specific URL, though: the one corresponding to the twitter API, the version of the API the we are using (1.1). The other components correspond to the data we would like to get out (here, the statuses on my home timeline, and I'm going to get it out as json--the only data supported by Twitter). We pass it the authentication we used with oauth sign in. What we get out is the page that corresponds to the url, which is acutally a JSON data.

##########Converting the json object

json1 = content(homeTL)
json2 = jsonlite::fromJSON(toJSON(json1)) #Reformat the data as a Data Frame, so what we do is taking the JSON structure that we got out from the original command to content, and we convert that strucred R objects back into JSON. The result is a dataframe where each row corresponds to a tweet in the timline.
json2[1,1:4]

#Summirizing
if(!file.exists("./data")) {dir.create("./data")}
fileUrl <- "https://data.baltimorecity.gov/api/views/k5ry-ef3g/rows.csv?accessType=DOWNLOAD"
download.file(fileUrl, destfile="./data/restaurants.csv", method="curl")
restData <- read.csv("./data/restaurants.csv")

head(restData,n=3)
tail(restData,n=3) #last 3 rows
summary(restData) #Overall summary of the dataset, with some information per variable
str(restData) #N of varables and observations + information about each of the variable with examples
quantile(restData$councilDistrict, na.rm=TRUE)
quantile(restData$councilDistrict, probs=c(0.5,0.75,0.9)) #telling it to look at different probabilities (quantiles)

table(restData$zipCode,useNA="ifany") #table of a specific varaible. "useNA=ifany" adds a column if there are any NAs. By default "table" does not show the number of missing values.
table(restData$councilDistrict,restData$zipCode) #two-dimensional tables

sum(is.na(restData$councilDistrict)) #check for missing values. Sums the times that is.na is true
any(is.na(restData$councilDistrict)) #TRUE or FALSE, whether there are any NA

all(restData$zipCode > 0) #Check and see if every single value satisfy the condition (eg, here we would expect that all the zip codes are greater than 0)

colSums(is.na(restData)) #sums the NAs by column
all(colSums(is.na(restData))==0) #Checks whether all the sum of NAs across columns are equal to 0 (i.e. there are no missing values in the dataset)

table(restData$zipCode %in% c("21212","21213")) #to find all the zipcodes that are equal to 21212 or 21213. Is there any value of this variables that falls into this vector. Return the number of Trues and Falses
restData[restData$zipCode %in% c("21212","21213"),]#returns all the cases where the zip code falls into the vector

data(UCBAdmissions)
DF = as.data.frame(UCBAdmissions) #creates a dataframe from the database UCBAdmissions
summary(DF) # summary of the variables at datafram DF

xt <- xtabs(Freq ~ Gender + Admit,data=DF) #creates a crosstab displaying sums of Freq, broke down by Gender and Admit

warpbreaks$replicate <- rep(1:9, len = 54) #using warpbreaks dataset, adding another replicate variable
xt = xtabs(breaks~., data=warpbreaks) #creats a crosstab, displaying the sums of breaks, breaking that down by all the other variables in the dataset. To do that it creats mutiple two-dimensional tables
ftable(xt) #transforms xt into a one single flat table

fakeData = rnorm(1e5)
object.size(fakeData) #returns the size in memory of the object
print(object.size(fakeData), units="MB") #same, expressed in MB rather than bytes

#Creating variables (example using restData dataset downloaded in the previous set of commands)
if(!file.exists("./data")) {dir.create("./data")}
fileUrl <- "https://data.baltimorecity.gov/api/views/k5ry-ef3g/rows.csv?accessType=DOWNLOAD"
download.file(fileUrl, destfile="./data/restaurants.csv", method="curl")
restData <- read.csv("./data/restaurants.csv")

###### creating sequences (eg to index operations done on data)
s1 <- seq(1,10,by=2) #sequence from 1 to 10, distnace betweewn numbers is 2 (1 is the minimum value, 10 is the maximum value. Here there is no 10 because of the minimum and the "by 2" rule)
s2 <- seq(1,10, length=3) #same specifiying that 3 values must be created. It, by defalut, makes the distance between values equal
x <- c(1,2,8,25,100); seq(along = x) #generates a sequence along the elements of vector x. Useful for looping or subseting the elements of x.

####### subsetting variables
restData$nearMe = restData$neighborhood %in% c("Roland Park", "Homeland") #creates teh variable "nearMe" indicating whehter the neigbohood falls into the the vector
table(restData$nearMe)

####### creating binary variables
restData$zipWrong = ifelse(restData$zipCode < 0, TRUE, FALSE) #creates a variable "zipWrong" indicating whether zipCode is negative, assigninng value "TRUE" to zipCode values below 0 and value "FALSE" to zipCode values 0 or more
table(restData$zipWrong,restData$zipCode < 0) #table of zipWrong and the check of values of zipCode below 0

###### creating categorical variables
restData$zipGroups = cut(restData$zipCode, breaks=quantile(restData$zipCode)) #creates variable "zipGroups" which values are quantiles of zipCode. Alternatively, we can assign "breaks" to specific cut values.
table(restData$zipGroups)
table(restData$zipGroups,restData$zipCode)

library(Hmisc)
restData$zipGroups = cut2(restData$zipCode, g=4) #same, more easily done by asking to generate four groups along the quantiles (groups of roughly equal size). It creates them as a factor variable.
table(restData$zipGroups)

###### creating factor variables
restData$zcf <- factor(restData$zipCode) #just copying ZipCode and pasting it as a factor rather than numeric var
restData$zcf[1:10]
class(restData$zcf)

yesno <- sample(c("yes","no"), size=10, replace=TRUE) #creats a dummy vector which is not made of factors yet
yesnofac = factor(yesno, levels=c("yes","no")) #turns yesno into factor. By default it assign values to lables in alphabetical order (so, "no" would be 1 and "yes" would be 2). So we tell it to creat them in to the order especified in "levels" argument.
relevel(yesnofac,ref="yes") # makes the references class equal to "yes" (no sé para qué sirve esto)
as.numeric(yesnofac) #treats the factor as numeric.


###### mutate function
library(Hmisc); library(plyr)
restData2 = mutate(restData, zipGroups=cut2(zipCode,g=4)) #creates a new dataframe, adding the variable zipGroups
table(restData2$zipGroups)


# Reshaping data
library(reshape2)
head(mtcars)

######### Melting data frames
mtcars$carname <- rownames(mtcars) #generates a new variable "carnames" that is equal to the row names (we know that they are row names and not a variables because it does not have a variable lable)
carMelt <- melt(mtcars, id=c("carname","gear","cyl"), measure.vars=c("mpg","hp")) # we tell them which variables are id variables and which are measure variables. It then creats a bunch of id values, and then it melt all the other values (i.e., put all the measure variables into two: one identifying the variables and other its value, maintaining the values for all the other variables specified as id), so there is one row for every mpg and one row for every hp.
head(carMelt, n=3)
tail(carMelt, n=3)

######### Casting data frames (works with melted dataframes)
cylData <- dcast(carMelt, cyl ~ variable) # dcast recast the data set into a particular shape, into a particular data frame. Here, we want to see cyl broken down by the differen values of "variable" that where melted. It summarizes the dataset along the values of cyl (by default it returns the lenght, i.e., the number of measures of the values of "variable" for each value of "cyl").
cylData

cylData <- dcast(carMelt, cyl~ variable, mean) #same, but returns the mean of each value of "variable"
cylData

######### Averaging values
head(InsectSprays)
tapply(InsectSprays$count, InsectSprays$spray, sum) # Sums the "count" per each "spray" (tapply menas apply, along an index, a particular function. here, the sum of cout along spray)

spIns = split(InsectSprays$count, InsectSprays$spray) #create a list of  "count" per each "spray"
spIns
sprCount = lapply(spIns,sum) # applies the function "sum" to each of the elements in the list spIns
sprCount
unlist(sprCount) # transforms the list into a vector
sapply(spIns,sum) #a simpler way to do all the previous steps

library(plyr)
ddply(InsectSprays, .(spray), summarize, sum=sum(count)) #same. using plyr package. ".(spray)" are the variables we would like to summarize, the method of summarizing (sum) is by summing the variable "count"

########creating a new variable
spraySums <- ddply(InsectSprays, .(spray), summarize, sum=ave(count, FUN=sum)) # summarize the spray variable, hte method of summarizing is summing the counts. The ave function... Basically, it creates a new variables assigning to each value of "spray" the sum of "count" per the corresponding value of "spray" (e.g., for each spray=a, the new variable will be sum of count among spreay=a)
head(spraySums)



# dplyr
library(dplyr)

######## select
chicago <- readRDS("chicago.rds")
dim(chicago)
str(chicago)
names(chicago) #useful with dplyr because the we can access the columns using only the names

head(select(chicago, city:dptp)) #shows all the heads for the columns between city and dptp
head(select(chicago, -(city:dptp))) #the heads for all the columns, except the ones indicated in the range


###### filter (selecting rows based on values of the indicated variables)
chic.f <- filter(chicago, pm25tmean2 > 30)
head(chic.f)

chic.f <- filter(chicago, pm25tmean2 > 30 & tmpd > 80)
head(chic.f)


####### arrange (used to reorder the rows of a dataframe based on the values of a column)
chicago <- arrange(chicago, date)
head(chicago)
tail(chicago)

chicago <- arrange(chicago, desc(date)) #to arrange in descending order
head(chicago)
tail(chicago)


####### rename (to rename variables)
chicago <- rename(chicago, pm25 = pm25tmean2, dewpoint = dptp) # the new name goes to the left of the "=".
head(chicago)


######## mutate (to transform existing variables or create new ones)
chicago <- mutate(chicago, pm25detrend = pm25-mean(pm25, na.rm = TRUE)) #creates the variables that is equal to pm25 with the mean substracted off (i.e. centered)

head(select(chicago, pm25, pm25detrend))


####### group_by (split a dataframe behind the scene, according to certain categorical variables)
chicago <- mutate(chicago, tempcat = factor(1 * (tmpd > 80), label = c("cold", "hot")))
hotcold <- group_by(chicago, tempcat)
hotcold

summarize(hotcold, pm25 = mean(pm25, na.rm = TRUE), o3 = max(o3tmean2), no2 = median(no2tmean2)) #returns all theses summary statistics for both hot and cold categories. For each statistics, to the left of the "=" is the label that will be displayed in the output

chicago <- mutate(chicago, year =as.POSIXlt(date)$year + 1900) #generates a variable of the year (+ 1900, I assume, is because the reference point of posixlt is that year)
years <- group_by(chicago, year)
summarize(years, pm25 = mean(pm25, na.rm = TRUE), o3 = max(o3tmean2), no2 = median(no2tmean2))

###### to chain different operations together: %>% (pipeline operator)
chicago %>% mutate(month = as.POSIXlt(date)$mon + 1) %>% group_by(month) %>% summarize(pm25 = mean(pm25, na.rm = TRUE))


## Merging data
if(!file.exists("./dta")) {dir.create("./dta")}
fileUrl1 = "https://raw.githubusercontent.com/jtleek/dataanalysis/master/week2/007summarizingData/data/reviews.csv"
fileUrl2 = "https://raw.githubusercontent.com/jtleek/dataanalysis/master/week2/007summarizingData/data/solutions.csv"
download.file(fileUrl1, destfile="./data/reviews.csv", method="curl")
download.file(fileUrl2, destfile="./data/solutions.csv", method="curl")

reviews = read.csv("./data/reviews.csv"); solutions <- read.csv("./data/solutions.csv")
head(reviews,2) #variable "soution_id" correspond to the "id" variable in the "solutions" database
head(solutions,2)

names(reviews)
names(solutions)

mergedData = merge(reviews, solutions, by.x="solution_id", by.y="id", all=TRUE) #"reviews" is the "x" dataframe; "solutions" is the "y" dataframe
head(mergedData)

########## Merge using join in the plyr package

df1 = data.frame(id=sample(1:10), x=rnorm(10))
df2 = data.frame(id=sample(1:10), y=rnorm(10))
arrange(join(df1,df2),id) # easier, but can only do the merge on the bases of a common variable name between the two datasets. The "arrange" command to arrange them in increasing order by id

df3 = data.frame(id=sample(1:10), z=rnorm(10))
dfList =list(df1,df2,df3)
join_all(dfList) #creates a new dataset by joining all the datasets in the list


##Editing text data
if(!file.exists("./dta")) {dir.create("./dta")}
fileUrl <- "https://data.baltimorecity.gov/api/views/2weq-566u/rows.csv?accessType=DOWNLOAD"
download.file(fileUrl, destfile = "./data/cameras.csv", method = "curl") 
cameraData <- read.csv("./data/cameras.csv")
names(cameraData)

tolower(names(cameraData)) #Makes all letters lowercase (here, for the variable names). There is also a "toupper" command, in case we want to do the opposite.

############ Fixing Character vectors -strsplit() (good for automatically splitting variable names)
splitNames = strsplit(names(cameraData), "\\.") # splitting var names on periods (we have to use the "\\" because the period is a reserved character)
splitNames[[6]] # the name of this variable (Location.1) had two words separated with a dot, now it has two components and the dot has been removed. This has been sotored in the "splitNames" object

########### Quick aside - lists
mylist <- list(letters = c("A", "b", "c"), numbers =1:3, matrix(1:25, ncol = 5))
head(mylist)
mylist[1]
mylist$letters
mylist[[1]]

########## Fixing character vectors - sapply() (applias a function to each element in a vector or list)
splitNames [[6]] [1]
firstElement <- function(x) {x[1]}
sapply(splitNames, firstElement) #Returns the first element in each of the varnames (as specified in firstElement)


reviews = read.csv("./data/reviews.csv"); solutions <- read.csv("./data/solutions.csv")
names(reviews)
sub("_", "", names(reviews),) #substituting characters (here substituting "_" for "")

testName <- "this_is_a_test"
gsub("_", "", testName) # "gsub" to resplace multiple instances of a particular character (sub only removes the first instance in each name)

###### Finding values -grep(), grepl()
grep("Alameda", cameraData$intersection) # show all the cases where "Alameda" is found in the variable "intersection"
grep("Alameda", cameraData$intersection, value=TRUE) # returns the value of the variable intersection in the cases where "Alameda" appears instead of the number of the cases
grep("JeffStreet", cameraData$intersection, value=TRUE) #Looking for values that don't appear in the intersection variable (returns a 0)

table(grepl("Alameda", cameraData$intersection)) # Returns a table indicating the number of times where "Alameda" was present and absent in the variable "intersection"

cameraData2 <- cameraData[!grepl("Alameda",cameraData$intersection),] #subsets the database excluding all the cases where "Alameda" was present in the variable "intersection"

#########Other string functions
 library(stringr)
nchar("Jeffrey Leek") #counting characters
substr("Jeffrey Leek",1,7) #substract the first through the seventh letter (i.e. eliminates the rest)
paste("Jeffrey", "Leek") #concatenate strings (default separation is a space. Can be changed with the "sep" argument)
paste0("Jeffrey", "Leek") #paste strings without spaces in it
str_trim("Jeff     ") #eliminates any excess space at the beginning or the end of the string


# Working with dates
d1 = date() #returns current date and time. Stored as a character
d1
class(d1)

d2 = Sys.Date() #returns current date. Stored as date.
d2
class(d2)

format(d2,"%a %b %d") #%d = day as number (0-31); %a = abbreviated weekday;  %A = unabbreviated weekday; %m = month(00-12); %b = abbreviated month; %B = unabbreviated month; %y = 2 digit year; %Y = four digit year

########## Creating dates
x = c("1jan1960", "2jan1960", "31mar1960", "30jul1960"); z = as.Date(x, "%d%b%Y")
z
z[1] - z[2] #returns the day difference between the two dates
as.numeric(z[1]-z[2]) #same, returning the day difference in plain numbers

weekdays(d2)
months(d2)
julian(d2) #converts the date to Julian (n of days since the origin--in this case the origin is 1970-01-01)

######## Lubridate package (converts a number into a date, regardless of the format of the number)
library(lubridate); ymd("20140108")
mdy("08/04/2013")
dmy("03-04-2013")

ymd_hms("2011-08-03 10:15:03")
ymd_hms("2011-08-03 10:15:03", tz="Pacific/Auckland")

?Sys.timezone #help for the function Sys.timezone

x = dmy(c("1jan2013", "2jan2013", "31mar2013", "30jul2013"))
wday(x[1])
wday(x[1], label=TRUE)
