#install.packages('jsonlite')
library('jsonlite')
data <- fromJSON("https://api.apify.com/v2/datasets/LYeOfHQwsv7FsfdGV/items?format=json&clean=1")
data <- as.data.frame(data)
data