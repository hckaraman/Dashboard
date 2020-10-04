library(vroom)
library(rgdal)

data <-
  vroom(url("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv"))

# data <-
#   vroom('./Data/test.csv')
data$date <- as.Date(strptime(data$date, "%Y-%m-%d"))
cont <- unique(data$location)
states <- readOGR("./Data/countries_simplifed.geojson")

date <- Sys.Date() -1
df <- data[which(data$date == date),] %>% arrange(desc(total_cases))
