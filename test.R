library('ggplot2')

data <- read.csv(url("https://covid19.who.int/WHO-COVID-19-global-data.csv"))
data$Date <- as.Date(strptime(data$ï..Date_reported,"%Y-%m-%d"))
temp <- data[ which(data$Country_code=='TR'), ]

ggplot(temp) + 
  geom_line(aes(Date, as.numeric(New_cases)), group = 1,color="red", size=1) 



