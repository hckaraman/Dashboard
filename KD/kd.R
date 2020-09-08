library("readxl")
library(RSQLite)
library(stringr)
library(tidyverse)
library(tidyr)
require(zoo)
library(plotly)
library(ggplot2)

f1 <- '/home/cak/Desktop/Dashboard/Data/Results_1.xlsx'
f2 <- '/home/cak/Desktop/Dashboard/Data/Results_2.xlsx'
f_n <- '/home/cak/Desktop/Dashboard/Data/HES_Noktalari_06092020_mkd_head.xlsx'
db <- '/home/cak/Desktop/Dashboard/Data/data.db'

# df1 <- read_excel(f1)
# df2 <- read_excel(f2)
# df_n <- read_excel(f_n)
# 
# df <- rbind(df1,df2)
conn <- dbConnect(RSQLite::SQLite(),db)

# dbWriteTable(conn, "data", df)
# dbWriteTable(conn, "hes", df_n)
##########################

# hes <- df_n$PROJECT_NA


drno <- 'B_180'
model <- 'MPI-ESM-MR'
senaryo <-'RCP8.5'

query <- "SELECT DISTINCT Drenaj_No from data ;"
dr <- dbGetQuery(conn, query)


query = str_interp("SELECT * from data d join hes h on d.Drenaj_No = h.BID where d.Drenaj_No = '${drno}' and Model = '${model}' and Senaryo = '${senaryo}'")
data <- dbGetQuery(conn, query)
data <- tidyr::as_tibble(data)
data$Date <- as.yearmon(paste(data$Yil, data$Ay), "%Y %m")


fig <- plot_ly(data)
fig <- fig %>% add_trace(x = ~Date, y = ~Discharge, type = 'scatter',mode='lines', name = 'Discharge',
                         line = list(color = "#26828e"),
                         # line = list(color = '#45171D'),
                         
                         hoverinfo = "text",
                         text = ~paste(Discharge, ' m3'))

fig <- fig %>% layout(title = 'Total Discharge ',
                      xaxis = list(title = "Date"),
                      yaxis = list(side = 'left', title = 'Total Discharge , m3/s', showgrid = FALSE, zeroline = FALSE))

fig

