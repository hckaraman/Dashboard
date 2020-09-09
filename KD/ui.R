library(shiny)
library(RSQLite)
library(leaflet)
library(plotly)

db <- '/home/cak/Desktop/Dashboard/KD/Data/data.db'
conn <- dbConnect(RSQLite::SQLite(),db)
query <- "SELECT * from hes h "
dr <- dbGetQuery(conn, query)

shinyUI(pageWithSidebar(
  headerPanel("Climate Change Analysis"),
  sidebarPanel(selectInput("hpp", "Select HPP", dr$PROJECT_NA, selected = "Kalecik HPP", multiple = FALSE),
               selectInput("model", "Select Model", c("MPI-ESM-MR","HadGEM2-ES","CNRM-CM5"), selected = "MPI-ESM-MR", multiple = FALSE),
               selectInput("senaryo", "Select Senario", c('RCP4.5','RCP8.5'), selected = 'RCP4.5', multiple = FALSE),
               downloadButton("downloadData", "Download Results")),
  mainPanel(
    
    tabsetPanel(type = "tabs",
                tabPanel("Data Summary",
                         fluidRow(
                           plotlyOutput("plot",width = "100%",height = "600px"),
                           leafletOutput("mymap",width = "100%",height = "600px")))
                #plotOutput("main_plot", height = "800px")
                
    )
  )))