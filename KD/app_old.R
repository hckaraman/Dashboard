library("readxl")
library(RSQLite)
library(stringr)
library(tidyverse)
library(tidyr)
require(zoo)
library(plotly)
library(ggplot2)
library(shiny)
library(leaflet.extras)
require("shinydashboard")




db <- '/home/cak/Desktop/Dashboard/Data/data.db'
conn <- dbConnect(RSQLite::SQLite(),db)
query <- "SELECT * from hes h "
dr <- dbGetQuery(conn, query)


header <- dashboardHeader(title = "Covid-19 Dashboard")

sidebar <- dashboardSidebar(sidebarMenu(
  menuItem(
    "Dashboard",
    tabName = "dashboard",
    icon = icon("dashboard")
  )
  #menuItem("Widgets", tabName = "widgets", icon = icon("th"))
))


body <-   dashboardBody(tabItems(
  # First tab content
  tabItem(tabName = "dashboard",
          
          fluidRow(box(
            width = 2,
            selectInput(
              "cont",
              "Select country",
              dr$PROJECT_NA,
              selected = 'Turkey',
              multiple = FALSE
            )
          ),
          
          fluidRow(
            box(
              title = "New Cases",
              width = 10,
              solidHeader = TRUE,
              collapsible = TRUE,
              plotlyOutput("plot1", height = 500)
            )),
          
          fluidRow(
            box(
              title = "Daily Cases",
              width = 10,
              solidHeader = TRUE,
              collapsible = TRUE,
              leafletOutput("map")  
            ),
            
            
          )),
  
  # # Second tab content
  # tabItem(tabName = "widgets",
  #         h2("Widgets tab content"))
)))





ui <- dashboardPage(header,
                    sidebar,
                    body)


server <- function(input, output) {
  
  
  query = str_interp("SELECT * from data d join hes h on d.Drenaj_No = h.BID where d.Drenaj_No = '${drno}' and Model = '${model}' and Senaryo = '${senaryo}'")
  
  
  output$plot1 <- renderPlotly({
    
    
  
    
      
  })
  

}

shinyApp(ui, server)

