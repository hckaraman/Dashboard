library(shiny)
library(RSQLite)
library(leaflet)
library(plotly)
library(shinyWidgets)


files <- Sys.glob(file.path('./Data', "*.csv"))
print(files)
shinyUI(pageWithSidebar(
  headerPanel("Drought Analysis"),
  sidebarPanel(h3("Upload Data"),fluidRow( column(6,numericInput("area1", "Area:", 100, min = 10, max = 1000))),
               h3("NAM Parameters"),fluidRow(
    column(6,numericInput("area", "Area:", 100, min = 10, max = 1000)),
    column(6,selectInput("basin", "Select Basin", c("Alihoca","Cakit","Darbogaz"), selected = "Alihoca", multiple = FALSE)),
    column(6,materialSwitch(inputId = "cal", label = "Calibration", inline = T,right = FALSE,status = "danger")),
    column(6,selectInput("objective", "Objective Function:",
                         c("NSE","KGE","RMSE","MAE","Volume Error"), selected = "NSE", multiple = FALSE)),
    column(6,sliderInput("umax", "Umax:",
                         min = 0, max = 50, value = 25)),
    column(6,sliderInput("lmax", "Lmax:",
                         min = 0, max = 1000, value = 500)),
    column(6,sliderInput("cqof", "Cqof:",
                           min = 0, max = 1, value = 0.5)),
    column(6,sliderInput("ckif", "Ckif:",
                           min = 200, max = 1000, value = 600)),
    column(6,sliderInput("ck12", "Ck12:",
                           min = 10, max = 50, value = 30)),
    column(6,sliderInput("tof", "Tof:",
                           min = 0, max = 1, value = 0.5
               )),
    column(6,sliderInput("tif", "Tif:",
                           min = 0, max = 1, value = 0.5
               )),
    column(6,sliderInput("tg", "Tg:",
                           min = 0, max = 1, value = 0.5
               )),
               
    column(6,sliderInput("ckbf", "Ckbf:",
                           min = 500, max =5000, value = 2500
               )),
    column(6,sliderInput("csnow", "Csnow:",
                           min = 0, max = 4, value = 2,step = 0.25
               )),
               downloadButton("downloadData", "Download Results"))),
  mainPanel(
    
    tabsetPanel(type = "tabs",
                tabPanel("Plot", plotlyOutput(outputId ="data_plot", height = "800px")),
                tabPanel("Data", DT::dataTableOutput("ysummary")),
                tabPanel("Model Performance", DT::dataTableOutput("stats",width = "75%")),
                tabPanel("Model Parameters", DT::dataTableOutput("parameters",width = "75%")),
                tabPanel("Flow Duration", plotlyOutput(outputId ="flowdur", height = "500px")),
                tabPanel("Map", leafletOutput("mymap",height = "800px"))
                
                
                #plotOutput("main_plot", height = "800px")
                
    )
  )))