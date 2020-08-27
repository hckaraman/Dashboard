## app.R ##
library(plotly)
library(shinydashboard)
library('ggplot2')


data <-
  read.csv(url("https://covid19.who.int/WHO-COVID-19-global-data.csv"))
data$Date <- as.Date(strptime(data$ï..Date_reported, "%Y-%m-%d"))
cont <- unique(data$Country)


header <- dashboardHeader(title = "Covid-19 Dashboard")

sidebar <- dashboardSidebar(sidebarMenu(
  menuItem(
    "Dashboard",
    tabName = "dashboard",
    icon = icon("dashboard")
  ),
  menuItem("Widgets", tabName = "widgets", icon = icon("th"))
))


body <-   dashboardBody(tabItems(
  # First tab content
  tabItem(tabName = "dashboard",
          
          fluidRow(box(
            width = 2,
            selectInput(
              "cont",
              "Select country",
              cont,
              selected = 'Turkey',
              multiple = FALSE
            )
          ),
          valueBoxOutput("totalCase", width = 2),
          valueBoxOutput("totalDeath", width = 2),
          valueBoxOutput("newCase", width = 2),
          valueBoxOutput("newDeath", width = 2),),
          
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
                title = "New Cases",
                width = 10,
                solidHeader = TRUE,
                collapsible = TRUE,
                leafletOutput("map")  
              ),
            
            
          )),
  
  # Second tab content
  tabItem(tabName = "widgets",
          h2("Widgets tab content"))
))





ui <- dashboardPage(header,
                    sidebar,
                    body)



server <- function(input, output) {
  output$plot1 <- renderPlotly({
    temp <- data[which(data$Country == input$cont), ]
    p <- plot_ly(
      temp,
      x = ~ Date,
      y = ~ New_cases,
      type = 'scatter',
      mode = 'markers+lines',
      size = ~ New_deaths,
      colors = 'silver',
      hoverinfo = 'x+text',
      text = ~ paste(
        "New Cases: ",
        New_cases,
        "<br>",
        "New Deaths: ",
        New_deaths,
        "<br>"
      ),
      
      #Choosing the range of the bubbles' sizes:
      sizes = c(5, 30),
      marker = list(
        opacity = 0.5,
        sizemode = 'diameter',
        color = ~ New_deaths,
        line = list(color = ~ New_deaths,
                    width = 1),
        colorbar = list(title = 'New Deaths')
      )
    ) %>%
      layout(
        title = paste('New Cases in ' ,input$cont),
        xaxis = list(title = "Date", showgrid = TRUE),
        yaxis = list(title = "New Cases", showgrid = TRUE),
        showlegend = FALSE
      )
    
    p
    
  })
  
  output$totalCase <- renderValueBox({
    temp <- data[which(data$Country == input$cont), ]
    valueBox(
      paste0(tail(temp$Cumulative_cases, n=1)), "Total Cases", icon = icon("list"),
      color = "purple"
    )
  })
  
  output$totalCase <- renderValueBox({
    temp <- data[which(data$Country == input$cont), ]
    valueBox(
      paste0(tail(temp$Cumulative_cases, n=1)), "Total Cases", icon = icon("list"),
      color = "aqua"
    )
  })
  
  output$totalDeath <- renderValueBox({
    temp <- data[which(data$Country == input$cont), ]
    valueBox(
      paste0(tail(temp$Cumulative_deaths, n=1)), "Total Deaths", icon = icon("list"),
      color = "red"
    )
  })
  
  output$newCase <- renderValueBox({
    temp <- data[which(data$Country == input$cont), ]
    valueBox(
      paste0(tail(temp$New_cases, n=1)), "New Case", icon = icon("list"),
      color = "purple"
    )
  })
  
  output$newDeath <- renderValueBox({
    temp <- data[which(data$Country == input$cont), ]
    valueBox(
      paste0(tail(temp$New_deaths, n=1)), "New Deaths", icon = icon("list"),
      color = "yellow"
    )
  })
  
  output$map <- renderLeaflet({
    # reactive expression code required here to connect with ui selection?
    leaflet() %>% addProviderTiles("Esri.OceanBasemap") %>% 
      fitBounds(160, -30, 185, -50)
  })

  
  
}

shinyApp(ui, server)