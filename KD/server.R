library(shiny)
library(RSQLite)
library(stringr)
library(ggplot2)
library('zoo')
# library(DT)
library(sf)



db <- '/home/cak/Desktop/Dashboard/Data/data.db'



shinyServer(function(input, output) {
  
  output$plot <- renderPlotly({
    query = str_interp("SELECT * from data d join hes h on d.Drenaj_No = h.BID where Model = '${input$model}' and Senaryo = '${input$senaryo}' and h.PROJECT_NA = '${input$hpp}'")
    print(query)
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
    
    
    
    
  })
  
  
  output$mymap <- renderLeaflet({
    # reactive expression code required here to connect with ui selection?
    
    query <- "SELECT * from hes h "
    data <- dbGetQuery(conn, query)
    
    m <- leaflet() %>%
      setView(37, 37.8, 2) %>%
      addProviderTiles("Esri.OceanBasemap")
    
   
    m
    
    
    
    
    
  })

  
})

# coordinates(data) <- cbind(data$Longtitude , data$Latitude)
# 
# spdf <- SpatialPointsDataFrame(coords = xy, data = data,
#                                proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
