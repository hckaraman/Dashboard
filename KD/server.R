library(shiny)
library(RSQLite)
library(stringr)
library(ggplot2)
library('zoo')
# library(DT)
library(sf)
require("rgdal")
library(tidyverse)





db <- '/home/cak/Desktop/Dashboard/KD/Data/data.db'
conn <- dbConnect(RSQLite::SQLite(),db)
basins <- readOGR("/home/cak/Desktop/Dashboard/KD/Data/HES_b.geojson")
hes <- readOGR("/home/cak/Desktop/Dashboard/KD/Data/HES_p.geojson")


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
      setView(37, 37.8, 4) %>%
      addProviderTiles("Esri.OceanBasemap")
    
    # labels <- sprintf(
    #   "<strong>%s</strong><br/>%g daily cases",
    #   df$BID
    # ) %>% lapply(htmltools::HTML)
    
    # m <- leaflet(df) %>%
    #   setView(37, 37.8, 2) %>%
    #   addProviderTiles("Esri.OceanBasemap")
    
    
   
    m
    
    
    
    
    
  })
  
  observeEvent(input$hpp, {
    
    
    query = str_interp("SELECT DISTINCT h.BID from data d join hes h on d.Drenaj_No = h.BID where h.PROJECT_NA  = '${input$hpp}'")
    bid <- dbGetQuery(conn, query)
    temp <- basins[which(basins$BID == as.character(bid)), ]
    temp_hes <- hes[which(hes$BID == as.character(bid)), ]
    temp_hes$popup <- str_c(temp_hes$COMPANY_NA,
                            temp_hes$PROJECT_NA,
                            temp_hes$Type,
                            temp_hes$City,
                            temp_hes$Capacity_U,
                            temp_hes$Service_Ca,
                           sep = "<br/>")
    
    lat <- (ymax(temp)+ymin(temp))/2
    lon <- (xmax(temp)+xmin(temp))/2
    leafletProxy("mymap") %>%
      setView(lng = lon,
              lat = lat,
              zoom = 8)
    
    leafletProxy("mymap",data = temp) %>% addPolygons(
      color = "#c9554d",
      weight = 2,
      opacity = 1,
      dashArray = "3",
      fillOpacity = 0.7,
      layerId = "foo",
      highlight = highlightOptions(
        weight = 5,
        color = "#666",
        dashArray = "",
        fillOpacity = 0.7,
        bringToFront = TRUE),
      # label = labels,
      # labelOptions = labelOptions(
      # style = list("font-weight" = "normal", padding = "3px 8px"),
      # textsize = "15px",
      # direction = "auto"))
    ) 
    
    leafletProxy("mymap",data = temp_hes) %>% addMarkers(layerId = "foop") %>% addMarkers(label = ~popup)
    
  })
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$hpp, ".csv", sep = "")
    },
    content = function(file) {
      
      query = str_interp("SELECT * from data d join hes h on d.Drenaj_No = h.BID where Model = '${input$model}' and Senaryo = '${input$senaryo}' and h.PROJECT_NA = '${input$hpp}'")
      data <- dbGetQuery(conn, query)
      
      # spei1 <- v(input$Freq)$spei1
      # res <- data.frame(as.matrix(spei1$fitted), date=time(spei1$fitted))
      # res$year <- trunc(res$date)
      # res$month <- (res$date - res$year) * 12 + 1
      
      write.csv(data, file, row.names = FALSE)
    }
  )

  
})

# coordinates(data) <- cbind(data$Longtitude , data$Latitude)
# 
# spdf <- SpatialPointsDataFrame(coords = xy, data = data,
#                                proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# 
# homicides %>%
#   mutate(popup = str_c(Date,
#                        Block,
#                        str_c("Location type:", `Location Description`,
#                              sep = " "),
#                        sep = "<br/>")) %>%
#   leaflet() %>%
#   addTiles() %>%
#   addMarkers(popup = ~popup) %>%
#   frameWidget()
