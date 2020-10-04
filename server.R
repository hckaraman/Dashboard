## app.R
require('ggplot2')
require("sf")
require("plotly")
require("rgdal")
require("leaflet")
require("dplyr")
library(leaflet.extras)
library(vroom)
library(raster)


server <- function(input, output) {
  output$plot1 <- renderPlotly({
    temp <- data[which(data$location == input$cont), ]
    p <- plot_ly(
      temp,
      x = ~ date,
      y = ~ new_cases,
      type = 'scatter',
      mode = 'markers+lines',
      size = ~ new_deaths,
      colors = 'silver',
      hoverinfo = 'x+text',
      text = ~ paste(
        "New Cases: ",
        new_cases,
        "<br>",
        "New Deaths: ",
        new_deaths,
        "<br>"
      ),
      
      #Choosing the range of the bubbles' sizes:
      sizes = c(5, 30),
      marker = list(
        opacity = 0.5,
        sizemode = 'diameter',
        color = ~ new_deaths,
        line = list(color = ~ new_deaths,
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
    temp <- data[which(data$location == input$cont), ]
    valueBox(
      paste0(tail(temp$total_cases, n=1)), "Total Cases", icon = icon("list"),
      color = "purple"
    )
  })
  # 
  # output$totalCase <- renderValueBox({
  #   temp <- data[which(data$Country == input$cont), ]
  #   valueBox(
  #     paste0(tail(temp$Cumulative_cases, n=1)), "Total Cases", icon = icon("list"),
  #     color = "aqua"
  #   )
  # })
  
  output$totalDeath <- renderValueBox({
    temp <- data[which(data$location == input$cont), ]
    valueBox(
      paste0(tail(temp$total_deaths, n=1)), "Total Deaths", icon = icon("list"),
      color = "red"
    )
  })
  
  output$newCase <- renderValueBox({
    temp <- data[which(data$location == input$cont), ]
    valueBox(
      paste0(tail(temp$new_cases, n=1)), "New Case", icon = icon("list"),
      color = "purple"
    )
  })
  
  output$newDeath <- renderValueBox({
    temp <- data[which(data$location == input$cont), ]
    valueBox(
      paste0(tail(temp$new_deaths, n=1)), "New Deaths", icon = icon("list"),
      color = "yellow"
    )
  })
  
  output$map <- renderLeaflet({
    # reactive expression code required here to connect with ui selection?
    
   
    temp <- data %>%
      group_by(location) %>%
      slice(c(n())) %>%
      ungroup()
    
    df <- merge(states, temp, by.x="location", by.y="location")
    
    
    bins <- c(0, 50, 100, 500, 1000, 5000, 10000, 100000, Inf)
    pal <- colorBin("YlOrRd", domain = df$new_cases, bins = bins)
    
    labels <- sprintf(
      "<strong>%s</strong><br/>%g daily cases",
      df$location, df$new_cases
    ) %>% lapply(htmltools::HTML)
    
    m <- leaflet(df) %>%
      setView(37, 37.8, 2) %>%
      addProviderTiles("Esri.OceanBasemap")
    
    m <- m %>% addPolygons(
      fillColor = ~pal(df$new_cases),
      weight = 2,
      opacity = 1,
      color = "white",
      dashArray = "3",
      fillOpacity = 0.7,
      highlight = highlightOptions(
        weight = 5,
        color = "#666",
        dashArray = "",
        fillOpacity = 0.7,
        bringToFront = TRUE),
      label = labels,
      labelOptions = labelOptions(
        style = list("font-weight" = "normal", padding = "3px 8px"),
        textsize = "15px",
        direction = "auto"))
    m
    
    
    
    
    
  })
  
  observeEvent(input$cont, {
    temp <- states[which(states$location == input$cont), ]
    
    lat <- (ymax(temp)+ymin(temp))/2
    lon <- (xmax(temp)+xmin(temp))/2
    leafletProxy("map") %>%
      setView(lng = lon,
              lat = lat,
              zoom = 5)
  })
  
  
}

# shinyApp(ui, server)