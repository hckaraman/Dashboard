library(shiny)
library(leaflet)
library(leaflet.extras)
require("rgdal")
library(sf)
library(rgeos)
library(raster)
library(whitebox)

latitude <-
  c(35.94077, 35.83770, 35.84545, 35.81584, 35.79387, 36.05600)
longitude <-
  c(-78.58010,
    -78.78084,
    -78.72444,
    -78.62568,
    -78.64262,
    -78.67600)
radius <- c(15, 12, 12, 12, 12, 15)
ids <- c("a", "b", "c", "d", "e", "f")

shinyApp(
  ui = fluidPage(titlePanel("Column Plot"),
                 tabsetPanel(
                   tabPanel(
                     "Upload File",
                     titlePanel("Uploading Files"),
                     sidebarLayout(
                       sidebarPanel(
                         fileInput('myFile', 'Choose JSON File',
                                   accept = c(".json")),
                         
                         # added interface for uploading data from
                         # http://shiny.rstudio.com/gallery/file-upload.html
                         tags$br(),
                         checkboxInput('header', 'Header', TRUE),
                         radioButtons('sep', 'Separator',
                                      c(
                                        Comma = ',',
                                        Semicolon = ';',
                                        Tab = '\t'
                                      ),
                                      ','),
                         radioButtons(
                           'quote',
                           'Quote',
                           c(
                             None = '',
                             'Double Quote' = '"',
                             'Single Quote' = "'"
                           ),
                           '"'
                         )
                         
                       ),
                       mainPanel(leafletOutput(
                         "mymap", width = "100%", height = "800px"
                       ))
                     )
                   )
                 )),
  
  
  # fluidRow(
  #   leafletOutput("mymap",height=400)),
  #   # leafletMap(
  #   #   "map", "100%", 400,s
  #   #   initialTileLayer = "//{s}.tiles.mapbox.com/v3/jcheng.map-5ebohr46/{z}/{x}/{y}.png",
  #   #   initialTileLayerAttribution = HTML('Maps by <a href="http://www.mapbox.com/">Mapbox</a>'),
  #   #   options=list(
  #   #     center = c(37.45, -93.85),
  #   #     zoom = 4,
  #   #     maxBounds = list(list(17, -180), list(59, 180))))),
  # fluidRow(verbatimTextOutput("text"))),
  
  
  server = function(input, output, session) {
    value <- reactiveValues()
    SpatialPolygonsDataFrame(SpatialPolygons(list()),
                             data = data.frame (notes = character(0), stringsAsFactors = F)) -> value$drawnPoly
    
    
    data <- reactive({
      req(input$file1)
      
      inFile <- input$file1
      if (is.null(inFile))
        return(NULL)
      return(inFile)
    })
    
    output$mymap <- renderLeaflet(
      leaflet() %>%
        addProviderTiles("Esri.OceanBasemap", group = "Ocean Basemap") %>%
        setView(lng = 37, lat = 38.0, zoom = 2) %>%
        addDrawToolbar(
          targetGroup = 'draw',
          editOptions = editToolbarOptions(selectedPathOptions = selectedPathOptions())
        )
      
    )
    
    observeEvent(input$mymap_draw_new_feature, {
      feature <- input$mymap_draw_new_feature
      coor <-
        unlist(input$mymap_draw_new_feature$geometry$coordinates)
      Longitude <- coor[seq(1, length(coor), 2)]
      Latitude <- coor[seq(2, length(coor), 2)]
      
      
      if (feature$properties$feature_type == "marker")
      {
        pnt_attrib = data.frame(                           # data.frame object
          name = "point"
        )
        lnd_sf <- SpatialPointsDataFrame(cbind(Longitude, Latitude), pnt_attrib)    # sf object
        proj4string(lnd_sf) <-
          "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
        writeOGR(
          lnd_sf,
          dsn = "./Data/point.shp",
          layer = "shpExport",
          driver = "ESRI Shapefile",
          overwrite_layer = TRUE
        )
      } 
      
      poly <- Polygon(cbind(Longitude, Latitude))
      polys <-
        Polygons(list(poly),
                 ID = input$mymap_draw_new_feature$properties$`_leaflet_id`)
      spPolys <- SpatialPolygons(list(polys))
      
      df <-
        SpatialPolygonsDataFrame(spPolys, data = data.frame(notes = NA, row.names =
                                                              row.names(spPolys)))
      
      
      proj4string(df) <-
        "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
      writeOGR(
        df,
        dsn = "shpExport.shp",
        layer = "shpExport",
        driver = "ESRI Shapefile",
        overwrite_layer = TRUE
      )
      
      v()
      # print(feature)
      
    })
    
    observeEvent(input$mymap_click, {
      feature <- input$mymap_click
      
      print(feature)
      print(getwd())
    })
    
    
    observeEvent(input$myFile, {
      inFile <- input$myFile
      if (is.null(inFile))
        return()
      
      nycounties <- rgdal::readOGR(inFile$datapath)
      # file.copy(inFile$datapath, file.path("c:/temp", inFile$name) )
      
      sf_cent <- gCentroid(nycounties)
      x <- data.frame(sf_cent)
      lon <- as.numeric(x[1])
      lat <- as.numeric(x[2])
      
      leafletProxy("mymap", data = nycounties) %>%
        clearShapes() %>%
        addPolygons(stroke = FALSE,
                    smoothFactor = 0.2,
                    fillOpacity = 0.6) %>%
        setView(lng = lon,
                lat = lat,
                zoom = 7)
    })
    
    v<- function() {
      
      wbt_snap_pour_points(
          pour_pts = "./Data/point.shp",
          flow_accum = "./Data/out_acc.tif",
          output = "./Data/snapped.shp",
          snap_dist = 0.02
        )
      
      wbt_watershed(d8_pntr = "./Data/flow_dir.tif",
                                  pour_pts = "./Data/snapped.shp",
                                  output = "./Data/Watershed.tif")
      
      r1 <- raster("./Data/Watershed.tif")
      r.to.poly<-rasterToPolygons(r1, dissolve = T)
      writeOGR(
          r.to.poly,
          dsn = "./Data/Polygon.json",
          layer = "shpExport",
          # driver = "ESRI Shapefile",
          driver = "GeoJSON",
          overwrite_layer = TRUE
        )
      
      nycounties <- rgdal::readOGR("./Data/Polygon.json")
      sf_cent <- gCentroid(nycounties)
      x <- data.frame(sf_cent)
      lon <- as.numeric(x[1])
      lat <- as.numeric(x[2])
      
      leafletProxy("mymap", data = nycounties) %>%
        clearShapes() %>%
        addPolygons(stroke = FALSE,
                    smoothFactor = 0.2,
                    fillOpacity = 0.6) %>%
        setView(lng = lon,
                lat = lat,
                zoom = 7)
    }
    
    
    
    
    
    # })
    
    
    
    
  }
)






# wbt_flow_accumulation_full_workflow(dem = dem,out_dem = "./Test/Data/out_dem.tif",out_pntr = "./Test/Data/flow_dir.tif",out_accum =  "./Test/Data/out_acc.tif" )
# wbt_extract_streams(flow_accum = "./Test/Data/out_acc.tif",
#                     output = "Test/Data/streams.tif",
#                     threshold = 3)
# wbt_snap_pour_points(
#   pour_pts = "./Test/Data/point.shp",
#   flow_accum = "./Test/Data/out_acc.tif",
#   output = "./Test/Data/snapped.shp",
#   snap_dist = 0.02
# )
# wbt_watershed(d8_pntr = "./Test/Data/flow_dir.tif",
#               pour_pts = "./Test/Data/snapped.shp",
#               output = "./Test/Data/Watershed.tif")
# 
# r1 <- raster("./Test/Data/Watershed.tif")
# r.to.poly<-rasterToPolygons(r1, dissolve = T)
# writeOGR(
#   r.to.poly,
#   dsn = "./Test/Data/Polygon.json",
#   layer = "shpExport",
#   # driver = "ESRI Shapefile",
#   driver = "GeoJSON",
#   overwrite_layer = TRUE
# )

