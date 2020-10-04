require("shiny")
require("shinydashboard")
require("plotly")
require("leaflet")





header <- dashboardHeader(title = "Covid-19 Dashboard")

sidebar <- dashboardSidebar(sidebarMenu(
  menuItem(
    "Dashboard",
    tabName = "dashboard",
    icon = icon("dashboard")
  )
  #menuItem("Widgets", tabName = "widgets", icon = icon("th"))
)
)


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
              title = "Daily Cases",
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