## app.R ##
library(shinydashboard)
library('ggplot2')


data <- read.csv(url("https://covid19.who.int/WHO-COVID-19-global-data.csv"))
data$Date <- as.Date(strptime(data$ï..Date_reported, "%Y-%m-%d"))
cont <- unique(data$Country)


header <- dashboardHeader(title = "Basic dashboard")

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
          fluidRow(
            box(
              title = "Histogram",
              solidHeader = TRUE,
              collapsible = TRUE,
              plotOutput("plot1", height = 250)
            ),
            
            box(
              title = "Controls",
              selectInput(
                "cont",
                "Select Index",
                cont,
                selected = 'Turkey',
                multiple = FALSE
              )
            )
          )),
  
  # Second tab content
  tabItem(tabName = "widgets",
          h2("Widgets tab content"))
))





ui <- dashboardPage(header,
                    sidebar,
                    body)



server <- function(input, output) {
  
  
  
  output$plot1 <- renderPlot({
    temp <- data[which(data$Country == input$cont),]
    ggplot(temp) +
      geom_line(
        aes(Date, as.numeric(New_cases)),
        group = 1,
        color = "red",
        size = 1
      )
  })
}

shinyApp(ui, server)