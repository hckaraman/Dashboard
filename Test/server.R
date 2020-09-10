library(shiny)
library(RSQLite)
library(stringr)
library(ggplot2)
library('zoo')
library(DT)
library(sf)
require("rgdal")
library(tidyverse)
library(raster)
library(reticulate)
library(plotly)
use_python(python = "C:\\Users\\cagri\\AppData\\Local\\Programs\\Python\\Python38\\python.exe", required = TRUE)
# py_config()
source_python("./nam.py")
# print(getwd())

shinyServer(function(input, output, session) {
  
  newData <- reactive({
    params <- c(input$umax, input$lmax, input$cqof, input$ckif
                ,input$ck12, input$tof, input$tif, input$tg
                , input$ckbf,input$csnow, 2)
    print(input$objective)
    temp <- run(input$area,params,input$cal,input$basin,input$objective) 
    temp[[1]]$date <- as.Date(row.names(temp[[1]]), "%Y-%m-%d")
    temp
  })
  
  output$data_plot <- renderPlotly({
    
    df <- newData()[[1]]
    
    pt <- ggplot(df, aes(x=date)) + 
      geom_line(aes(y = Q), color = "darkred") + 
      geom_line(aes(y = Qsim), color="steelblue", linetype="twodash") +
      xlab("Date") + ylab("Discharge , m3/s") 
    
    fig <- ggplotly(pt)
    fig
    
  })
  
  output$ysummary = DT::renderDataTable({
    df <- newData()[[1]]
    rownames(df) <- df$date
    df
  })
  
  output$parameters = DT::renderDataTable({
    df <- newData()[[2]]
    par <- c("Umax","Lmax","CQOF","CKIF","CK12","TOF","TIF","TG","CKBF","Csnow","Snowtemp")
    par_v <- c("Max W.C in the surface storage","Max W.C in root zone storage","Overland flow runoff coefficient","Time Constant for Interflow","Time constant for routing overland flow","Root zone treshold value for overland flow","Root zone treshold value for inter flow","Root zone treshold value for groundwater recharge","Time constant for routing base flow","Degree day coefficient","Snow melt temp")
    df <- data.frame(par,par_v,df)
    names(df) <- c("Parameter","Description","Value")
    df
  })
  
  output$stats = DT::renderDataTable({
    df <- newData()[[3]]
    metric <- c()
    value <- c()
    for(i in 1:length(df)){
      
      metric <- c(metric, df[[i]][[1]])
      value <- c(value, df[[i]][[2]])
      
    }
    
    par <- c("Agreement Index (d) developed by Willmott (1981)","BIAS","Correlation Coefficient","Covariance ","Decomposed MSE developed by Kobayashi and Salam (2000","Kling-Gupta Efficiency","Logarithmic probability distribution","Log Nash-Sutcliffe model efficiency","Mean Absolute Error","Mean Squared Error","Nash-Sutcliffe model efficinecy","Procentual Bias","Root Mean Squared Error","Relative Root Mean Squared Error","Coefficient of Determination","RMSE-observations standard deviation ratio","Volume Error (Ve)")
    df_stat <- data.frame(par, value) 
    names(df_stat) <- c("Metric","Value")
    
    df_stat
  })
  
  output$flowdur <- renderPlotly({
    
    df <- newData()[[4]]
    
    pt <- ggplot(df) + 
      geom_line(aes(x=Qsim_x,y = Qsim_y), color = "darkred") + 
      geom_line(aes(x=Qobs_x,y = Qobs_y), color="steelblue", linetype="twodash") +
      xlab("Percentage Exceedence (%)") + ylab("Discharge m$^3$/s") 
    
    fig <- ggplotly(pt)
    fig
    
  })
  
  # observe({
  #   x <- input$cal
  #   if (x == T ){
  #     value <- newData()[[2]]
  #     updateSliderInput(session, "umax",value = value[1])
  #     updateSliderInput(session, "lmax",value = value[2])
  #     updateSliderInput(session, "cqof",value = value[3])
  #     updateSliderInput(session, "ckif",value = value[4])
  #     updateSliderInput(session, "ck12",value = value[5])
  #   }
  # })
  
  
})



