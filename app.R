# BEMoDA shiny v1.0 - BiowaivEr aid for Model Dependent-Independent Approach application for in-vitro dissolution profile comparison
# 
# Model Dependent-Independent Approach application for in-vitro dissolution profile comparison as proposed Tsong et al. in 1996
# (Tsong Y, Hammerstrom T, Sathe P, Shah VP. (1996) Statistical Assessment of Mean Differences between Two Dissolution Data Sets, Drug Info. J. 30:1105-1112), and by Sathe et al. in 1996
# (Sathe PM, Tsong Y, Shah VP. In-vitro dissolution profile comparison: statistics and analysis, model dependent approach. Pharm Res. 1996 Dec;13(12):1799-803).
# 
# Copyright (C) 2017 Jakub Szlęk, Aleksander Mendyk
# 
# Authors: 
# Jakub Szlęk, Aleksander Mendyk
# 
# Affiliation: 
# Jagiellonian University Medical College,
# Faculty of Pharmacy,
# Department of Pharmaceucial Technology and Biopharmaceutics,
# Medyczna 9 st.,
# 30-688 Kraków
# Poland
# 
# Bugs, issues, please e-mail to maintainer
# Jakub Szlęk: j.szlek@uj.edu.pl
# 
# Copyright (C) 2017 Jakub Szlęk, Aleksander Mendyk
# 
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General 
# Public License as published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program.
# If not, see <http://www.gnu.org/licenses/>.
# 

require(shiny)
require(DT)
require(shinyBS)
require(dplyr)
require(optimx)
require(ggplot2)
require(nloptr)
require(GenSA)
require(MASS)
require(reshape)
require(gridExtra)

source("BEMoDA_auxiliary_functions.R")
source("BEMoDA_shiny_InDep.R")
source("BEMoDA_shiny_Dep.R")
source("BEMoDA_shiny_MANOVA.R")

#' busyIndicator START
#' 
#' busyIndicator from https://github.com/AnalytixWare/ShinySky/blob/master/R/busy-indicator.r by  xiaodaigh
#' 
#' A busy indicator
#' 
#' @param text The text to show
#' @param img An anitmated gif
#' @param wait The amount of time to wait before showing the busy indicator. The
#'   default is 1000 which is 1 second.
#'   
#' @export

busyIndicator <- function(text = "Calculation in progress..",img = "ajax-loader.gif", wait=1000) {
  tagList(
    singleton(tags$head(
      tags$link(rel="stylesheet", type="text/css",href="busyIndicator.css")
    ))
    ,div(class="shinysky-busy-indicator",p(text),img(src=img))
    ,tags$script(sprintf(
      "	setInterval(function(){
      if ($('html').hasClass('shiny-busy')) {
      setTimeout(function() {
      if ($('html').hasClass('shiny-busy')) {
      $('div.shinysky-busy-indicator').show()
      }
      }, %d)  		    
      } else {
      $('div.shinysky-busy-indicator').hide()
      }
},100)
      ",wait)
    )
  )	
}

#'
#'
#' busyIndicator END
#' 
#' 


# Define UI for application that draws a histogram
ui <- navbarPage(
  
  # Application title
  "BEMoDA_shiny v0.1",
  tabPanel("Home",
             titlePanel(h4("Introduction")),
             mainPanel(
               tags$br(),
               includeHTML("Help.html")
               )
             ),
  
  
             tabPanel("Reference",
                      
                    sidebarLayout(
                      sidebarPanel(titlePanel(h4(
                        "Data set for reference product")),
                        tags$hr(),
                        # Input: Select a file ----
                        fileInput("ref_file", "Choose CSV File",
                                  multiple = TRUE,
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),
                        
                        # Horizontal line ----
                        tags$hr(),
                        
                        # Input: Checkbox if file has header ----
                        checkboxInput("ref_header", "Header", TRUE),
                        
                        # Input: Select separator ----
                        radioButtons("ref_sep", "Separator",
                                     choices = c(Comma = ",",
                                                 Semicolon = ";",
                                                 Tab = "\t"),
                                     selected = ","),
                        
                        # Input: Select quotes ----
                        radioButtons("ref_quote", "Quote",
                                     choices = c(None = "",
                                                 "Double Quote" = '"',
                                                 "Single Quote" = "'"),
                                     selected = '"')
                        
                        
                      ),
                      
                      mainPanel(
                       DT::dataTableOutput("ref_file_info"),
                       plotOutput("ref_file_plot", click = "plot_click1")
                        )
                      )
                    ),
             tabPanel("Test",
                      
                    sidebarLayout(  
                      sidebarPanel(titlePanel(h4(
                        "Data set for test product")),
                        tags$hr(),
                        # Input: Select a file ----
                        fileInput("test_file", "Choose CSV File",
                                  multiple = TRUE,
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),
                        
                        # Horizontal line ----
                        tags$hr(),
                        
                        # Input: Checkbox if file has header ----
                        checkboxInput("test_header", "Header", TRUE),
                        
                        # Input: Select separator ----
                        radioButtons("test_sep", "Separator",
                                     choices = c(Comma = ",",
                                                 Semicolon = ";",
                                                 Tab = "\t"),
                                     selected = ","),
                        
                        # Input: Select quotes ----
                        radioButtons("test_quote", "Quote",
                                     choices = c(None = "",
                                                 "Double Quote" = '"',
                                                 "Single Quote" = "'"),
                                     selected = '"')
                        
                        
                        ),
                      
                        mainPanel(
                          DT::dataTableOutput("test_file_info"),
                          plotOutput("test_file_plot", click = "plot_click2")
                        )
                      
                      )

                  ),
             tabPanel("Standard",
                    
                    sidebarLayout(
                      sidebarPanel(titlePanel(h4(
                        "Data sets for standard batches")),
                        tags$hr(),
                        # Input: Select a file ----
                        h4("Standard batch 1"),
                        fileInput("std_file1", "Choose CSV File",
                                  multiple = TRUE,
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),
                        
                        # Input: Checkbox if file has header ----
                        checkboxInput("std_header1", "Header", TRUE),
                        
                        # Input: Select separator ----
                        radioButtons("std_sep1", "Separator",
                                     choices = c(Comma = ",",
                                                 Semicolon = ";",
                                                 Tab = "\t"),
                                     selected = ","),
                        
                        # Input: Select quotes ----
                        radioButtons("std_quote1", "Quote",
                                     choices = c(None = "",
                                                 "Double Quote" = '"',
                                                 "Single Quote" = "'"),
                                     selected = '"'),
                        
                        tags$hr(),
                        
                        # Input: Select a file ----
                        h4("Standard batch 2"),
                        fileInput("std_file2", "Choose CSV File",
                                  multiple = TRUE,
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),
                        
                        # Input: Checkbox if file has header ----
                        checkboxInput("std_header2", "Header", TRUE),
                        
                        # Input: Select separator ----
                        radioButtons("std_sep2", "Separator",
                                     choices = c(Comma = ",",
                                                 Semicolon = ";",
                                                 Tab = "\t"),
                                     selected = ","),
                        
                        # Input: Select quotes ----
                        radioButtons("std_quote2", "Quote",
                                     choices = c(None = "",
                                                 "Double Quote" = '"',
                                                 "Single Quote" = "'"),
                                     selected = '"'),
                        
                        tags$hr(),
                        # Input: Select a file ----
                        h4("Standard batch 3"),
                        fileInput("std_file3", "Choose CSV File",
                                  multiple = TRUE,
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),
                        
                        # Input: Checkbox if file has header ----
                        checkboxInput("std_header3", "Header", TRUE),
                        
                        # Input: Select separator ----
                        radioButtons("std_sep3", "Separator",
                                     choices = c(Comma = ",",
                                                 Semicolon = ";",
                                                 Tab = "\t"),
                                     selected = ","),
                        
                        # Input: Select quotes ----
                        radioButtons("std_quote3", "Quote",
                                     choices = c(None = "",
                                                 "Double Quote" = '"',
                                                 "Single Quote" = "'"),
                                     selected = '"')
                        
                      ),
                      
                      mainPanel(
                        DT::dataTableOutput("std_file_info1"),
                        plotOutput("std_file_plot1", click = "plot_click3"),
                        DT::dataTableOutput("std_file_info2"),
                        plotOutput("std_file_plot2", click = "plot_click4"),
                        DT::dataTableOutput("std_file_info3"),
                        plotOutput("std_file_plot3", click = "plot_click5")
                        )
                      )
                  ),
             
            navbarMenu("Profile comparison",
                  
                  tabPanel("Model dependent", fluid=TRUE,
                             
                               titlePanel(h4("Settings for model dependent approach")),
                                            
                               mainPanel(
                                 fluidRow(
                                   column(4,
                                          
                                          radioButtons("model_choose","Please choose mechanistic model which will be used to fit the dissolution profiles: ",
                                                       choices=c("Zero-order with Tlag", "Zero-order with F0",
                                                        "First-order with Tlag", "First-order with Fmax","Gompertz",
                                                        "Higuchi with Tlag", "Higuchi with F0", 
                                                        "Hixon-Crowell with Tlag", "Hopfenberg", "Korsmeyer-Peppas", "Logistic",
                                                        "Peppas-Sahlin", "Quadratic", "Weibull"
                                                        ),selected="Weibull"),
                                          
                                          tags$hr(),
                                          
                                          radioButtons("model_optim","Model parameters optimization method",choices=c("optimx","nloptr","genSA","nls"),selected="nloptr"),
                                          bsTooltip(id="model_optim", title=" optimx - function for optimization in R, which implements BFGS or Nelder-Mead, nloptr - nonlinear optimization, providing a common interface for a number of different free optimization routines, genSA - Generalized Simulated Annealing, nls - Nonlinear Least Square",
                                                    placement = "right", options = list(container = "body")),
                                          uiOutput("optimx_method"),
                                          uiOutput("optimx_iter"),
                                          uiOutput("optimx_toler"),
                                          uiOutput("nloptr_iter"),
                                          uiOutput("nloptr_toler"),
                                          uiOutput("genSA_iter"),
                                          uiOutput("nls_iter")
                                          
                                          ),
                                   
                                   column(3,
                                          sliderInput("sr_level","Similarity region (SR) confidence interval",min=0,max=1.0,step=0.01,value=0.99),
                                          sliderInput("cr_level","Critical region (CR) confidence interval",min=0,max=1.0,step=0.01,value=0.9),
                                          sliderInput("sr_pts","No of points detected at boundary of SR",min=0,max=500,step=5,value=50),
                                          bsTooltip(id="sr_pts", title="According to selected SR the algorithm detects points(x,y) at the boundary of SR",
                                                    placement = "right", options = list(container = "body")),
                                          sliderInput("cr_pts","No of points detected at boundary of CR",min=0,max=500,step=5,value=50),
                                          bsTooltip(id="cr_pts", title="According to selected CR the algorithm detects points(x,y) at the boundary of CR",
                                                    placement = "right", options = list(container = "body")),
                                          radioButtons("sr_region", "Draw SR as:",
                                                       choices = c("Rectangle",
                                                                   "Ellipse"),
                                                       selected = 'Ellipse'),
                                          checkboxInput(inputId = "log_transform", label = "Log transform parameters before CR an SR calculations", TRUE)
                                          ),
                                   
                                   column(4,
                                          h5(tags$b("Settings for nloptr to search for critical and similarity regions.")),
                                          numericInput("nloptr_reg_iter","No of steps in optimization of points(x,y) coordinates",max=100000,min = 0,value=4000),
                                          numericInput("nloptr_reg_toler","Relative tolerance of nloptr",max=1e-14,min = 1e-24,value=format(1e-20,scientific=TRUE))
                                          ),
                                   
                                   column(1,
                                          actionButton("BEMoDA_Dep","Run calculations"),
                                          busyIndicator()
                                          )
                                   
                                   
                                 )
                                 
                               )
                             
                           ),
                  
                  
                             
                    tabPanel("Mahalanobis distance", fluid=TRUE,
                             sidebarLayout(
                               sidebarPanel(titlePanel(h4("Settings for model independent approach.")),
                                            h5("Multivariate statistical distance method."),
                                            
                                            tags$hr(),
                                            numericInput("modindg", label="mean % of difference", value=10, min = 0, max = 100, step = 1,
                                                         width = 200),
                                            bsTooltip(id="modindg", title="dg is mean percent of difference in dissolution between each time points. typically 10. The assumption is that the SD is normally distributed across the variables",
                                                      placement = "right", options = list(container = "body")),
                                            tags$hr(),
                                            numericInput("modinprob", label="vector of probabilities", value=0.90, min = 0, max = 1, step = 0.01,
                                                         width = 200),
                                            bsTooltip(id="modinprob", title=" prob is the vector of probabilities calculated via qf() function in R which is density, distribution function, quantile function and random generation for the F distribution with df1() and df2() degrees of freedom and optional non-centrality parameter ncp",
                                                      placement = "right", options = list(container = "body")),
                                            actionButton("BEMoDA_InDep","Run calculations")
                                            
                                            ),
                               mainPanel(
                                 busyIndicator(),
                                 DT::dataTableOutput("BEMoDA_shiny_InDep_output")
                               )
                             )
                           ),
                  
                    tabPanel("MANOVA", fluid=TRUE,
                             sidebarLayout(
                               sidebarPanel(titlePanel(h4("Settings for model independent MANOVA repeated measures.")),
                                            h5("One way MANOVA"),
                                            
                                            tags$hr(),
                                            selectInput("manova_test", label="Test statistic", choices=c("Pillai",  "Wilks", "Hotelling-Lawley", "Roy"), selected = "Pillai",
                                                         width = 200),
                                            
                                            bsTooltip(id="manova_test", title="The name of the test statistic to be used.",
                                                      placement = "right", options = list(container = "body")),
                                            
                                            tags$hr(),
                                            checkboxInput("manova_intercept", label="Intercept", value=FALSE,
                                                         width = 200),
                                            bsTooltip(id="manova_intercept", title="If ‘TRUE’, the intercept term is included in the table.",
                                                      placement = "right", options = list(container = "body")),
                                            
                                            tags$hr(),
                                            numericInput("manova_tol", label="Tolerance", value = 1e-7, min = 0, max = 1, width = 200),
                                            bsTooltip(id="manova_intercept", title="If ‘TRUE’, the intercept term is included in the table.",
                                                      placement = "right", options = list(container = "body")),
                                            
                                            tags$hr(),
                                            actionButton("BEMoDA_MANOVA","Run calculations")
                                            
                               ),
                               mainPanel(
                                 busyIndicator(),
                                 verbatimTextOutput("BEMoDA_shiny_MANOVA_output")
                               )
                             )
                             )
                      
                  ),
  
            navbarMenu("Results",
                       tabPanel("Model fit",
                                mainPanel(
                                  tags$h4("Chosen model"),
                                  tags$h5(htmlOutput("eq_print")),
                                  tags$h4("Reference plots"),
                                  busyIndicator(),
                                  plotOutput("res_plot_ref", height = 840, width = 750),
                                  DT::dataTableOutput("res_ref_par"),
                                  tags$br(),
                                  tags$h4("Test plots."),
                                  busyIndicator(),
                                  plotOutput("res_plot_test", height = 840, width = 750),
                                  DT::dataTableOutput("res_test_par"),
                                  tags$br(),
                                  tags$h4("Standard batch no 1 plots."),
                                  busyIndicator(),
                                  tryCatch(plotOutput("res_plot_std1", height = 840, width = 750), error=function(e){print("Data wasn't loaded, please check the input files")}),
                                  tryCatch(DT::dataTableOutput("res_std_par1"), error=function(e){print("Data wasn't loaded, please check the input files")}),
                                  tags$br(),
                                  tags$h4("Standard batch no 2 plots."),
                                  busyIndicator(),
                                  tryCatch(plotOutput("res_plot_std2", height = 840, width = 750), error=function(e) print("Data wasn't loaded, please check the input files")),
                                  tryCatch(DT::dataTableOutput("res_std_par2"), error=function(e) print("Data wasn't loaded, please check the input files")),
                                  tags$br(),
                                  tags$h4("Standrad batch no 3 plots."),
                                  busyIndicator(),
                                  tryCatch(plotOutput("res_plot_std3", height = 840, width = 750), error=function(e) print("Data wasn't loaded, please check the input files")),
                                  tryCatch(DT::dataTableOutput("res_std_par3"), error=function(e) print("Data wasn't loaded, please check the input files")),
                                  tags$br()
                                  )
                                ),
                       tabPanel("SR & CR plots",
                                mainPanel(
                                  tags$h4("Similarity region (SR) and critical region (CR) plot"),
                                  busyIndicator(),
                                  plotOutput("res_plot_sr", height = 840, width = 750),
                                  htmlOutput("res_M_dist"),
                                  tags$br(),
                                  DT::dataTableOutput("res_s_cr"),
                                  tags$br(),
                                  DT::dataTableOutput("res_s_sr"),
                                  tags$br(),
                                  DT::dataTableOutput("res_log_mt"),
                                  tags$br(),
                                  DT::dataTableOutput("res_log_mr"),
                                  tags$br(),
                                  DT::dataTableOutput("res_log_ms"),
                                  tags$br(),
                                  DT::dataTableOutput("res_log_diff"),
                                  tags$br(),
                                  htmlOutput("res_min_max_sr_cr")
                                  
                                )
                                )
                          
                       ),
  
            tabPanel("About",
                       pre(includeText("README.TXT"))
            ),
  
            # tabPanel("License",
            #            h3("License, BEMoDA shiny v0.1"),
            #            br(""),
            #            p("BEMoDA shiny v0.1, shiny port for Bioequivalence Model Dependent-Independent Approach (BEMoDA v1.0) script", align="justified"),
            #            p("Copyright (C) 2018  Jakub Szlęk, Aleksander Mendyk, GPLv3",align="justified"),
            #            p("This program is free software: you can redistribute it and/or modify
            #                 it under the terms of the GNU General Public License as published by
            #                 the Free Software Foundation, either version 3 of the License, or
            #                 (at your option) any later version.",align="justified"),
            #           p("This program is distributed in the hope that it will be useful,
            #                 but WITHOUT ANY WARRANTY; without even the implied warranty of
            #                 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
            #                 GNU General Public License for more details.",align="justified"),
            #           p("You should have received a copy of the GNU General Public License
            #                 along with this program.  If not, see ", a("<https://www.gnu.org/licenses/gpl-3.0.html>", 
            #                                                            href = "http://www.gnu.org/licenses/gpl-3.0.html"),align="justified")
            # 
            #   ),

      
      
      # main panel
      mainPanel(
        # tableOutput(outputId = "loaded.ref.data")
        

      )
  
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  
  results <- reactiveValues()
  
  results$Dep <- NULL
  
  
  
  # the data read in by read.csv will now be accessible throughout the app
  # by calling this reactive, e.g. in_ref_data().
  in_ref_data <- reactive({
    shiny::validate(
      need(input$ref_file, "Select a file!")
    )
    
    tbl_ref <- read.csv(input$ref_file$datapath, header = input$ref_header,
             sep = input$ref_sep, quote = input$ref_quote, row.names = 1)
    
      return(tbl_ref)

  })
  
  
  
  # the data read in by read.csv will now be accessible throughout the app
  # by calling this reactive, e.g. in_ref_data().
  in_test_data <- reactive({
    shiny::validate(
      need(input$test_file, "Select a file!")
    )
    
    tbl_test <- read.csv(input$test_file$datapath, header = input$test_header,
                        sep = input$test_sep, quote = input$test_quote, row.names = 1)

      return(tbl_test)
    
  })
  
  
  # the data read in by read.csv will now be accessible throughout the app
  # by calling this reactive, e.g. in_ref_data().
  in_std_data1 <- reactive({
    shiny::validate(
      need(input$std_file1, "Select a file!")
    )
    
    tbl_std1 <- read.csv(input$std_file1$datapath, header = input$std_header1,
                        sep = input$std_sep1, quote = input$std_quote1, row.names = 1)

      return(tbl_std1)

  })
  
  in_std_data2 <- reactive({
    shiny::validate(
      need(input$std_file2, "Select a file!")
    )
    
    tbl_std2 <- read.csv(input$std_file2$datapath, header = input$std_header2,
                         sep = input$std_sep2, quote = input$std_quote2, row.names = 1)

      return(tbl_std2)
    
  })
  
  in_std_data3 <- reactive({
    shiny::validate(
      need(input$std_file3, "Select a file!")
    )
    
    tbl_std3 <- read.csv(input$std_file3$datapath, header = input$std_header3,
                         sep = input$std_sep3, quote = input$std_quote3, row.names = 1)
    
      return(tbl_std3)
  
  })
  
  
  # the reactive to access the data loaded
  output$ref_file_info <- DT::renderDataTable(datatable(in_ref_data(), colnames=c("Product"=1), caption = "Table: Reference product dissolution profiles.", options = list(
    searching = FALSE, autoWitdth=TRUE, bLengthChange=0, bInfo=0)) %>% formatRound(c(1:12), 2))
  output$test_file_info <- DT::renderDataTable(datatable(in_test_data(), colnames=c("Product"=1), caption = "Table: Test product dissolution profiles.", options = list(
    searching = FALSE, autoWitdth=TRUE, bLengthChange=0, bInfo=0)) %>% formatRound(c(1:12), 2))
  output$std_file_info1 <- DT::renderDataTable(datatable(in_std_data1(), colnames=c("Product"=1), caption = "Table: Standard batch 1 product dissolution profiles.", options = list(
    searching = FALSE, autoWitdth=TRUE, bLengthChange=0, bInfo=0)) %>% formatRound(c(1:12), 2))
  output$std_file_info2 <- DT::renderDataTable(datatable(in_std_data2(), colnames=c("Product"=1), caption = "Table: Standard batch 2 product dissolution profiles.", options = list(
    searching = FALSE, autoWitdth=TRUE, bLengthChange=0, bInfo=0)) %>% formatRound(c(1:12), 2))
  output$std_file_info3 <- DT::renderDataTable(datatable(in_std_data3(), colnames=c("Product"=1), caption = "Table: Standard batch 3 product dissolution profiles.", options = list(
    searching = FALSE, autoWitdth=TRUE, bLengthChange=0, bInfo=0)) %>% formatRound(c(1:12), 2))
  
              # BEMoDA_InDep run button
               BEMoDA_InDep_df <-eventReactive(input$BEMoDA_InDep, {
                 BEMoDA_shiny_InDep(input$ref_file$datapath,input$test_file$datapath,input$modindg,input$modinprob)
               })
               
               output$BEMoDA_shiny_InDep_output <- DT::renderDataTable({datatable(BEMoDA_InDep_df(), rownames=FALSE, caption= "Table: Results of model independent approach.", options = list(
                 searching = FALSE, autoWitdth=TRUE, bLengthChange=0, bInfo=0)) %>% formatRound(c(1:3,5,6), 2)})
               
               #BEMoDA_MANOVA run button
               
               BEMoDA_MANOVA_output <- eventReactive(input$BEMoDA_MANOVA, {
                 BEMoDA_MANOVA(input$ref_file$datapath, input$test_file$datapath, input$manova_test, input$manova_intercept,
                               input$manova_tol)
               })
               
               output$BEMoDA_shiny_MANOVA_output <- renderPrint({
                 print(BEMoDA_MANOVA_output())
                 })
               

               
               # model dependent settings for optimx
               
               output$optimx_method <- renderUI({
                 if(input$model_optim == "optimx"){
                     selectInput("optimx_method", "Optimx method",
                                 list("BFGS", "Nelder-Mead"))
                 }
               })
               
               output$optimx_iter <- renderUI({
                 if(input$model_optim == "optimx"){
                   numericInput("optimx_iter", "Number of maximum iteriations",value=10000,min=0,max=100000)
                 }
               })
               
               output$optimx_toler <- renderUI({
                 if(input$model_optim == "optimx"){
                   numericInput("optimx_toler", "Optimization relative tolerance",value=format(1e-20, scientific = TRUE), min=1e-24,max=1e-14)
                 }
               })
               
               # model dependent settings for nloptr
               
               output$nloptr_iter <- renderUI({
                 if(input$model_optim == "nloptr"){
                   numericInput("nloptr_iter", "Number of maximum iteriations",value=10000,min=0,max=100000)
                 }
               })
               
               output$nloptr_toler <- renderUI({
                 if(input$model_optim == "nloptr"){
                   numericInput("nloptr_toler", "Optimization relative tolerance",value=format(1e-20, scientific = TRUE), min=1e-24)
                 }
               })
               
               
               # model dependent settings for genSA
               
               output$genSA_iter <- renderUI({
                 if(input$model_optim == "genSA"){
                   numericInput("genSA_iter", "Number of maximum iteriations",value=5000,min=0,max=50000)
                 }
               })
               
               
               # model dependent settings for nloptr
               
               output$nls_iter <- renderUI({
                 if(input$model_optim == "nls"){
                   numericInput("nls_iter", "Number of maximum iteriations",value=5000,min=0,max=10000)
                 }
               })
               
               
               # Plot loaded data
               output$ref_file_plot <- renderPlot({
                 data.tmp <- in_ref_data()
                 data.time <- extract.timepoints(data.tmp)
                 colnames(data.tmp) <- data.time
                 data.tmp[ "Prod" ] <- rownames(data.tmp)
                 df.molten <- melt(data.tmp, id.vars="Prod")
                 tmp.plot <- ggplot(df.molten, aes(x = variable, y = value)) + geom_line(aes(color = Prod, group = Prod))
                 return(tmp.plot)
               })
               
               output$test_file_plot <- renderPlot({
                 data.tmp <- in_test_data()
                 data.time <- extract.timepoints(data.tmp)
                 colnames(data.tmp) <- data.time
                 data.tmp[ "Prod" ] <- rownames(data.tmp)
                 df.molten <- melt(data.tmp, id.vars="Prod")
                 tmp.plot <- ggplot(df.molten, aes(x = variable, y = value)) + geom_line(aes(color = Prod, group = Prod))
                 return(tmp.plot)
               })
               
               output$std_file_plot1 <- renderPlot({
                 data.tmp <- in_std_data1()
                 data.time <- extract.timepoints(data.tmp)
                 colnames(data.tmp) <- data.time
                 data.tmp[ "Prod" ] <- rownames(data.tmp)
                 df.molten <- melt(data.tmp, id.vars="Prod")
                 tmp.plot <- ggplot(df.molten, aes(x = variable, y = value)) + geom_line(aes(color = Prod, group = Prod))
                 return(tmp.plot)
               })
               
               output$std_file_plot2 <- renderPlot({
                 data.tmp <- in_std_data2()
                 data.time <- extract.timepoints(data.tmp)
                 colnames(data.tmp) <- data.time
                 data.tmp[ "Prod" ] <- rownames(data.tmp)
                 df.molten <- melt(data.tmp, id.vars="Prod")
                 tmp.plot <- ggplot(df.molten, aes(x = variable, y = value)) + geom_line(aes(color = Prod, group = Prod))
                 return(tmp.plot)
               })
               
               output$std_file_plot3 <- renderPlot({
                 data.tmp <- in_std_data3()
                 data.time <- extract.timepoints(data.tmp)
                 colnames(data.tmp) <- data.time
                 data.tmp[ "Prod" ] <- rownames(data.tmp)
                 df.molten <- melt(data.tmp, id.vars="Prod")
                 tmp.plot <- ggplot(df.molten, aes(x = variable, y = value)) + geom_line(aes(color = Prod, group = Prod))
                 return(tmp.plot)
               })
               
               
               
               observeEvent(input$BEMoDA_Dep, {
                 
                 ref <- in_ref_data()
                 test <- in_test_data()
                 std1 <- try(in_std_data1())
                 std2 <- try(in_std_data2())
                 std3 <- try(in_std_data3())
                 std.list <- list(std1, std2, std3)
                 std.list <- std.list[sapply(std.list, function(x) !inherits(x, "try-error"))]
                 
                 dissol.model <- input$model_choose
                 
                 if(input$model_optim == "optimx"){
                 
                   optim.model.params <- list(method=input$model_optim, optimx_method=input$optimx_method, maxit=input$optimx_iter, toler=input$optimx_toler)
                 
                 } else if (input$model_optim == "nloptr"){
                   
                   optim.model.params <- list(method=input$model_optim, maxit=input$nloptr_iter, toler=input$nloptr_toler)
                   
                 } else if (input$model_optim == "genSA"){
                   
                   optim.model.params <- list(method=input$model_optim, maxit=input$genSA_iter)
                   
                 } else if (input$model_optim == "nls"){
                   
                   optim.model.params <- list(method=input$model_optim, maxit=input$nls_iter)
                   
                 }
                 
                 optim.ellipse.params <-list(sr_level=input$sr_level, cr_level=input$cr_level, sr_pts=input$sr_pts, cr_pts=input$cr_pts,
                                             iter=input$nloptr_reg_iter, toler=input$nloptr_reg_toler, log_trans=input$log_transform) 
                 
                 results$Dep <- BEMoDA_shiny_Dep(ref, test, std.list, optim.model.params, dissol.model, optim.ellipse.params)
                 
               })
               
               
               output$res_plot_ref <- renderPlot({
                 
                 plots <- results$Dep[[1]]
                 numPlots <-length(plots)
                 nCol <- floor(sqrt(numPlots))
                 do.call("grid.arrange", c(plots, ncol=nCol))

                })
               
               output$res_plot_test <- renderPlot({
                 
                 plots <- results$Dep[[2]]
                 numPlots <-length(plots)
                 nCol <- floor(sqrt(numPlots))
                 do.call("grid.arrange", c(plots, ncol=nCol))
                 
               })
               
               output$res_plot_std1 <- tryCatch(renderPlot({
                 
                 plots <- results$Dep[[3]][[1]]
                 numPlots <-length(plots)
                 nCol <- floor(sqrt(numPlots))
                 do.call("grid.arrange", c(plots, ncol=nCol))
                 
               }), error=function(e) print("Data wasn't provided!"))
               
               output$res_plot_std2 <- tryCatch(renderPlot({
                 
                 plots <- try(results$Dep[[3]][[2]])
                 numPlots <-try(length(plots))
                 nCol <- try(floor(sqrt(numPlots)))
                 try(do.call("grid.arrange", c(plots, ncol=nCol)))
                 
               }), error=function(e) print("Data wasn't provided!"))
               
               output$res_plot_std3 <- tryCatch(renderPlot({
                 
                 plots <- try(results$Dep[[3]][[3]])
                 numPlots <-try(length(plots))
                 nCol <- try(floor(sqrt(numPlots)))
                 try(do.call("grid.arrange", c(plots, ncol=nCol)))
                 
               }), error=function(e) print("Data wasn't provided!"))
               
               # access to the results
               output$res_std_par1 <- try(DT::renderDataTable(datatable(results$Dep[[12]][[1]], colnames=c("Tablet"=1), caption = "Table: Standard batch 1 parameters.", options = list(
                 searching = FALSE, autoWitdth=TRUE, bLengthChange=0, bInfo=0)) %>% formatRound(c(1:4), 3)))
               output$res_std_par2 <- try(DT::renderDataTable(datatable(results$Dep[[12]][[2]], colnames=c("Tablet"=1), caption = "Table: Standard batch 2 parameters.", options = list(
                 searching = FALSE, autoWitdth=TRUE, bLengthChange=0, bInfo=0)) %>% formatRound(c(1:4), 3)))
               output$res_std_par3 <- try(DT::renderDataTable(datatable(results$Dep[[12]][[3]], colnames=c("Tablet"=1), caption = "Table: Standard batch 3 parameters.", options = list(
                 searching = FALSE, autoWitdth=TRUE, bLengthChange=0, bInfo=0)) %>% formatRound(c(1:4), 3)))
               output$res_ref_par <- DT::renderDataTable(datatable(results$Dep[[10]], colnames=c("Tablet"=1), caption = "Table: Reference product parameters.", options = list(
                 searching = FALSE, autoWitdth=TRUE, bLengthChange=0, bInfo=0)) %>% formatRound(c(2:5), 3))
               output$res_test_par <- DT::renderDataTable(datatable(results$Dep[[11]], colnames=c("Tablet"=1), caption = "Table: Test product paramters.", options = list(
                 searching = FALSE, autoWitdth=TRUE, bLengthChange=0, bInfo=0)) %>% formatRound(c(2:5), 3))
               
               # Plot SR and CR
               output$res_plot_sr <- renderPlot({
                 
                 how.draw.sr.region <- input$sr_region
                 
                 ####################################################################                            
                 #####     PLOT ELLIPSE BASED ON THE POINTS GENERETED BY NLOPTR  ####
                 ####################################################################
                 
                 ellipse.cr <- results$Dep[[7]]
                 ellipse.sr <- results$Dep[[8]]
                 rect.sr1 <- results$Dep[[9]][[1]]
                 rect.sr2 <- results$Dep[[9]][[2]]
                 rect.sr3 <- results$Dep[[9]][[3]]
                 
                 ellipse1 <- draw.ellipse(ellipse.cr[,1], ellipse.cr[,2])
                 ellipse2 <- draw.ellipse(ellipse.sr[,1], ellipse.sr[,2])
                 ellipse1.pts <- contourLines(ellipse1$u, ellipse1$v, ellipse1$z, levels=0)
                 ellipse2.pts <- contourLines(ellipse2$u, ellipse2$v, ellipse2$z, levels=0)
                 
                 if(how.draw.sr.region == "Ellipse"){
                 
                   xlim <- c(min=min(ellipse.cr[,1],ellipse.sr[,1]), max=max(ellipse.cr[,1],ellipse.sr[,1]))
                   ylim <- c(min=min(ellipse.cr[,2],ellipse.sr[,2]), max=max(ellipse.cr[,2],ellipse.sr[,2]))
                 
                 } else if(how.draw.sr.region == "Rectangle"){
                   
                   xlim <- c(min=min(ellipse.cr[,1],ellipse.sr[,1],rect.sr3$x2), max=max(ellipse.cr[,1],ellipse.sr[,1], rect.sr3$x1))
                   ylim <- c(min=min(ellipse.cr[,2],ellipse.sr[,2], rect.sr3$y2), max=max(ellipse.cr[,2],ellipse.sr[,2],rect.sr3$y1))
                   
                 }
                   
                 # Plot and save resulting graph!
                 
                 # png("output.png", width=1800, height=1800, res=300)
                 if(input$log_transform == TRUE){
                    
                   xlabel <- paste("Log-normal (param_A)")
                   ylabel <- paste("Log-normal (param_B)")
                 
                 } else if (input$log_transform == FALSE){
                   
                   xlabel <- paste("param_A")
                   ylabel <- paste("param_B")
                   
                 }
                 
                 contour(ellipse1$u, ellipse1$v, ellipse1$z, levels=0, lwd=2, xlab=xlabel, ylab=ylabel, asp=1,col="red", xlim=c(xlim[1]*1.10,xlim[2]*1.10),ylim=c(ylim[1]*1.10,ylim[2]*1.10), label="test-ref CR", labcex=1.1)
                 
                 
                 if(how.draw.sr.region == "Ellipse"){
                   contour(ellipse2$u, ellipse2$v, ellipse2$z, levels=0, lwd=2, xlab=xlabel, ylab=ylabel, asp=1,col="blue",add=TRUE, label=paste("SR, ", input$sr_level), labcex=1.1)
                 }
                 
                 
                 if(how.draw.sr.region == "Rectangle"){
                   rect(xleft=rect.sr1$x2, ybottom=rect.sr1$y2, xright=rect.sr1$x1, ytop=rect.sr1$y1, border = "black", lty=1) # 1 SD
                   rect(xleft=rect.sr2$x2, ybottom=rect.sr2$y2, xright=rect.sr2$x1, ytop=rect.sr2$y1, border = "black", lty=2) # 2 SD
                   rect(xleft=rect.sr3$x2, ybottom=rect.sr3$y2, xright=rect.sr3$x1, ytop=rect.sr3$y1, border = "black", lty=3) # 3 SD
                   
                 }
                 
                 # dev.off()

               })
               
               output$res_M_dist <- renderUI({
                 M.dist.cr <- results$Dep[[15]]
                 M.dist.cr.scaled <- results$Dep[[16]]
                 str1 <- paste("MAHALANOBIS DISTANCE (D^2) FOR TEST-REFERENCE:")
                 str2 <- paste(M.dist.cr)
                 str3 <- paste("")
                 str4 <- paste("SCALED MAHALANOBIS DISTANCE (HOTTELING'S T^2) FOR TEST-REFERENCE:")
                 str5 <- paste(M.dist.cr.scaled)
                 
                 HTML(paste(str1, str2, str3, str4, str5, sep = '<br/>'))
                 
               })
               
               # Print-out model equation, units and parameters
               output$eq_print <- renderUI({
                 
                 
                 if(input$model_choose == "Weibull"){
                   
                   eq_1 <- paste("Q = 100*(1-exp(-param_A*(time^(param_B))))")
                   eq_2 <- paste("Q = [%]")
                   eq_3 <- paste("time = [hrs]")
                   eq_4 <- paste("param_A = alpha")
                   eq_5 <- paste("param_B = beta")
                   eq_6 <- paste("Reference: Sathe PM, Tsong Y, Shah VP. In-vitro dissolution profile comparison: statistics and analysis, model dependent approach. Pharm Res. 1996;13:1799–803. doi: 10.1023/A:1016020822093.")
                   # y=100*(1-e^(-alpha*((t)^beta)))
                   # Sathe PM, Tsong Y, Shah VP. In-vitro dissolution profile comparison: statistics and analysis, model dependent approach. Pharm Res. 1996;13:1799–803. doi: 10.1023/A:1016020822093.
                   
                 } else if (input$model_choose == "Korsmeyer-Peppas"){
                   
                   eq_1 <- paste("Q = 100*param_A*time^param_B")
                   eq_2 <- paste("Q = [%]")
                   eq_3 <- paste("time = [hrs]")
                   eq_4 <- paste("param_A = k")
                   eq_5 <- paste("param_B = n")
                   eq_6 <- paste("Reference: Korsmeyer RW, Gurny R, Doelker E, Buri P, Peppas NA. Mechanisms of solute release from porous hydrophilic polymers. Int J Pharm. 1983;15:25–35. doi: 10.1016/0378-5173(83)90064-9")
                   # y=k*t^n
                   # Korsmeyer RW, Gurny R, Doelker E, Buri P, Peppas NA. Mechanisms of solute release from porous hydrophilic polymers. Int J Pharm. 1983;15:25–35. doi: 10.1016/0378-5173(83)90064-9
                   
                 } else if(input$model_choose == "Peppas-Sahlin"){
                   
                   eq_1 <- paste("Q = (alpha*time^(0.5)+beta*t)")
                   eq_2 <- paste("Q = [%]")
                   eq_3 <- paste("time = [min]")
                   eq_4 <- paste("param_A = k1")
                   eq_5 <- paste("param_B = k2")
                   eq_6 <- paste("Reference: Peppas NA, Sahlin JJ. A simple equation for the description of solute release III. Coupling of diffusion and relaxation. Int J Pharm. 1989;57:169–72. doi: 10.1016/0378-5173(89)90306-2.")
                   
                   # y=k1*t^(0.5)+k2*t
                   # Peppas NA, Sahlin JJ. A simple equation for the description of solute release III. Coupling of diffusion and relaxation. Int J Pharm. 1989;57:169–72. doi: 10.1016/0378-5173(89)90306-2.
                   
                 } else if (input$model_choose=="Quadratic"){
                   
                   eq_1 <- paste("Q = 100*(param_A*time^2+param_B*time)")
                   eq_2 <- paste("Q = [%]")
                   eq_3 <- paste("time = [min]")
                   eq_4 <- paste("param_A = k1")
                   eq_5 <- paste("param_B = k2")
                   eq_6 <- paste("Reference: Costa P, Sousa Lobo JM. Modeling and comparison of dissolution profiles. Eur J Pharm Sci. 2001;13:123–33. doi: 10.1016/S0928-0987(01)00095-1.")
                   
                   # y=100*(k1*t^2+k2*t)
                   # Costa P, Sousa Lobo JM. Modeling and comparison of dissolution profiles. Eur J Pharm Sci. 2001;13:123–33. doi: 10.1016/S0928-0987(01)00095-1.
                   
                 } else if(input$model_choose == "Logistic"){
                   
                   eq_1 <- paste("Q = 100*((exp(param_A + param_B*log10(time)))/(1+exp(param_A + param_B*log10(time))))")
                   eq_2 <- paste("Q = [%]")
                   eq_3 <- paste("time = [min]")
                   eq_4 <- paste("param_A = alpha")
                   eq_5 <- paste("param_B = beta")
                   eq_6 <- paste("Reference: Sathe PM, Tsong Y, Shah VP. In-vitro dissolution profile comparison: statistics and analysis, model dependent approach. Pharm Res. 1996;13:1799–803. doi: 10.1023/A:1016020822093.")
                   
                   # y=100*((e^(alpha+beta*log(t)))/(1+e^(alpha+beta*log(t))))
                   # Sathe PM, Tsong Y, Shah VP. In-vitro dissolution profile comparison: statistics and analysis, model dependent approach. Pharm Res. 1996;13:1799–803. doi: 10.1023/A:1016020822093.
                   
                 } else if(input$model_choose == "Gompertz"){
                   
                   eq_1 <- paste("Q = 100*exp(-param_A*exp(-param_B*log10(time)))")
                   eq_2 <- paste("Q = [%]")
                   eq_3 <- paste("time = [min]")
                   eq_4 <- paste("param_A = alpha")
                   eq_5 <- paste("param_B = beta")
                   eq_6 <- paste("Reference: Tsong Y, Hammerstrom T, Chen JJ. Multipoint dissolution specification and acceptance sampling rule based on profile modeling and principal component analysis. J Biopharm Stat. 1997;7:423–39.")
                   
                   # y=100*e^(-alpha*e^(-beta*log(t)))
                   # Tsong Y, Hammerstrom T, Chen JJ. Multipoint dissolution specification and acceptance sampling rule based on profile modeling and principal component analysis. J Biopharm Stat. 1997;7:423–39.
                   
                   # } else if (which.model=="Probit"){ 
                   #   fi <- pnorm(time,mean=mean(time),sd=sd(time)) # should be normal distribution function of time, but it seems that smth is wrong here =========> switch-off
                   #   diss <- 100*fi*(alpha+beta*log(t))
                   #   # y = 100*fi*(alpha+beta*log(t))
                   #   # Tsong Y, Hammerstrom T, Chen JJ. Multipoint dissolution specification and acceptance sampling rule based on profile modeling and principal component analysis. J Biopharm Stat. 1997;7:423–39.
                 }  else if (input$model_choose == "Zero-order with Tlag"){
                   
                   eq_1 <- paste("Q = param_A*(time - param_B)")
                   eq_2 <- paste("Q = [%]")
                   eq_3 <- paste("time = [min]")
                   eq_4 <- paste("param_A = k0")
                   eq_5 <- paste("param_B = T_lag")
                   eq_6 <- paste("Reference: Borodkin S, Tucker FE. Linear drug release from laminated hydroxypropyl cellulose-polyvinyl acetate films. J Pharm Sci. 1975;64:1289–94.")
                   
                   # y = alpha*(t - beta)
                   # Borodkin S, Tucker FE. Linear drug release from laminated hydroxypropyl cellulose-polyvinyl acetate films. J Pharm Sci. 1975;64:1289–94.
                 } else if(input$model_choose == "Zero-order with F0"){
                   
                   eq_1 <- paste("Q = param_A + param_B*time")
                   eq_2 <- paste("Q = [%]")
                   eq_3 <- paste("time = [min]")
                   eq_4 <- paste("param_A = F0")
                   eq_5 <- paste("param_B = k0")
                   eq_6 <- paste("Reference: Costa P, Sousa Lobo JM. Modeling and comparison of dissolution profiles. Eur J Pharm Sci. 2001;13:123–33.")
                   
                   # y = alpha + beta*t
                   # Costa P, Sousa Lobo JM. Modeling and comparison of dissolution profiles. Eur J Pharm Sci. 2001;13:123–33.
                 } else if(input$model_choose == "First-order with Tlag"){
                   
                   eq_1 <- paste("Q = 100*(1-exp(-param_A*(time - param_B)))")
                   eq_2 <- paste("Q = [%]")
                   eq_3 <- paste("time = [min]")
                   eq_4 <- paste("param_A = k1")
                   eq_5 <- paste("param_B = T_lag")
                   eq_6 <- paste("Reference: Phaechamud T, Pitaksantayothin K, Kositwattanakoon P, Seehapong P, Jungvivatanavong S. Sustainable release of propranolol hydrochloride tablet using chitin as press-coating material. Silpakorn Univ Int J. 2002;2:147–59.")
                   
                   # y = 100*(1-exp(-alpha*(t-beta)))
                   # Phaechamud T, Pitaksantayothin K, Kositwattanakoon P, Seehapong P, Jungvivatanavong S. Sustainable release of propranolol hydrochloride tablet using chitin as press-coating material. Silpakorn Univ Int J. 2002;2:147–59.
                 } else if (input$model_choose == "First-order with Fmax"){
                   
                   eq_1 <- paste("Q = param_A*(1-exp(-param_B*time))")
                   eq_2 <- paste("Q = [%]")
                   eq_3 <- paste("time = [min]")
                   eq_4 <- paste("param_A = F_max")
                   eq_5 <- paste("param_B = k1")
                   eq_6 <- paste("Reference: Tsong Y, Hammerstrom T, Chen JJ. Multipoint dissolution specification and acceptance sampling rule based on profile modeling and principal component analysis. J Biopharm Stat. 1997;7:423–39.")
                   
                   # y = alpha*(1-exp(-beta*t))
                   # Tsong Y, Hammerstrom T, Chen JJ. Multipoint dissolution specification and acceptance sampling rule based on profile modeling and principal component analysis. J Biopharm Stat. 1997;7:423–39.
                 } else if(input$model_choose == "Higuchi with Tlag"){
                   
                   eq_1 <- paste("Q = param_A*(time - param_B)^0.5")
                   eq_2 <- paste("Q = [%]")
                   eq_3 <- paste("time = [min]")
                   eq_4 <- paste("param_A = k_H")
                   eq_5 <- paste("param_B = T_lag")
                   eq_6 <- paste("Reference: Tarvainen M, Peltonen S, Mikkonen H, Elovaara M, Tuunainen M, Paronen P et al. Aqueous starch acetate dispersion as a novel coating material for controlled release products. J Control Release. 2004;96:179–91.")
                   
                   # y = alpha*(t-beta)^0.5
                   # Tarvainen M, Peltonen S, Mikkonen H, Elovaara M, Tuunainen M, Paronen P et al. Aqueous starch acetate dispersion as a novel coating material for controlled release products. J Control Release. 2004;96:179–91.
                 } else if(input$model_choose == "Higuchi with F0"){
                   
                   eq_1 <- paste("Q = param_A + param_B*(time)^0.5")
                   eq_2 <- paste("Q = [%]")
                   eq_3 <- paste("time = [min]")
                   eq_4 <- paste("param_A = F0")
                   eq_5 <- paste("param_B = k_H")
                   eq_6 <- paste("Reference: Ford JL, Mitchell K, Rowe P, Armstrong DJ, Elliott PNC, Rostron C et al. Mathematical modelling of drug release from hydroxypropylmethylcellulose matrices: effect of temperature. Int J Pharm. 1991;71:95–104.")
                   
                   # y = alpha + beta*t^0.5
                   # Ford JL, Mitchell K, Rowe P, Armstrong DJ, Elliott PNC, Rostron C et al. Mathematical modelling of drug release from hydroxypropylmethylcellulose matrices: effect of temperature. Int J Pharm. 1991;71:95–104.
                 } else if(input$model_choose == "Hixon-Crowell with Tlag"){
                   
                   eq_1 <- paste("Q = 100*(1 - (1 - param_A*(time - param_B))^3)")
                   eq_2 <- paste("Q = [%]")
                   eq_3 <- paste("time = [min]")
                   eq_4 <- paste("param_A = k_HC")
                   eq_5 <- paste("param_B = T_lag")
                   eq_6 <- paste("Reference: Mollo AR, Corrigan OI. An investigation of the mechanism of release of the amphoteric drug amoxycillin from poly( D , L -lactide-co-glycolide) matrices. Pharm Dev Technol. 2002;7:333–43.")
                   
                   # y = 100*(1-(1-alpha*(t-beta))^3)
                   # Mollo AR, Corrigan OI. An investigation of the mechanism of release of the amphoteric drug amoxycillin from poly( D , L -lactide-co-glycolide) matrices. Pharm Dev Technol. 2002;7:333–43.
                 } else if(input$model_choose == "Hopfenberg"){
                   
                   eq_1 <- paste("Q = 100*(1 - (1 - param_A*time)^param_B)")
                   eq_2 <- paste("Q = [%]")
                   eq_3 <- paste("time = [min]")
                   eq_4 <- paste("param_A = k_HB")
                   eq_5 <- paste("param_B = n")
                   eq_6 <- paste("Reference: Enscore DJ, Hopfenberg HB, Stannett VT. Effect of particle size on the mechanism controlling n-hexane sorption in glassy polystyrene microspheres. Polymer. 1977;18:793–800.")
                   
                   # y = 100*(1-(1-alpha*t)^beta)
                   # Enscore DJ, Hopfenberg HB, Stannett VT. Effect of particle size on the mechanism controlling n-hexane sorption in glassy polystyrene microspheres. Polymer. 1977;18:793–800.
                 }
                 
                 HTML(paste(eq_1, eq_2, eq_3, eq_4, eq_5, eq_6, sep = '<br/>'))
                 
                 
               })

               output$res_s_cr <- try(DT::renderDataTable(datatable(results$Dep[[13]],
                                                                    caption = "VARIANCE/COVARIANCE MATRIX FOR CRITICAL REGION CALCULATIONS:",
                                                                    options = list(searching = FALSE, autoWitdth=TRUE, bLengthChange=0, bInfo=0))
                                                          %>% formatRound(c(1:2), 8)))
               output$res_s_sr <- try(DT::renderDataTable(datatable(results$Dep[[14]],
                                                                    caption = "VARIANCE/COVARIANCE MATRIX FOR SIMILARITY REGION CALCULATIONS:",
                                                                    options = list(searching = FALSE, autoWitdth=TRUE, bLengthChange=0, bInfo=0))
                                                          %>% formatRound(c(1:2), 8)))
               output$res_log_mt <- try(DT::renderDataTable(datatable(results$Dep[[17]],
                                                                    caption = if(input$log_transform){
                                                                      paste("MEAN LOG-NORMAL VALUES FOR TEST:")
                                                                      } else {
                                                                        paste("MEAN VALUES FOR TEST:")
                                                                      },
                                                                    options = list(searching = FALSE, autoWitdth=TRUE, bLengthChange=0, bInfo=0))
                                                          %>% formatRound(c(1:2), 8)))
               output$res_log_mr <- try(DT::renderDataTable(datatable(results$Dep[[18]],
                                                                      caption = if(input$log_transform){
                                                                        paste("MEAN LOG-NORMAL VALUES FOR REFERENCE:")
                                                                      } else {
                                                                        paste("MEAN VALUES FOR REFERENCE:")
                                                                      },
                                                                      options = list(searching = FALSE, autoWitdth=TRUE, bLengthChange=0, bInfo=0))
                                                            %>% formatRound(c(1:2), 8)))
               output$res_log_ms <- try(DT::renderDataTable(datatable(results$Dep[[19]],
                                                                      caption = if(input$log_transform){
                                                                        paste("MEAN LOG-NORMAL VALUES FOR STANDARD BATCHES:")
                                                                      } else {
                                                                        paste("MEAN VALUES FOR STANDARD BATCHES:")
                                                                      },
                                                                      options = list(searching = FALSE, autoWitdth=TRUE, bLengthChange=0, bInfo=0))
                                                            %>% formatRound(c(2:3), 8)))
               output$res_log_diff <- try(DT::renderDataTable(datatable(results$Dep[[20]],
                                                                      caption = if(input$log_transform){
                                                                        paste("MEAN LOG-NORMAL DIFFERENCE (TEST-REF):")
                                                                      } else {
                                                                        paste("MEAN DIFFERENCE (TEST-REF):")
                                                                      },
                                                                      options = list(searching = FALSE, autoWitdth=TRUE, bLengthChange=0, bInfo=0))
                                                            %>% formatRound(c(1:2), 8)))

               output$res_min_max_sr_cr <- renderUI({
                 ellipse.cr <- results$Dep[[7]]
                 ellipse.sr <- results$Dep[[8]]
                 rect.sr1 <- results$Dep[[9]][[1]]
                 rect.sr2 <- results$Dep[[9]][[2]]
                 rect.sr3 <- results$Dep[[9]][[3]]

                 ellipse1 <- draw.ellipse(ellipse.cr[,1], ellipse.cr[,2])
                 ellipse2 <- draw.ellipse(ellipse.sr[,1], ellipse.sr[,2])
                 ellipse1.pts <- contourLines(ellipse1$u, ellipse1$v, ellipse1$z, levels=0)
                 ellipse2.pts <- contourLines(ellipse2$u, ellipse2$v, ellipse2$z, levels=0)
                 
                 str1 <- paste("CRITICAL REGION MIN-MAX")
                 str2 <- paste("LOG-NORMAL PARAM_A")
                 str3 <- paste("(",min(ellipse1.pts[[1]]$x),", ", max(ellipse1.pts[[1]]$x) ,")")
                 str4 <- paste("LOG-NORMAL PARAM_B")
                 str5 <- paste("(",min(ellipse1.pts[[1]]$y),", ", max(ellipse1.pts[[1]]$y) ,")")
                 str6 <- paste("")
                 str7 <- paste("SIMILARITY REGION MIN-MAX")
                 str8 <- paste("LOG-NORMAL PARAM_A")
                 str9 <- paste("(",min(ellipse2.pts[[1]]$x),", ", max(ellipse2.pts[[1]]$x) ,")")
                 str10 <- paste("LOG-NORMAL PARAM_B")
                 str11 <- paste("(",min(ellipse2.pts[[1]]$y),", ", max(ellipse2.pts[[1]]$y) ,")")
                 
                 if (input$log_transform == FALSE){
                   
                   str2 <- paste("PARAM_A")
                   str4 <- paste("PARAM_B")
                   str8 <- paste("PARAM_A")
                   str10 <- paste("PARAM_B")
                   
                 }
                 
                 
                 HTML(paste(str1, str2, str3, str4, str5, str6, str7, str8, str9, str10, str11, sep = '<br/>'))
                 
               })
               

}

# Run the application 
shinyApp(ui = ui, server = server)

