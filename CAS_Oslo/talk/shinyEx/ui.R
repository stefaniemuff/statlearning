library(shiny)

# setwd("/home/steffi/Shiny/MEApp/")
# runApp(MEC_ChooseL)

shinyUI(fluidPage(
  titlePanel("Classical measurement error in linear, logit and Poisson regression"),
  
  sidebarLayout(
    sidebarPanel(
#       radioButtons("regMod", 
#                          label = "Select regression model:", 
#                          choices = list("Linear" = 1, 
#                                         "Logistic" = 2, "Poisson" = 3),
#                          selected = 1),
      selectInput("regMod", label = "Select regression model:", 
                  choices = list("Linear" = 1, "Logistic" = 2,
                                 "Poisson (log link)" = 3), 
                  selected = 1),
      
      p("You can now check the effect of classical measurement error in a covariate of linear, logistic and Poisson regression. The linear predictor of the model is given as",
      #uiOutput("tex1"),
      withMathJax(('$$ \\eta = \\beta_0 + \\beta_1 \\cdot x + \\epsilon$$')),
      " but covariate", em("x"), "is not directly observable. Instead, a substitute"),
      "$$w = x + u$$",
      "is observed, assuming that" ,
      "$$u \\sim N(0,\\sigma_u^2) .$$",
     "To check what happens when the error increases, simply move the slider below to increase the error variance in the x covariate.",
      br(),


  sliderInput("errVar", label = ("$$ \\text{Select an error variance } \\sigma_u^2 :$$"),
                                 min = 0, max = 8, value = 1,step=0.5)
  
 # textInput("text", label = "Choose a title for your plot",  value = "Error in linear regression")
  ), 
    
    mainPanel(
      plotOutput("plot",width="400pt",height="400pt"),
      textOutput("slope"),
      textOutput("resVar")
    )

  )
))

# how to publish app on the server:
# library(rsconnect)
# rsconnect::deployApp('/home/steffi/Shiny/MEApp/MEClassical/')

#       selectInput("var", 
#                   label = "Choose a variable to display",
#                   choices = c("Percent White", "Percent Black",
#                               "Percent Hispanic", "Percent Asian"),
#                   selected = "Percent White"),


#       numericInput("num", 
#                    label = h3("Select an error variance for w"), 
#                    value = 1)

# 
# 
# fluidRow(
#   
#   column(3, # width of column
#          h3("Buttons"),
#          actionButton("action", label = "Action!!"),
#          br(),
#          br(), 
#          submitButton("Submit")),
#   
#   column(3,
#          h3("Single checkbox"),
#          checkboxInput("checkbox", label = "Choice A", value = TRUE)),
#   
#   column(3, 
#          checkboxGroupInput("checkGroup", 
#                             label = h3("Checkbox group"), 
#                             choices = list("Choice 1" = 1, 
#                                            "Choice 2" = 2, "Choice 3" = 3),
#                             selected = 1)),
#   
#   column(3, 
#          dateInput("date", 
#                    label = h3("Date input"), 
#                    value = "2014-01-01"))   
# ),
# 
# fluidRow(
#   
#   column(3,
#          dateRangeInput("dates", label = h3("Date range"))),
#   
#   column(3,
#          fileInput("file", label = h3("File input"))),
#   
#   column(3, 
#          h3("Help text"),
#          helpText("Note: help text isn't a true widget,", 
#                   "but it provides an easy way to add text to",
#                   "accompany other widgets.")),
#   
#   column(3, 
#          numericInput("num", 
#                       label = h3("Numeric input"), 
#                       value = 1))   
# ),
# 
# fluidRow(
#   
#   column(3,
#          radioButtons("radio", label = h3("Radio buttons"),
#                       choices = list("Choice 1" = 1, "Choice 2" = 2,
#                                      "Choice 3" = 3),selected = 1)),
#   
#   column(3,
#          selectInput("select", label = h3("Select box"), 
#                      choices = list("Choice 1" = 1, "Choice 2" = 2,
#                                     "Choice 3" = 3), selected = 1)),
#   
#   column(3, 
#          sliderInput("slider1", label = h3("Sliders"),
#                      min = 0, max = 100, value = 50),
#          sliderInput("slider2", "",
#                      min = 0, max = 100, value = c(25, 75)),
#          sliderInput("slider3", "",
#                      min = 0, max = 100, value = c(50,75))
#   ),
#   
#   column(3, 
#          textInput("text", label = h3("Text input"), 
#                    value = "Enter text..."))   
# )