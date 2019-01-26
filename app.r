library(ggplot2)
library(shiny)
library(reshape2)
library(wesanderson)
library(shinythemes)
library(DT)


stat_smooth_func_with_pval <- function(mapping = NULL, data = NULL,
                                       geom = "smooth", position = "identity",
                                       ...,
                                       method = "auto",
                                       formula = y ~ x,
                                       se = TRUE,
                                       n = 80,
                                       span = 0.75,
                                       fullrange = FALSE,
                                       level = 0.95,
                                       method.args = list(),
                                       na.rm = FALSE,
                                       show.legend = NA,
                                       inherit.aes = TRUE,
                                       xpos = NULL,
                                       ypos = NULL,
                                       xpos2 = NULL,
                                       ypos2 = NULL) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatSmoothFunc,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      method = method,
      formula = formula,
      se = se,
      n = n,
      fullrange = fullrange,
      level = level,
      na.rm = na.rm,
      method.args = method.args,
      span = span,
      xpos = xpos,
      ypos = ypos,
      xpos2 = xpos2,
      ypos2 = ypos2,
      ...
    )
  )
}

StatSmoothFunc <- ggproto("StatSmooth", Stat,
                          
                          setup_params = function(data, params) {
                            # Figure out what type of smoothing to do: loess for small datasets,
                            # gam with a cubic regression basis for large data
                            # This is based on the size of the _largest_ group.
                            if (identical(params$method, "auto")) {
                              max_group <- max(table(data$group))
                              
                              if (max_group < 1000) {
                                params$method <- "loess"
                              } else {
                                params$method <- "gam"
                                params$formula <- y ~ s(x, bs = "cs")
                              }
                            }
                            if (identical(params$method, "gam")) {
                              params$method <- mgcv::gam
                            }
                            
                            params
                          },
                          
                          compute_group = function(data, scales, method = "auto", formula = y~x,
                                                   se = TRUE, n = 80, span = 0.75, fullrange = FALSE,
                                                   xseq = NULL, level = 0.95, method.args = list(),
                                                   na.rm = FALSE, xpos=NULL, ypos=NULL, 
                                                   xpos2=NULL, ypos2=NULL) {
                            if (length(unique(data$x)) < 2) {
                              # Not enough data to perform fit
                              return(data.frame())
                            }
                            
                            if (is.null(data$weight)) data$weight <- 1
                            
                            if (is.null(xseq)) {
                              if (is.integer(data$x)) {
                                if (fullrange) {
                                  xseq <- scales$x$dimension()
                                } else {
                                  xseq <- sort(unique(data$x))
                                }
                              } else {
                                if (fullrange) {
                                  range <- scales$x$dimension()
                                } else {
                                  range <- range(data$x, na.rm = TRUE)
                                }
                                xseq <- seq(range[1], range[2], length.out = n)
                              }
                            }
                            # Special case span because it's the most commonly used model argument
                            if (identical(method, "loess")) {
                              method.args$span <- span
                            }
                            
                            if (is.character(method)) method <- match.fun(method)
                            
                            base.args <- list(quote(formula), data = quote(data), weights = quote(weight))
                            model <- do.call(method, c(base.args, method.args))
                            
                            m = model
                            eq1 <- substitute(italic(y) == a + b %.% italic(x), 
                                              list(a = format(coef(m)[[1]], digits = 3), 
                                                   b = format(coef(m)[[2]], digits = 3)))
                            
                            eq2 <- substitute(italic(r)^2~"="~r2*","~~italic(p)~"="~pval, 
                                              list(r2 = format(summary(m)$r.squared, digits = 3),
                                                   pval = format(summary(m)$coef[2,4], digits = 3)))
                            
                            func_string1 = as.character(as.expression(eq1))
                            func_string2 = as.character(as.expression(eq2))
                            
                            if(is.null(xpos)) xpos = min(data$x)
                            if(is.null(ypos)) ypos = max(data$y)
                            if(is.null(xpos2)) xpos2 = xpos
                            if(is.null(ypos2)) ypos2 = ypos - (ypos*0.02)
                            
                            data.frame(x = rbind(xpos, xpos2), 
                                       y = rbind(ypos, ypos2), 
                                       label = rbind(func_string1, func_string2))
                            
                          },
                          
                          required_aes = c("x", "y")
)

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# Read CSV into R, it can take a while.
# Change csv paths to file paths if nessesery, but remeber that you need to use "\" separators instead of "/". 
# You can download csv files for yourself on: https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/
WhiteWines <- read.csv(file=url("https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-white.csv"), header=TRUE, sep=";")
RedWines <- read.csv(file=url("https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-red.csv"), header=TRUE, sep=";")

ui <- fluidPage(
  theme = shinytheme("sandstone"),
  
  navbarPage("WineData",
             tabPanel("Plots created from data",
                      h1("Wine data Plots", 
                         align="center"),
                      
                      sidebarLayout(
                        
                        position = "right",
                        
                        sidebarPanel(
                          
                          checkboxInput("checkboxTheme", 
                                        label = "Dark theme", 
                                        value = FALSE),
                          
                          radioButtons("radioPlotStyle", 
                                       inline = TRUE,
                                       label = NA,
                                       choices = list("Jiter", 
                                                      "Point",
                                                      "Smooth",
                                                      "Column",
                                                      "Violin", 
                                                      "Boxplot"), 
                                       selected = "Jiter"),
                          
                          radioButtons("radioPlot", 
                                       inline = TRUE, 
                                       label = NA,
                                       choices = list("Normal plot", 
                                                      "Plot for each Quality")),
                          
                          selectInput("selectX", 
                                      label = h3("Select X and Y for plot drawing"), 
                                      choices = list("Fixed acidity",
                                                     "Volatile acidity",
                                                     "Citric acid",
                                                     "Residual sugar",
                                                     "Chlorides",
                                                     "Free sulfur dioxide",
                                                     "Total sulfur dioxide",
                                                     "Density", 
                                                     "pH",
                                                     "Sulphates",
                                                     "Alcohol",
                                                     "Quality"), 
                                      selected = "Quality"),
                          
                          selectInput("selectY", 
                                      label = NA, 
                                      choices = list("Fixed acidity",
                                                     "Volatile acidity",
                                                     "Citric acid",
                                                     "Residual sugar",
                                                     "Chlorides",
                                                     "Free sulfur dioxide",
                                                     "Total sulfur dioxide",
                                                     "Density", 
                                                     "pH",
                                                     "Sulphates",
                                                     "Alcohol",
                                                     "Quality")),
                          
                          p("To zoom, please select a zone you want to zoom and double click it.", 
                            align = "center"),
                          p("Double click without a zone to return to normal stage.", 
                            align = "center"),
                          
                          br(),
                          p("You can download data for yourself ", a("here.", href="https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/"), 
                            align = "center")
                          
                        ),
                        
                        mainPanel(
                          tabsetPanel(
                            tabPanel("White Wines", 
                                     plotOutput(outputId = "whitePlot", 
                                                height = "900px",
                                                dblclick = "PlotDblclick",
                                                brush = brushOpts(
                                                id = "PlotBrush",
                                                resetOnNew = TRUE)
                                                ), 
                                     h4("Selected points", 
                                        align = "center"),
                                     verbatimTextOutput("plotWhiteBrushInfo"), 
                                     p("Note that only in Point mode, selecting is precise.", 
                                       align = "center"),
                                     tags$head(tags$style(HTML("
                                    #plotWhiteBrushInfo {
                                                               font-size: 11.5px;
                                                               }
                                                               ")))),
                            tabPanel("Red Wines", 
                                     plotOutput(outputId = "redPlot", 
                                                height = "900px",
                                                dblclick = "PlotDblclick",
                                                brush = brushOpts(
                                                id = "PlotBrush",
                                                resetOnNew = TRUE)
                                                ),
                                     h4("Selected points", 
                                        align = "center"),
                                     verbatimTextOutput("plotRedBrushInfo"),
                                     p("Note that only in Point mode, selecting is precise.", 
                                       align = "center"),
                                    tags$head(tags$style(HTML("
                                    #plotRedBrushInfo {
                                                              font-size: 11.5px;
                                                              }
                                                              ")))), 
                            tabPanel("Correlation plot for White wines", 
                                     plotOutput(outputId = "colWhitePlot")),
                            tabPanel("Correlation plot for Red wines", 
                                     plotOutput(outputId = "colRedPlot"))
                          )
                        )
                      )),
             
             tabPanel("Info about data", 
                      h1("Wine data Info", 
                         align="center"),
                      h3("Informations are from ", 
                         a("winequality.names.", 
                           href="https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality.names"), 
                         align = "center"),
                      htmlOutput("modelParameters")
                      ),
             
             tabPanel("Data for view", 
                      tabsetPanel(
                        tabPanel("White Wines", 
                                 br(),
                                 DT::dataTableOutput("dataToShowWhite")),
                        tabPanel("Red Wines", 
                                 br(),
                                 DT::dataTableOutput("dataToShowRed"))
                                )
                      )
             )
  
 
  )


server <- function(input, output, session) {
  
  observe({
    choicesRadioPlotStyle = if(input$selectX == "Quality"){
      list("Jiter", 
           "Point",
           "Column",
           "Violin", 
           "Boxplot")
    }
    else{
      list("Jiter", 
           "Point",
           "Column",
           "Smooth")
    }
    
    
    updateRadioButtons(session, "radioPlotStyle",
                       choices = choicesRadioPlotStyle,
                       inline = TRUE)
    
  })
  
  observe({
    choices = list("Fixed acidity",
                   "Volatile acidity",
                   "Citric acid",
                   "Residual sugar",
                   "Chlorides",
                   "Free sulfur dioxide",
                   "Total sulfur dioxide",
                   "Density", 
                   "pH",
                   "Sulphates",
                   "Alcohol",
                   "Quality")

    updateSelectInput(session, "selectY",
                      choices = choices[choices != input$selectX],
                      selected = input$selectY)
    
    updateSelectInput(session, "selectX",
                      choices = choices[choices != input$selectY],
                      selected = input$selectX)
  })
  
  output$dataToShowWhite = DT::renderDataTable({WhiteWines})
  
  output$dataToShowRed = DT::renderDataTable({RedWines})
  
  output$plotWhiteBrushInfo <- renderPrint({
    
    fixedNamesY <- switch(input$selectY, 
                          "Fixed acidity" = "fixed.acidity", 
                          "Volatile acidity" = "volatile.acidity",
                          "Citric acid" = "citric.acid",
                          "Residual sugar" = "residual.sugar", 
                          "Chlorides" = "chlorides",
                          "Free sulfur dioxide" = "free.sulfur.dioxide", 
                          "Total sulfur dioxide" = "total.sulfur.dioxide",
                          "Density" = "density", 
                          "pH" = "pH",
                          "Sulphates" = "sulphates",
                          "Alcohol" = "alcohol",
                          "Quality" = "quality")
    
    fixedNamesX <- switch(input$selectX, 
                          "Fixed acidity" = "fixed.acidity", 
                          "Volatile acidity" = "volatile.acidity",
                          "Citric acid" = "citric.acid",
                          "Residual sugar" = "residual.sugar", 
                          "Chlorides" = "chlorides",
                          "Free sulfur dioxide" = "free.sulfur.dioxide", 
                          "Total sulfur dioxide" = "total.sulfur.dioxide",
                          "Density" = "density", 
                          "pH" = "pH",
                          "Sulphates" = "sulphates",
                          "Alcohol" = "alcohol",
                          "Quality" = "quality")
    
    brushedPoints(WhiteWines, input$PlotBrush, fixedNamesX, fixedNamesY)
  })
  
  output$plotRedBrushInfo <- renderPrint({
    
    fixedNamesY <- switch(input$selectY, 
                          "Fixed acidity" = "fixed.acidity", 
                          "Volatile acidity" = "volatile.acidity",
                          "Citric acid" = "citric.acid",
                          "Residual sugar" = "residual.sugar", 
                          "Chlorides" = "chlorides",
                          "Free sulfur dioxide" = "free.sulfur.dioxide", 
                          "Total sulfur dioxide" = "total.sulfur.dioxide",
                          "Density" = "density", 
                          "pH" = "pH",
                          "Sulphates" = "sulphates",
                          "Alcohol" = "alcohol",
                          "Quality" = "quality")
    
    fixedNamesX <- switch(input$selectX, 
                          "Fixed acidity" = "fixed.acidity", 
                          "Volatile acidity" = "volatile.acidity",
                          "Citric acid" = "citric.acid",
                          "Residual sugar" = "residual.sugar", 
                          "Chlorides" = "chlorides",
                          "Free sulfur dioxide" = "free.sulfur.dioxide", 
                          "Total sulfur dioxide" = "total.sulfur.dioxide",
                          "Density" = "density", 
                          "pH" = "pH",
                          "Sulphates" = "sulphates",
                          "Alcohol" = "alcohol",
                          "Quality" = "quality")
    
    brushedPoints(WhiteWines, input$PlotBrush, fixedNamesX, fixedNamesY)
  })
  
  ranges <- reactiveValues(x = NULL, y = NULL)

  observeEvent(input$PlotDblclick, {
    brush <- input$PlotBrush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  
  output$brush_info <- renderPrint({
    brushedPoints(mtcars2, input$plot1_brush)
  })
  
  output$modelParameters <- renderUI(
    HTML(
      paste(
        c("<pre>", "Citation Request:
                        This dataset is public available for research. The details are described in [Cortez et al., 2009]. 
          Please include this citation if you plan to use this database:
          
          P. Cortez, A. Cerdeira, F. Almeida, T. Matos and J. Reis. 
          Modeling wine preferences by data mining from physicochemical properties.
          In Decision Support Systems, Elsevier, 47(4):547-553. ISSN: 0167-9236.
          
          Available at: [@Elsevier] http://dx.doi.org/10.1016/j.dss.2009.05.016
          [Pre-press (pdf)] http://www3.dsi.uminho.pt/pcortez/winequality09.pdf
          [bib] http://www3.dsi.uminho.pt/pcortez/dss09.bib
          
          1. Title: Wine Quality 
          
          2. Sources
          Created by: Paulo Cortez (Univ. Minho), Antonio Cerdeira, Fernando Almeida, Telmo Matos and Jose Reis (CVRVV) @ 2009
          
          3. Past Usage:
          
          P. Cortez, A. Cerdeira, F. Almeida, T. Matos and J. Reis. 
          Modeling wine preferences by data mining from physicochemical properties.
          In Decision Support Systems, Elsevier, 47(4):547-553. ISSN: 0167-9236.
          
          In the above reference, two datasets were created, using red and white wine samples.
          The inputs include objective tests (e.g. PH values) and the output is based on sensory data
          (median of at least 3 evaluations made by wine experts). Each expert graded the wine quality 
          between 0 (very bad) and 10 (very excellent). Several data mining methods were applied to model
          these datasets under a regression approach. The support vector machine model achieved the
          best results. Several metrics were computed: MAD, confusion matrix for a fixed error tolerance (T),
          etc. Also, we plot the relative importances of the input variables (as measured by a sensitivity
          analysis procedure).
          
          4. Relevant Information:
          
          The two datasets are related to red and white variants of the Portuguese Vinho Verde wine.
          For more details, consult: http://www.vinhoverde.pt/en/ or the reference [Cortez et al., 2009].
          Due to privacy and logistic issues, only physicochemical (inputs) and sensory (the output) variables 
          are available (e.g. there is no data about grape types, wine brand, wine selling price, etc.).
          
          These datasets can be viewed as classification or regression tasks.
          The classes are ordered and not balanced (e.g. there are munch more normal wines than
          excellent or poor ones). Outlier detection algorithms could be used to detect the few excellent
          or poor wines. Also, we are not sure if all input variables are relevant. So
          it could be interesting to test feature selection methods. 
          
          5. Number of Instances: red wine - 1599; white wine - 4898. 
          
          6. Number of Attributes: 11 + output attribute
          
          Note: several of the attributes may be correlated, thus it makes sense to apply some sort of
          feature selection.
          
          7. Attribute information:
          
          For more information, read [Cortez et al., 2009].
          
          Input variables (based on physicochemical tests):
          1 - fixed acidity
          2 - volatile acidity
          3 - citric acid
          4 - residual sugar
          5 - chlorides
          6 - free sulfur dioxide
          7 - total sulfur dioxide
          8 - density
          9 - pH
          10 - sulphates
          11 - alcohol
          Output variable (based on sensory data): 
          12 - quality (score between 0 and 10)
          
          8. Missing Attribute Values: None", "</pre>"),
        collapse = "<br>"
      )
    )
  )
  
  output$whitePlot <- renderPlot({
    
    selectY <- switch(input$selectY, 
                     "Fixed acidity" = WhiteWines$fixed.acidity, 
                     "Volatile acidity" = WhiteWines$volatile.acidity,
                     "Citric acid" = WhiteWines$citric.acid,
                     "Residual sugar" = WhiteWines$residual.sugar, 
                     "Chlorides" = WhiteWines$chlorides,
                     "Free sulfur dioxide" = WhiteWines$free.sulfur.dioxide, 
                     "Total sulfur dioxide" = WhiteWines$total.sulfur.dioxide,
                     "Density" = WhiteWines$density, 
                     "pH" = WhiteWines$pH,
                     "Sulphates" = WhiteWines$sulphates,
                     "Alcohol" = WhiteWines$alcohol,
                     "Quality" = WhiteWines$quality)
    
    selectX <- switch(input$selectX, 
                     "Fixed acidity" = WhiteWines$fixed.acidity, 
                     "Volatile acidity" = WhiteWines$volatile.acidity,
                     "Citric acid" = WhiteWines$citric.acid,
                     "Residual sugar" = WhiteWines$residual.sugar, 
                     "Chlorides" = WhiteWines$chlorides,
                     "Free sulfur dioxide" = WhiteWines$free.sulfur.dioxide, 
                     "Total sulfur dioxide" = WhiteWines$total.sulfur.dioxide,
                     "Density" = WhiteWines$density, 
                     "pH" = WhiteWines$pH,
                     "Sulphates" = WhiteWines$sulphates,
                     "Alcohol" = WhiteWines$alcohol,
                     "Quality" = WhiteWines$quality)
    
    
    fixedNamesY <- switch(input$selectY, 
                         "Fixed acidity" = "fixed.acidity", 
                         "Volatile acidity" = "volatile.acidity",
                         "Citric acid" = "citric.acid",
                         "Residual sugar" = "residual.sugar", 
                         "Chlorides" = "chlorides",
                         "Free sulfur dioxide" = "free.sulfur.dioxide", 
                         "Total sulfur dioxide" = "total.sulfur.dioxide",
                         "Density" = "density", 
                         "pH" = "pH",
                         "Sulphates" = "sulphates",
                         "Alcohol" = "alcohol",
                         "Quality" = "quality")
    
    fixedNamesX <- switch(input$selectX, 
                         "Fixed acidity" = "fixed.acidity", 
                         "Volatile acidity" = "volatile.acidity",
                         "Citric acid" = "citric.acid",
                         "Residual sugar" = "residual.sugar", 
                         "Chlorides" = "chlorides",
                         "Free sulfur dioxide" = "free.sulfur.dioxide", 
                         "Total sulfur dioxide" = "total.sulfur.dioxide",
                         "Density" = "density", 
                         "pH" = "pH",
                         "Sulphates" = "sulphates",
                         "Alcohol" = "alcohol",
                         "Quality" = "quality")
  
    normalPlot <- ggplot(WhiteWines, aes(x = selectX, color = quality,  y = selectY)) +
      scale_color_gradientn(colours = wes_palette("Zissou1", 100, type = "continuous"))+
      labs(x=input$selectX, y=input$selectY)+
      coord_cartesian(xlim = ranges$x, ylim = ranges$y)
    
    radioSwitchPlot <- switch(input$radioPlot,
                              "Normal plot" = normalPlot,
                              "Plot for each Quality" = normalPlot + facet_grid(rows = vars(quality))
    )
    
    if(input$checkboxTheme == TRUE){
      radioSwitchPlot <- radioSwitchPlot + theme_dark()
    }  
    
    radioSwitchPlot <- switch(input$radioPlotStyle, 
                              "Jiter" = radioSwitchPlot + geom_jitter(height = min(selectX)*0.05, width = min(selectX)*0.05) + stat_smooth_func_with_pval(geom="text",method="lm",hjust=0,parse=TRUE) + geom_smooth(method = "lm", se = FALSE), 
                              "Point" = radioSwitchPlot + geom_point() + stat_smooth_func_with_pval(geom="text",method="lm",hjust=0,parse=TRUE) + geom_smooth(method = "lm", se = FALSE),
                              "Smooth" = radioSwitchPlot + geom_smooth(),
                              "Column" = radioSwitchPlot + geom_col(),
                              "Violin" = radioSwitchPlot + geom_violin(scale = "area", aes(group = WhiteWines$quality)), 
                              "Boxplot" = radioSwitchPlot + geom_boxplot(aes(group = WhiteWines$quality)))

   radioSwitchPlot
    
  })
  
  output$redPlot <- renderPlot({
    
    selectY <- switch(input$selectY, 
                      "Fixed acidity" = RedWines$fixed.acidity, 
                      "Volatile acidity" = RedWines$volatile.acidity,
                      "Citric acid" = RedWines$citric.acid,
                      "Residual sugar" = RedWines$residual.sugar, 
                      "Chlorides" = RedWines$chlorides,
                      "Free sulfur dioxide" = RedWines$free.sulfur.dioxide, 
                      "Total sulfur dioxide" = RedWines$total.sulfur.dioxide,
                      "Density" = RedWines$density, 
                      "pH" = RedWines$pH,
                      "Sulphates" = RedWines$sulphates,
                      "Alcohol" = RedWines$alcohol,
                      "Quality" = RedWines$quality)
    
    selectX <- switch(input$selectX, 
                      "Fixed acidity" = RedWines$fixed.acidity, 
                      "Volatile acidity" = RedWines$volatile.acidity,
                      "Citric acid" = RedWines$citric.acid,
                      "Residual sugar" = RedWines$residual.sugar, 
                      "Chlorides" = RedWines$chlorides,
                      "Free sulfur dioxide" = RedWines$free.sulfur.dioxide, 
                      "Total sulfur dioxide" = RedWines$total.sulfur.dioxide,
                      "Density" = RedWines$density, 
                      "pH" = RedWines$pH,
                      "Sulphates" = RedWines$sulphates,
                      "Alcohol" = RedWines$alcohol,
                      "Quality" = RedWines$quality)
    
    
    fixedNamesY <- switch(input$selectY, 
                          "Fixed acidity" = "fixed.acidity", 
                          "Volatile acidity" = "volatile.acidity",
                          "Citric acid" = "citric.acid",
                          "Residual sugar" = "residual.sugar", 
                          "Chlorides" = "chlorides",
                          "Free sulfur dioxide" = "free.sulfur.dioxide", 
                          "Total sulfur dioxide" = "total.sulfur.dioxide",
                          "Density" = "density", 
                          "pH" = "pH",
                          "Sulphates" = "sulphates",
                          "Alcohol" = "alcohol",
                          "Quality" = "quality")
    
    fixedNamesX <- switch(input$selectX, 
                          "Fixed acidity" = "fixed.acidity", 
                          "Volatile acidity" = "volatile.acidity",
                          "Citric acid" = "citric.acid",
                          "Residual sugar" = "residual.sugar", 
                          "Chlorides" = "chlorides",
                          "Free sulfur dioxide" = "free.sulfur.dioxide", 
                          "Total sulfur dioxide" = "total.sulfur.dioxide",
                          "Density" = "density", 
                          "pH" = "pH",
                          "Sulphates" = "sulphates",
                          "Alcohol" = "alcohol",
                          "Quality" = "quality")
    
    normalPlot <- ggplot(RedWines, aes(x = selectX, color = quality,  y = selectY)) +
      scale_color_gradientn(colours = wes_palette("Zissou1", 100, type = "continuous"))+
      labs(x=input$selectX, y=input$selectY)+
      coord_cartesian(xlim = ranges$x, ylim = ranges$y)
    
    radioSwitchPlot <- switch(input$radioPlot,
                              "Normal plot" = normalPlot,
                              "Plot for each Quality" = normalPlot + facet_grid(rows = vars(quality))
    )
    
    if(input$checkboxTheme == TRUE){
      radioSwitchPlot <- radioSwitchPlot + theme_dark()
    }  
    
    radioSwitchPlot <- switch(input$radioPlotStyle, 
                              "Jiter" = radioSwitchPlot + geom_jitter(height = min(selectX)*0.05, width = min(selectX)*0.05) + stat_smooth_func_with_pval(geom="text",method="lm",hjust=0,parse=TRUE) + geom_smooth(method = "lm", se = FALSE), 
                              "Point" = radioSwitchPlot + geom_point() + stat_smooth_func_with_pval(geom="text",method="lm",hjust=0,parse=TRUE) + geom_smooth(method = "lm", se = FALSE),
                              "Smooth" = radioSwitchPlot + geom_smooth(),
                              "Column" = radioSwitchPlot + geom_col(),
                              "Violin" = radioSwitchPlot + geom_violin(scale = "area", aes(group = RedWines$quality)), 
                              "Boxplot" = radioSwitchPlot + geom_boxplot(aes(group = RedWines$quality)))
    
    radioSwitchPlot
    
  })
  
  output$colWhitePlot <- renderPlot({
    
    cormat <- round(cor(WhiteWines),2)
    cormat <- reorder_cormat(cormat)
    upper_tri <- get_upper_tri(cormat)
    melted_cormat <- melt(upper_tri, na.rm = TRUE)
    ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
      geom_tile(color = "white") +
      scale_fill_gradientn(colours = wes_palette("Zissou1", 100, type = "continuous"), limit = c(-1,1), space = "Lab", 
                           name="Pearson\nCorrelation") +
      theme_minimal() + 
      geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.6, 0.7),
        legend.direction = "horizontal") +
      guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                   title.position = "top", title.hjust = 0.5)) +
      coord_fixed()
    
    print(ggheatmap)
  })
  
  output$colRedPlot <- renderPlot({
    
    cormat <- round(cor(RedWines),2)
    cormat <- reorder_cormat(cormat)
    upper_tri <- get_upper_tri(cormat)
    melted_cormat <- melt(upper_tri, na.rm = TRUE)
    
    ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
      geom_tile(color = "white") +
      scale_fill_gradientn(colours = wes_palette("Zissou1", 100, type = "continuous"), limit = c(-1,1), space = "Lab", 
                           name="Pearson\nCorrelation") +
      theme_minimal() + 
      geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.6, 0.7),
        legend.direction = "horizontal") +
      guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                   title.position = "top", title.hjust = 0.5)) +
      coord_fixed()
    
    print(ggheatmap)
  })
  
}

shinyApp(ui = ui, server = server)